/*
Code used in Melo, HPM, et al. 
"Combining experiments and in silico modeling to infer the role of 
adhesion and proliferation on the collective dynamics of cells." bioRxiv (2021).

Please, consider to cite this paper in case you find this code useful. 

01-2021
Hygor Piaget
*/

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <string>
# include <sstream>
# include <ctime>
# include <vector>
# include <omp.h>
# include "boostrng.h"
# include "boostrng.cpp"
#define PI  3.141592653589793 


double PI2 = 2. * PI;

using namespace std;


double lambda;    //proliferation rate
double D=5000.0;         //diffusion coefficient
double raio=10.0;        //nuclei radius
double tau;         // adhesion intensity

int seed=797797515;



int Lx=1388; //system size X
int Ly=1040; //system size Y
int N_max=300; //this is the number of cells where you want to end the simulation and print the K-function //largest value of the number of cells

double r=0.0;
double dt=0.0001;
int t_final=100000000;



vector<double> x(N_max,0);
vector<double> y(N_max,0);
vector<double> generation(N_max,0);

double Dist[600][600]; //This matrix will keep the distance between pair of cells. Use a pair of number larger than N_max


void initial_conditions(int N, BoostRNG *rand);
void diffusion_past(int N, BoostRNG *rand);
int reproduction(int N, double lambda , int index, BoostRNG *rand);
int superposition(int N,double xx, double yy, vector <double> x_old, vector <double> y_old, int part_index);

double Ripley_K(int N, int Lx, int Ly, double r, double Dist[][600]);
double area_inside_box(double x0,double y0,double r,int Lx,int Ly);

int main(){



	ifstream infofile;
	infofile.open("data.txt", ios::in);

	char temp[100];
	infofile.getline(temp,100);
	lambda=strtof(temp,0);
	infofile.getline(temp,100);
	tau=strtof(temp,0);
	infofile.getline(temp,100);
	seed=strtof(temp,0);
	infofile.close();


	BoostRNG rand (seed);
	char filename[91]; 
	double K;



// The next two loops are only necessary in case you want to simulate for a set of  parameters values at the same run, 
// without having to call this code multiple times and using the data.txt!!!


	for(int b = 0; b < 1; b++){ 
		for(int l = 0; l < 1; l++){


			for(int sample = 0; sample < 10; sample++){

				snprintf(filename, sizeof(filename), "K_%f_%f_%d.dat", tau,lambda,sample);  // at this file you will find the K-function given tau, lambda and N_max
				std::ofstream outfile(filename);

				//N_max=109; // this is the number of cells where you want to end the simulation and print the K-function

				int N=50;   // this is your initial number of cells at the sample
				int N_old=50;

				initial_conditions(N,&rand);

				for(int t = 0; t < 20000000; t++){
					if(N>=N_max){
						for(int i=0; i<N; i++){
							for(int j=0; j<N; j++){ 
								Dist[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
							}	
						}
						for(int i=0; i<50; i++){
							r=10.+i*(100.-10.)/49.;
							K=Ripley_K(N, Lx, Ly, r, Dist);			
							outfile<<r<<" "<<K/(r*r)<<endl;
						}

						break;
					}
					
					diffusion_past(N,&rand);



					N_old=N;
					for(int i = 0; i < N_old; i++){
						if(N<N_max){if(reproduction(N,lambda,i,&rand)==1){N=N+1;}}
					}

				}
			}
		}
	}
}


void initial_conditions(int N, BoostRNG *rand){
// The function creates the initial random distribution of cells

	for(int i = 0; i < N; i++){
		int bloqueado=1;
		while(bloqueado==1){
			double r1 = rand->get_number();
			double r2 = rand->get_number();
			double xx = Lx*r1;
			double yy = Ly*r2;
			if(superposition(N,xx,yy,x,y,-1)==0){
				x[i]=xx;
				y[i]=yy;
				generation[i]=0;
				bloqueado=0;
			}
		}
	}
}



void diffusion_past(int N, BoostRNG *rand){
// this function updates the position of all cells in the system using Brownian Dynamics and taking into account the adhesion interaction.
	double xx;
	double yy;

	vector<double> x_dt(N,0);
	vector<double> y_dt(N,0);

	for(int i = 0; i < N; i++){
		double rr1 = rand->get_number();
		double rr2 = rand->get_number();
		double rx;
		double rr3 = rand->get_number();
		double rr4 = rand->get_number();
		double ry;

	 	rx = sqrt(-2 * log(rr1)) * cos(PI2 * rr2); //Box Muller method to generate Gaussian numbers
	 	ry = sqrt(-2 * log(rr3)) * cos(PI2 * rr4);



		if(superposition(N,x[i],y[i],x,y,i)==0){
			x_dt[i]=x[i]+sqrt(2*D*dt)*rx;
			y_dt[i]=y[i]+sqrt(2*D*dt)*ry;
		}
		else{
			x_dt[i]=x[i]+sqrt(2*tau*D*dt)*rx; 
			y_dt[i]=y[i]+sqrt(2*tau*D*dt)*ry; 

		}		

		if(x_dt[i] >= Lx){x_dt[i] -= Lx;}		
		if(x_dt[i] < 0. ){x_dt[i] += Lx;}
		if(y_dt[i] >= Ly){y_dt[i] -= Ly;}
		if(y_dt[i] < 0. ){y_dt[i] += Ly;}

	}

	for(int i = 0; i < N; i++){
		x[i]=x_dt[i];
		y[i]=y_dt[i];
	}


}



int reproduction(int N, double lambda, int index, BoostRNG *rand){

// this function adds a new cell by the process of proliferation controlled by the parameter lambda

	double ang = rand->get_number();
	double r = rand->get_number();
	double xx,yy;
	double theta = ang*2*PI;


	if(r<lambda){
		xx = x[index]+2*raio*cos(theta);
		yy = y[index]+2*raio*sin(theta);

		if(xx >= Lx){xx -= Lx;}		
		if(xx < 0. ){xx += Lx;}
		if(yy >= Ly){yy -= Ly;}
		if(yy < 0. ){yy += Ly;}

		if(superposition(N,xx,yy,x,y,-1)==0){
			x[N]=xx;
			y[N]=yy;
			generation[N]=1;		


			return 1;
		}else{return 0;}
	}
	else{return 0;} 

}


int superposition(int N,double xx, double yy, vector <double> x_old, vector <double> y_old, int part_index){

// this function checks if two cells are overlapping  

	double d;


	for(int i = 0; i < N; i++){
		if(i!=part_index){
			d=sqrt((xx-x[i])*(xx-x_old[i])+(yy-y_old[i])*(yy-y[i]));
			if(d<=2*raio){
				return 1;
			}
		}
	}
	return 0;

}


double area_inside_box(double x0,double y0,double r,int Lx,int Ly){

// this function is used to calculate the K-function (spatial correlation). 
// To estimate the K-function, it is important to correct the border effects. 
// This function is used to calculate the area of a given circle inside a rectangular box.


	double dtheta=0.05; //how precise do you whant to be?
	double area = 0.0;
	double theta_var,aux_x,aux_y,d,Distancia;
    
	if ((r<x0)and(r<y0)and(r<(Lx-x0))and(r<(Ly-y0))){  //checking if the radius dont cross the edge
		area = PI*r*r;}
	else{
		theta_var=0.0;
		for(int i = 0; i < int(PI2/dtheta); i++){
			theta_var = theta_var+dtheta;

			aux_x=x0+r*cos(theta_var);
			aux_y=y0+r*sin(theta_var);
        
			if ((aux_x>=0.0 and aux_x<=Lx) and (aux_y>=0.0 and aux_y<=Ly)){
				area += dtheta*r*r/2; } //inside the box, then sum dtheta
			else{
				Distancia=sqrt((aux_x-x0)*(aux_x-x0)+(aux_y-y0)*(aux_y-y0));
				if ((aux_x>=0.0 and aux_x<=Lx) and (aux_y>0.0)){ d=Distancia*(Ly-y0)/(aux_y-y0);}
				if ((aux_x>=0.0 and aux_x<=Lx) and (aux_y<0.0)){ d=Distancia*(y0)/(y0-aux_y);}
				if (aux_x>Lx){d=Distancia*(Lx-x0)/(aux_x-x0);}
				if (aux_x<0){d=Distancia*x0/(x0-aux_x);}
				area += dtheta*d*d/2;
				}
		}
			  
	}

	return area;

}

double Ripley_K(int N, int Lx, int Ly, double r, double Dist[][600]){

// this function calculates the value of K-function given a distance r


	double density = N/(1.*Lx*Ly);
	double K=0.0;
	double S=0.0;

	for(int i = 0; i < N; i++){
		S=0.0;
		for(int j = 0; j < N; j++){
			if(Dist[i][j]<=r){ S=S+1;}
		}
		K=K+(S-1)*(PI*r*r)/area_inside_box(x[i],y[i],r,Lx,Ly);
    }        
	K = K/(density*N);



	return K;
}














