# Combining experiments and in silico modeling to infer the role of adhesion and proliferation on the collective dynamics of cells


### This is the repository for all data and codes used in the paper [Combining experiments and in silico modeling to infer the role of adhesion and proliferation on the collective dynamics of cells, H. P. M. Melo, F. R. Maia, A. S. Nunes, R. L. Reis, J. M. Oliveira, N. A. M. Araújo. Scientific Reports 11, 19894 (2021)](https://www.nature.com/articles/s41598-021-99390-x)

## If you use any data or code, please considerer to cite. 

1. The folder **\images** contains all experimental images of GBM cells proliferation on a flat substrate. The file name contains the time points (2h, 4h, 6h, 8h, 12h and 24h) and the sample number;
1. The file **Nucleus_Segmentation.ipynb** is a Python code that reads all images at the folder **\images** and finds the location of cell nuclei;
2. The file **Mol_Dynamics_Model.cpp** is a C++ code used to simulate the model described at the [paper](https://www.nature.com/articles/s41598-021-99390-x);
3. In order to run the Molecular Dynamics simulation you first need to give a set of parameters using the file **data.txt**:  
    - First line: proliferation rate <img src="https://render.githubusercontent.com/render/math?math=\lambda">;
    - Second line: intensity of adhesion <img src="https://render.githubusercontent.com/render/math?math=\tau">; 
    - Third line: seed used for the random number generator. 


