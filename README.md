# Combining experiments and in silico modeling to infer the role of adhesion and proliferation on the collective dynamics of cells


### This is the repository for all data and codes used in the paper [Melo, HPM, Maia, FR, Nunes, AS, Reis, RL, Oliveira, JM, and Araújo, NA.  "Combining experiments and in silico modeling to infer the role of adhesion and proliferation on the collective dynamics of cells." bioRxiv (2021)](https://www.biorxiv.org/content/10.1101/2021.03.29.437400v1.abstract)

* If you use any data or code, please considerer to cite. 

1. The file **Nucleus_Segmentation.ipynb** is a Python code that reads all images at the folder **\images** and finds the location of cells nuclei.
2. The file **Mol_Dynamics_Model.cpp** is a C++ code used to simulate the model described at the [paper](https://www.biorxiv.org/content/10.1101/2021.03.29.437400v1.abstract).
3. In order to run the Molecular Dynamics simulation you first need to give a set of parameters using the file **data.txt**:  
    - First line: proliferation rate <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">
    - Second line: tau, intensity of adhesion
    - Third line: seed, seed used for the random number generator 


