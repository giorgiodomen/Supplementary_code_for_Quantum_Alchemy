The folder "APDFT 14e diatomics" contains a Dataset Supplementary to the article "Effects of perturbation order and basis set on alchemical predictions" by Giorgio Domenichini, Guido Falk von Rudorff and O.Anatole von Lilienfeld.
The data stored in  '.\data:' as .csv files contain all the predictions within the 14 electron diatomics series, those are refered to the sections IIIA-IIIE of the article.
To every file correspond one of the eigth basis set used. 
The data are stored using consistently atomic units: Bohrs for length, electrons for charge and Hartrees for energy.
For a better visualization using Jupyter notebooks and Pandas library the notebook "Display_data.ipynb" already contains the code to open the data as pandas' data frames.

The folder "Alchemical CPHF perturbator" contains the subroutines for the CPHF evaluation of alchemical derivatives and alchemical forces.
The code content is made available to the public audience, free of charge under the MIT License.
The code works using PySCF v. 1.7.6 - accessed February 2021 (https://pyscf.org/), no support for later versions is guaranteed.
This project can be considered a research code supplementary to the paper "Alchemical predictions of relaxed geometries throughout chemical space" by Giorgio Domenichini and O.Anatole von Lilienfeld.
