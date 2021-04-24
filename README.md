# Dscribe-for-Machine-Learning
#Eigen value for Coulomb Matrix
# prerequisite: ASE, Pymatgen, DScribe descriptors
# this python script reads data from h5py file.
# it reads the position coordinates and energies of an atomic cluster(55 atoms) at different temperature
# the data file have 40001 different configuration
# the script first reads the data file, then convert it into .xyz format
# then it calculates Coulomb Matrix (CM) descriptors from .xyz format
# Then it find eigenvalues of CMs at different configurations and writes in .csv format
#Then it is ready to  use in machine learning model.
