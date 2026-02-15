# Project_Dual_inhibitor
All Files Associated with the Dual Inhibitor Design for 14-3-/Synuclein
https://doi.org/10.1002/adts.202501455

#Supervised_ML_code  < br / >

Contains two python scripts namely:  < br / >

01_distance_matrix.py : it reads as input .prmtop and .nc files generates two numpy arrays, which constitutes the 1/distance matrix along with residue pair information.  < br / >

02_ML_prediction.ipynb : this uses the 4 np arrays for the 4 systems, viz., APO, ORT, ALO, and DUO, along with another np array for residue pair information and generates important residue pairs based on LR, RF, and MLP .< br / >

	
