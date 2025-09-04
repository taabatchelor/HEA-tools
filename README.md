Python classes for the linear regression model
The model is implemented using a collection of python classes found in the "/model/" folder.

The fingerprint generation,  linear regression algorithm and predicted binding energy calculations are implemented using the class found in "LinearRegression.py". Descriptions of the different functions are detailed in this python class.

Details of the training and test slabs are extracted from an atoms object using the class found in "Slab.py". Descriptions of the different functions are detailed in this python class.

Methods used to aid different calculations made in these classes are found in "helperMethods.py".

Functions to establish a matplotlib plot style are found in "plot.py"

The Database
Training and test set data can be found in the folder "/DFT_calculations/", as well as python scripts for running the slab relaxations.

Key-value pairs: Description
Key words for train.db database:

    type: 'OH','O','slab'

    slabID: '0','1','2','3'...N where N is the total number (-1) of each slab for each type; 

                  for 'OH' N = 876 (877 slabs in total)

                  for 'O' N = 999 (1000 slabs in total)

                  for 'slab' N = 999 (1000 slabs in total)

Key words for the test.db database:

    type: 'OH','O','slab'

    slabID:'0','1','2','3'...N where N is the total number (-1) of each slab for each type; 

                  for 'OH' N = 75 (76 slabs in total)

                  for 'O' N = 35 (36 slabs in total)

                  for 'slab' N = 75 (76 slabs in total)

The slab IDs exist for consistency, so that the structure with ID for an OH adsorbed slab corresponds with the same structure in O adsorbed or adsorbant-less slabs with the same ID.

Data structure stored in .csv files

*Oxygen adsorbed predicted adsorption energy csv files have been omitted as they take up a lot of space*
The data in the csv file containing predicted adsorption energies (found in the "/pred_histogram/" folder) have the following format:

'fingerprint, adsorption energy'

The fingerprint is set up as 15 (for *OH) or 55 (for *O) integers describing the number of each element found in different atomic zones composing the binding site and nearest neighbors (see paper linked above for details)

e.g. for *OH adsorption

0,0,0,1,0,0,1,2,1,2,1,0,2,0,0,1.123

Here the first 15 integers are the fingerprint ([0,0,0,1,0] = site, [0,1,2,1,2] = surface nearest neighbors, [1,0,2,0,0] = subsurface nearest neighbors, 1.123 = adsorption energy in eV. The [x,x,x,x,x] shown here corresponds to the number of of atoms from each element [Pt,Pd,Ir,Rh,Ru]). 

The same goes for *O adsorption but with 55 integers, with the first 35 integers describing the 3-atom binding site ensemble (35 possible combinations with 5 elements) and the remaining 20 integers describing the 4 zones of nearest neighbors, similar as for *OH.

This is better described in the paper.

![alt text](https://github.com/taabatchelor/HEA-tools/blob/main/DFT_histogram/OH_DFT_histogram.png "DFT calculated *OH adsorption energies on IrPdPtRhRu")

![alt text](https://github.com/taabatchelor/HEA-tools/blob/main/DFT_histogram/O_DFT_histogram.png "DFT calculated *O adsorption energies on IrPdPtRhRu")



Scripts for plotting the previous two histograms can be found in "/DFT_histogram/" under the name "X_DFT_histogram.py" with corresponding data in "X_train.csv" (where X is "OH" or "O").



![alt text](https://github.com/taabatchelor/HEA-tools/blob/main/regression_evaluation/OH_DFT_vs_pred.png "Predicted vs. DFT *OH adsorption energies on IrPdPtRhRu")



![alt text](https://github.com/taabatchelor/HEA-tools/blob/main/regression_evaluation/O_DFT_vs_pred.png "Predicted vs. DFT *O adsorption energies on IrPdPtRhRu")


Scripts for plotting the previous two scatter plots can be found in "/regression_evaluation/" under the name "X_DFT_vs_pred.py" (where X is "OH" or "O") with corresponding data in "X_train_pred.csv" and "X_test_pred.csv".



![alt text](https://github.com/taabatchelor/HEA-tools/blob/main/pred_histogram/OH_pred_histogram.png "Full span of predicted *OH adsorption energies on IrPdPtRhRu")


![alt text](https://github.com/taabatchelor/HEA-tools/blob/main/pred_histogram/O_pred_histogram.png "Full span of predicted *O adsorption energies on IrPdPtRhRu")


Scripts for plotting the previous two histograms can be found in pred_histogram under the name "X_pred_histogram.py" (where X is "OH" or "O") with corresponding *OH data in "OH_all_slabs.csv" and *O data in "O_all_slabs_Y.csv" (where Y is 0,1,2,3....34. The *O data had to be split into separate csv files).

