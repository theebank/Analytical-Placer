# To run

Unzip the folder

CD into the root folder of the program and type the following command:

`make ./visualizer`

Once the executable file is created, run it like so:

`./visualizer <test circuit Num> <weight factor between fixed node and movables> <0 for Parts 1/2, 1 for Part 3.1, 2 for Part 3.2>`

Note: the window may take a bit to open depending on the chosen test file and the weight for part 2

## To adjust the Psi value

in analyticalPlacer.cpp:

line 333: `int iter = 2;` is the starting value for the iteration

line 335: `double factor = 1;` is the factor which to multiply the iteration function by

line 343: `psi = factor * (iter * iter);` is the formula to calculate psi

## To adjust the weight between the anchors and each cell

in analyticalPlacer.cpp:

line 707: `matrix_UMF matrixToSolve = createAnchoredMatrix(0);` adjust the parameter for createAnchoredMatrix to what you wish the weight to be. 


### If the project does not compile
Try deleting Suitesparse then following the instructions to compile again as described in: https://janders.eecg.utoronto.ca/1387_2024/ECE1387_lin_solver.html