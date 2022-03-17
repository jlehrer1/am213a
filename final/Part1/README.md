# Final Project: How to Run

This folder contains the programming work for the AM213A final with Dongwook Lee. 

`LinAl.f90` contains all numerical routines. `parta.f90` contains the code for the image compression via the SVD composition, while `visualize.py` generates the images from the raw matrices as well as the error plot. `partb_1.f90` contains the Gauss-Jacobi and Gauss-Seidel methods for solving `Ax=b`, and is run with the matrix given in the project spec. Unlike the project spec, I do not request which iterative method should be run, as Ian May said this is okay to change. Finally, `partb_2.f90` contains the driver code for solving `Ax=b` via the Conjugate Gradient algorithm.

To build and run the three files:
Part a: `make parta && ./parta`  
Part b.1: `make bone && ./bone`  
Part b.2: `make btwo && ./btwo` 

To generate the visualizations for all three parts, run the `make`s and then run `python visualize.py`, or `python visualize.py --force` to regenerate the plots even if they exist.  
