# Finite Element Analysis for Modal and Transient Analysis of Beam in Matlab

The Matlab Partial Differential Equations (PDE) toolbox provides a powerful
and flexible setup for Finite Element Analysis (FEA) in multi-physics problems.
However, two factors limit its application specially in FEA models with large
Degrees of Freedom (DOFs).

The first comes from its memory storage scheme as it stores highly sparse global matrices such as mass, damping and stiffness in full square 2D arrays. Regarding the fact that indirect (iterative) solvers could be more efficient in solving large scale linear systems, the full matrix storage model needs many unnecessary Floating Point Operations (FLOPs) on zero components that increase computational cost.

The second factor is its memory usage in time domain (transient) analysis where it stores time-dependent data for all DOFs at all time steps. In practice, however, the data at some predefined points are interested. Accordingly, even a medium-sized FEA model can easily go out of memory for small time steps.

In this work FEA is implemented from scratch for modal and transient analysis for a beam. Although the full matrices are used for storing data, it allows to easily to replace them with Matlab's built-in sparse matrices. Moreover, the Newmark-beta method is used for numerical integartion of equation in transient analysis for taking full control over time steps at different time intervals. 

The test example is here a beam under point load. However, every geometric shape may be chosen as input. It only needed to store the geometry in STL file and specify suitbale boundary/initial conditions and different load combinations. The Matlab's built-in meshing algorithms are used for generating tetrahedons. The other parts of the FEA are coded that where it can be optimized by techniques like parallelizing matrix operations.

