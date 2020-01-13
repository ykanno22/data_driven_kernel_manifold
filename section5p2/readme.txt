Copyright (C) Yoshihiro Kanno, 2020. All rights reserved.
===============================================================================
Y. Kanno: "A kernel method for learning constitutive relation in data-driven computational elasticity." Submitted to Japan Journal of Industrial and Applied Mathematics.

The codes can run on Matlab ver.9.7 with Optimization Toolbox. 
Also, the following Matlab codes available at 
   <http://au.mathworks.com/support/books/book118765.html>
are required for FEM. 
  - elem_q4_fun.m
  - fmlin.m
  - formbee.m
  - formdsig.m
  - gauss.m
  - Q4_mesh_fun.m
These FEM codes are contained in the following book:
  - A. Khennane:
    Introduction to Finite Element Analysis Using MATLAB and Abaqus.
    CRC Press, Boca Raton (2013).

To reproduce the result presented in section 5.2, run the following sequentially: 
  - nonlin_6d_eig.m
  - loop_fem_nonlin_6d_eig.m
Part of these codes and the functions used by them was based on [Khennane, 2013] cited above.
