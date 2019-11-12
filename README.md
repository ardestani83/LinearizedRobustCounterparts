This repository includes data sets and codes used in the following paper:

Ardestani-Jaafari, Amir, and Erick Delage. Linearized robust counterparts of two-stage robust optimization problems with applications in operations management. GERAD, École des Hautes études commerciales, 2016.

It contains three folders: "Section5" includes the Matlab code used for the numerical study of Section 5 focusing on the robust location-transportation problem. "Section6" includes the Matlab code used for the numerical study of Section 6 focusing on a robust newsvendor problem. Finally, "Appendices" includes the  Matlab code used for the numerical study of all examples presented in the Appendix section of the paper.

Implementing the repository requires the following software:

Matlab: it is the main platform of this repository (see https://www.mathworks.com/products/matlab.html)

YALMIP: it is used as a mathematical programming language (more detail: https://yalmip.github.io/)

ROME : it is used to implement some of the robust optimization models (see: https://robustopt.com/)

*: There is an error in using ROME in the recent versions of matlab. It can be fixed as replacing rome_var.m file in ROME package with the availabel rome_var.m file in the repository.

CPLEX: it is required for solving linear programming and mixed-integer linear programming (see http://cplex.com)

MOSEK: it is required for solving for semi-definite programming (see https://www.mosek.com/)
