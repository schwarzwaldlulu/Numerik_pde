# Numerik_pde
code from the exercise in the course Numerik für patielle DIfferentialgleichungen ( numerical methodes for pde )

the whole project consists of several header files and main file to solve different pdes based on a given grid, where the mesh file contains coordinates of each
vertex on one line and indices of the corners of the triangles

cg.hh: methode to solve linear equation system, the algorithm is described in wikipedia article

globalrefine.hh: refine the grid i.e. add the middle points of every edge as vertices and cut every triangle to four new ones

assemble_linear_equation.hh: methodes to fill the stiffness matrix, the mass matrix and the right-hand side of the linear equation to find the approximate to the 
weak solution of pde

grid.hh: define the struct Vertex, Triangle and class Grid with all necessary data memeber and member functions

matrix.hh & sparsematrix.hh & position.hh & vector.hh: user-defined data structures to falicitate the calculation 
