# NFoldIP technical documentation
This program enables creating an instance of integer n-fold programming problem. In the initialization function the properties for solver can be chosen (see the description of the __init__ function below). Then the instance can be solved in O(n3) time using dynamic programming. And that is not the only way the problem could be solved. There is also an option to use the standard GLPK solver (see the function solve below) or an experimental solver which uses MILP and runs faster than the solver using dynamic programming. 

All of the important implemented functions are described below.

The code was written in Python/Cython as a SageMath's module using the package 4ti2 -- a software package for algebraic, geometric and combinatorial problems on linear spaces. Available at www.4ti2.de. 


## IMPLEMENTATION - functions:
The __init__ function - takes twelve arguments:

 1. self 
 2. A - diagonal matrix (vectors of machines when scheduling)
 3. D - upper matrix (matrix of jobs when scheduling) 
 4. n - (number of jobs when scheduling)
 5. b -right hand-side vector 
 6. l - lower bound 
 7. u - upper bound 
 8. w - objective function 
 9. verbose - the identifier of logging, standard system of
    loggers - notset, debug, info, warning, error, critical (see
    https://docs.python.org/2/library/logging.html ), not obligatory,
    default is set to error level
 10. graver_complexity - possible values are: “exact” - counts the value of graver complexity from the definition, “approximate” - returns an approximate value own integer value 
 11. current_solution - not obligatory, if no current solution is given, the algorithm counts the initial solution from an auxiliary, program, 
 12. experimental - boolean value, the false option uses a dynamic program, true for an experimental algorithm which uses MILP


> **Implementation:** 

> -- checking the validity of the given data by a function check_validity_of_data, following must hold:
 >> general size of matrices: 
 >![enter image description here](https://lh3.googleusercontent.com/qvX3mVPqOhcma_LnSqK4uZ-4olQnBYPfmI9X3TP7JBMmFyCNAHzjqJRqt7i7LzPGLi1u8XvFFXZd "general size of matrices")
>>together with the sizes of l/u bounds, x and b:![enter image description here](https://lh3.googleusercontent.com/dAUwdfJjT0pZwB1_Hk_XngEUHih6wT4i1BBg5oWTiX5pJvAOXeFQMfFApVp-OaGmCUhViqbvmMvc "picture of sizes")
 -- lower and upper bounds can’t be infinity
 -- and for every number l in the lower bound vector and for every number u in the lower bound vector must hold that l<=u (for l,u in the same position)

>-- the logging level + format, self arguments as self.A, self.D, self.t, self.r, self.s etc. are being initialized
-- if any initial feasible solution was given, it is being set to self.current_solution; if there isn’t any, the initial feasible solution will be computed by the find_init_feasible_solution function later
-- graver complexity is  being computed and set to self.graver_complexity by a function approximate_graver_complexity, or exact_graver_complexity, or only the value is being applied (according to the graver_complexity option)
-- ZE to self.ZE is being set through a function _construct_ZE (only if experimental == false)
-- checking the validity of the data (whether the sizes of matrices, l/u bounds etc. make sense

Two functions for computing the graver complexity of given matrix (at most one of these is chosen according to the tenth (graver_complexity) argument in the __init__ function):

**exact_graver_complexity** returns an exact value of the Graver complexity of the given data from the definition
>**Implementation**:

>  takes the maximum entry of a matrix which  is product of multiplication of D by the graver basis of the matrix A ( the graver base has been computed using 4ti2/set as an integer in the init function)

**approximate_graver_complexity** return an approximate value of the Graver complexity of the given matrices
>**Implementation**:


>-- computes graver complexity according to the following formula [2]
$$g(A)\leq p(r*|| D*GA||_{\infty})^{r} $$
where:
>> -- g(A) is the graver  complexity of matrix A
-- p is the number of elements in the graver basis of A
-- GA is the graver basis of A

Function for computing Z(E) - **construct_ZE**. Z(E) is the sum of at most Graver complexity elements of the matrix D. This function is also called from the __init__ function.

>Implementation:

etc...


## SOURCES:

[1] R. Hammecke, S. Onn and L. Romanchuk, “N-fold integer programming in cubic time,” Mathematical Programming, pp. 1-17, 2013. 

[2] S. Onn, Nonlinear discrete optimization, Zurich Lectures in Advanced Mathematics, European Mathematical Society, 2010. 

<!--stackedit_data:
eyJoaXN0b3J5IjpbNDIwODU3NjkxXX0=
-->