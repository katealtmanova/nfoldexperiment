# NFoldIP technical documentation
We provide a class to handle instances of n-fold integer programming.
Most configuration is done via the arguments of the ``__init__`` function (see below).
The majority of computational work is done in the ``_find_good_step`` function, which searches for a good augmenting step of a specified length.
The algorithm computes such a step for each step-length from a set of step-lengths Gamma and then uses the best of these to augment the current solution; it terminates when no further augmenting steps can be found.
The ``_find_good_step`` procedure is implemented using the dynamic programming approach of [1]; because of its poor performance it can also be solved as an ILP subproblem, which is solved using a MILP solver such as GLPK, Coin-OR or Gurobi; see the function ``solve``. 

All of the important implemented functions are described below.

The code was written in Python/Cython as a SageMath's module using the package 4ti2 -- a software package for algebraic, geometric and combinatorial problems on linear spaces. Available at www.4ti2.de. 


## IMPLEMENTATION - functions:
The ``__init__`` function - takes these arguments:

 1. ``self`` 
 2. ``A`` - diagonal matrix (vectors of machines when scheduling)
 3. ``D`` - upper matrix (matrix of jobs when scheduling) 
 4. ``n`` - (number of jobs when scheduling)
 5. ``b`` -right hand-side vector 
 6. ``l`` - lower bound 
 7. ``u`` - upper bound 
 8. ``w`` - objective function 
 9. ``verbose`` - the identifier of logging, standard system of
    loggers - ``notset``, ``debug``, ``info``, ``warning``, ``error`` (default), ``critical`` (see
    https://docs.python.org/2/library/logging.html ), optional.
 10. ``graver_complexity`` - possible values are: ``“exact”`` - computes the value of Graver complexity from the definition, ``“approximate”`` (default) - computes an upper bound using a formula, integer - user-input; turns the algorithm into a heuristic
 11. ``current_solution`` - optional, if no current solution is given, the algorithm computes the initial solution using an auxiliary instance, 
 12. ``experimental`` - optional, ``False`` (default) uses dynamic programming (slow), ``True`` uses a MILP solver to solve a subinstance equivalent to the DP, ``"ng1"`` uses a MILP solver to find augmenting steps with 1-norm bounded by ``graver_complexity``, ``"nginfty"`` uses a MILP solver to find augmenting steps with \infty-norm bounded by ``graver_complexity``.
 13. ``instancename`` - optional, a string which will be used for output file logging, defaults to ``"instancename"``
 14. ``gamma`` - optional, a string, ``"best"`` uses the "best step" augmentation strategy, ``"logarithmic"`` (default) uses the "approximate best step" augmentation strategy, ``"unit"`` uses the "any step" augmentation strategy,
 15. ``solver`` - optional, a string, determines which MILP solver to use if ``experimental is not False``, ``"GLPK"`` (default), ``"Coin"``, ``"Gurobi"``.

> **Implementation:** 

> 
> - checking the validity of the given data by a function check_validity_of_data, following must hold:
 >> general size of matrices: 
 >![enter image description here](https://lh3.googleusercontent.com/qvX3mVPqOhcma_LnSqK4uZ-4olQnBYPfmI9X3TP7JBMmFyCNAHzjqJRqt7i7LzPGLi1u8XvFFXZd "general size of matrices")
>>together with the sizes of l/u bounds, x and b:![enter image description here](https://lh3.googleusercontent.com/dAUwdfJjT0pZwB1_Hk_XngEUHih6wT4i1BBg5oWTiX5pJvAOXeFQMfFApVp-OaGmCUhViqbvmMvc "picture of sizes")

> - lower and upper bounds can’t be infinity
> - and for every number l in the lower bound vector and for every number u in the upper bound vector must hold that l<=u (for l,u in the same position)
>- the logging level + format, self arguments such as self.A, self.D, self.t, self.r, self.s etc. are being initialized
>- if any initial feasible solution was given, it is being set to self.current_solution; if there isn’t any, the initial feasible solution will be computed by the find_init_feasible_solution function later
>- graver complexity is  being computed and set to self.graver_complexity by a function approximate_graver_complexity, or exact_graver_complexity, or only the value is being applied (according to the graver_complexity option)
>- ZE to self.ZE is being set through the function ``_construct_ZE`` (only if ``experimental == false``)
>- checking the validity of the data (whether the sizes of matrices, l/u bounds etc. make sense

Two functions for computing the graver complexity of given matrix (at most one of these is chosen according to the tenth (graver_complexity) argument in the __init__ function):

**exact_graver_complexity** returns an exact value of the Graver complexity of the given data from the definition
>**Implementation**:

>  takes the maximum entry of a matrix which  is product of multiplication of D by the graver basis of the matrix A ( the graver base has been computed using 4ti2/set as an integer in the init function)

**approximate_graver_complexity** return an approximate value of the Graver complexity of the given matrices
>**Implementation**:

>
>- computes graver complexity according to the following formula [2]
$$g(A)\leq p(r*|| D*GA||_{\infty})^{r} $$
where:
>> 
>>- g(A) is the graver  complexity of matrix A
>>- p is the number of elements in the graver basis of A
>>- GA is the graver basis of A

Function for computing Z(E) - **construct_ZE**. Z(E) is the sum of at most Graver complexity elements of the matrix A. This function is also called from the ``__init__`` function.
>**Implementation**
>- at first creates vector of zeros (it is definitely in ZE) 
>- then in a three inner for cycles happen following:
>>- 1st cycle:   graver complexity times new empty set is created
>>- 2nd cycle:   depending on the size of yet computed unique elements of ZE
>>- 3rd cycle:   two vectors are computed (graver complexity times) – it’s a sum/difference of one vector from yet computed ZE with a vector from graver basis of ``A``


Finding feasible solution - ``find_init_feasible_solution``. It computes the initial solution if it has not been given in the ``__init__`` function. It consists of two methods - **create_auxiliary_program**, which creates the instance of an auxiliary program. Then there is a method for solving the aux instance - **solve_auxiliary_program**.


>**Implementation  (**``create_auxiliary_program``**)**:

>- at first it constructs two matrices:
>>- $$D_{2} = [D Z_{sr} Z_{sr} I_{s} -I_{s}]$$
>>- $$A_{2} = [A I_{r} -I_{r} Z_{rs} Z_{rs}]$$
>>>- where $$M_{ab}$$ means a matrix of $$a$$ rows and $$b$$ columns, $$I$$ is an identity matrix, $$Z$$ is a matrix full of zeros
>- then then it computes new lower and upper bounds:
>>- lower bound vector consists of $$(t+2*r+2*s)*n$$ zeros
>>- upper bound vector consist of $$(t+2*r+2*s)*n$$ times the max value in the ``self.b`` vector
>- creates new objective function:
>>- it consist of a vector of $$t$$ times zero and $$2*r+2*s$$ times one which is n-times copied
>- makes the initial feasible solution of the auxiliary program
>>- the vector of the initial feasible solution consists of 2 types of vectors:
>>>- the first type has $$t+2r+2s$$ numbers and each number is a value from the lower vector (on the corresponding position) if it’s not minus infinity, elif it’s a value from the upper bound vector if it’s not an infinity, elif it is zero; then is this vector filled with positive numbers form the $$b[0]$$ vector (or zeros when not positive), then the negative numbers from the $$b[0]$$ vector (or zeros when not negative) and the same with the $$b[1]$$ vector (positive and negative part)
>>>- the second type of the vector of the init feasible solution is actually the same as the first type but the first part is different, there are only zeros
>>>-note: the vector ``b`` consists of two parts - one corresponds to the ``D`` matrix, the second part corresponds to the ``A`` matrix
>- finally it returns an instance of NFoldIP


>**Implementation (**``solve_auxiliary_program``**):**
>- its argument is an auxiliary instance of NFoldIP with its initial solution
>- uses the algorithm to minimize the auxiliary variables in order to have the initial solution for the main program
>- then it checks whether the auxiliary vars are zero -- if yes, we have an initial feasible solution (returns the init feasible solution), otherwise the main program has no feasible solution (returns None)


If the initial feasible solution exists we are searching for the augmenting steps by the function **find_graverbest_step**.

>**Implementation:**
>- at the beginning the generator of Gamma of gammas has to be computed (see the function ``construct_Gamma below``)
>- then for the gamma in Gamma:
>> while the dot product of vector ``w`` with the ``good_step`` (which was computed from the ``_find_good_step(gamma)`` function) is not greater or equal to zero:
>>>- try to take bigger gamma in order to prolong the good step
>- at the end the function returns the maximum step, which is ``gamma*good_step`` with the best dot product with corresponding ``w``

Now there is a short look into the construction of Gamma - **construct_Gamma**. It is an iterator of gammas for possible extension of the lengths of the feasible steps. It is used in the ``find_graverbest_step`` function as it has been mentioned before.

>**Implementation:**

>- it is an iterator of logarithmic values
>- it makes a vector of the difference from upper-lower bounds, chooses the biggest element of the difference, the max value which this iterator returns is the number ``n``(floor) for which holds following:
>>- $$n^{r} = max_value$$
>- the first value is 1, every next value is two times the previous value if it is lower than the ``max_value``, otherwise ``StopIteration`` is raised


There is the option not to use ``find_good_step(gamma)``, but its experimental version - **find_good_step_ecperimental(gamma, time limit)**.
>**Implementation:**
>- TODO

And finally the crucial function for solving a given n-fold program: ``solve(“solver”)``. It takes one argument which says which solver should be used. 
- if ``"native"``, it uses the ``find_graverbest_step`` to compute the solution 
- if ``“GLPK”`` it creates a MILP instance and solves it with GLPK 

If we choose the native solver the function **native_solve** is called. It solves the given problem with the implemented algorithm.
>**Implementation:**
>- it finds and applies the graver best steps by the function ``find_graverbest_step`` until there are no more graver best steps

The second option for solving calls the function **glpk_solve**. This function builds MILP in a standard form and uses GLPK for solving the problem.
>**Implementation:**
>- sets ``maximization = False``, add new nonnegative variable ``d``
>- builds the form of ``Ad=b``, $$d≥0$$
>- and minimizes ``wd``
>- it does not respect any initial solution if given at the beginning 


## SOURCES:

[1] R. Hemmecke, S. Onn and L. Romanchuk, “N-fold integer programming in cubic time,” Mathematical Programming, pp. 1-17, 2013. 

[2] S. Onn, Nonlinear discrete optimization, Zurich Lectures in Advanced Mathematics, European Mathematical Society, 2010. 

<!--stackedit_data:
eyJoaXN0b3J5IjpbNDIwODU3NjkxXX0=
-->
