from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.rings.infinity import Infinity
from sage.matrix.special import block_matrix
from sage.matrix.special import diagonal_matrix
from sage.matrix.special import zero_matrix
from sage.matrix.special import identity_matrix
from sage.functions.other import floor, ceil
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.structure.sage_object import SageObject
from sage.interfaces.four_ti_2 import four_ti_2
from sage.numerical.mip import MIPSolverException

from itertools import combinations_with_replacement
from math import log

import itertools
import logging
import time
import sys

class construct_Gamma_iterator:
    '''iterator that yields numbers in Gamma'''

    def __init__(self,u,l,n):
        
        #make the difference upper-lower bounds
        difference = []
        for i in range(n):
            dif = (u[i]-l[i])
            difference.append(dif)
        
        #choose the biggest element of the difference
        max_val = float('-inf')
        for dif in difference:
            max_entry = dif.norm(Infinity)
            if max_entry>max_val:
                max_val = max_entry
                
        self.logarithm = floor(log(max_val,2)) #floor nebo ceil?
        
        #counters and last item
        self.num_of_iterations = 0
        self.last = 1
        

    def __iter__(self):
        return self

    def next(self):
        last_gamma = self.last
        self.num_of_iterations = self.num_of_iterations +1
        
        if self.num_of_iterations > self.logarithm:
            raise StopIteration #pokud bude problematicke pracovat s exception, klidne vratit tady -1 nebo jiny podobny marker
        self.last = last_gamma * 2
        return last_gamma


class NFoldIP(SageObject):
    
    def __init__(self, A, D, n, b, l, u, w, verbose = logging.ERROR, graver_complexity = "exact", current_solution = None, experimental=False, instancename="instancename", gamma="logarithmic",solver="GLPK"):  
        
        #SET LOGGING
        self.NFoldLogging = logging.getLogger(__name__)
        if not getattr(self.NFoldLogging, 'handler_set', None):
            formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d  %(message)s',datefmt='%H:%M:%S') #formates time for logging
            handler1 = logging.StreamHandler() #stdout
            handler1.setFormatter(formatter)
            if(verbose!="None"):
                handler1.setLevel(verbose)
                self.NFoldLogging.setLevel(verbose)
            self.NFoldLogging.addHandler(handler1)
            self.NFoldLogging.handler_set = True
            if(verbose=="None"):
                self.NFoldLogging.disabled = True
        self.instancename = instancename
        self.filelog = open(self.instancename+'.log', 'w')
        self.solver = solver
        
        #CHECK THE VALIDITY OF GIVEN DATA
        self._check_validity_of_data(A,D,n,b,l,u,w,current_solution)
        
        #MATRICES AND VECTORS OF THE GIVEN PROBLEM
        self.A = A #matrix A
        self.D = D #matrix D
        self.n = n
        self.b = b
        self.l = l #lower bound
        self.u = u #upper bound
        self.w = w #objective function
        self.w_vector = vector(sum((ww.list() for ww in w), []))
        self.t = D.ncols() #number of columns of A,D
        self.s = A.nrows() #number of rows of A
        self.r = D.nrows() #number of rows of D
       
        #INITAL FEASIBLE SOLUTION 
        self.current_solution = current_solution #can be None
        
        if experimental == True:
            self._find_good_step = self._find_good_step_experimental
            self.experimental = True
            self.GA = four_ti_2.graver(self.A) #Graver basis of A
        elif experimental == "ng1":
            self.experimental = "ng1"
            self._find_good_step = self._find_good_step_experimental_ng
        elif experimental == "nginfty":
            self.experimental = "nginfty"
            self._find_good_step = self._find_good_step_experimental_ng

        self.average_good_step_time = -1
        self.max_good_step_time = 0
        self.feasible_good_steps = 0
        
        self.gamma = gamma
        
        #COMPUTE/SET THE GRAVER COMPLEXITY
        if not experimental:
            self.GA = four_ti_2.graver(self.A) #Graver basis of A
        if graver_complexity == "exact":
            self.graver_complexity = self.exact_graver_complexity()
        elif graver_complexity =="approximate":
            self.graver_complexity = self.approximate_graver_complexity()
        elif type(graver_complexity) == int:
            self.graver_complexity = graver_complexity
        else:
            print ("not a valid argument for the graver_complexity keyword")
        


        if not experimental:
            #CONSTRUCT ZE
            self.ZE = self._construct_ZE() #construct ZE 
            self.experimental = False
            

        #LOGGING  PROPERTIES
        self.graver_best_counter = 0		
        self.all_find_good_step_counter = 0		
        self.find_good_step_counter = 0		
        self.set_of_time_find_good_step = []
        
        #VARS FOR SOLUTION
        self.native_solution = None
        self.glpk_solution = None
        
        
    def _check_validity_of_data(self,A,D,n,b,l,u,w,current_solution):
        "Check if the sizes and bounds are valid."
        
        self.NFoldLogging.info('Checking the validity of input data.')    
        
        #sizes of A,D
        if D.ncols() != A.ncols():
            raise ValueError("Matrices A,D don't have the same number of columns! The instance was not initialized.")
        
        #vectors of lists of l/u bounds and current solution (if any)
        if len(l) != n:
            print("zde-----------")
            print l
            print n
            raise ValueError("The length of the lower bound vector must have the same number of vectors as the input n! The instance was not initialized.")
        if len(u) != n:
            raise ValueError("The length of the upper bound vector must have the same number of vectors as the input n! The instance was not initialized.")
        if current_solution != None:
            if len(current_solution) != n:
                raise ValueError("The length of the current solution vector must have the same number of vectors as the input n! The instance was not initialized.")
                
        #sizes of vectors in l/u bounds and current solution (if any)                 
        for i in range(n):    
            if len(l[i]) != D.ncols():
                raise ValueError("The length of the first vector in lower bounds must have the same lenght as the number of columns in D! The instance was not initialized.")
            if len(u[i]) != D.ncols():
                raise ValueError("The length of the first vector in upper bounds must have the same lenght as the number of columns in D! The instance was not initialized.")
            if current_solution != None:
                if len(current_solution[i]) != D.ncols():
                    raise ValueError("The length of the first vector in the current solution must have the same lenght as the number of columns in D! The instance was not initialized.")
            
        #upper/lower can't be infinity!
        for v in l:
            for vi in v:
                if v==float('inf'):
                    raise ValueError("Lower bound can't be infinity! The instance was not initialized.")
                if v==float('-inf'):
                    raise ValueError("Lower bound can't be -infinity! The instance was not initialized.")
        for v in u:
            for vi in v:
                if u==float('inf'):
                    raise ValueError("Upper bound can't be infinity! The instance was not initialized.")
                if u==float('-inf'):
                    raise ValueError("Upper bound can't be -infinity! The instance was not initialized.")
        
        #if upper<=lower in the first vector according to current solution (if any)
        for i in range(n):
            for j in range(D.ncols()):
                if l[i][j]> u[i][j]:
                     raise ValueError("Upper bounds must be greater or equal to the lower bounds! The instance was not initialized.")
                if current_solution != None:
                    if l[i][j] > current_solution[i][j]:
                        raise ValueError("Given current solution does not satisfy lower bounds! The instance was not initialized.")
                    if current_solution[i][j]>u[i][j]:
                        raise ValueError("Given current solution does not satisfy upper bounds! The instance was not initialized.")
                         
        #size of b
        if len(b) != n+1:
            raise ValueError("The length of b vector must have the same number of vectors as the input n + 1! The instance was not initialized.")
        if len(b[0]) != D.nrows():
            raise ValueError("The first vector of b must have the same size as the number of rows in D! The instance was not initialized.")
        for i in range(n):
            if len(b[i+1])!=A.nrows():
                raise ValueError("The second and next vectors of b must have the same size as the number of rows in A! The instance was not initialized.")          
            
        #size of w    
        if len(w) != n:
            raise ValueError("The length of the objective function w must have the same number of vectors as the input n! The instance was not initialized.")
        for i in range(n):
            if len(w[i]) != D.ncols():
                raise ValueError("The vectors of w must have the same size as the number of columns in D/A! The instance was not initialized.")
        
        self.NFoldLogging.info('Input data checked.')      
        
    def exact_graver_complexity(self):
        """
        Return an exact value of the Graver complexity of the given data from the definition.
        """
        value = four_ti_2.graver((self.D)* self.GA.transpose()).norm(Infinity)
        self.NFoldLogging.info('Exact Graver complexity: ',value)
        return value
      
        
    def approximate_graver_complexity(self):
        """
        Return an approximate value of the Graver complexity of the given matrices.
        """
        p = self.GA.nrows() * 2
        G = self.GA.transpose()
        mult = self.D*G #ZDE BYL BUG! Bylo tu self.A*G
        biggestValue = mult.height() 
        value = ceil(p*( (self.r* biggestValue)^self.r)) # A TADY BYLO (sqrt(self.r)) uvnitr
        self.NFoldLogging.info('Approximate Graver complexity: ',value)     
        return value
        
        
    def _construct_ZE_slower(self):
        """
        Construct ZE - the sum of at most Graver complexity elements of the matrix D using itertools.
        """
        #start = time.time() 
        #nejprve se v GA museji zdvojit vysledky (prenasobit minus jednickou vsechno)
        GA2 = []
        for i in self.GA:
            GA2.append(i)
            GA2.append(-i)
        
        #pridat nulovy vektor
        GA2.append(vector((0,)*self.t))
        
        R = set()
        for x in combinations_with_replacement(GA2,self.graver_complexity):
            vv=sum(vector(v) for v in x)
            vv.set_immutable()
            R.add(vv)
            
        end = time.time()
        #print end-start    
        
        self.NFoldLogging.info('COMPUTING ZE - ZE has {} unique vectors'.format(len(R)))
                
        return R
    
    def _construct_ZE(self):
        """
        Construct ZE - the sum of at most Graver complexity elements of the matrix A.
        """
        self.NFoldLogging.debug('COMPUTING ZE:')
        #self.NFoldLogging.fairytale('I HAVE FOUND OUT FOLLOWING VECTORS (not every is unique): ')
        vector0 = vector((0,)*self.t)
        vector0.set_immutable()
        count = 0
        R = {vector0}
        for i in range(self.graver_complexity):
            R_new = set()    
            for x in R:
                count = count+len(R)*self.GA.nrows() #counter of iterations in this inner forcycle
                for y in self.GA:
                    z1 = x+y
                    z2 = x-y
                    z1.set_immutable()
                    z2.set_immutable()
                    #self.NFoldLogging.fairytale('vector z1: {}, vector z2:{}'.format(z1,z2))
                    if z1 not in R:
                        self.NFoldLogging.debug('new unique vector: {}'.format(z1))
                    R_new.add(z1)
                    if z2 not in R:
                        self.NFoldLogging.debug('new unique vector: {}'.format(z2))
                    R_new.add(z2)
            R = R.union(R_new)
        R.add(vector0)              
        self.NFoldLogging.info('COMPUTING ZE - ZE has {} unique vectors, but {} not necessarily unique vectors were computed'.format(len(R),self.graver_complexity*count))
        return R
                
        
    def _find_good_step(self,gamma):
        """
        Dynamic program runs exactly here. Finding the best path in a weighted acyclic diagraph of layers from ZE elements.
        
        """

        from copy import copy
        current = [(vector((0,)*self.t), [], 0)]
            
        for i in range(self.n):
            last = current
            current = []
                
            for h in self.ZE:
                best_val = float('inf')  # int??
                best_path = []
                if (i == self.n-1) and (self.D*h != vector((0,)*self.r)):
                    continue
                for g, g_path, g_val in last:
                    hg = vector(h-g)
                    hg.set_immutable()
                    if not (hg in self.ZE):
                        continue
                        
                    help = self.current_solution[i] + (hg * gamma)
                    lbub = all((ll <= hh for ll, hh in zip(self.l[i],help) )) and all((uu >= hh for hh,uu in zip(help, self.u[i]) )) #when lower and 
                    if lbub:
                        new_value = (hg * self.w[i]) + g_val
                        if new_value < best_val:
                            best_val = new_value
                            best_pred = g
                            best_path = copy(g_path)
    
                v = vector(h)
                v.set_immutable()
                #pokud h neni dosazitelny, pak best_path neni inicializovane
                if best_val != float('inf'):
                    best_path.append(h)
                    current.append((v, best_path, best_val)) #(tuple([v, best_path, best_val]))
            
        last = current
        
        min = float('inf')
        path = []
        for (g, g_path, g_val) in last:
            if g_val<min:
                min = g_val
                path = g_path
        
        expanded_path = [path[0]] + [path[i] - path[i-1] for i in range(1,self.n)]
         
        return expanded_path
            
        
    def find_graverbest_step(self):
        """
        Generate generator for Gamma, which is a set of gammas.
        Then for each gamma in Gamma count good steps (function find_good_step(gamma)) and then choose the one step with the biggest dot product with the subjective function w.
        """
        self.NFoldLogging.debug('FINDING THE GRAVER BEST STEP: ')
        # MAX run time hardcoded to be one hour.
        previous_time = 3600
        current_min = 0
        min_step = vector((0,)*(self.n*self.t))
        
        if self.gamma == "unit":
            Gamma = [1]
        elif self.gamma == "best" and (self.experimental in ["ng1", "nginfty"]):
            Gamma = self.construct_best_Gamma()
        else:
            # logarithmic Gamma = default
            Gamma = construct_Gamma_iterator(self.u, self.l, self.n)        
        
        for gamma in Gamma:
            self.NFoldLogging.warning('Finding good step for GAMMA: {}'.format(gamma))
            good_step = self._find_good_step(gamma)
            start_time = time.time()
            if (self.average_good_step_time == -1) or (gamma==1):
                good_step = self._find_good_step(gamma)
            else:
                good_step = self._find_good_step(gamma, timelimit=min((3*self.max_good_step_time, 5)))
                #good_step = self._find_good_step(gamma)
            end_time = time.time()
            previous_time = end_time - start_time
            self.NFoldLogging.warning('Finding good step took {}'.format(previous_time))
            self.NFoldLogging.debug('good_step from find graverbest: {}'.format(good_step))
            
            if self.w_vector.dot_product(good_step) >= 0:
                break
            else:
                self.feasible_good_steps += 1
                fgs = self.feasible_good_steps
                self.max_good_step_time = max((self.max_good_step_time, previous_time))
                if self.average_good_step_time == -1:
                    self.average_good_step_time = previous_time
                else:
                    self.average_good_step_time = max(( self.average_good_step_time*( (fgs-1) / float(fgs)) + previous_time*(1/float(fgs)), previous_time))
                self.NFoldLogging.warning('Updated average time {}'.format(self.average_good_step_time))
                self.NFoldLogging.warning('Updated max time {}'.format(self.max_good_step_time))
            
            #look if it's possible to take bigger gamma for prolonging the good_step
            new_gamma = float('inf')
            
           
            #print("current_solution = {}".format(self.current_solution))
            #print("Gamma = {}".format(gamma))
            #print("good_step = {}".format(good_step))            
            #print("lower bound = {}".format(self.l))            
            #print("upper bound = {}".format(self.u)) 


            for i in range(len(good_step)):
                if good_step[i]!=0:
                    
                    ii = i/self.t
                    jj = i%self.t
                    
                    li = -(self.current_solution[ii][jj] - self.l[ii][jj])/good_step[i]
                    ui = (self.u[ii][jj]-self.current_solution[ii][jj])/good_step[i]
                    #take the bigger one
                    if ui>li:
                        bi = floor(ui)
                    else:
                        bi = floor(li)
                    #if new_gamma is too long, make it shorter
                    if new_gamma > bi:
                        new_gamma = bi
            
            if new_gamma != float('inf'):
                if gamma<new_gamma:
                    gamma = new_gamma         
                    self.NFoldLogging.warning('GAMMA CHANGED -- Finding good step for gamma: {}'.format(gamma)) 
                      
            
            self.NFoldLogging.debug('good_step with prolonging from find graverbest: {}'.format(good_step))
            self.find_good_step_counter += 1
            dot_product = self.w_vector.dot_product(gamma*good_step)
            self.NFoldLogging.info('for step length {} the best step has value {}'.format(gamma,dot_product))
            self.filelog.write(str(self.current_obj + dot_product) + " ")              
            
            if dot_product < current_min:
                self.NFoldLogging.warning('For gamma {} the improvement is {}'.format(gamma,dot_product))
                current_min =  dot_product
                min_step = gamma*good_step
                self.NFoldLogging.debug('New better dot product: current_min = {}, min_step = {}'.format(current_min,min_step))
            
        
        self.NFoldLogging.info('GRAVER BEST STEP has value {} and is = {}'.format(current_min,min_step))
        self.NFoldLogging.warning('in this call of find_graverbest_step the function _find good step was called {} times'.format(self.find_good_step_counter))                                
        self.all_find_good_step_counter = self.all_find_good_step_counter + self.find_good_step_counter		
        # We set it to one initially because on the last iteration, we do not add 1
        self.find_good_step_counter = 1
                   
        return min_step
        
    def construct_best_Gamma(self):
        Gamma = set()
        gc = self.graver_complexity
        for i in range(self.n):
            for j in range(self.t):
                uu = self.u[i][j]
                ll = self.l[i][j]
                xx = self.current_solution[i][j]
                #print "ll:",ll,"uu:",uu,"xx:",xx
                for k in range(1,gc+1):
                    gamma_1 = floor((uu-xx)/k)
                    gamma_2 = floor((xx-ll)/k)
                    Gamma.add(gamma_1)
                    Gamma.add(gamma_2)
        Gamma = list(Gamma)
        Gamma.sort()
        Gamma.remove(0)
        return Gamma

                
        
    def _find_good_step_experimental(self,gamma,timelimit=None):
        GA = [v for v in self.GA] 
        GA += [-v for v in GA] 
        GA += [vector((0,)*self.t)]
        DGA = [self.D*v for v in GA]
        p = MixedIntegerLinearProgram(maximization=False, solver = "GLPK")
        # indexed as x[layer, sublayer, index of element in GA]
        # NEW indexed as x[layer, index of element in GA]
        # number of layers is n
        # OLD number of sublayers is graver_complexity
        x = p.new_variable(integer=True, nonnegative=True)
        
        
        # at most graver_complexity elements selected across all layers
        p.add_constraint(sum((x[i,k] for i in range(self.n) for k in range(len(GA)))) <= self.graver_complexity)
        #p.add_constraint(-self.graver_complexity <=  sum((x[i,k] for i in range(self.n) for k in range(len(GA)))))
        
        # the sum of selected elements is in Ker(D)
        p.add_constraint(sum((x[i,k]*DGA[k] for i in range(self.n) for k in range(len(GA)))) == 0)
        
        # lower and upper bounds
        for i in range(self.n):
            p.add_constraint(sum((gamma*x[i,k]*GA[k] for k in range(len(GA)))) + self.current_solution[i] <= self.u[i])
            p.add_constraint(sum((gamma*x[i,k]*GA[k] for k in range(len(GA)))) + self.current_solution[i] >= self.l[i])
        
        # minimize the objective
        f = 0
        for i in range(self.n):
            for k in range(len(GA)):
                f += x[i,k]*self.w[i].dot_product(GA[k])
        p.set_objective(f)
        
        if timelimit:
            p.solver_parameter("timelimit", timelimit)
        
        try:
            p.solve()
        except MIPSolverException:
            return vector((0,)*(self.n*self.t))
        if p.get_objective_value() > -1:
            return vector((0,)*(self.n*self.t))
        
        # Construct the augmenting step
        h = []
        for i in range(self.n):
            h.append(sum((int(p.get_values(x[i,k]))*GA[k] for k in range(len(GA)))))
        hh = []
        for hhh in h:
            hh += hhh
        hh = vector(hh)
    
        return hh
        
    def _find_good_step_experimental_ng(self,gamma,timelimit=None):
        # norm can be "infty" or "1"; we constraint the solution to be bounded by self.graver_complexity of this norm
        # l_\infty
        #p = MixedIntegerLinearProgram(maximization=False, solver = solver)
        p = MixedIntegerLinearProgram(maximization=False, solver = self.solver)
        x = p.new_variable(integer=True, nonnegative=False)
        for i in range(self.n):
            for j in range(self.t):
                ub = floor( min(self.graver_complexity*gamma, self.u[i][j] - self.current_solution[i][j]) / gamma)
                lb = ceil( max(-self.graver_complexity*gamma, self.l[i][j] - self.current_solution[i][j]) / gamma)
                # print "setting lb/ub i,j", i,j, "to", lb, ub
                p.set_max(x[i,j],ub)
                p.set_min(x[i,j],lb)
        
        An = self.create_An_matrix()
        
        for row in An:
            p.add_constraint(sum((row[i*self.t + j]*x[i,j] for j in range(self.t) for i in range(self.n))) == 0)
        
        if self.experimental=="ng1":
            pos = p.new_variable(integer=True, nonnegative=True)
            neg = p.new_variable(integer=True, nonnegative=True)
            for i in range(self.n):
                for j in range(self.t):
                    p.add_constraint(x[i,j] == pos[i,j] - neg[i,j])
            p.add_constraint( sum(pos[i,j] + neg[i,j] for i in range(self.n) for j in range(self.t)) <= self.graver_complexity )
        
        
        # minimize the objective
        f = 0
        for i in range(self.n):
            for j in range(self.t):
                f += self.w[i][j]*x[i,j]
        p.set_objective(f)
        
        #return p
        
        if timelimit and self.solver != "Coin":
            p.solver_parameter("timelimit", timelimit)
        
        try:
            p.solve()
        except MIPSolverException:
            return vector((0,)*(self.n*self.t))
        if p.get_objective_value() > -1:
            return vector((0,)*(self.n*self.t))
        
        #for i, v in p.get_values(x).iteritems():
        #    print i, v
        # Construct the augmenting step
        h = []
        for i in range(self.n):
            for j in range(self.t):
                h.append(int(p.get_values(x[i,j])))
           
        hh = vector(h)
        return hh
        
        
    def _bricks_to_vector(self,bricks):
        v = []
        for i in range(self.n):           
            v+=bricks[i]
        v = vector(v)
        return v
        
        
    def _vector_to_bricks(self,vect):
        bricks = []
        for i in range(self.n):
            b = vector((0,)*self.t)
            for j in range(self.t):
                b[j]=vect[i*self.t + j]
            bricks.append(b)
        return bricks       
        
    
    def native_solve(self):
        """
        Solves the given problem with the implemented algorithm.
        """
        
        ### TODO: performance idea -- if next solve takes > 1.1*previous_solve, kill the solver and return the previous best step
        # Requires the following parameter. Not sure what p returns if it runs over the limit.
        # p = MixedIntegerLinearProgram(solver = "GLPK")
        # p.solver_parameter("timelimit", 60)
        
        
        self.NFoldLogging.info('NATIVE SOLVER: ')
        
        #initial feasible solution
        if self.current_solution == None:
            self.current_solution = self.find_init_feasible_solution()
        current_obj = self.w_vector.dot_product(self._bricks_to_vector(self.current_solution))
        self.current_obj = current_obj
        self.filelog.write(str(current_obj)+"\n")
        
        
        if self.current_solution != None:
            while True:
                start = time.time()
                step = self.find_graverbest_step()
                end = time.time()
                self.NFoldLogging.warning('time of find_graver_best_step: {}'.format(end-start))
                if (step == 0) or (self.w_vector.dot_product(step) > -1):
                    self.filelog.write(str(self.w_vector.dot_product(self._bricks_to_vector(self.current_solution))))
                    break
                self.NFoldLogging.info('Previous solution: {}'.format(self.current_solution))
                previous_obj = self.w_vector.dot_product(self._bricks_to_vector(self.current_solution))
                self.current_solution = self._vector_to_bricks(self._bricks_to_vector(self.current_solution) + step)
                current_obj = self.w_vector.dot_product(self._bricks_to_vector(self.current_solution))
                self.current_obj = current_obj
                #filelog.write(str(previous_obj)+"\n")
                self.filelog.write("\n")
                self.NFoldLogging.info('Previous objective: {}, New objective:{}'.format(previous_obj, current_obj))
                self.NFoldLogging.info('New solution: {}'.format(self.current_solution))
            #filelog.write(str(current_obj)+"\n")
            self.filelog.write("\n")
            self.filelog.close()
            self.NFoldLogging.info('find_good_step called {} times'.format(self.all_find_good_step_counter))
            self.NFoldLogging.info("The result is: {}".format(self.current_solution))     
            
            m = 0
            for i in range(self.n):
                for j in range(self.t):
                    m+=self.current_solution[i][j]*self.w[i][j]
            self.NFoldLogging.warning("Objective value from the Native solver: {}".format(m))
            self.native_solution = m
            
            
        else:
            self.NFoldLogging.warning('Initial feasible point not found.')
        
            
    def create_An_matrix(self):   
        #initalization of the new matrix A
        A_matrix = [[0 for x in range(self.n*self.t)] for y in range(self.r+self.n*self.s)] 
        
        #n-times copy the matrix D = the 1st self.r lines
        for row in range(self.r):
            for col in range(self.n*self.t):
                A_matrix[row][col] = self.D[row][col%self.t]
        
        #n-times copy the matrix A, but only to the main diagonal
        start_index_row = -self.t
        for row in range(self.n*self.s):
                  
            #change the start_index_row 
            if (row)%self.s == 0:
                start_index_row+= self.t
            
            #through columns    
            for col in range(self.n*self.t):
                if start_index_row <= col < (start_index_row+self.t):
                    #copy
                    A_matrix[row+self.r][col] = self.A[row%self.s][col%self.t]  
        
        #return the new A matrix
        
        #for v in A_matrix:
        #    print v
            
        return A_matrix
   
    def milp_instance(self):
        """
        Prepares and returns instance for GLPK.
        """
        
        #new A^(n) matrix 
        new_A_matrix = self.create_An_matrix()
        #print
        #for v in new_A_matrix:
        #    print v
        #print
        #print('len of A_matrix= {}'.format(len(new_A_matrix)))
        
        #initialization of MILP
        nameOfMilp = "milp"
        milp = MixedIntegerLinearProgram(maximization=False, solver = "GLPK")
        d = milp.new_variable(integer=True, nonnegative=True)
        
        #new format of self.b (only one long vector)
        bb = vector(itertools.chain.from_iterable(self.b))
            
        #adding constraints
        for row in range(len(new_A_matrix)):
            milp.add_constraint(sum(new_A_matrix[row][i]*d[i] for i in range(self.n*self.t)), max = bb[row], min = bb[row])
            
        #adding lower & upper bounds
        ll = vector(itertools.chain.from_iterable(self.l))
        uu = vector(itertools.chain.from_iterable(self.u))
        for i in range (self.r+self.n*self.s):
            milp.add_constraint(d[i], max = uu[i], min = ll[i]) 
        
        #set the objective function, which is w*d
        ww = vector(itertools.chain.from_iterable(self.w))  
        milp.set_objective(sum(ww[i]*d[i] for i in range(self.n*self.t)))
            
        #how does the MILP looks like
        #milp.show()
        
        return milp

    def glpk_solve_instance(self,milp):
        """
        Solves given MILP instance.
        """    
        self.NFoldLogging.info('I am going to solve the MILP instance by glpk.')
        #print milp.solve()
        #milp.solve()
        self.glpk_solution = milp.solve()
        self.NFoldLogging.warning('Objective Value from GLPK: {}'.format(self.glpk_solution))
                   
    def glpk_solve(self):
        """
        Creates and solves glpk instance
        """
        try:
            self.glpk_solve_instance(self.milp_instance())
        except MIPSolverException:
            print "GLPK has not found any solution"
    
    
    def solve(self,solver):
        """
        Search Graver best steps until there are no more Graver best steps.
        """
        if solver == "Native":
            self.native_solve()
                
        elif solver == "GLPK":
            self.glpk_solve()
            
        else:
            print ("Only Native and GLPK solver aviable.")
    
    def _give_D2(self):
        """
        Build D2 matrix for finding initial feasible solution.
        """
        Ir1 = identity_matrix(self.r)
        Ir2 = diagonal_matrix(self.r,[-1 for i in range(self.r)])
        
        Zrs = zero_matrix(self.r,self.s)
        
        D2 = matrix( block_matrix([self.D, Ir1, Ir2, Zrs, Zrs], nrows =1,ncols =5))
        return D2
        
        
    def _give_A2(self):
        """
        Build A2 matrix for finding initial feasible solution.
        """
                
        Is1 = identity_matrix(self.s)
        Is2 = diagonal_matrix(self.s,[-1 for i in range(self.s)])
        Zsr = zero_matrix(self.s,self.r)
                
        A2 = matrix( block_matrix([self.A,Zsr,Zsr,Is1,Is2],nrows =1,ncols =5))
        return A2       
        
    @staticmethod
    def _firs_part_vect_init_feas_sol(p):
        if p[0] != float('-inf'):
            return p[0]
        elif p[1]!= float('inf'):
            return p[1]
        else:
            return 0
            
    def pos(self,x):
        """
        When positive, return it. If not, return 0.
        """
        if x>=0:
            return x
        else:
            return 0

    def neg(self,x):
        """
        When negative, return it. If not, return 0.
        """
        if x>=0:
            return 0
        else:
            return x
    
    def _neg_or_pos_to_z(self,v):
        """
        Return feasible solution.
        """
        v1 = vector([self.pos(v[i]) for i in range(len(v))])
        v2 = vector([-self.neg(v[i]) for i in range(len(v))])
    
        return [v1,v2]
        
    def find_init_feasible_solution(self):
        """
        Creates and solves auxiliary_program.
        """
        self.solve_auxiliary_program(self.create_auxiliary_program())
            
    def create_auxiliary_program(self):
        """
        Creates aux program, returns istance for finding the initial feasible solution if no init solution has been given with the instance.
        """
        self.NFoldLogging.info('In create_aux_program')
        
        #new matrices
        D2 = self._give_D2()
        A2 = self._give_A2()
        self.NFoldLogging.info('Has new D = {}, and new A = {}'.format(D2,A2))
        
        #bounds
        l2 = [vector((0,)*(self.t + (2*self.s) + (2*self.r)))]*self.n
        max_in_b = 0 #and b should not contain infinity
        for vect in self.b:
            for v in vect:
                if max_in_b<abs(v):
                    max_in_b = abs(v)
        u2 = [vector((max_in_b,)*(self.t + (2*self.s) + (2*self.r)))]*self.n    
        self.NFoldLogging.info('Has new l,u bounds. l = {}, u= {}'.format(l2,u2)) 
        
        #objective function
        w2 = [vector((0,)*self.t + (1,)*((2*self.s) + (2*self.r)))]*self.n 
        self.NFoldLogging.info('Has new objective function: {}'.format(w2))
        
        #first vector in the feasible solution       
        first_part_z = [self._firs_part_vect_init_feas_sol(p) for p in zip (self.l[0], self.u[0])]
        D2_pos = [self.pos(self.b[0][i]) for i in range(len(self.b[0]))]
        D2_neg = [-self.neg(self.b[0][i]) for i in range(len(self.b[0]))]
        A2_pos = [self.pos(self.b[1][i]) for i in range(len(self.b[1]))]
        A2_neg = [-self.neg(self.b[1][i]) for i in range(len(self.b[1]))]        
        z = [vector(list(first_part_z)+ list(D2_pos) + list(D2_neg) + list(A2_pos) +list(A2_neg))]     

        #other vectors in the feasible solution
        for i in range(2,len(self.b)):
            zeros = [0 for y in range(2*self.r+self.t)]
            A2_pos = [self.pos(self.b[i][j]) for j in range(len(self.b[i]))]
            A2_neg = [-self.neg(self.b[i][j]) for j in range(len(self.b[i]))]
            z.append(vector((list(zeros)+list(A2_pos)+list(A2_neg))))
        
        self.NFoldLogging.info('Has the feasible point for the aux. program: {}'.format(z))
            
        #solving the auxiliary program
        auxiliary_program = NFoldIP(A2, D2, self.n, self.b, l2, u2, w2, verbose = logging.DEBUG, graver_complexity= self.graver_complexity, current_solution = z, experimental = self.experimental) 
        
        return auxiliary_program
        
        
    def solve_auxiliary_program(self,auxiliary_instance):      
        """
        Solves the given auxiliary instance.  
        """    
        self.NFoldLogging.info('Starting to solve the aux program')  
        auxiliary_instance.solve("Native")
        
        self.NFoldLogging.info('Aux program was solved. Its solution: {}'.format(auxiliary_instance.current_solution))
        
        #check wheather the auxiliary vars are zeros               
        bool_break = 0 
        for j in range(len(auxiliary_instance.current_solution)):
            if bool_break == 1:
                break
            for i in range(self.t,self.t-1+2*self.r+self.s*2):
                if bool_break == 1:
                    break
                if auxiliary_instance.current_solution[j][i] !=0:
                    self.NFoldLogging.warning('No feasible solution found!')
                    bool_break =1
                    break
        
        if bool_break == 0:        
            init_feas_solution = auxiliary_instance.current_solution
            print ("init_feas_solution: {}".format(init_feas_solution))
            return init_feas_solution    
        
        
        else:
            return None
