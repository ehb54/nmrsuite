#!/opt/miniconda3/bin/python
import numpy as np
import scipy.optimize
import math
import os

from scipy.optimize import Bounds, LinearConstraint

import cyipopt as ipopt

def matlab_range (start, stepSize, stop):
  # In Matlab, one can type something like "[-0.6:0.2:7.0]",
  # which would return '-0.6 -0.4 -0.2 ... 6.8 7.0',
  # but Python has no reliable way of doing this.
  # Therefore, I made this function to emulate it.
  stepDifference = stop - start
  numSteps = round(stepDifference/stepSize) + 1
  return np.linspace(start, stop, numSteps)

def run_max_entropy (A_filename, y_filename, lambda_lower, lambda_step, lambda_upper, run_directory):
  A = np.loadtxt(A_filename)
  y = np.loadtxt(y_filename)
  y = y[:, 0]

  lambda_ = np.power(2, matlab_range(lambda_lower, lambda_step, lambda_upper)) # 'lambda' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
  best = scipy.optimize.nnls(A, y)[0] # PRODUCES X VECTOR; Python function: Nonnegative least squares
  
  x0 = np.ones(np.shape(A)[1]) / np.shape(A)[1]

  w = np.ones(np.shape(A)[1]) * np.shape(A)[1] # Prior distribution vector

  """ y = np.round(y, 5)
  A = np.round(A, 5)
  w = np.round(w, 5) """

  xsol = mymaxent(A, y, lambda_, w, x0+1.0e-6)

  #print (xsol)
  
  xsol[:, 0] = best
  lambda_ = np.insert(lambda_, 0, 0)
  res = np.array([])
  lambdares = np.array([])

  for i in range (np.shape(xsol)[1]):
    res = np.append(res, np.linalg.norm(y - np.matmul(A, xsol[:, i]))/np.linalg.norm(y))
    with np.errstate(divide='ignore', invalid='ignore'): # Ignores division by zero or invalid value errors
      lambdares = np.append(lambdares, np.sum(np.multiply(xsol[:, i], np.log(np.multiply(w, xsol[:, i]))))+np.sum(np.size(xsol[:,i])))
  
  xnorm = np.array([])
  for i in range (np.shape(xsol)[1]):
    xnorm = np.append(xnorm, xsol[:,i]/np.sum(xsol[:,i]))

  if (os.path.exists(run_directory) == False):
    os.mkdir(run_directory)
  
  np.savetxt(os.path.join(run_directory, "A.txt"), A)#, np.round(A, 5),  fmt='%f')
  np.savetxt(os.path.join(run_directory, "lambda.txt"), lambda_)#, np.round(lambda_, 5),  fmt='%f')
  np.savetxt(os.path.join(run_directory, "x.txt"), xsol)#, np.round(xsol, 5),  fmt='%f')
  np.savetxt(os.path.join(run_directory, "y.txt"), y)#, np.round(y, 5),  fmt='%f')

  print ("Done.")

def mymaxent(A, y, lambda_, w, x0):
  Ac = np.ones(int(np.shape(A)[1]))
  bc = 1
  At = A.T
  """ print (Ac)
  print (bc)
  print (At) """

  xout = np.empty((np.shape(A)[1], 1))

  for i in range (len(lambda_)):
    print ("A")
    #print('iter %d of %d\n' % (i, len(lambda_)))

    if (np.shape(xout)[1] != 1):
      x0 = xout[:, -1]

    #options = optimset(options,'Hessian', {'lbfgs',3})
    # @(x)myfun(x,A,At,y,w,lambda(iter)) is used because myfun has more than one variable, but we want to look at x.
    # x0 = x0
    # A = []
    # b = []
    # Aeq = Ac
    # Beq = bc
    # lb = zeros(length(x0),1)
    # ub = ones(length(x0),1)
    # nonlcon = []
    # options = options
    
    # Aeq and Beq means Aeq*x = beq and A*x ≤ b
    # lb and ub means solution x always between A and B (inclusive), so it can be anywhere from 0 all the way to 1

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    # For lb and ub use bounds 
    # For Aeq and Beq use constraints
    # Aeq and Beq mean 

    # options = optimset('Display','off','Algorithm','interior-point','TolCon',1.0e-6,'TolFun',0,'TolX',1.0e-8,'GradObj','on','MaxIter',50000,...
    #  'DerivativeCheck','off','MaxFunEvals',50000,...
    #  'Diagnostics','off','FinDiffType','central','ScaleProblem','obj-and-constr','PlotFcns',[]);
    
    args = (A, At, y, w, lambda_[i])
    lambda_value = lambda_[i]
    #bounds = Bounds(0, 1)
    bounds = np.full((100,2),np.array([0, 1]))
    constraints = LinearConstraint(Ac, bc, bc)

    method = "L-BFGS-B"

    options = dict(
        maxcor=3,  # Number of previous gradients used to approximate the Hessian
        ftol=1e-6,
        gtol=1e-6,
        maxfun=500,  
        maxiter=500,
        iprint=99
    )

    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

    length = np.shape(x0)[0]

    lb = np.zeros(length)
    ub = np.ones(length)

    """ cl = np.ones((length, 1))
    cu = np.ones((length, 1)) """

    """ cl = np.ones((length, 1))
    cu = np.ones((length, 1)) """
    """ cl = np.full((length, 1), 0.9)
    cu = np.full((length, 1), 1.1) """

    cl = np.array([1])
    cu = np.array([1])

    Ac = np.ones(length)

    class hs(object):
        def __init__(self):
            pass

        def objective(self, x):
            #
            # The callback for calculating the objective
            #
            Ax = np.matmul(A, x)
            r = Ax-y #for chi^2
            logwx = np.log(np.multiply(w, x)) #for entropy, x is the probability distribution
            f = sum(np.square(r)) + (lambda_value ** 2) * np.dot(x, logwx) #eq 8 of jacs paper
            return f

        def gradient(self, x):
            #
            # The callback for calculating the gradient
            #
            Ax = np.matmul(A, x)
            r = Ax-y #for chi^2
            logwx = np.log(np.multiply(w, x)) #for entropy, x is the probability distribution
            g  = 2*np.dot(At, r) + lambda_value ** 2 * (1 + logwx) #what is that?
            #g = g.reshape((-1, 1))
            return g

        def constraints(self, x):
          #
          # The callback for calculating the constraints
          #
          return np.array([np.dot(Ac.flatten(), x)])

        def jacobian(self, x):
          return Ac#.flatten()

        def intermediate(
            self,
            alg_mod,
            iter_count,
            obj_value,
            inf_pr,
            inf_du,
            mu,
            d_norm,
            regularization_size,
            alpha_du,
            alpha_pr,
            ls_trials
            ):

            #
            # Example for the use of the intermediate callback.
            #
            print()
            #print ("Objective value at iteration #%d is - %g" % (iter_count, obj_value))

    nlp = ipopt.Problem(
        n=len(x0),
        m=1,
        problem_obj=hs(),
        lb=lb,
        ub=ub,
        cl=cl,
        cu=cu
    )

    nlp.add_option(b'limited_memory_update_type', b'bfgs')
    nlp.add_option(b'hessian_approximation', b'limited-memory')
    nlp.add_option(b'limited_memory_max_history', 3)
    nlp.add_option(b'constr_viol_tol', 1e-6)
    nlp.add_option(b'tol', 1e-8)
    nlp.add_option(b'max_iter', 50000)

    x, info = nlp.solve(x0)

    print (x)
    xout = np.column_stack((xout, x))
    #print (info)


    print ("Y")
    #[x, res] = fmincon(@(x)myfun(x,A,At,y,w,lambda(iter)),x0,[],[],Ac,bc,zeros(length(x0),1),ones(length(x0),1),[], options)
    #res = scipy.optimize.minimize(fun=myfun, x0=x0, args=args, method=method, bounds=bounds, constraints=constraints, options=options)

    """ args = (A, At, y, w, lambda_[i])
    bounds = Bounds(0, 1)
    constraints = LinearConstraint(Ac, bc, bc)

    method = "L-BFGS-B"

    options = dict(
        maxcor=3,  # Number of previous gradients used to approximate the Hessian
        ftol=1e-6,
        gtol=1e-6,
        maxfun=50000,  
        maxiter=50000,
    )
    
    res = scipy.optimize.fmin_l_bfgs_b(func=myfun, x0=x0, args=args, approx_grad=True, bounds=bounds,
    m=3, factr=1e7, pgtol=1e-5, maxfun=50000, maxiter=50000) """

    #print (res)
    #print (res)

    '''
    scipy.optimize.minimize
    x0 = x0
    
    '''
  print (xout)
  return (xout)

run_max_entropy("./Matrix_Ub2A_MNODES.txt", "./y_Ub2A.txt", -1, 0.2, 1, "result_of_raquel_test")
#run_max_entropy("./MatrixA_test.txt", "./y_test.txt", -0.6, 0.2, 7, "result_of_maxent")
