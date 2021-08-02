import numpy as np
import scipy.optimize
import math

from scipy.optimize import Bounds, LinearConstraint

def matlabRange (start, stepSize, stop):
  # In Matlab, one can type something like "[-.6:.2:7.0]",
  # but Python has no reliable way of doing this.
  # Therefore, I made this function.
  stepDifference = stop - start
  numSteps = round(stepDifference/stepSize) + 1
  return np.linspace(start, stop, numSteps)

def main ():
  A = np.loadtxt("MatrixA_test.txt")
  y = np.loadtxt("y_test.txt")
  y = y[:, 0]

  lambda_input = np.power(2, matlabRange(-0.6, 0.2, 7))

  best = scipy.optimize.nnls(A, y)[0]
  
  x0 = np.ones(np.shape(A)[1]) / np.shape(A)[1]

  w = np.ones(np.shape(A)[1]) * np.shape(A)[1]

  xsol = mymaxent(A, y, lambda_input, w, x0+1.0e-6)

  #w = ones(size(A,2),1)*size(A,2); #prior distribution
  
  """ xsol = [best xsol];
  lambda = [0 lambda];
    
  res = [];
  for iter=1:size(xsol,2)
    res(iter) = norm(y-A*xsol(:,iter))/norm(y);
    lambdares(iter) = sum(xsol(:,iter).*log(w.*xsol(:,iter)))+sum(length(xsol(:,iter)));
  end;
    
  xnorm = [];
  for iter=1:size(xsol,2)
    xnorm(:,iter) = xsol(:,iter)/sum(xsol(:,iter));
  end;

  result_of_maxent.x = xsol;
  result_of_maxent.lambda = lambda;
  result_of_maxent.A = A;
  result_of_maxent.y = y;

  save result_of_maxent.mat result_of_maxent
 """
def myfun(x, A, At, y, w, lambda_input): #
    Ax = np.matmul(A, x)
    r = Ax-y #for chi^2
    logwx = np.log(np.multiply(w, x)) #for entropy, x is the probability distribution
    f = sum(np.square(r)) + (lambda_input ** 2) * np.dot(x, logwx) #eq 8 of jacs paper

    g  = 2*np.dot(At, r) + np.multiply(np.square(lambda_input), (1 + logwx)) #what is that?
    g = g.reshape((-1, 1))
    #return np.array([f, g])#[f, g]
    return [f, g]

def mymaxent(A, y, lambda_input, w, x0):
  Ac = np.ones(int(np.shape(A)[1]))
  bc = 1
  At = A.T
  """ print (Ac)
  print (bc)
  print (At) """

  xout = np.array([])

  for i in range (len(lambda_input)):
    print ("A")
    #print('iter %d of %d\n' % (i, len(lambda_input)))

    if (len(xout) != 0):
      x0 = xout[:]

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

    
    args = (A, At, y, w, lambda_input[i])
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


    print ("Y")
    #[x, res] = fmincon(@(x)myfun(x,A,At,y,w,lambda(iter)),x0,[],[],Ac,bc,zeros(length(x0),1),ones(length(x0),1),[], options)
    res = scipy.optimize.minimize(fun=myfun, x0=x0, args=args, method=method, bounds=bounds, constraints=constraints, options=options)
   

    """ args = (A, At, y, w, lambda_input[i])
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

    print (res)
    #print (res)

    '''
    scipy.optimize.minimize
    x0 = x0
    
    '''

  return (x0)

main()
