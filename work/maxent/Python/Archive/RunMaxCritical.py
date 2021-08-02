import numpy as np
 
# This is the code run in the for loop
# I tried to provide the right amount of comments but please...
# ... let me know if you have any questions or clarifications

if (len(xout) != 0): # If 'xout' is not empty
    x0 = xout[:] # Set x0

length = np.shape(x0)[0] # Number of elements in x0

# The following two lines are bounds for the weights (elements of x)
lb = np.zeros((length, 1)) # Lb = Lower bound, this creates an array of dimensions Length x 1 with all 0s
ub = np.ones((length, 1)) # Ub = Upper bound, this creates an array of dimensions Length x 1 with all 1s


""" cl = np.ones((length, 1))
cu = np.ones((length, 1)) """

# The below numbers are the constraint lower and upper bounds
""" cl = np.full((length, 1), 0.9)
cu = np.full((length, 1), 1.1) """
cl = 0.9
cu = 1.1

Ac = np.ones((length, 1))

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
        return np.dot(Ac.flatten(), x)

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
        print ("Objective value at iteration #%d is - %g" % (iter_count, obj_value))

nlp = ipopt.problem(
    n=len(x0),
    m=len(cl),
    problem_obj=hs(),
    lb=lb,
    ub=ub,
    cl=cl,
    cu=cu
)

nlp.addOption(b'limited_memory_update_type', b'bfgs')
nlp.addOption(b'hessian_approximation', b'limited-memory')
nlp.addOption(b'limited_memory_max_history', 3)
#nlp.addOption(b'constr_viol_tol', 1e-6)
nlp.addOption(b'constr_viol_tol', 0.1)
nlp.addOption(b'tol', 1e-8)
nlp.addOption(b'max_iter', 50000)

x, info = nlp.solve(x0)