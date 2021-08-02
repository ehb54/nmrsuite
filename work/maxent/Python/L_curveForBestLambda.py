#This script tells you which of the values of lambda is the best lambda for you to use.
#You can have that information by looking at the variable index68 for
#example. This scrip plots the L-curve got by Log(chi^2) x Log(entropy).  

import json
import numpy as np
from scipy.interpolate import UnivariateSpline
from csaps import CubicSmoothingSpline

A = np.loadtxt("result_of_maxent/A.txt")
lambda_ = x = np.loadtxt("result_of_maxent/lambda.txt") #'lambda' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
x = np.loadtxt("result_of_maxent/x.txt")
y = np.loadtxt("result_of_maxent/y.txt")

y_prime = np.array([])

def index_max(array): #Return the (first) index of the maximum value and the maximum value in a given array
    _max = np.max(array) #'max' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
    index = np.where(array == _max)[-1][0]
    return ([index, _max])

def lcurve_points(y, A, x, lambda_):

    #x = xsol(:,2:end)
    Ax = np.matmul(A, x)

    w = np.ones((np.shape(A)[1], 1)) * np.shape(A)[1] #vector filled with the number of maximum size(A,1) of A

    r = np.empty((np.shape(Ax)[1], np.shape(Ax)[0]))
    rn = np.empty(np.shape(Ax)[1])
    xval = np.empty(np.shape(Ax)[1])
    f = np.empty(np.shape(Ax)[1])

    for i in range (np.shape(Ax)[1]): #Returns Chi^2 for each lambda.
        r[i, :] = Ax[:,i] - y
        rn[i] = np.linalg.norm(r[i,:])/np.linalg.norm(y)

        xval[i] = np.sum(np.multiply(x[:,[i]], np.log(np.multiply(w, x[:,[i]])))) #Entropy computation for each lambda
        f[i] = np.square(rn[i]) #Chi^2 for each lambda (and for best solution i = 1)

    #pp = csaps(x,y) returns the cubic smoothing spline interpolation to the given data (x,y) in ppform. 
    #The value of spline f at data site x(j) approximates the data value y(:,j) for j = 1:length(x).
    # csaps
    # x = log(xval(1:end))
    # y = log(f(1:end))
    # p = .995 (smooth par)
    """ x = np.log(xval[1:])
    y = np.log(f[1:])
    smoothing_factor = .995 """
    x_par = np.flip(np.log(xval)) #x has to be strictly increasing
    y_par = np.flip(np.log(f))

    """ pp = UnivariateSpline(x=x_par, y=y_par, s=0.995)
    y_prime = pp(x=np.log(xval), nu=2) """

    pp = CubicSmoothingSpline(x_par, y_par, smooth=0.995).spline #
    p_first = pp.derivative(nu=1)
    p_second = pp.derivative(nu=2)
    y_prime = p_second(np.log(xval))

    y_prime = np.flip(y_prime)

    [index,v] = index_max(y_prime)

    full_x = np.linspace(np.min(np.log(xval)),np.max(np.log(xval)), 200)
    ydf2_full = p_second(full_x)

    return [index, xval, f]

x = x[:, 1:] #excluding the results with lambda = 0
lambda68 = lambda_[1:] #from results_of_maxent

[index68,x68,y68] = lcurve_points(y, A, x, lambda68)

with open('index.txt', 'w') as f: f.write(str(index68))

""" figure(1)
plot(log(x68),log(y68),'m.')
xlabel('log(Entropy)')
ylabel('log(\chi^2)')
#legend(strcat('pH 4.5 best \lambda= ', num2str(lambda45(index45))),strcat('pH 6.8 best \lambda= ', num2str(lambda68(index68))),strcat('pH 7.6 best \lambda= ', num2str(lambda76(index76))),'location', 'best')
set(gca, 'FontName', 'Arial', 'FontSize', 14)   
#ylim([-6 0]) """