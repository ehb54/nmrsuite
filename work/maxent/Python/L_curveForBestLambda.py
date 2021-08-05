#This script tells you which of the values of lambda is the best lambda for you to use.
#You can have that information by looking at the variable index68 for
#example. This scrip plots the L-curve got by Log(chi^2) x Log(entropy).  

import json
import numpy as np
from scipy.interpolate import UnivariateSpline
from csaps import CubicSmoothingSpline

import plotly
from plotly.subplots import make_subplots
import plotly.graph_objects as go


A = np.loadtxt("result_of_maxent/A.txt")
lambda_ = np.loadtxt("result_of_maxent/lambda.txt") #'lambda' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
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
    #y_par = np.flip(np.log(f))
    y_par = np.log(f)

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

    return [index, xval, f, pp]

x = x[:, 1:] #excluding the results with lambda = 0
lambda68 = lambda_[1:] #from results_of_maxent

[index68,x68,y68, pp] = lcurve_points(y, A, x, lambda68)

with open('index.txt', 'w') as f: f.write(str(index68))

fig = make_subplots(rows=1, cols=2, subplot_titles=("Plot 1", "Plot 2"))

# If figure is not being saved correctly, switch "append_trace" to "add_trace"
fig.append_trace(go.Scatter(x=np.log(x68), y=np.log(y68), mode="markers", marker=dict(color="Black"), name=r"$All other \lambda$"), row=1, col=1) #, line=dict(color="Blue", dash="dash")
fig.append_trace(go.Scatter(x=np.array(np.log(x68[index68])), y=np.array(np.log(y68[index68])), mode="markers", marker=dict(color="Red"), name=r"$\lambda$ chosen"), row=1, col=1) #, line=dict(color="Blue", dash="dash")

'''
xi = np.linspace(np.log(x68)[0], np.log(x68)[-1], 300) #Plotly has no method to plot a function, so I just took 300 evenly spaced xs (300 is just a random large number I picked).
yi = pp(xi)
fig.append_trace(go.Scatter(x=xi, y=yi, line=dict(color="Blue", dash="dash"), name="Spline fit"), row=1, col=1)
'''
fig.append_trace(go.Scatter(x=np.log(x68), y=np.log(y68), line=dict(color="Blue", dash="dash"), name="Spline fit"), row=1, col=1)


fig.append_trace(go.Histogram(x=x[:,index68]*100, nbinsx=30), row=1, col=2)

fig.update_xaxes(title_text="log(Entropy)", row=1, col=1)
fig.update_yaxes(title_text=r"$\log{\chi^2}$", row=1, col=1)

fig.update_xaxes(title_text="Weight(%)", row=1, col=2)
fig.update_yaxes(title_text="# of Structures", row=1, col=2)

#fig.update_yaxes(title_text="Log Relative Error", type="log", row=1, col=2)

fig.update_layout(title_text="L-curve Plots", font=dict(family="Arial", size=14)) #height=600, width=1000, 

fig.show()
""" 
figure(1)
plot(log(x68),log(y68),'m.')
xlabel('log(Entropy)')
ylabel('log(\chi^2)') """
#legend(strcat('pH 4.5 best \lambda= ', num2str(lambda45(index45))),strcat('pH 6.8 best \lambda= ', num2str(lambda68(index68))),strcat('pH 7.6 best \lambda= ', num2str(lambda76(index76))),'location', 'best')
#set(gca, 'FontName', 'Arial', 'FontSize', 14)   
#ylim([-6 0])


"""

figure(1)
plot(log(x68),log(y68),'k.','MarkerSize', 12);
hold on
fnplt(pp,'b--');
hold on
plot(log(x68(index68)),log(y68(index68)),'r.','MarkerSize', 12);
xlabel('log(Entropy)');
ylabel('log(\chi^2)');
legend('All other \lambda ','Spline fit',strcat('\lambda chosen= ', num2str(lambda68(index68))),'location', 'best');
set(gca, 'FontName', 'Arial', 'FontSize', 14);   
%xlim([-8 1]);
title('L-curve for best \lambda');

avg = mean(x(:,index68)*100);
devi = std(x(:,index68)*100);

figure(2)
hist(x(:,index68)*100,30);
%xlabel('Entropy');
xlabel('Weight(%)');
ylabel('# of structures');
%legend(strcat('Mean = ', num2str(avg)),'location', 'best');
txt = strcat('Mean = ', num2str(sprintf('%.3f',avg)));
text(5*10^0,10,txt);
txt1 = strcat('Standard Deviation = ', num2str(sprintf('%.3f',devi)));
text(5*10^0,11,txt1);
set(gca, 'FontName', 'Arial', 'FontSize', 14); 
ylim([0 20]);
title('Distribution of weights');
"""