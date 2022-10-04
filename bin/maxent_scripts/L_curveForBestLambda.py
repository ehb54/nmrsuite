#!/opt/miniconda3/bin/python

#This script tells you which of the values of lambda is the best lambda for you to use.
#You can have that information by looking at the variable index68 for
#example. This scrip plots the L-curve got by Log(chi^2) x Log(entropy).

import os, json
import numpy as np
from scipy.interpolate import UnivariateSpline
from csaps import CubicSmoothingSpline

#import plotly
#from plotly.subplots import make_subplots
import plotly.graph_objects as go


def index_max(array): #Return the (first) index of the maximum value and the maximum value in a given array
    _max = np.max(array) #'max' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
    index = np.where(array == _max)[-1][0]
    return ([index, _max])

def lcurve_points(y, A, x, lambda_):

    #x = xsol(:,2:end)
    Ax = np.matmul(A, x)

    w = np.ones((np.shape(A)[1], 1)) * np.shape(A)[1] #vector filled with the number of maximum size(A,1) of A

    length = np.shape(Ax)[1]

    r = np.empty((length, np.shape(Ax)[0]))
    rn = np.empty(length)
    xval = np.empty(length)
    f = np.empty(length)

    for i in range (length): #Returns Chi^2 for each lambda.
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
    
    x_par = np.log(xval) #x has to be strictly increasing
    #y_par = np.flip(np.log(f))
    y_par = np.log(f)

    indices = np.argsort(x_par)
    x_par = x_par[indices]
    y_par = y_par[indices]

    """ pp = UnivariateSpline(x=x_par, y=y_par, s=0.995)
    y_prime = pp(x=np.log(xval), nu=2) """

    xval = xval[indices]
    f = f[indices]

    pp = CubicSmoothingSpline(x_par, y_par, smooth=0.995).spline
    p_second = pp.derivative(nu=2)
    y_prime = p_second(np.log(xval))

    y_prime = np.flip(y_prime)

    [index,v] = index_max(y_prime)
    index += 1

    full_x = np.linspace(np.min(np.log(xval)),np.max(np.log(xval)), 200)
    ydf2_full = p_second(full_x)

    xval = np.flip(xval)
    f = np.flip(f)

    return [index, xval, f, pp]


def L_curve_for_best_lambda (run_directory, sigma):
    A = np.loadtxt(os.path.join(run_directory, "A.txt"))
    lambda_ = np.loadtxt(os.path.join(run_directory, "lambda.txt")) #'lambda' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
    x = np.loadtxt(os.path.join(run_directory, "weights_for_all_lambdas.txt"))
    y = np.loadtxt(os.path.join(run_directory, "data.txt"))

    y_prime = np.array([])

    x = x[:, 1:] #excluding the results with lambda = 0
    lambda68 = lambda_#[1:] #from results_of_maxent
    
    #return lcurve_points(y, A, x, lambda68)
    [index,x68,y68, pp] = lcurve_points(y, A, x, lambda68)

    lambda_ = lambda68 #The lambda will have the size shortened, but we want to still have the full lambda68
    #Now the curve will be recalculated several times until no outlier is
    #present anymore.
    while True:
        maxf = np.max(pp(np.log(x68)))
        minf = np.min(pp(np.log(x68)))
        devif = np.abs(maxf-minf) #It takes into account the range of the curve
        sigma = sigma*devif #By doing it like that, this will work for different systems
        outliers = np.where(np.abs(np.log(y68)-pp(np.log(x68))) > sigma)[0]
        #print (outliers)
        
        np.delete(x, outliers, axis=1)
        np.delete(lambda_, outliers)

        oldx68 = x68
        x68 = []
        y68 = []
        index = []

        [index,x68,y68,pp]=lcurve_points(y, A, x, lambda_)  #Redo the curve now with shorter x and lambda
        
        #If the length of x68 and old68 remains the same after a new curve is
        #calculated, it means the curve is not changing anymore.
        if np.size(x68) == np.size(oldx68):
            break

    index68 = np.where(lambda68 == lambda_[index])[0]
    index68 += 1

    chosen_lambda = lambda_[index]

    with open(os.path.join(run_directory, "index.txt"), "w") as f: f.write(str(index68[0]))

    # X: y
    # Y: A*x[:,index68]

    line_plot = {
        "data": [
            {
                "x": np.log(x68).tolist(), # converting to list because the genapp plotly does not accept numpy arrays
                "y": np.log(y68).tolist(),
                "mode": "markers",
                "marker": {
                    "color": "Black"
                },
                "name": "All other &#955;" #"name": r"$\text{All other }\lambda$" # rstring to use latex characters
            },
            {
                "x": np.log(x68).tolist(),
                "y": np.log(y68).tolist(),
                "line": {
                    "color": "Blue",
                    "dash": "dash"
                },
                "name": "Spline fit"
            },
                        {
                "x": [np.log(x68[index-1])],
                "y": [np.log(y68[index-1])],
                "mode": "markers",
                "line": {
                    "color": "Red"
                },
                "name": "&#955; chosen = {}".format(round(chosen_lambda, 5))
                # the name variable was weird to work with because I needed to use a latex string but also add a variable in, so I had to use the slightly archaic .format and add an escape second bracket to each of the latex brackets
            }
        ],
        "layout": {
                "title": "L-curve for best &#955;",
                "xaxis": {
                    "title": "log(Entropy)"
                },
                "yaxis": {
                    "title": "log(&#967;<sup>2</sup>)"
                }
        }
    }


    histogram = {
            "data": [
                {
                    "type": "histogram",
                    "nbinsx": 30, # number of bins
                    "x": (x[:,index]*100).tolist()
                }
            ],
            "layout": {
                    "title": "Distribution of Weights",
                    "xaxis": {
                        "title": "Weight(%)"
                    },
                    "yaxis": {
                        "title": "# of Structures",
                        "range": [0, 20]
                    }
            }
        }

    scatter_x = y.tolist()
    scatter_y = np.reshape(np.matmul(A, x[:,index68]), -1).tolist()

    # The below code finds the lowest value between x and y

    minimum_x = min(scatter_x)
    minimum_y = min(scatter_y)
    first_point = min(minimum_x, minimum_y)

    # Same thing as above but for maximum

    maximum_x = max(scatter_x)
    maximum_y = max(scatter_y)
    second_point = max(maximum_x, maximum_y)


    scatter_plot = {
        "data": [
            {
                "x": [first_point, second_point],
                "y": [first_point, second_point],
                "mode": "lines",
                "marker": {
                    "color": "Black"
                },
                "name": "Agreement"
            },
            {
                "x": scatter_x,
                "y": scatter_y,
                "mode": "markers",
                "marker": {
                    "color": "Green"
                },
                "name": "Data"          
            }
        ],
        "layout": {
            "title": "Agreement between Experimental and Predicted Data",
            "xaxis": {
                "title": "Experimental Data"
            },
            "yaxis": {
                "title": "Predicted Data",
            }
        }
    }



    line_plot_figure = go.Figure(line_plot)
    line_plot_figure.write_image(os.path.join(run_directory, "l_curve.png"))


    """ fig = make_subplots(rows=1, cols=2, subplot_titles=("Plot 1", "Plot 2"))

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

    fig.show() """

    return line_plot, histogram, scatter_plot

#L_curve_for_best_lambda("result_of_maxent")

"""
figure(1)
plot(log(x68),log(y68),'m.')
xlabel('log(Entropy)')
ylabel('log(\chi^2)') """
#legend(strcat('pH 4.5 best \lambda= ', num2str(lambda45(index45))),strcat('pH 6.8 best \lambda= ', num2str(lambda68(index68))),strcat('pH 7.6 best \lambda= ', num2str(lambda76(index76))),'location', 'best')
#set(gca, 'FontName', 'Arial', 'FontSize', 14)
#ylim([-6 0])
