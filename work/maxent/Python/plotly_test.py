#from maxent_scripts.RunMaxEntropy import summation

import numpy as np

x68 = np.array([0.0168346992, 0.0110934191, 0.00713835199, 0.0044982689, 0.00278410992, 0.00169726938, 0.00102181541, 0.000608922197, 0.000359908896, 0.000211349185, 0.000123478586, 7.18547561e-05, 4.16852348e-05, 2.41256455e-05, 1.39375103e-05, 8.04056733e-06, 4.63368518e-06, 2.66817644e-06, 1.53544699e-06, 8.83184719e-07, 5.07824484e-07, 2.9191628e-07, 1.6776984e-07, 9.64055114e-08, 5.53909181e-08, 3.18226533e-08, 1.8281204e-08, 1.05014861e-08, 6.03225647e-09, 3.46494228e-09, 1.99022633e-09, 1.1431457e-09, 6.56591233e-10, 3.77124226e-10, 2.16606032e-10, 1.24409699e-10, 7.14555963e-11, 4.10408433e-11, 2.35719113e-11])

y68 = np.array([0.13197019, 0.14221298, 0.15151641, 0.15970526, 0.16671678, 0.17257981, 0.1773858, 0.18126095, 0.1843439, 0.18677029, 0.18866368, 0.19013125, 0.19126283, 0.19213181, 0.19279704, 0.19330507, 0.19369235, 0.19398716, 0.19421134, 0.19438167, 0.19451101, 0.19460917, 0.19468365, 0.19474014, 0.19478299, 0.19481547, 0.19484009, 0.19485876, 0.19487291, 0.19488364, 0.19489177, 0.19489793, 0.1949026, 0.19490614, 0.19490883, 0.19491086, 0.1949124, 0.19491357, 0.19491445])

index68 = 0

import os
run_directory = "result_of_maxent"
x = np.loadtxt(os.path.join(run_directory, "x.txt"))

lambda_ = np.loadtxt(os.path.join(run_directory, "lambda.txt"))
lambda68 = lambda_[1:] #from results_of_maxent

x = x[:, 1:]

chosen_lambda = lambda68[index68]

line_plot = {
        "data": [
            {
                "x": np.log(x68).tolist(), # converting to list because the genapp plotly does not accept numpy arrays
                "y": np.log(y68).tolist(),
                "mode": "markers",
                "marker": {
                    "color": "Black"
                },
                "name": "All other λ" # rstring to use latex characters
            },
            {
                "x": [np.log(x68[index68])],
                "y": [np.log(y68[index68])],
                "mode": "markers",
                "line": {
                    "color": "Red"
                },
                "name": r"$\lambda\text{{ chosen = }}{}$".format(round(chosen_lambda, 5))
                # the name variable was weird to work with because I needed to use a latex string but also add a variable in, so I had to use the slightly archaic .format and add an escape second bracket to each of the latex brackets
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
        ],
        "layout": {
                "title": r"$\text{L-curve for best }\lambda$",
                "xaxis": {
                    "title": "log(Entropy)"
                },
                "yaxis": {
                    "title": r"$\log{\chi^2}$"
                }
        }
    }


histogram = {
        "data": [
            {
                "type": "histogram",
                "nbinsx": 30, # number of bins
                "x": (x[:,index68]*100).tolist()
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

import plotly.io as pio

pio.show(line_plot)
pio.show(histogram)

