#!/opt/miniconda3/bin/python
import numpy as np
from matplotlib import pyplot as plt

A = np.loadtxt("result_of_maxent/A.txt")
lambda_ = x = np.loadtxt("result_of_maxent/lambda.txt") # 'lambda' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
x = np.loadtxt("result_of_maxent/x.txt")
y = np.loadtxt("result_of_maxent/y.txt")

with open('index.txt', 'r') as f: index68 = int(f.read())

xsolplot = x[:,index68]*100

low_limit = 1.5 # pick a low limit weight appropriate for your situation
a = 1

x = np.arange(1, len(xsolplot) + 1)
y = xsolplot

fig = plt.figure()

plt.plot(x, y)

plt.title('Max Entropy Result', fontsize=15)
plt.xlabel('Structure', fontsize=15)
plt.ylabel('Weight (%)', fontsize=15)

plt.ticklabel_format(useOffset=False)

plt.tight_layout()

plt.show()
