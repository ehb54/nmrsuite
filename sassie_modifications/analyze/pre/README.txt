Inputs:

1. 1DEZ_f.pdb (file)
This is a PDB (Protein Data Bank) file that contains the protein's atomic coordinates and secondary structure.

2. ratio_4_SLfit_test_Q49.txt (file)
Label: 
Help:

3. T2dia (float)
Label:
Unit:
Help:
Default:
Min:
Max:

4. Htime (float)
Label:
Unit:
Help:
Default:
Min:
Max:

5. freq (float)
Label:
Unit:
Help:
Default:
Min:
Max:

6. TAUc (float)
Label:
Unit:
Help:
Default:
Min:
Max:

Outputs:

1. Chi2 (array containing a float)
The sum of the squared errors of the model generated that maps residue numbers to ratios.

2. Position of the Spin Label (array containing three floats)
A prediction of the X, Y, and Z coordinates of the spin label.

3. Figure 1 (graph)
A plot of residue numbers vs. ratios, in which experimental values are plotted in a bar graph and prediction/actual fit values are plotted as a line graph with filled in points for the former and empty ones for the latter.

4. Figure 2 (graph)
A plot of distances (to the spin label?) vs. ratios. The relationship between distances and experimental values is shown with a scatterplot, while the relationship between distances and predicted values is shown with a curve of best fit.

5. Figure 3 (graph)
A plot of experimental vs. predicted values, along with a diagonal line showing the locations of perfect fit. The qR-factor and correlation coefficient are labeled clearly with the title at the top.


Test Cases:
I have only tested this program with the ratio_4_SLfit_test_Q49.txt and 1DEZ_f.pdb files, both of which are contained within this folder.
