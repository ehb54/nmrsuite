<META http-equiv="Content-Style-Type" content="text/css">
<LINK href="../sassie_style.css" rel="stylesheet" type="text/css"> 

###Diffusion Tensor-Guided Docking (ELMDOCK)

This module performs a rigid-body docking of macromolecular complexes guided by experimental rotational diffusion tensors using the program ELMDOCK [7,8].
Caution: this part is currently designed only for 2-domain/component systems. 
Please use the examples “diUb-AB-relax.txt” and “diUb-AB.pdb” to go through the process.

To select this module from the general ROTDIF-Web page (shown below), click the button “Diffusion Tensor-Guided Docking” 
![inline image](images/elmdock.fld/image001.png) 

The interface will look like:
![inline image](images/elmdock.fld/image002.png) 

*	Run name: The name of the folder to store of the current session
*	Relaxation Data File: Upload a text file containing spin relaxation data. The format of the data file is the same as spelled above in the ROTDIF & Dynamics Analysis, except that it has to contain data for two molecules (chainID), and the chainID letters must match those in the PDB file. 
*	PDB Coordinates File: Upload the starting coordinates file to be used for the docking. The file must contain coordinates for both molecules (chainID) listed in the relaxation data file.
*	PDB File model: If the coordinates file contains more than one structural model, please specify the model you want to analyze. The default model number is 1. 
*	Temperature: The desired temperature (in K) of the aqueous medium. The value of the temperature could vary from complex to complex and might need to be adjusted to avoid steric clashes. For the diUb_AB example provided we suggest 301 K 
*	Advanced model: As in the ELM Predictor, these options are designed for advanced adjustment of the hydration layer thickness and the water radius. It is recommended that you use the default values. 

Prior to running ELMDOCK, please make sure the project name is specified. If it is not specified or you would like to change the project name, click on the “head” icon ![inline image](images/elmdock.fld/image003.png) in the top right corner. The following selection window will appear:
![inline image](images/elmdock.fld/image004.png)

You can select an existing project from the pull-down list or click “New project” and create a new project by providing the project name and the description, for example:
![inline image](images/elmdock.fld/image005.png)

Then click the ![inline image](images/elmdock.fld/image006.png) button. The “Status” information will read like the following: ![inline image](images/elmdock.fld/image007.png)

Now you can run Diffusion Tensor-Guided Docking (ELMDOCK). Here is an example of input parameters; you can also try other values. ![inline image](images/elmdock.fld/image008.png)

Fill out the required fields, shown above, then click the Submit button. After few minutes, the page will look like this (wait till the molecular drawing appears on the bottom of the output screen) ![inline image](images/elmdock.fld/image009.png)

Output files accessible/downloadable directly from the GUI:

*	Detailed Results File: elmdock_log.out. This file contains the records of all the details
*	ELMDOCK results: ELMDOCK. This file contains the following information on the found solutions:
![inline image](images/elmdock.fld/image010.png)
 
*	ELMDOCK PDB file: out_dock.pdb. After the docking is completed. ELMDOCK will start Jmol viewer that allows you to visualize the results. Alternatively, you can download this coordinates file and view it in any molecular viewing program. Here is the view from Jmol. The temperature was set to 301K.
![inline image](images/elmdock.fld/image011.png)
For more options for display of the molecular image, you can right click on the molecule (see below). Typically, ELMDOCK returns multiple docked solutions reflecting the degeneracy associated with the relative positioning and orientation of the two proteins in a complex, reflecting the inherent symmetry of the rotational diffusion tensor, see [13]. By default, model 1 will be shown. You can select which solution/model (or all models) to display using the popup menu (as shown below). 
![inline image](images/elmdock.fld/image012.png)
As a comparison, here are the positions and orientations of the two proteins used in this example before the docking:
![inline image](images/elmdock.fld/image013.png)

###Exercise: 

This exercise includes diffusion-guided rigid-body docking of two ubiquitin molecules using 15N relaxation data for K48-linked di-ubiquitin [7, 9]. Use files “diUb-AB.pdb” (starting coordinates) and “diUb-AB-relax.txt” (relaxation data) to run ELMDOCK. (Note that diUb-AB.pdb contains two arbitrarily positioned ubiquitin molecules. The atom coordinates for each ubiquitin molecule are from PDB ID 1D3Z) 


###Reference(s) and Citations
1. O. Walker, R. Varadan, D. Fushman,[”Efficient and accurate determination of the overall rotational diffusion tensor of a molecule from 15N relaxation data using computer program ROTDIF,”](https://www.ncbi.nlm.nih.gov/pubmed/15140445) J. Magn. Reson. (2004) 168, 336-345.
2. K. Berlin, A. Longhini, T. K. Dayie, D. Fushman, [“Deriving Quantitative Dynamics Information for Proteins and RNAs using ROTDIF with a Graphical User Interface”,](https://www.ncbi.nlm.nih.gov/pubmed/24170368) J Biomol NMR (2013) 57, 333-352.
3. D. Fushman, S. Cahill, D. Cowburn, [“The main chain dynamics of the dynamin Pleckstrin Homology (PH) domain in solution: Analysis of 15N relaxation with monomer/dimer equilibration,”](https://www.ncbi.nlm.nih.gov/pubmed/9054979) J. Mol. Biol. 266 (1997) 173-194. 
4. J. B. Hall, and D. Fushman, [“Characterization of the overall and local dynamics of a protein with intermediate rotational anisotropy: Differentiating between conformational exchange and anisotropic diffusion in the B3 domain of protein G,“](https://www.ncbi.nlm.nih.gov/pubmed/12975584) J. Biomol. NMR (2003) 27, 261-275. 
5. J. B. Hall, D. Fushman, [“Variability of the 15N chemical shielding tensors in the B3 domain of protein G from 15N relaxation measurements at several fields. Implications for backbone order parameters,”](https://pubs.acs.org/doi/abs/10.1021/ja060406x) J. Am. Chem. Soc. (2006) 128, 7855-70.
6. Y. Ryabov, C. Geraghty, A. Varshney, D. Fushman, [“An efficient computational method for predicting rotational diffusion tensors of globular proteins using an ellipsoid representation,”](https://pubs.acs.org/doi/abs/10.1021/ja062715t) J. Am. Chem. Soc. (2006) 128, 15432-15444.
7. Y. Ryabov, D. Fushman, [“Structural assembly of multidomain proteins and protein complexes guided by the overall rotational diffusion tensor,”](https://www.ncbi.nlm.nih.gov/pubmed/17550252) J. Am. Chem. Soc. (2007) 129, 7894-7902.
8. K. Berlin, D. P. O’Leary, D. Fushman, [“Fast Approximations of the Rotational Diffusion Tensor and their Application to Structural Assembly of Molecular Complexes”,](https://www.ncbi.nlm.nih.gov/pubmed/21604302) Proteins (2011) 79, 2268-2281.
9. R. Varadan, O. Walker, C. Pickart, D. Fushman, [“Structural properties of polyubiquitin chains in solution,”](https://www.ncbi.nlm.nih.gov/pubmed/12460567) J. Mol. Biol. (2002) 324, 637-647
10. D. Fushman, R. Ghose, D. Cowburn, [“The effect of finite sampling on the determination of orientational properties: A theoretical treatment with application to interatomic vectors in proteins,”](https://pubs.acs.org/doi/abs/10.1021/ja001128j) J. Am. Chem. Soc. 122 (2000) 10640-9
11. D. Fushman, [“Determining protein dynamics from 15N relaxation data by using DYNAMICS”,](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4361738/) in Protein NMR Techniques, Third Edition, Eds. A. Shekhtman, D. S. Burz; Methods in Molecular Biology, 2012, Volume 831, 485-511, Springer Science, DOI: 10.1007/978-1-61779-480-3_24. PubMed: PMC4361738
12. R. C. Weast, [Handbook of Chemistry and Physics](http://hbcponline.com/faces/contents/ContentsSearch.xhtml), 59th ed.; CRC Press: West Palm Beach, FL, 1978.
13. D. Fushman, R. Varadan, M. Assfalg, O.Walker,  [“Determining domain orientation in macromolecules by using spin-relaxation and residual dipolar coupling measurements,”](https://s3.amazonaws.com/academia.edu.documents/112727/progress.pdfAWSAccessKeyId=AKIAIWOWYYGZ2Y53UL3A&Expires=1556145392&Signature=hJCC7jXX39HvdpwzfOW6hmA1nCs%3D&response-contentdisposition=inline%3B%20filename%3DDetermining_domain_orientation_in_macrom.pdf) Progress in NMR Spectroscopy, (2004) 44, 189-214.

--
####[Return to Introduction](./rotdif_intro.html)
####[Tutorial of ROTDIF & DYNAMICS](./rotdif_all.html)
####[Tutorial of ELM](./elm.html) 

<a href=#>Return to top</a>

<footer>
  <ul>
  Supported via NSF grants OAC-1739549 (to D.F.) OAC-1740087 (to E.B.) and Science Gateways Community Institute (summer internship to Y.C.)
  </ul>
</footer> 
