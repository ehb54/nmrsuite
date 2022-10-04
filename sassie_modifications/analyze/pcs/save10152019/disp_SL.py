#Python script to display SL in PyMol
from pymol import cmd
#paramagnetic center parameters
resnum = 999
coor = [63.763839205291006,-91.63064830827406,-11.426006879530368]
atsize = 10
chainID = ""
atlabel = ""
atcolor = "orange"
atelement = "None"
#create pseudoatom and show
SL_atom = "SL"
cmd.pseudoatom(SL_atom, resi = resnum, chain = chainID, b = atsize, color = atcolor, elem = atelement, label = atlabel)
cmd.show("spheres",SL_atom)