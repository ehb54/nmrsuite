import numpy as np

def write_pymol(resnum = None, atcoord = None, atsize = 10, atcolor = 'orange', atlabel = '', chainID = 'A', param_type = 0):
    #TODO, Cheol: ChainID should come from pdb file
    if param_type == 0:
        atelement = 'None'

    output = ""
    output += "#Python script to display SL in PyMol"
    output += "\n"
    output += "from pymol import cmd"
    output += "\n"
    output += "#paramagnetic center parameters"
    output += "\n"
    output += ("resnum = " + str(resnum))
    output += "\n"
    output += ("coor = [" + str(atcoord[0]) + "," + str(atcoord[1]) + "," + str(atcoord[2]) +"]")
    output += "\n"
    output += ("atsize = " + str(atsize))
    output += "\n"
    output += ("chainID = \"" + str(chainID) + "\"")
    output += "\n"
    output += ("atlabel = \"" + str(atlabel) + "\"")
    output += "\n"
    output += ("atcolor = \"" + atcolor + "\"")
    output += "\n"
    output += ("atelement = \"" + atelement + "\"")
    output += "\n"
    output += "#create pseudoatom and show"
    output += "\n"
    output += "SL_atom = \"SL\""
    output += "\n"
    output += "cmd.pseudoatom(SL_atom, resi = resnum, chain = chainID, b = atsize, color = atcolor, pos = coor, elem = atelement, label = atlabel)"
    output += "\n"
    output += "cmd.show(\"spheres\",SL_atom)"
    
    return output
