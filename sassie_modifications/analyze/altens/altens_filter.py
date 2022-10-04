'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os, sys, string, locale, math
#import input_filter
import sassie.interface.input_filter as input_filter
import sasmol.sasmol as sasmol

def check_altens(variables,**kwargs):

    runname = variables['runname'][0]
    pdbfile = variables['pdbfile'][0]
    dcdfile = variables['dcdfile'][0]

    rdc_input_file = variables['rdc_input_file'][0]

    residue_list_file_flag = variables['residue_list_file_flag'][0]
    residue_list_file = variables['residue_list_file'][0]

    use_monte_carlo_flag = variables['use_monte_carlo_flag'][0]
    number_of_monte_carlo_steps = variables['number_of_monte_carlo_steps'][0]

    seed = variables['seed'][0]

    error=[]
    error = input_filter.check_name(runname)

    if(error!=[]):
        return error

    error=input_filter.check_file_exists(pdbfile)
    if(len(error) != 0):
        #error.append('input pdb file, '+pdbfile+', does not exist')
        return error
    ev,value=input_filter.check_pdb_dcd(pdbfile,'pdb')
    if(ev == 0):
        error.append('check input pdb file: '+pdbfile)
        return error
    if(value == 0):
        error.append( 'input pdb file, '+pdbfile+', is not a valid pdb file')
        return error
    try:
        m1 = sasmol.SasMol(0)
        m1.read_pdb(pdbfile)
        number_of_frames = m1.number_of_frames()
        print '> found '+str(number_of_frames)+' frames in PDB file'
    except:
        error.append('could not open PDB file '+pdbfile+' to check number of frames')
        return error
    if(number_of_frames < 1):
        error.append('PDB file has no frames : '+pdbfile)
        return error

    error=input_filter.check_file_exists(dcdfile)
    if(len(error) != 0):
        #error.append('input dcd file, '+dcdfile+', does not exist')
        return error

    error=input_filter.check_file_exists(rdc_input_file)
    if(len(error) != 0):
        return error

    if (residue_list_file_flag):
        error=input_filter.check_file_exists(residue_list_file)
        if(len(error) != 0):
            return error

    '''
    ev,value=input_filter.check_pdb_dcd(dcdfile,'dcd')
    if(ev == 0):
        error.append('check input dcd filename : '+dcdfile)
        return error
    elif(value == 0):
        ev,value=input_filter.check_pdb_dcd(dcdfile,'pdb')
        if(value ==0):
            error.append( 'input dcd file, '+dcdfile+', is not a valid pdb file')
            return error
        infile_type = 'pdb'
    else:
        infile_type = 'dcd'

    if(infile_type == 'dcd'):
        value=input_filter.certify_pdb_dcd(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and dcd file '+dcdfile+', are not compatiable')
            return error
    elif(infile_type == 'pdb'):
        value=input_filter.certify_pdb_pdb(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and dcd file '+dcdfile+', are not compatiable')
            return error
    '''
    if dcdfile[-3:] == 'dcd':
        infile_type = 'dcd'
        ev,value=input_filter.check_pdb_dcd(dcdfile,'dcd')
        if(ev == 0):
            error.append('check input dcd filename : '+dcdfile)
            return error
        elif(value == 0):
            error.append( 'input dcd file, '+dcdfile+', is not a valid dcd file')
            return error

        value=input_filter.certify_pdb_dcd(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and dcd file '+dcdfile+', are not compatible')
            return error

    elif dcdfile[-3:] == 'pdb':
        infile_type = 'pdb'
        ev,value=input_filter.check_pdb_dcd(dcdfile,'pdb')
        if(ev == 0):
            error.append('check input dcd filename : '+dcdfile)
            return error
        elif(value == 0):
            error.append( 'input dcd file, '+dcdfile+', is not a valid pdb file')
            return error

        value=input_filter.certify_pdb_pdb(pdbfile,dcdfile)
        if(value == 0):
            error.append('input pdb file '+pdbfile+' and pdb file '+dcdfile+', are not compatible')
            return error
    try:
        m1 = sasmol.SasMol(0)
        if(infile_type == 'dcd'):
            dcdinputfile = m1.open_dcd_read(dcdfile)
            number_of_frames = dcdinputfile[2]
        elif(infile_type == 'pdb'):
            m1.read_pdb(dcdfile)
            number_of_frames = m1.number_of_frames()
    except:
        error.append('could not open dcd file '+dcdfile+' to check number of frames')
        return error
    if(number_of_frames < 1):
        error.append('dcd file has no frames : '+dcdfile)
        return error

    msg,rdc_type,rdc_resinfo = check_altens_file(rdc_input_file)
    if(len(msg) != 0):
        error.append(msg)
        return error

    msg = check_rdc_type_supported(rdc_type,rdc_input_file)
    if(len(msg) != 0):
        error.append(msg)
        return error

    msg = check_rdc_data_column(rdc_resinfo,rdc_input_file)
    if(len(msg) != 0):
        error.append(msg)
        return error

    msg = check_rdc_pdb(rdc_resinfo,rdc_input_file,pdbfile)
    if(len(msg) != 0):
        error.append(msg)
        return error

    if (residue_list_file_flag):
        msg,exclusion_type,exclusion_resinfo = check_altens_file(residue_list_file)
        if (len(msg) != 0):
            error.append(msg)
            return error
        else: 
            msg = check_exclusion_type(rdc_type, exclusion_type, rdc_input_file, residue_list_file)       
            if(len(msg) != 0):
                error.append(msg)
                return error
            else:
                msg = check_exclusion_data(exclusion_resinfo, rdc_resinfo, residue_list_file, rdc_input_file)
                if(len(msg) != 0):
                    error.append(msg)
                    return error

    if (use_monte_carlo_flag):
        if (number_of_monte_carlo_steps <= 0):
            error.append("Number of Monte Carlo steps is smaller than or equal to 0")
            return error

    return error

def check_altens_file(altens_infile):
    '''
    Check the format of headers and data entries to define RDC data 
    '''

    msg=""

    inlines = open(altens_infile,'rU').readlines()
#    inlines = open(altens_infile,'r').readlines()

    rdc_type = []
    rdc_resinfo = []
    count_line = 0
    count_header = 0
    count_data = 0

# Check if the number of line is zero
    if (len(inlines) == 0 ):
        msg += "Input error: no data found in " + str(altens_infile)
        return msg, rdc_type, rdc_resinfo

    for line in inlines:
        lin = string.split(line)
        
# Check if blank line exists
        if (len(lin) == 0):
            msg += "Input error: found blank line in " + str(altens_infile)
            return msg, rdc_type, rdc_resinfo
        else: # not blank line

# Check if the first line is header
            if (count_line == 0 ):
                if (lin[0][0] != "#"): 
                    msg += "Input error: the first line is not header starting # in " + str(altens_infile)
                    return msg, rdc_type, rdc_resinfo

            count_line += 1

            if (lin[0][0] == "#" and len(lin[0]) > 1):
                lin[0] = lin[0].replace("#","")
                if (lin[0] == ""):
                    lin[0] = "#"
                else:
                    lin.insert(0,"#") 

            if (lin[0] == "#"):
                count_header += 1
                
# Check if headers have enough info. to define chain, atom1, and atom2  
                if ( len(lin) < 4 ):
                    msg += "Input error: check missing information in header in " + str(count_line) + "th line of " + str(altens_infile)
                    return msg, rdc_type, rdc_resinfo
                else:
                    chain_id, atom_1, atom_2 = lin[1:4]
                    rdc_type.append([chain_id, atom_1, atom_2])
            else:  # not header
                count_data += 1

# Check if the data is number....
                for k in range(len(lin)):
                    try: 
                        locale.atof(lin[k])
                        #test_lin = float(lin[k])
                        #if ( math.isnan(test_lin) or math.isinf(test_lin)):
                        #    msg += "Input error: invalid data found " +str(lin) + " of " + str(altens_infile)
                        #    return msg, rdc_type, rdc_resinfo
                    except:
                        msg += "Input error: non-numeric data found at line " +str(count_line) + " of " + str(altens_infile)
                        return msg, rdc_type, rdc_resinfo
                    else:
                        test_lin = float(lin[k])
                        if ( math.isnan(test_lin) or math.isinf(test_lin)):
                            msg += "Input error: non-numeric data found at line " +str(count_line) + " of " + str(altens_infile) 
                            return msg, rdc_type, rdc_resinfo

                rdc_resinfo.append([ len(lin), chain_id, int(locale.atof(lin[0])), atom_1, atom_2])

# Check if header exits 
    if ( count_header == 0 ):
        msg += "Input error: could not find any header in " + str(altens_infile)
        return msg,rdc_type,rdc_resinfo

# Check if duplicated header exists 
    if ( count_header > 1 ):
        for i in range(count_header-1):
            for j in range(i+1, count_header):    
                if (rdc_type[i][0] == rdc_type[j][0] and rdc_type[i][1] == rdc_type[j][1] and rdc_type[i][2] == rdc_type[j][2]):
                    msg += "Input error: found duplicated headers in" + str(altens_infile)
                    return msg, rdc_type,rdc_resinfo
                if (rdc_type[i][0] == rdc_type[j][0] and rdc_type[i][1] == rdc_type[j][2] and rdc_type[i][2] == rdc_type[j][1]):
                    msg += "Input error: found duplicated headers in" + str(altens_infile) 
                    return msg, rdc_type,rdc_resinfo

# Check if RDC data entry exists
    if (count_data == 0):
        msg += "Input error: could not find any residue information " + str(altens_infile)

    return msg, rdc_type, rdc_resinfo

def check_rdc_type_supported(rdc_type,rdc_file):

    '''
    Check if the atom pairs in the header of rdc_input_file supported 
    '''

    msg = ""
    altens_valid_pairs = ["HN", "NH", "HC", "CH", "CN", "NC", "CC"]

    for i in range(len(rdc_type)):
        atom_pair = str(rdc_type[i][1][0]) + str(rdc_type[i][2][0])
        if ( atom_pair not in altens_valid_pairs ):
            msg += "Input error: invalid atom pair " + str(rdc_type[i][1:3]) + " for chain " + str(rdc_type[i][0]) + " was chosen in " + str(rdc_file)   
        return msg

    return msg

def check_rdc_data_column(rdc_resinfo,rdc_file):

    '''
    Check if the number of data columns in experimental RDC file is greater than 1
    and Check if the number of input values (residid, RDC, error) remains same for all residues
    '''
    msg=""
    num_column = 0

    for i in range(len(rdc_resinfo)):

# Check if RDC data is given for all resiudes in the input rdc file
        if (rdc_resinfo[i][0] < 2):
            msg += "Input error: could not find experimental RDC value in " + str(rdc_file)
            return msg

# Check if any unnecessary value is given in the input rdc file
        elif (rdc_resinfo[i][0] > 3):
            msg += "Input error: invalid data found in " + str(rdc_file) 
            return msg
        else: 
            if (i == 0):
                num_column = rdc_resinfo[i][0]
            else:

# Check if RDC errors are defined for all or not 
                if (num_column != rdc_resinfo[i][0]):
                    msg += "Input error: inconsistent number of data entry found "
                    return msg  

# Check if duplicated resiudes with same vector type exist
    for i in range(len(rdc_resinfo)-1):
        for j in range(i+1,len(rdc_resinfo)):
            if (rdc_resinfo[i][1:5] == rdc_resinfo[j][1:5]):
                msg += "Input error: found duplicated RDC entry for resiude " + str(rdc_resinfo[i][1:5]) 
                return msg
    return msg

def check_rdc_pdb(rdc_resinfo,rdc_file,pdbfile):

    '''
    Check if the number of data columns in experimental RDC file is greater than 1
    and Check if the number of input values (residid, RDC, error) remains same for all residues
    '''
    msg=""
    num_column = 0

    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdbfile)

    seg_pdb = m1.chains()
    res_pdb = m1.resids()
    natoms = m1.natoms()

    for i in range(len(rdc_resinfo)):

# Check if segname in input rdc file exists in pdb file
        if (rdc_resinfo[i][1] not in seg_pdb):
            msg += "Input error: could not find chain " +str(rdc_resinfo[i][1])+ " in input pdb file due to wrong chain name or inappropriate usage of header section."
            return msg
        else: 
# Check if resid in input rdc file exists in pdb file
            if (rdc_resinfo[i][2] not in res_pdb):
                msg += "Input error: could not find residue " +str(rdc_resinfo[i][1:3])+ " in input pdb file" 
                return msg
# Check if atom names defined in the header exists in pdb file
    atom1_exist = [["",0] for x in range(len(rdc_resinfo))]
    atom2_exist = [["",0] for x in range(len(rdc_resinfo))]
    
    for i in range(natoms):
        atom_pdb = m1.name()[i]
        for j in range(len(rdc_resinfo)):
            chain_id, res_id, atom_1, atom_2 = rdc_resinfo[j][1:5]
            if (m1.chain()[i] == chain_id and m1.resid()[i] == res_id):
                if (atom_pdb == atom_1 ):
                    atom1_exist[j][0] = atom_pdb
                    atom1_exist[j][1] = 1
                if (atom_pdb == atom_2 ):
                    atom2_exist[j][0] = atom_pdb
                    atom2_exist[j][1] = 1
                if ( m1.resname()[i] == "GLY" and (atom_1 == "CA" or atom_2 == "CA")):
                    if ( atom_pdb == "CA" ):
                        atom1_exist[j][0] = atom_pdb
                        atom1_exist[j][1] = 1
                    if (atom_pdb == "1HA" or atom_pdb == "HA1"):
                        atom2_exist[j][0] = atom_pdb
                        atom2_exist[j][1] = 1

    for i in range(len(rdc_resinfo)):
        if ( atom1_exist[i][1] == 0 ):
            msg += "Input error: could not find atom " + str(rdc_resinfo[i][3]) + " within residue " + str(rdc_resinfo[i][2]) + " of chain_id = "+ str(rdc_resinfo[i][1])+ " in input pdb"
            return msg
        if ( atom2_exist[i][1] == 0 ):
            msg += "Input error: could not find atom " + str(rdc_resinfo[i][4]) + " within residue " + str(rdc_resinfo[i][2]) + " of chain_id = "+ str(rdc_resinfo[i][1]) + " in input pdb"  
            return msg
    return msg

def check_exclusion_type(rdc_type, exclusion_type, rdc_file, exclusion_file):

    '''
    Check if the headers in residue_list_file matches headers in rdc_input_file. Order of headers doesn't matter.
    '''
    msg = ""

    for i in range(len(exclusion_type)):
        type_missing = 0 
        for j in range(len(rdc_type)):
            if (exclusion_type[i][0] == rdc_type[j][0]):
                if (exclusion_type[i][1] == rdc_type[j][1] and exclusion_type[i][2] == rdc_type[j][2]):
                    type_missing = 1
                elif (exclusion_type[i][1] == rdc_type[j][2] and exclusion_type[i][2] == rdc_type[j][1]):
                    type_missing = 1
        if (type_missing == 0):
            msg += "Exclusion input error: could not match residue" +str(exclusion_type[0:2]) +" with " + str(rdc_file)
            return msg

    return msg

def check_exclusion_data(exclusion_resinfo, rdc_resinfo, exclusion_file, rdc_file):

    '''
    Check if residues in residue_list_file for exclusion exist in the rdc_input_file
    '''
    msg=""
    num_column = 0
    exclusion_missing = 0

# Check if too many residues were excluded

    if (len(exclusion_resinfo) >= len(rdc_resinfo)):
        msg += "Exclusion input error: all RDC data will be excluded for ALTENS analysis"
        return msg
    number_res = len(rdc_resinfo) - len(exclusion_resinfo) 
    if ( number_res < 4 ):
        msg += "Exclusion input error: After exclusion of residues, there left only " + str(number_res) + " that istoo small for analysis"   
        return msg 

    for i in range(len(exclusion_resinfo)):

# Check if only resid is given in data entry
        if (exclusion_resinfo[i][0] > 1):
            msg += "Exclusion input error: invalid data found in " + str(exclusion_file) 
            return msg
        else:

# Check if the residue in exclusion list can be found in input rdc file
            exclusion_missing = 0
            for j in range(len(rdc_resinfo)):
                if (exclusion_resinfo[i][1:3] == rdc_resinfo[j][1:3]):
                    if (exclusion_resinfo[i][3] == rdc_resinfo[j][3] and exclusion_resinfo[i][4] == rdc_resinfo[j][4]):
                        exclusion_missing = 1
                    elif (exclusion_resinfo[i][3] == rdc_resinfo[j][4] and exclusion_resinfo[i][4] == rdc_resinfo[j][3]):
                        exclusion_missing = 1
            if (exclusion_missing == 0):
                msg += "Exclusion input error: could not find a residue" + str(exclusion_resinfo[i][1:3]) + " in " + str(rdc_file)
                return msg
 
# Check if duplicated resiudes with same vector type exist
    for i in range(len(exclusion_resinfo)-1):
        for j in range(i+1,len(exclusion_resinfo)):
            if (exclusion_resinfo[i][1:5] == exclusion_resinfo[j][1:5]):
                msg += "Exclusion Input error: found duplicated RDC entry for resiude " + str(exclusion_resinfo[i][1:5])
                return msg

    return msg

