import sys


pdb_file = sys.argv[1]
suspects = [["MET"," HB3","HB1"], \
            ["MET"," HG3","HG1"], \
            ["GLY"," HA3","HA1"], \
            ["TYR"," HB3","HB1"], \
            ["GLU"," HB3","HB1"], \
            ["GLU"," HG3","HG1"], \
            ["ASP"," HB3","HB1"], \
            ["ARG"," HB3","HB1"], \
            ["ARG"," HG3","HG1"], \
            ["ARG"," HD3","HD1"], \
            ["PHE"," HB3","HB1"], \
            ["SER"," HB3","HB1"], \
            ["SER"," HG ","HG1"], \
            ["GLN"," HB3","HB1"], \
            ["TRP"," HB3","HB1"], \
            ["GLN"," HG3","HG1"], \
            ["LEU"," HB3","HB1"],\
            ["PRO"," HD3","HD1"],\
            ["PRO"," HB3","HB1"],\
            ["PRO"," HG3","HG1"],\
            ["ASN"," HB3","HB1"],\
            ["ILE","HG13","HG11"],\
            ["ILE"," CD1","CD"],\
            ["ILE","HD11","HD1"],\
            ["ILE","HD12","HD2"],\
            ["ILE","HD13","HD3"],\
             ["HIS"," HB3","HB1"],\
            ["LYS"," HB3","HB1"],\
            ["LYS"," HG3","HG1"], \
            ["LYS"," HE3","HE1"], \
            ["LYS"," HD3","HD1"]]
# Read in file - Only model 1 and ATOM informations are read
 #COLUMNS        DATA  TYPE    FIELD        DEFINITION
# -------------------------------------------------------------------------------------
#  1 -  6        Record name   "ATOM  "
#  7 - 11        Integer       serial       Atom  serial number.
#  13 - 16        Atom          name         Atom name.
#  17             Character     altLoc       Alternate location indicator.
#  18 - 20        Residue name  resName      Residue name.
#  22             Character     chainID      Chain identifier.
#  23 - 26        Integer       resSeq       Residue sequence number.
#  27             AChar         iCode        Code for insertion of residues.
#  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#  39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#  47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#  55 - 60        Real(6.2)     occupancy    Occupancy.
#  61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#  77 - 78        LString(2)    element      Element symbol, right-justified.
#  79 - 80        LString(2)    charge       Charge  on the atom.
#

def check_me(atom,res):
    
    for item in suspects:
        if(res == item[0] and atom == item[1]):
            print "# Renaming ", item[1], 'in', res, 'to', item[2]
            return item[2]
    else:
        return atom
    
def fix_n_term(atom):
    if (atom == ' HT1'):
        return ' H1 '
    if (atom == ' HT2'):
        return ' H2 '
    if (atom == ' HT3'):
        return ' H3 '
    return atom

def fix_c_term(atom):
    if(atom == ' OT1'):
        return ' O'
    if(atom ==" OT2"):
        return "OXT"
    return atom

fh = open(pdb_file,'r')
data = []
# Read it once and get first/last residue
first = 0
last = 0
count = 0
for line in fh:
    split = line.split()
    if(len(split) != 0):
        if(split[0] == "ATOM"):
            res_nr = int(line[22:26])
            if(count ==0):
                first = res_nr 
            last = res_nr
            count += 1
fh.close()
# read it again
fh = open(pdb_file,'r')
for line in fh:
    split = line.split()
    if(len(split) != 0):
        if(split[0] == "ATOM"):
            atom_nr = int(line[6:11])
            atom_name = line[12:16]
            res_name = line[17:20]
            res_nr = int(line[22:26])
            ## Stupid check - Gromacs defaults to HSE
            if(res_name == "HIS"):
                res_name = "HSE"
                if(atom_name == ' HD1'):
                    continue
            atom_name = check_me(atom_name,res_name)
            # check N/C term
            if(res_nr == first):
                atom_name = fix_n_term(atom_name)
            if(res_nr == last):
                atom_name = fix_c_term(atom_name)

            # Handle BB hydrogens
            if(atom_name==" HN "):
                atom_name = "H"
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            v = [atom_nr,atom_name,res_name, res_nr, x, y, z]
            data.append(v)
fh.close()

new_file = sys.argv[2]
fh = open(new_file,'w')
for j in range (0,len(data)):
    
    string = "ATOM  %5i %4s %3s  %4i    %8.3f%8.3f%8.3f \n" % ( data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5], data[j][6])
    fh.write(string)

fh.close()

