'''
A collection of functions to parse information from Gaussian09.

Programmer: Lo-Yun Lee (josephleee@gmail.com), 2016
'''

import numpy as np

def File2List(filename):
    f = open(filename)
    f_lst = f.readlines()
    f.close()
    return f_lst

#### Boolean ####
def contain_hessian(g09file):
    f_lst = File2List(g09file)
    hess_exist = False
    for line in f_lst:
        if '#p' in line.lower() or 'iop(7/33=1)' in line.lower():
            hess_exist = True
            break
    return hess_exist

def using_hpmode(g09file):
    f_lst = File2List(g09file)
    hpmode = False
    keyw = ['#', 'freq', 'hpmode']
    for i, line in enumerate(f_lst):
        if all( [k in line.lower() for k in keyw] ):
            hpmode = True
            break
    return hpmode
#################

def get_hessian_iop(filename_log, num_coords):
    if contain_hessian(filename_log):
        f_lst = File2List(filename_log)
        hessian = np.zeros([num_coords, num_coords])
        for i, line in enumerate(f_lst):
            if 'Z-matrix is all fixed cartesians, so copy forces.' in line: i_start = i+2
            if 'Final forces over variables, Energy' in line: i_end = i
        #end loop
        for line in f_lst[i_start:i_end]:
            if '.' in line and 'D' in line:   ## At a line of data
                word = line.split()
                row = int(word[0])-1
                data = [ float(tmp.replace('D', 'e')) for tmp in word[1:] ]
                for ii, elem in enumerate(data):
                    i, j = ( row, col_line[ii] )
                    hessian[i, j] = elem
                    hessian[j, i] = elem
            else:  ## At a line of indices
                col_line = [int(tmp)-1 for tmp in line.split()]
        #end loop
        return hessian
    else:
        print "No Hessian in file \'%s\'. Please check!!"%(filename_log)
        return None

def hessian(filename, NAtom):
    '''
    Prerequisite: G09-job output with '#p' as its keyword header
    
    Arguments
    ---------
    filename: str
        File name of Gaussian09 output.
    NAtom: int
        Number of atoms.

    Returns
    -------
    hess: np.ndarray
        2-D array using [hartree] unit system.
        This Hessian is not mass-weighted, and it use Cartesian coordinates as basis.
        Format of basis: [X1, Y1, Z1, X2, Y2, Z2, ... ]
    '''
    if contain_hessian(filename):
        f_lst = File2List(filename)
        dim = 3*NAtom
        hess = np.zeros([dim, dim])
        ## Set pointers:
        i_go, i_stop, NBunch, NLine = (0,0,0,0)
        TxtFlag1 = 'Z-matrix is all fixed cartesians, so copy forces.'
        TxtFlag2 = 'Force constants in Cartesian coordinates:'
        for i, line in enumerate(f_lst):
            if TxtFlag2 in line:# or TxtFlag2 in f_lst[i+1]:
                i_go = i+1#2
                #break
        NBunch = ( dim/5 + int(dim%5 != 0) )
        NLine = NBunch * ( 2*dim - 5*(NBunch - 1) ) / 2   # arithmetic series
        i_stop = i_go + NBunch + NLine
        ## Read Hessian:
        for line in f_lst[i_go: i_stop]:
            word = line.replace('\n','').split()
            if word[0].isdigit() and word[-1].isdigit():
                col_line = [ int(tmp)-1 for tmp in word ]
            else:
                row = int(word[0])-1
                data = [ float(tmp.replace('D', 'e')) for tmp in word[1:] ]
                hess[ row, col_line[:len(data)] ] = np.asarray(data)
                hess[ col_line[:len(data)], row ] = np.asarray(data)
        return hess
    else:
        print "No Hessian in file \'%s\'. Please check!!"%(filename)
        return None

####################################################
####################################################
####################################################
def get_num_atoms(filename_log):
    # Get number of atoms in the system
    # Returns number of atoms as an integer
    f_lst = File2List(filename_log)
    for line in f:
        if 'NAtoms=' in line:
            strs = line.lstrip().rstrip().split()  #lee# Simplify: strs=line.split()
            break
    return int(strs[1])

def get_atomic_masses(filename_log, num_atoms):
    f_lst = File2List(filename_log)
    TxtFlag = ['Atom', 'has atomic number', 'and mass']
    ms = np.zeros(num_atoms)  #lee# creating an array, set the dimension of it to be: num_atoms
    for line in f_lst:
        if all([ tf in line for tf in TxtFlag]):
            # set masses
            strs = line.lstrip().rstrip().split()
            #lee# get the labeling number of an atom. (from G09)
            ms[ int(strs[1]) - 1 ] = float(strs[-1])
            #lee# the value of atomic mass is at the end of a line.
        elif 'Principal axes and moments of inertia in atomic units:' in line:
            break
    return ms

def get_atom_info(filename_log, num_atoms, n=1):
    # Return atomic coordinates and
    # atomic numbers of nth standard configuration
    f = open(filename_log)
    coords_atoms = []
    Zs_atoms = []
    count = 0
    strings = ('Input orientation:', 'Z-Matrix orientation:')
    for line in f:
    #lee# Locating the line contains 'Input orientation:' 
    #lee# or 'Z-Matrix orientation:' in the file.
        if any(s in line for s in strings):  #lee# Simplify: if line in strings:
            count += 1
            if count == n:
                break
    #-------------#
    for _ in range(4):    #lee# while found the line, slip four lines bellow.
        line = next(f)
    for i in range(num_atoms):    #lee# Starting reading data from the file.
        nextline = next(f)        #lee# Each line describes an atom.
        strs = nextline.lstrip().rstrip().split()
        Z = float(strs[1])    #lee# atomic number
        x, y, z = ( float(strs[i]) for i in [3,4,5] )  #lee# coordinates
        coords_atoms.append(np.array([x, y, z]))    
        #lee# storing coordinates in array-of-array 'coords_atoms'
        #lee# coords_atoms ==> [array([x1,y1,z1]), array([x2,y2,z2]), ...]
        Zs_atoms.append(Z)    #lee# storing atomic number in array 'Zs_atoms'
    f.close()
    return coords_atoms, Zs_atoms
##================================================##
##def atom_info(filename):
##    f_lst = File2List(filename_log)
##    atom = {}
##    ## Set Text Flags
##    TxtFlag_mass = ['Atom', 'has atomic number', 'and mass']
##
##    atom['mass'] = np.zeros(num_atoms)
##    for line in f_lst:
##        if 'NAtoms=' in line:
##            atom['NAtom'] = line.split()[1]
##        elif all([ tf in line for tf in TxtFlag_mass]):
##            strs = line.split()
##            atom['mass'][ int(strs[1]) - 1 ] = float(strs[-1])
##        elif
        
####################################################
####################################################
####################################################

def NorModes_hpmode(g09file, NAtom, NVib):
    '''
    Arguments
    ---------
    g09file: str
        File name of Gaussian09 output.
    NAtom: int
        Number of atoms.
    NVib:  int
        Number of vibrational modes.

    Returns
    -------
    A two-dimensional np.ndarray with all elements as floats.
    To access the i-th normal mode: NMode[i,:] or NMode[i] for short.

    Note
    ----
    G09 keyword 'freq=hpmode' will give you high precision format (to five figures).
    For more information, please check: http://www.gaussian.com/g_tech/g_ur/k_freq.htm
    '''
    print 'using hpmode'
    f_lst = File2List(g09file)
    text_sign = 'Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering'
    bool_data = False
    MLstPtr = []  ## pointer for reading mode lists
    VctPtr = []    ## pointer for reading normal-mode vectors
    NMode = np.zeros([NVib, 3*NAtom])
    ## Set pointers:
    for i, line in enumerate(f_lst):
        if not bool_data and text_sign in line:
            bool_data = True
        elif  bool_data  and text_sign in line:
            bool_data = False
            break
        #endif
        word = line.replace('\n','').split()
        if bool_data and word[0].isdigit() and word[-1].isdigit():  MLstPtr.append(i)
        if bool_data and 'Coord Atom Element:' in line:  VctPtr.append(i+1)
    #print MLstPtr, '\n', VctPtr
    ## Input normal modes:
    for ibunch, i in enumerate(VctPtr):
        MList = [ int(ii)-1 for ii in f_lst[MLstPtr[ibunch]].split() ]
        for j, line in enumerate( f_lst[i: i+3*NAtom] ):
            tmp= [ float(ii) for ii in line.split()[3:] ]
            NMode[ MList[0]: MList[-1]+1, j ] = np.asarray(tmp)
    return NMode

def NorModes(g09file, NAtom, NVib):
    '''
    Arguments
    ---------
    g09file: str
        File name of Gaussian09 output.
    NAtom: int
        Number of atoms.
    NVib:  int
        Number of vibrational modes.

    Returns
    -------
    A two-dimensional np.ndarray with all elements as floats.
    To access the i-th normal mode: NMode[i,:] or NMode[i] for short.
    '''
    if using_hpmode(g09file):
        return NorModes_hpmode(g09file, NAtom, NVib)
    else:
        f_lst = File2List(g09file)
        go_sign = 'Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering'
        stop_sign = '- Thermochemistry -'
        bool_data = False
        MLstPtr = []  ## pointer for reading mode lists
        VctPtr = []   ## pointer for reading normal-mode vectors
        NMode = np.zeros([NVib, 3*NAtom])
        ## Set pointers:
        for i, line in enumerate(f_lst):
            if not bool_data and go_sign in line:
                bool_data = True
            elif bool_data and ( stop_sign in line or line.isspace() ):  # change to line.isspace() ??
                bool_data = False
                break
            #endif
            word = line.replace('\n','').split()
            if bool_data and word[0].isdigit() and word[-1].isdigit():  MLstPtr.append(i)
            if bool_data and 'AtomANXYZ' in line.replace(' ',''):  VctPtr.append(i+1)
        #print MLstPtr, '\n', VctPtr
        ## Input normal modes:
        for ibunch, i in enumerate(VctPtr):
            MList = [ int(ii)-1 for ii in f_lst[MLstPtr[ibunch]].split() ]
            for j, line in enumerate( f_lst[i: i+NAtom] ):
                data = [ float(ii) for ii in line.split()[2:] ]
                tmpx = [ xx for i, xx in enumerate(data) if i % 3 == 0]
                tmpy = [ xx for i, xx in enumerate(data) if i % 3 == 1]
                tmpz = [ xx for i, xx in enumerate(data) if i % 3 == 2]
                NMode[ MList[0]: MList[-1]+1, j*3+0 ] = np.asarray(tmpx)
                NMode[ MList[0]: MList[-1]+1, j*3+1 ] = np.asarray(tmpy)
                NMode[ MList[0]: MList[-1]+1, j*3+2 ] = np.asarray(tmpz)
        return NMode

