
def intro():
    introString = '''
********************************************************************************

ADPcomp: Finds atomic mean-square displacement amplitudes along arbitrary 
    directions.
    
Daniel W. Paley, 2016.

Usage:

adpcomp -i: start in interactive mode.
adpcomp aipname [CIFname]: run in file mode with instructions from aipname.aip.
    Optionally, read structure from CIFname; otherwise structure is given on 
    FILE line in aip file. Output is written to aipname.aop.
adpcomp -h: print instructions.

********************************************************************************
    '''
    print(introString)


def intro_short():
    introShortString = '''
********************************************************************************

ADPcomp: Finds atomic mean-square displacement amplitudes along arbitrary 
    directions.
    
Daniel W. Paley, 2016.

********************************************************************************
'''
    print(introShortString)

def help():
    helpString = '''
********************************************************************************


ADPcomp calculates the mean-square deviation amplitude (MSDA) of atomic 
displacement parameters along arbitrary direct- or reciprocal-space vectors. 
The program is controlled interactively or by an input file but works the same
in either mode. 

The interactive mode is called by ADPcomp -i. The user is prompted for the cif
file name and then enters all instructions line-by-line at a prompt.

The file mode is called by ADPcomp name [CIFname]. The current directory is 
checked for name.aip. If CIFname is not given, the aip file is searched for a
line 'file CIFname'. Coordinates and ADPs are read from the cif file, then all 
other lines in the aip file are processed. A file name.aop is written with a
measured MSDA (with esd) appended to each line. 

An example aip file is shown here:

----------------------------------------
file perovtest.cif

eqiv $1 3/2-x, 1/2-y, -1/2+z

ip br02 pb01 pb01_$1
oop br02 pb01 pb01_$1

bond br02 pb01
axis br02 pb01 pb01_$1

uvw br02 0 0 1
hkl br02 2 0 1
----------------------------------------


INSTRUCTIONS:

FILE: ADPcomp will read the structure file given in the FILE instruction if a 
structure file was not given on the command line. If a structure file was given,
FILE is ignored.

EQIV $n x', y', z': identical to ShelXL-format EQIV instructions. This creates a
$n flag that may be used on atoms. 

IP atom1 atom2 atom3: Calculates the MSDA of atom1 for a vector lying in the 
atom1-atom2-atom3 plane and between the directions atom1-atom2 and atom1-atom3.

OOP atom1 atom2 atom3: MSDA of atom1 perpendicular to the atom1-
atom2-atom3 plane.

BOND atom1 atom2: MSDA of atom1 in the direction toward atom2.

AXIS atom1 atom2 atom3: MSDA of atom1 along the direction atom2-atom3. 

UVW atom u v w; HKL atom h k l: MSDA along a direct or reciprocal vector.


The input file will be echoed to <name>.aop with calculated MSDA appended to
each line.



INFORMATION:

ADPcomp calculates the mean-square displacement amplitude <u^2> of an atom along
an arbitrary vector. Note (as always for ADPs) that static and dynamic 
displacements cannot be distinguished and that displacement parameters are often 
a catch-all for other deficiencies in the model. The results are probably most
appropriate for comparison of closely related structures, e.g. multi-temperature
data sets.

The formula for <u^2> is in Busing and Levy, Acta Cryst. 1964 (17), 142-146:

<u^2> = [h(T) Beta(ij) h] / [2 pi^2 h(T) g^-1 h]

where h is an arbitrary hkl vector, (T) indicates the transpose, Beta is the 
unitless ADP tensor and g^-1 is the reciprocal metric tensor.

Note this formula was subsequently given incorrectly in Trueblood et al., Acta
Cryst. 1996 (A52), 770-781:

<u^2> = n(T) U n

where n is a unit vector referred to the unit vectors parallel to the reciprocal
axes. This formula is correct only when the axes are orthogonal. The correction
is 

1 / (n(T) [cosinesRS] n)

where cosinesRS (cosines in reciprocal space) is the symmetric matrix:

    [ 1 cos(ga*) cos(be*) ]
    [ cos(ga*) 1 cos(al*) ]
    [ cos(be*) cos(al*) 1 ]

This is the implementation in ADPcomp since it does not require the Uij matrix
to be transformed.

Estimated standard deviation for the MSDA is calculated approximately from the
esd's of the Uij components. Esd's of cell parameters and covariance between Uij
components are neglected.

Note that the cif parsing is rather inflexible. The atom positions loop must
contain at least the following lines:

loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z

The ADPs must be given separately in Uij form in the following order:

loop_
 _atom_site_aniso_label
 _atom_site_aniso_U_11
 _atom_site_aniso_U_22
 _atom_site_aniso_U_33
 _atom_site_aniso_U_23
 _atom_site_aniso_U_13
 _atom_site_aniso_U_12

This is the standard output of ShelXL but may cause problems for older data or
macromolecular structures (which probably should not be analyzed with this
program anyway!)


********************************************************************************
'''
    print(helpString)
            
