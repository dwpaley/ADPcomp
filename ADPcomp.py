#!/usr/bin/env python

import sys
import numpy as np
import sympy as sp
import read_cif, help_files
import atoms as atom_lib



def make_proc_atom(atomDict, eqivDict, xyzDict, gij, UijDict, sUijDict):
    '''This is a closure to lock in all the dict names so they don't have
    to be passed as arguments every time you process an atom.'''

    def proc_atom(atName):
        '''This checks the atomDict for the atom, creates it if necessary,
        and returns it.'''
        if atName not in atomDict.keys():
            atomDict[atName] = atom_lib.Atom(atName, eqivDict, xyzDict, gij,
                    UijDict, sUijDict)
        return atomDict[atName]


    return proc_atom



def process_eqiv(body, eqivDict):
    '''Takes the body of an EQIV instruction. Constructs a Seitz operator as
    a 4x4 matrix (to transform [x, y, z, 1]). Puts the Seitz matrix in a 
    dictionary.
    Notes: upper or lower case ok. fractions and decimals ok. spaces ignored.
    Uses sympy imported as sp.
    '''
    trList = ''.join(body[1:]).split(',')

    x, y, z = sp.Symbol('x'), sp.Symbol('y'), sp.Symbol('z')
    eqns = [sp.sympify(eq.lower()) for eq in trList]
    seitz = np.array([
            [sp.diff(eqn, var) for var in [x, y, z]]
            for eqn in eqns])
    seitz = np.append(seitz,
            [[eqn.evalf(subs={x:0,y:0,z:0})] for eqn in eqns],
            axis=1)
    seitz = np.append(seitz, [[0, 0, 0, 1]], axis=0)

    eqivDict[body[0]] = seitz.astype(float)


def get_file_name(aip):
    ''' finds a line beginning with 'file' and returns the first following word
    '''
    with open(aip, 'r') as inFile:
        for line in inFile.readlines():
            if line.lower().startswith('file'): return line.split()[1]
        return None




def main(intRun=False):

    if intRun:
        print('ADPcomp: Interactive mode.')
        cifName = input('Enter the structure file name (CIF format): ') 
        print('\nEnter instructions: eqiv, bond, axis, ip, oop, uvw, hkl')
        print('Type q to exit at any time.\n')
    else:
        aipName = sys.argv[1]
        outFile = open(aipName+'.aop', 'w')
        aipName += '.aip'
        cifName = sys.argv[2] if len(sys.argv)>2 else get_file_name(aipName)
        with open(aipName, 'r') as aipFile:
            aipFileLines = aipFile.readlines()
        maxLen = len(max(aipFileLines, key=len))
        aipIter = iter(aipFileLines)

    #make and populate the run-level data objects and the proc_atom function
    gij = np.array([[0,0,0],[0,0,0],[0,0,0]])
    xyzDict, UijDict, sUijDict, atomDict, eqivDict = [dict() for i in range(5)] 
    read_cif.read_cif(cifName, gij, xyzDict, UijDict, sUijDict)
    proc_atom = make_proc_atom(atomDict, eqivDict, xyzDict, gij, 
            UijDict, sUijDict)


    while 1:
        #Get a line. Check for exit status. Ignore short lines.
        try: line = input('>> ') if intRun else next(aipIter)
        except StopIteration: break
        calcFlag=False #Are we going to calculate an MSDA in this loop?  
        if intRun and line=='q': break
        if len(line.split()) < 2: 
            if not intRun: outFile.write(line)
            continue


        #Process atoms or eqiv instructions in the line
        head = line.lower().split()[0] 
        body = line.lower().split()[1:] 
        if head == 'eqiv': 
            process_eqiv(body, eqivDict)
        elif head in ['bond', 'axis', 'ip', 'oop']:
            calcFlag = True
            atoms = [proc_atom(atName) for atName in body]
        elif head in ['hkl', 'uvw']:
            calcFlag = True
            atoms = [proc_atom(body[0])]
            vector = np.array([float(x) for x in body[1:4]])
        else: pass


        #Construct the vector to measure along
        if head == 'bond':
            uvw = atoms[1].xyz - atoms[0].xyz
        elif head == 'axis':
            uvw = atoms[2].xyz - atoms[1].xyz
        elif head == 'ip':
            t1 = atoms[1].xyz - atoms[0].xyz
            t2 = atoms[2].xyz - atoms[0].xyz
            uvw = (t1/np.linalg.norm(t1)) + (t2/np.linalg.norm(t2))
        elif head == 'oop':
            uvw = np.cross( np.dot(gij, (atoms[1].xyz-atoms[0].xyz)),
                            np.dot(gij, (atoms[2].xyz-atoms[0].xyz)) )
        elif head == 'hkl':
            uvw = np.dot(np.linalg.inv(gij), vector)
        elif head == 'uvw':
            uvw = vector
        else: pass


        #Calculate and output results 
        if not intRun: print(line.rstrip('\n'))
        if calcFlag:
            uvw /= max(uvw, key=abs)
            print('MSDA for {0} along uvw=[{2:.3f} {3:.3f} {4:.3f}] = {1} A^2\n'
                    .format(atoms[0].atName.capitalize(),
                    atoms[0].calc_msda2(uvw), *uvw))
            if not intRun:
                outFile.write(line.rstrip('\n') + (maxLen + 5 -len(line))*' ')
                outFile.write(atoms[0].calc_msda2(uvw) + '\n')
        elif not intRun:
            outFile.write(line + '\n')
        else: 
            pass


    if not intRun: 
        print('Output written to {}\n'.format(outFile.name))
        outFile.close()



if len(sys.argv) == 1: help_files.intro()
elif '-h' in sys.argv: help_files.help()
else: 
    help_files.intro_short()
    main('-i' in sys.argv)



#Some useful references: Acta Cryst 17, 142 (This contains the right formula for
#<u^2>); Acta Cryst A32, 239; Acta A71, 59; Acta C44, 775; Acta B46, 683
#
#Copyright 2016, Daniel W. Paley
#Contact: dwp2111@columbia.edu
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


