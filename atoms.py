import numpy as np
from math import sqrt, pi, log10, ceil

class Atom(object):
    def __init__(self, atName, eqivDict, xyzDict, gij, UijDict, sUijDict):

        self.atName = atName
        self.process_name(atName, eqivDict)
        self.xyz = np.array(xyzDict[self.baseName])
        if self.trFlag: self.transform_xyz()

        #it's helpful for the atom to have some metric pars of the structure
        self.gij = gij 
        gInv = self.gInv = np.linalg.inv(gij)
        #Nij is the diagonal matrix with values a*, b* and c*.
        self.Nij = np.array([[sqrt(gInv[0][0]), 0, 0],
                             [0, sqrt(gInv[1][1]), 0],
                             [0, 0, sqrt(gInv[2][2])]])
        NInv = self.NInv = np.linalg.inv(self.Nij)
        #cosinesRS (cosines in reciprocal space) are the square symmetric 
        #matrix with diagonal terms all 1 and off-diagonal terms cos(al*), etc.
        self.cosinesRS = np.linalg.multi_dot([NInv, gInv, NInv]) 

        #set Uij (will be None if not available)
        self.set_Uij(UijDict, sUijDict)
        if self.trFlag and (self.Uij is not None): self.transform_Uij()

    def process_name(self, atName, eqivDict):
        '''Checks whether the atom ends with _$n and sets attributes baseName,
        seitz (a 4x4 numpy array from the eqivDict), and trFlag (shows whether
        the Seitz operator should be applied to xyz and Uij).
        '''
        if '$' in atName:
            eqNum = atName[atName.find('$'):]
            self.baseName = atName[:(atName.find('$')-1)]
            self.seitz = eqivDict[eqNum]
            self.trFlag = True
        else: 
            self.baseName = atName
            self.seitz = None
            self.trFlag = False

    def transform_xyz(self):
        '''no arguments bc self.xyz and self.seitz are already available'''
        self.xyz = np.dot(self.seitz, np.append(self.xyz, 1))[0:3]

    def set_Uij(self, UijDict, sUijDict):
        if self.baseName in UijDict.keys():
            U11, U22, U33, U23, U13, U12 = UijDict[self.baseName]
            sU11, sU22, sU33, sU23, sU13, sU12 = sUijDict[self.baseName]
            self.Uij = np.array([
                [U11, U12, U13],
                [U12, U22, U23],
                [U13, U23, U33]])
            self.sUij = np.array([
                [sU11, sU12, sU13],
                [sU12, sU22, sU23],
                [sU13, sU23, sU33]])
        else:
            self.Uij = self.sUij = None



    def transform_Uij(self):
        '''The transformation of Uij is given as U*' = R U* R(t), where U* has 
        fractional elements referred to the reciprocal basis, in Grosse- 
        Kunstleve and Adams, J. Appl. Cryst. (2002), 35, 477-480. This form is
        convenient for ADP transformations, but Uij in the CIF definition is
        referred to unit vectors in reciprocal space, giving a matrix 
        representation with units of A^2. The necessary transformations are
        Uij = N-1 U* N-1(t); U* = N Uij N(t); where N (as above) is the 
        diagonal matrix with elements a*, b*, c*. The matrix product in this
        method converts to U*, transforms U* by the rotational part of the Seitz
        operator, and converts back to Uij.
        '''
        R = np.array([row[0:3] for row in self.seitz[0:3]])
        N = self.Nij
        NInv = self.NInv
        self.Uij = np.linalg.multi_dot(
                [NInv, R, N, self.Uij, N, R.transpose(), NInv])
        self.sUij = np.linalg.multi_dot(
                [NInv, R, N, self.sUij, N, R.transpose(), NInv])


    def calc_msda2(self, vector, type='uvw'):
        '''returns the msda and esd as a string with an esd given between
        2 and 19.'''

        def _format_esd(value, esd):
            exp = ceil(log10(1.949) - log10(esd))
            esd = int(round(esd * 10**exp))
            return '{0:.{1}f}({2})'.format(value, exp, esd)

        hkl = vector if type=='hkl' else np.dot(self.gij, vector)
        nUvrs = np.dot(self.Nij, hkl) #uvrs = Unit Vectors in Reciprocal Space
        msda = (np.linalg.multi_dot([nUvrs, self.Uij, nUvrs]) / 
                np.linalg.multi_dot([nUvrs, self.cosinesRS, nUvrs])) 
        sMsda = (np.linalg.multi_dot([nUvrs, self.sUij, nUvrs]) / 
                np.linalg.multi_dot([nUvrs, self.cosinesRS, nUvrs]))    
        return _format_esd(msda, sMsda)
