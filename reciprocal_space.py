# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 14:50:49 2019

@author: p.zhilyaev
"""
import numpy as np
from input_parameters import InputParameters
from crystal_structure import CrystalStructure


class ReciprocalSpace():
    """
    """    
    def __init__(self, energyCutoff, latticeTranslationVectors,
                 latticeReciprocalVectors):
        """
        """
        self.energyCutoff = energyCutoff
        self.latticeTranslationVectors = latticeTranslationVectors
        self.latticeReciprocalVectors = latticeReciprocalVectors      
        self.millerIndexes = self.calculateMillerIndexes()
        self.gVectors = self.calculateGVectors()
        self.gVectorsLengths = self.calculateGVectorsLengths()
        self.numberOfGVectors = self.calculateNumberOfGVectors()
        # Size of Reciprocal grid
        self.ngx, self.ngy, self.ngz = self.calculateSizeOfReciprocalGrid()
        # Name of variables is way to long
        self.mapMillerIndexesToGVectorIndex = self.calculateMapMillerIndexesToGVectorIndex()
                                                          

    def calculateMillerIndexes(self):
        """
        Calculate number of G-vectors such that (\hbar^2 / 2m) * Gmax^2 < 4 * EnergyCutoff
        |Gmax| - maximum lentgh of G-vector that is 2 * sqrt(EnergyCutoff)
        Keep in mind that in this formula G-vector has a dimension Angstrom^-1
        ijkMax[3] = [imax, jmax, kmax] are estimate of max index i, j, k used in
        the generation of G-vectors: G(i, j, k) =  i * b1 + j * b2 + k * b3, where
        b1, b2, b3 - lattice reciprocal vectors 
        Since i = (G, a1) / (b1, a1) then abs(imax) <= |Gmax| * |a1| / (a1, b1), where
        a1, a2, a3 - lattice translation vectors
        """   
        millerIndexes = []
        ijkMax = np.zeros(3, dtype='int')     
        for index in range(3):
            ijkMax[index] = np.rint(2 * np.sqrt(self.energyCutoff) * np.linalg.norm(self.latticeTranslationVectors[index, :]) /
                                    np.dot(self.latticeTranslationVectors[index, :], self.latticeReciprocalVectors[index, :])
                                    + 0.5)
        for i in range(-ijkMax[0], ijkMax[0] + 1):
            for j in range(-ijkMax[1], ijkMax[1] + 1):
                for k in range(-ijkMax[2], ijkMax[2] + 1):
                    gVector = ( i * self.latticeReciprocalVectors[0, :] + 
                                j * self.latticeReciprocalVectors[1, :] +
                                k * self.latticeReciprocalVectors[2, :] )
                                      
                    if (np.linalg.norm(gVector) <= 2 * np.sqrt(self.energyCutoff)): 
                        millerIndexes.append((i, j, k))
        return np.array(millerIndexes, dtype='int')


    def calculateGVectors(self):
        """
        """
        gVectors = []
        for ijk in self.millerIndexes:
            gvector = ( ijk[0] * self.latticeReciprocalVectors[0, :] + 
                        ijk[1] * self.latticeReciprocalVectors[1, :] +
                        ijk[2] * self.latticeReciprocalVectors[2, :] )
            gVectors.append(gvector)
        return np.array(gVectors, dtype='float')

    def calculateGVectorsLengths(self):
        gVectorsLengths = []
        for gVector in self.gVectors:
            gVectorsLengths.append(np.linalg.norm(gVector))
        return np.array(gVectorsLengths, dtype='float')


    def calculateNumberOfGVectors(self):
        """
        """
        return self.gVectors.shape[0]


    def calculateSizeOfReciprocalGrid(self):
        """
        """
        ngx =  2 * max(self.millerIndexes[:, 0]) + 1
        ngy =  2 * max(self.millerIndexes[:, 1]) + 1
        ngz =  2 * max(self.millerIndexes[:, 2]) + 1
        return ngx, ngy, ngz


    def calculateMapMillerIndexesToGVectorIndex(self):
        """
        """
        mapMillerIndexesToGVectorIndex = np.zeros( (self.ngx, self.ngy, self.ngz), dtype='int')
        for gVectorIndex in range(self.numberOfGVectors):
            i, j, k = self.millerIndexes[gVectorIndex, :] 
            mapMillerIndexesToGVectorIndex[i, j, k] = gVectorIndex            
        return mapMillerIndexesToGVectorIndex


    def printAttributes(self):
        """
        """
        print('Number Of G-vectors: ', self.numberOfGVectors)
        print('Reciprocal-space grid ngx, ngy, ngz: ', self.ngx, self.ngy, self.ngz)

if __name__ == '__main__':
    inputParameters = InputParameters(20.0, 8, 0.5, 1E-4, 100, 4)
    #### Defining crystalStructure
    systemName = 'Si diamond'
    scaleFactor = 10.26
    latticeTranslationVectors = np.zeros((3, 3))
    latticeTranslationVectors[0, :] = np.array([0.5, 0.5, 0.0])
    latticeTranslationVectors[1, :] = np.array([0.5, 0.0, 0.5])
    latticeTranslationVectors[2, :] = np.array([0.0, 0.5, 0.5])
    numberOfAtoms = 2
    basisVectors = np.zeros((numberOfAtoms, 3))
    basisVectors[0, :] = np.array([0.125, 0.125, 0.125])
    basisVectors[1, :] = np.array([-0.125, -0.125, -0.125])    
    crystalStructure = CrystalStructure(systemName, numberOfAtoms, scaleFactor,
                                        latticeTranslationVectors, basisVectors)
    reciprocalSpace = ReciprocalSpace(inputParameters.getEnergyCutoff(), 
        crystalStructure.latticeTranslationVectors,
        crystalStructure.latticeReciprocalVectors)
    print reciprocalSpace.millerIndexes[22]
    print reciprocalSpace.gVectors[22]
    print reciprocalSpace.gVectorsLengths[22]
    print reciprocalSpace.ngx, reciprocalSpace.ngy, reciprocalSpace.ngz
    print reciprocalSpace.mapMillerIndexesToGVectorIndex[1, 2, 3]
    reciprocalSpace.printAttributes()