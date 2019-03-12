# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 18:57:21 2017

@author: p.zhilyaev
"""
import numpy as np
from input_parameters import InputParameters
from crystal_structure import CrystalStructure
from reciprocal_space import ReciprocalSpace

class Kpoints():
    """
    """
    def __init__(self, kPoints, weights, gVectors, energyCutoff):
        """
        """
        self.kPoints = kPoints
        self.numberOfKPoints = len(kPoints)
        self.weights = weights
        self.gVectors = gVectors
        self.energyCutoff = energyCutoff
        self.numberOfPlainWaves = self.getNumberOfPlainWaves()
        self.indexesOfGVectors = self.getIndexesOfGVectors()
        
        
    def getNumberOfPlainWaves(self):
        """
        """
        numberOfPlainWaves = np.zeros(self.numberOfKPoints, dtype='int')
        for kPointIndex, kPoint in enumerate(self.kPoints):
            plainWavesTotal = 0
            for gVector in self.gVectors:
                if np.linalg.norm(kPoint + gVector)**2 < self.energyCutoff:
                    plainWavesTotal += 1
            numberOfPlainWaves[kPointIndex] = plainWavesTotal
        return numberOfPlainWaves
     
     
    def getIndexesOfGVectors(self):
        """
        """
        indexesOfGVectors = np.zeros((self.numberOfKPoints, 
                                      max(self.numberOfPlainWaves)), dtype='int')
        for kPointIndex, kPoint in enumerate(self.kPoints):
            plainWavesTotal = 0
            for gVectorIndex, gVector in enumerate(self.gVectors):
                if np.linalg.norm(kPoint + gVector)**2 < self.energyCutoff:                    
                    indexesOfGVectors[kPointIndex, plainWavesTotal] = gVectorIndex
                    plainWavesTotal += 1
            
        return indexesOfGVectors
        
    def printAttributes(self):
        """
        """
        print('Number of plain waves for first k-point = ', self.numberOfPlainWaves[0])
        
        
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
    
    kPoints = [2 * np.pi / scaleFactor * np.array([0.6223, 0.2953, 0.0], dtype='float')]
    weights = [1.0]
    
    kpoints = Kpoints(kPoints, weights, reciprocalSpace.gVectors, inputParameters.energyCutoff)   
    kpoints.printAttributes()
    