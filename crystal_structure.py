# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:18:20 2019

@author: p.zhilyaev
"""
import numpy as np


class CrystalStructure:
    """
    """
    def __init__(self, systemName, numberOfAtoms, scaleFactor,
                 latticeTranslationVectors, basisVectors):
        self.systemName = systemName
        self.numberOfAtoms = numberOfAtoms
        self.scaleFactor = scaleFactor
        self.latticeTranslationVectors = scaleFactor * latticeTranslationVectors
        self.basisVectors = scaleFactor * basisVectors
        self.unitCellVolume = self.calculateUnitCellVolume()
        self.latticeReciprocalVectors = self.calculateLatticeReciprocalVectors()

        
    def calculateUnitCellVolume(self):
        """
        """
        return np.abs(np.linalg.det(self.latticeTranslationVectors))
    
    
    def calculateLatticeReciprocalVectors(self):
        """
        """
        latticeReciprocalVectors = np.zeros((3, 3))
        latticeReciprocalVectors[0, :] = np.cross(self.latticeTranslationVectors[1, :],
                                                  self.latticeTranslationVectors[2, :])
        latticeReciprocalVectors[1, :] = np.cross(self.latticeTranslationVectors[2, :],
                                                  self.latticeTranslationVectors[0, :])                                                
        latticeReciprocalVectors[2, :] = np.cross(self.latticeTranslationVectors[0, :],
                                                  self.latticeTranslationVectors[1, :])
        # Note here we can get volume < 0 and we use abs in calculateUnitCellVolume
        # because we want volume to be > 0                                                  
        volume = np.linalg.det(self.latticeTranslationVectors)
        return 2 * np.pi / volume * latticeReciprocalVectors 

    