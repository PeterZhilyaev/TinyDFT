# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 12:50:37 2019

@author: p.zhilyaev
"""
import numpy as np
from input_parameters import InputParameters
from crystal_structure import CrystalStructure
from tiny_DFT import TinyDFT


if __name__ == "__main__":
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
    
    tinyDFT = TinyDFT(inputParameters, crystalStructure)
    tinyDFT.run()