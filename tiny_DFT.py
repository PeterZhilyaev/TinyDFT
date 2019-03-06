# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 12:52:10 2019

@author: p.zhilyaev
"""

class TinyDFT:
    """
    """
    def __init__(self, inputParameter, crystalStructure):
        self.inputParameter = inputParameter
        self.crystalStructure = crystalStructure
        
        
    def run(self):
        print(self.inputParameter.getEnergyCutoff())
        print(self.inputParameter.getDeltaDensityTolerance())
        print(self.crystalStructure.latticeTranslationVectors)
        print(self.crystalStructure.calculateLatticeReciprocalVectors())
        