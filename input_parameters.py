# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 12:56:23 2019

@author: p.zhilyaev
"""

class InputParameters:
    """
    """
    def __init__(self, energyCutoff, numberOfElectrons, mixingParameter,
                 deltaDensityTolerance, maximumNumberOfIterations,
                 numberOfBands):
        """
        """
        # energy cutoff in Rynbergs
        self.energyCutoff = energyCutoff
        self.numberOfElectrons = numberOfElectrons
        self.mixingParameter = mixingParameter
        self.deltaDensityTolerance = deltaDensityTolerance
        self.maximumNumberOfIterations = maximumNumberOfIterations
        self.numberOfBands = numberOfBands
        
        
    def getEnergyCutoff(self):
        """
        """
        return self.energyCutoff
        
    
    def getNumberOfElectrons(self):
        """
        """
        return self.numberOfElectrons
    
    
    def getMixingParamter(self):
        """
        """
        return self.mixingParameter
    
    
    def getDeltaDensityTolerance(self):
        """
        """
        return self.deltaDensityTolerance
    
    
    def getMaximumNuberOfInterations(self):
        """
        """
        return self.maximumNumberOfIterations
    
        
    def getNumberOfBands(self):
        """
        """
        return self.numberOfBands