#===================================Start import and load packages==================================================
import time
import matplotlib
from neuron import h ,crxd as rxd, gui
import sys, math, itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import  cm
import scipy.interpolate as interp
rxd.options.enable.extracellular = True
h.load_file("import3d.hoc")
h.load_file("stdgui.hoc")
h.load_file("stdrun.hoc")
#===================================End import and load packages==================================================
h.celsius = 34

#===================================Start cell class==================================================
class HL23PV:

    '''
    Stated in Migliore 2005, the basic set of active dendritic properties included 
    sodium, DR- and A-type potassium conductances (Na, K_DR, KA, respectively), and I_h current
    '''

    '''
    NeuroMorph definition: Axon = axon, Soma = soma, Apical dendrite = apic, Basilar dendrite = dend
    '''

    '''
    Standard method of setting a specific parameter of a mechanim:
    soma(0.5).hh.gkbar
    '''

    def __init__(self, morphology):
        self.morphology = morphology
        self.load_morphology()
        self.discretize()
        # self.general_mech()
        # self.soma_mech()
        # self.axon_mech()   

    def __str__(self):
        return self.morphology

    def load_morphology(self):
        cell = h.Import3d_SWC_read()
        cell.input(f"{self.morphology}.swc")
        i3d = h.Import3d_GUI(cell, False)
        i3d.instantiate(self)
    
    def discretize(self):
        """
        dlamb rule : see NEURON :  a tool for neuroscientist, 2001
        """
        freq, d_lam = 100, 0.1
        for sec in self.all:
            sec.nseg = math.ceil( (sec.L/(d_lam * h.lambda_f(freq) ))/2.0 )*4 +1 

#===================================End cell class==================================================

#===================================Start cell instantiation==================================================
c = HL23PV("HL23PV")

ps = h.PlotShape(False)
ps.variable("v")
ps.plot(plt)
plt.show()

input('') # Ensure the closure of the program running, press enter to terminate the run. 


#===================================End cell instantiation==================================================