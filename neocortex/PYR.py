import time
tic_0 = time.perf_counter() #script runtime calculation value
import os
from os.path import join
import sys
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy import stats as st
import neuron
from neuron import h, gui
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from neuron.units import mV, ms

# Simulation hyperparameters
dt = 0.025 #for both cell and network
tstart = 0.
tstop = 4500.
celsius = 34.
v_init = -80. #for both cell and network




class PYR:
    def __init__(self):

        self._setup_morphology()
        self._setup_biophysics()
        
    
    def _setup_morphology(self):
        h.load_file('stdrun.hoc')
        h.nrn_load_dll('nrnmech.dll')
        h.load_file("import3d.hoc")

        h.load_file('NeuronTemplate.hoc')

        self.cell = h.NeuronTemplate("HL23PYR.swc")
    
    def _setup_biophysics(self):
        h.load_file('biophys_HL23PYR.hoc')
        h.biophys_HL23PYR(self.cell)
        
        


# pyr = PYR()
        
# h.load_file('stdrun.hoc')
# h.nrn_load_dll('nrnmech.dll')
# h.load_file("import3d.hoc")

# h.load_file('PYRtemplate.hoc')
# h.load_file('biophys_HL23PYR.hoc')

# neu = h.PYRtemplate()

# h.biophys_HL23PYR(neu)

# for sec in neu.all:
            
#     # Insertion of mechanism
#     h.pas.insert(sec)
#     h.Ih.insert(sec)

            
#     # Biophysical parameters
#     sec.Ra = 100
#     sec.cm = 1
#     sec.e_pas = -80
#     sec.g_pas = 0.0000954

# for sec in neu.somatic:
            
#     # Insertion of mechanism
#     h.NaTg.insert(sec)
#     h.K_P.insert(sec)
#     h.K_T.insert(sec)
#     h.Kv3_1.insert(sec)
#     h.Im.insert(sec)
#     h.SK.insert(sec)
#     h.Ca_HVA.insert(sec)
#     h.Ca_LVA.insert(sec)
#     h.CaDynamics.insert(sec)
            

#     sec.ek = -85
#     sec.ena = 50
#     sec.gamma_CaDynamics = 0.0005
#     sec.vshiftm_NaTg = 13
#     sec.vshifth_NaTg = 15
#     sec.slopem_NaTg = 7
#     sec.gbar_Ih = 0.000148

# for sec in neu.apical:

#     # Biophysical parameters
#     sec.cm = 2

# for sec in neu.basal:
            
#     # Biophysical parameters 
#     sec.cm = 2
            
#     sec.gbar_Ih = 0.000000709

# # h.distribute_channels(neu.apical, "apic","gbar_Ih",2,-0.8696,3.6161,0.0,2.0870,neu.soma.gbar_Ih)
# neu.distribute_channels("apic","gbar_Ih",2,-0.8696,3.6161,0.0,2.0870,0.000148)
    
# for sec in neu.axonal:
            
#     # Insertion of mechanism
#     h.SK.insert(sec)
#     h.Ca_HVA.insert(sec)
#     h.Ca_LVA.insert(sec)
#     h.K_T.insert(sec)
#     h.K_P.insert(sec)
#     h.Nap.insert(sec)
#     h.Kv3_1.insert(sec)
#     h.NaTg.insert(sec)
#     h.CaDynamics.insert(sec)
#     h.Im.insert(sec)
            

#     sec.ek = -85
#     sec.ena = 50
#     sec.vshiftm_NaTg = 10 
#     sec.slopem_NaTg = 9
#     sec.gamma_CaDynamics = 0.0005

# neu.distribute_channels("axon","decay_CaDynamics",0,1.000000,0.000000,0.000000,0.000000,226.0000000000)
# neu.distribute_channels("axon","gbar_SK",0,1.000000,0.000000,0.000000,0.000000,0.0145000000)
# neu.distribute_channels("axon","gbar_Ca_LVA",0,1.000000,0.000000,0.000000,0.000000,0.0439000000)
# neu.distribute_channels("axon","gbar_Ca_HVA",0,1.000000,0.000000,0.000000,0.000000,0.0003060000)
# neu.distribute_channels("axon","gbar_Kv3_1",0,1.000000,0.000000,0.000000,0.000000,0.9410000000)
# neu.distribute_channels("axon","gbar_K_T",0,1.000000,0.000000,0.000000,0.000000,0.0424000000)
# neu.distribute_channels("axon","gbar_K_P",0,1.000000,0.000000,0.000000,0.000000,0.3380000000)
# neu.distribute_channels("axon","gbar_Nap",0,1.000000,0.000000,0.000000,0.000000,0.0084200000)
# neu.distribute_channels("axon","gbar_NaTg",0,1.000000,0.000000,0.000000,0.000000,1.3800000000)
# neu.distribute_channels("axon","gbar_Im",0,1.000000,0.000000,0.000000,0.000000,0.0000000000)
# neu.distribute_channels("soma","decay_CaDynamics",0,1.000000,0.000000,0.000000,0.000000,20.0000000000)
# neu.distribute_channels("soma","gbar_Im",0,1.000000,0.000000,0.000000,0.000000,0.0003060000)
# neu.distribute_channels("soma","gbar_Ca_LVA",0,1.000000,0.000000,0.000000,0.000000,0.0029600000)
# neu.distribute_channels("soma","gbar_Ca_HVA",0,1.000000,0.000000,0.000000,0.000000,0.0015500000)
# neu.distribute_channels("soma","gbar_Kv3_1",0,1.000000,0.000000,0.000000,0.000000,0.0424000000)
# neu.distribute_channels("soma","gbar_SK",0,1.000000,0.000000,0.000000,0.000000,0.0008530000)
# neu.distribute_channels("soma","gbar_K_T",0,1.000000,0.000000,0.000000,0.000000,0.0605000000)
# neu.distribute_channels("soma","gbar_K_P",0,1.000000,0.000000,0.000000,0.000000,0.0002080000)
# neu.distribute_channels("soma","gbar_NaTg",0,1.000000,0.000000,0.000000,0.000000,0.2720000000)

# for sec in neu.myelin:
#     sec.Ra = 100
#     sec.cm = 0.02

# ps = h.PlotShape(False)
# ps.variable("v")
# ps.plot(plt, cmap=cm.jet)
# plt.show()

# # START the current clamp
# stim = h.IClamp(pyr.cell.soma[0](0.5))
# stim.delay = 100  # start of the current injection (ms)
# stim.dur = 600  # duration of the injection (ms)
# stim.amp = -0.05 # amplitude (nA)
# # END the current clamp

# # START the recording
# # Record time
# t = h.Vector()
# t.record(h._ref_t)

# # Record voltage from soma
# v = h.Vector()
# v.record(pyr.cell.soma[0](0.5)._ref_v)
# # END the recording


# h.finitialize(v_init * mV)
# h.continuerun(1000 * ms)

# plt.plot(t, v)
# plt.xlabel('Time (ms)')
# plt.ylabel('Membrane potential (mV)')
# plt.show()

# input('') # Ensure the closure of the program running, press enter to terminate the run. 
    