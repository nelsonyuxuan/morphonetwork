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

h.load_file('stdrun.hoc')
h.nrn_load_dll('nrnmech.dll')
h.load_file("import3d.hoc")

h.load_file('NeuronTemplate.hoc')
h.load_file('biophys_HL23PV.hoc')

neu = h.NeuronTemplate("HL23PV.swc")

h.biophys_HL23PV(neu)

# START the current clamp
stim = h.IClamp(neu.soma[0](0.5))
stim.delay = 100  # start of the current injection (ms)
stim.dur = 600  # duration of the injection (ms)
stim.amp = -0.05 # amplitude (nA)
# END the current clamp

# START the recording
# Record time
t = h.Vector()
t.record(h._ref_t)

# Record voltage from soma
v = h.Vector()
v.record(neu.soma[0](0.5)._ref_v)
# END the recording


h.finitialize(v_init * mV)
h.continuerun(1000 * ms)

plt.plot(t, v)
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.show()

input('') # Ensure the closure of the program running, press enter to terminate the run. 
    