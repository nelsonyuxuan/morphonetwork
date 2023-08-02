# Morpholog file and biophysical parameters

# Packages
from neuron import h
from neuron.units import mV, ms
import sys, math
from matplotlib import pyplot
h.load_file('stdrun.hoc')
h.nrn_load_dll('nrnmech.dll')
h.load_file("import3d.hoc")
import numpy as np
from matplotlib import cm

class HL23SST:
    def __init__(self, morphology):
        self.morphology = morphology
        self.load_morphology()
        self.discretize()

        # Global variable that stored the reversal potential for initialization 
        self.rest = -81.5

        # When importing the single neuron in the network, comment out below 
        # self.simply_passive()
     
        # When performing the effect of DC on single neuron, commment out the parts below for 
        self.general_mech()
        self.soma_mech()
        self.dend_mech()
        self.axon_mech()


        # Insert of extracellular mechanism for all sections
        # h.extracellular.insert([sec for sec in self.all])

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
    
    # This function is only applying passive membrane potential for the effect of DC 
    # for single neuron
    def simply_passive(self):
        h.pas.insert(self.all)
        for sec in self.all:
            sec.Ra = 100
            for seg in sec:
                seg.g_pas = 0.0000232
    
    def general_mech(self):
        
        for sec in self.all:

            # Insertion of mechanism
            h.pas.insert(self.all)
            h.Ih.insert(self.all)
            
            # Biophysical parameters
            sec.Ra = 100
            sec.cm = 1
            sec.e_pas = -81.5
            sec.g_pas = 0.0000232

    
    def soma_mech(self):
        
        for sec in self.soma:
            
            # Insertion of mechanism
            h.NaTg.insert(sec)
            h.K_P.insert(sec)
            h.K_T.insert(sec)
            h.Kv3_1.insert(sec)
            h.Im.insert(sec)
            h.SK.insert(sec)
            h.Ca_HVA.insert(sec)
            h.Ca_LVA.insert(sec)
            h.CaDynamics.insert(sec)
            
            for seg in sec:

                # Biophysical parameters
                seg.ek = -85
                seg.ena = 50
                seg.gamma_CaDynamics = 0.0005
                seg.vshiftm_NaTg = 13
                seg.vshifth_NaTg = 15
                seg.slopem_NaTg = 7
                seg.gbar_Ih = 0.0000431
    
    def dend_mech(self):
        
        for sec in self.dend:
            
            # Biophysical parameters 
            sec.cm = 1
            
            for seg in sec:
                seg.gbar_Ih = 0.0000949

    def axon_mech(self):
        
        for sec in self.axon:
            
            # Insertion of mechanism
            h.SK.insert(sec)
            h.Ca_HVA.insert(sec)
            h.Ca_LVA.insert(sec)
            h.K_T.insert(sec)
            h.K_P.insert(sec)
            h.Nap.insert(sec)
            h.Kv3_1.insert(sec)
            h.NaTg.insert(sec)
            h.CaDynamics.insert(sec)
            h.Im.insert(sec)
            
            for seg in sec:

                # Biophysical parameters
                seg.ek = -85
                seg.ena = 50
                seg.vshiftm_NaTg = 10 
                seg.slopem_NaTg = 9
                seg.gamma_CaDynamics = 0.0005

    def fromtodistance(self, origin, to):
        ''''
        Initilize the zero point used in the distance function and return the distance a
        certain section from the center of the soma
        '''

        h.distance(0, origin.x, sec = origin.sec)
        
        return h.distance(to.x, sec = to.sec)

    def returnSegmentCoordinates(self, section):
    # Get section 3d coordinates and put in numpy array
        n3d = section.n3d()
        x3d = np.empty(n3d)
        y3d = np.empty(n3d)
        z3d = np.empty(n3d)
        L = np.empty(n3d)
        for i in range(n3d):
            x3d[i]=section.x3d(i)
            y3d[i]=section.y3d(i)
            z3d[i]=section.z3d(i)

        # Compute length of each 3d segment
        for i in range(n3d):
            if i==0:
                L[i]=0
            else:
                L[i]=np.sqrt((x3d[i]-x3d[i-1])**2 + (y3d[i]-y3d[i-1])**2 + (z3d[i]-z3d[i-1])**2)

        # Get cumulative length of 3d segments
        cumLength = np.cumsum(L)

        N = section.nseg
        
        # Now upsample coordinates to segment locations
        xCoord = np.empty(N)
        yCoord = np.empty(N)
        zCoord = np.empty(N)
        dx = section.L / (N-1)
        for n in range(N):
            if n==N-1: 
                xCoord[n]=x3d[-1]
                yCoord[n]=y3d[-1]
                zCoord[n]=z3d[-1]
            else:
                cIdxStart = np.where(n*dx >= cumLength)[0][-1] # which idx of 3d segments are we starting at
                cDistFrom3dStart = n*dx - cumLength[cIdxStart] # how far along that segment is this upsampled coordinate
                cFraction3dLength = cDistFrom3dStart / L[cIdxStart+1] # what's the fractional distance along this 3d segment
                # compute x and y positions
                xCoord[n] = x3d[cIdxStart] + cFraction3dLength*(x3d[cIdxStart+1] - x3d[cIdxStart])
                yCoord[n] = y3d[cIdxStart] + cFraction3dLength*(y3d[cIdxStart+1] - y3d[cIdxStart])
                zCoord[n] = z3d[cIdxStart] + cFraction3dLength*(z3d[cIdxStart+1] - z3d[cIdxStart])
        return xCoord, yCoord, zCoord
    
    def extrema(self):
        """Give the bounding box that contains the cell"""
        # Units used in the extrema function are in microns um
        xlo = ylo = zlo = xhi = yhi = zhi = None
        for sec in self.all:
            n3d = sec.n3d()
            xs = [sec.x3d(i) for i in range(n3d)]
            ys = [sec.y3d(i) for i in range(n3d)]
            zs = [sec.z3d(i) for i in range(n3d)]
            my_xlo, my_ylo, my_zlo = min(xs), min(ys), min(zs)
            my_xhi, my_yhi, my_zhi = max(xs), max(ys), max(zs)
            if xlo is None:
                xlo, ylo, zlo = my_xlo, my_ylo, my_zlo
                xhi, yhi, zhi = my_xhi, my_yhi, my_zhi
            else:
                xlo, ylo, zlo = min(xlo, my_xlo), min(ylo, my_ylo), min(zlo, my_zlo)
                xhi, yhi, zhi = max(xhi, my_xhi), max(yhi, my_yhi), max(zhi, my_zhi)
        return (xlo, ylo, zlo, xhi, yhi, zhi)

    def apply_Eext_uniform_py(self, E0):
        """
        Trying to apply directly the vext via e_extracellular in python
        """


        for i,sec in enumerate(self.all):

            ## taking center coordinates
            # x_c = (sec.x3d(0) + sec.x3d(sec.n3d() -1))/2
            # y_c = (sec.y3d(0) + sec.y3d(sec.n3d() -1))/2
            # z_c = (sec.z3d(0) + sec.z3d(sec.n3d() -1))/2     
            # v_ext = -1e-3*(E0[0] * (x_c ) + \
            #               E0[1] * (y_c ) + \
            #               E0[2] * (z_c ) )

            xCoord, yCoord, zCoord = self.returnSegmentCoordinates(sec)

            for j,seg in enumerate(sec):
                
                # Obtaining the segment xyz coordinates from the primary function
                # v_ext = -1e-3*(E0[0] * self.x[i][j]  + \
                #                 E0[1] * self.y[i][j] + \
                #                 E0[2] * self.z[i][j] )

                # Obtaining the segment xyz coordinates from the second function
                v_ext = -1e-3*(E0[0] * xCoord[j]  + \
                                 E0[1] * yCoord[j] + \
                                 E0[2] * zCoord[j] )

                seg.extracellular.e = v_ext               
    # END apply_Eext_uniform_py()

    def range_depolarization(self):

        for i,sec in enumerate(self.all):
            if i==0: vmax = vmin = sec.v 
            if sec.v > vmax:  vmax = sec.v 
            if sec.v < vmin:  vmin = sec.v

        self.vmin = vmin
        self.vmax = vmax
    ## END range_depolarization()

    
            
# Building the cell
c = HL23SST("HL23SST")

# # Set up the electrical field 
# theta = 90 * np.pi/180
# E = [-np.sin(theta),-np.cos(theta),0]

# # Apply extracellular field 
# c.apply_Eext_uniform_py(E)


# h.finitialize(c.rest * mV)
# h.continuerun(200 * ms)

# c.range_depolarization()
# ps = h.PlotShape(False)
# ps.variable("v")
# # ps.scale(-85.5, -83.6)
# ps.scale(c.vmin, c.vmax)
# # ps.exec_menu("Shape Plot")
# ps.plot(pyplot, cmap=cm.jet)
# pyplot.show()

# START the current clamp
stim = h.IClamp(c.soma[0](0.5))
stim.delay = 100  # start of the current injection (ms)
stim.dur = 1200 # duration of the injection (ms)
stim.amp = 0.5 # amplitude (nA)
# END the current clamp

# START the recording
# Record time
t = h.Vector()
t.record(h._ref_t)

# Record voltage from soma
v = h.Vector()
v.record(c.soma[0](0.5)._ref_v)
# END the recording


h.finitialize(c.rest * mV)
h.continuerun(1500 * ms)

pyplot.plot(t, v)
pyplot.xlabel('Time (ms)')
pyplot.ylabel('Membrane potential (mV)')
pyplot.show()


input('') # Ensure the closure of the program running, press enter to terminate the run. 
    