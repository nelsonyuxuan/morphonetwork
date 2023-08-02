import neuron
from neuron import h
from netpyne import specs, sim

# Create a cell parameter dictionary
cellParams = specs.CellParams()

# Load your the necessary HOC file
h.load_file('stdrun.hoc')
h.nrn_load_dll('nrnmech.dll')
h.load_file("import3d.hoc")
h.load_file('NeuronTemplate.hoc')
h.load_file('biophys_HL23SST.hoc')
h.load_file('PYRtemplate.hoc')

netParams = specs.NetParams()  # object of class NetParams to store the network parameters

netParams.sizeX = 500.0 # x-dimension (horizontal length) size in um
netParams.sizeY = 500.0 # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = 950.0 # z-dimension (horizontal depth) size in um
netParams.propVelocity = 100.0 # propagation velocity (um/ms)
netParams.probLengthConst = 100.0 # length constant for conn probability (um)

# SUCCESSFUL import using the HOC template directly
# cellRule = netParams.importCellParams(
#         label='PYR',
#         conds={'cellType': 'PYR'},
#         fileName='PYRtemplate.hoc',
#         cellName='PYRtemplate',
#         importSynMechs=False)

#==========================LOAD CELL PARAMETERS (TEST COMPLETED)==============================================
# PYR cell loading
cellRule = netParams.importCellParams(
        label='PYR',
        conds={'cellType': 'PYR'},
        fileName='PYR_netpyne.py',
        cellName='PYR_netpyne',
        importSynMechs=False)

# SST cell loading
cellRule = netParams.importCellParams(
        label='SST',
        conds={'cellType': 'SST'},
        fileName='SST_netpyne.py',
        cellName='SST_netpyne',
        importSynMechs=False)

# VIP cell loading
# cellRule = netParams.importCellParams(
#         label='VIP',
#         conds={'cellType': 'VIP'},
#         fileName='VIP_netpyne.py',
#         cellName='VIP_netpyne',
#         importSynMechs=False)

# PV cell loading
# cellRule = netParams.importCellParams(
#         label='PV',
#         conds={'cellType': 'PV'},
#         fileName='PV_netpyne.py',
#         cellName='PV_netpyne',
#         importSynMechs=False)
#========================================================================

#==========================POPULATION PARAMETERS (TEST COMPLETED)==============================================
netParams.popParams['pre'] = {'cellType': 'PYR', 'numCells': 1}
netParams.popParams['sst'] = {'cellType': 'SST', 'numCells': 1}
netParams.popParams['post'] = {'cellType': 'PYR', 'numCells': 1}
#========================================================================

#==========================SYNAPTIC PARAMETERS (TEST COMPLETED)==============================================
# Assume that the synaptic mechanisms defined by Exp2Syn is AMPA, which mainly contributes to the generation of EPSPs
# netParams.synMechParams['excite_pyr_pre'] = {'mod': 'Exp2Syn', 'tau1': 0.3, 'tau2': 3, 'e': 0}  # excitatory synaptic mechanism
# netParams.synMechParams['inhibit_sst_pre'] = {'mod': 'Exp2Syn', 'tau1': 1, 'tau2': 10, 'e': -80}  # inbitary synaptic mechanism
netParams.synMechParams['NMDA'] = {'mod': 'Exp2Syn', 'tau1': 2, 'tau2': 65, 'e': 0}  # excitatory synaptic mechanism, NMDA
netParams.synMechParams['AMPA'] = {'mod': 'Exp2Syn', 'tau1': 0.3, 'tau2': 3, 'e': 0}  # excitatory synaptic mechanism, AMPA
netParams.synMechParams['GABA'] = {'mod': 'Exp2Syn', 'tau1': 1, 'tau2': 10, 'e': -80}  # inbitary synaptic mechanism, GABA
#========================================================================


#==========================STIMULATION PARAMETERS (TEST COMPLETED)==============================================
# Insert a current clamp to generate initial spike in the PYR cell.
netParams.stimSourceParams['CurrentClamp'] = {
    'type': 'IClamp',
    # 'amp': 1,  # amplitude of current, in nA
    # 'dur': 600,  # duration of current, in ms
    # 'delay': 50  # delay before current onset, in ms
    'amp': 0.5,  # amplitude of current, in nA
    'dur': 1000,  # duration of current, in ms
    'delay': 30  # delay before current onset, in ms
}

netParams.stimTargetParams['IClamp->PYR'] = {
    'source': 'CurrentClamp',
    'conds': {'pop': 'pre'},
    'sec': 'soma_0',
    'loc': 0.5  # location of the stimulation
}

#========================================================================


#==========================CONNECTION PARAMETERS==============================================
# Excitatary synapses (TEST COMPLETED)
netParams.connParams['PYR->SST'] = {
    'preConds': {'pop': 'pre'},
    'postConds': {'pop': 'sst'},
    # 'weight': 0.01,  # Not really sure about the weight parameter values used for the parameters, waiting for adjustment. It seems like increasing the weight value will decrease the amplitude of the EPSPs
    # 'delay': 10, # Millisecond delay between when the pre-synaptic neuron fires and when that signal affects the post-synaptic neuron  
    # 'synMech': 'excite_pyr_pre',
    'weight': 0.38e-4,  # Not really sure about the weight parameter values used for the parameters, waiting for adjustment. It seems like increasing the weight value will decrease the amplitude of the EPSPs
    # 'weight': '0.005*post_ynorm',   # proportional to the cortical depth of cell
    # 'delay': 0.1, # Millisecond delay between when the pre-synaptic neuron fires and when that signal affects the post-synaptic neuron  
    'delay': 'dist_3D/propVelocity',    # transmission delay (ms) = distance / propagation velocity
    'synMech': ['AMPA', 'NMDA'],
    'sec': 'dend', # section of the postsynaptic cell to connect to 
    # 'synsPerConn': 25  
    'synsPerConn': 8,  
    # 'probability': 0.19
}

# Inhibitory synapses (TEST COMPLETED)
netParams.connParams['SST->PYR'] = {
    'preConds': {'pop': 'sst'},
    'postConds': {'pop': 'post'},
    # 'weight': 0.01,  # Not really sure about the weight parameter values used for the parameters, waiting for adjustment
    # 'delay': 1, # Millisecond delay between when the pre-synaptic neuron fires and when that signal affects the post-synaptic neuron  
    # 'synMech': 'inhibit_sst_pre',
    'weight': 1.24e-4,  # Not really sure about the weight parameter values used for the parameters, waiting for adjustment
    # 'weight': '0.005*post_ynorm',   # proportional to the cortical depth of cell
    # 'delay': 10, # Millisecond delay between when the pre-synaptic neuron fires and when that signal affects the post-synaptic neuron  
    'delay': 'dist_3D/propVelocity',    # transmission delay (ms) = distance / propagation velocity
    'synMech': 'GABA',
    'sec': 'apical', # sections of the postsynaptic cell to connect to 
    # 'synsPerConn': 25 
    'synsPerConn': 12, 
    # 'probability': 0.19
}
#========================================================================

#==========================STIMULATION PARAMETERS==============================================
# Set up simulation configuration
simConfig = specs.SimConfig()  
simConfig.duration = 1.0*1e3  # Duration of the simulation, in ms
simConfig.dt = 0.025  # Internal integration timestep to use
simConfig.verbose = False  # Show detailed messages

# Record traces
# simConfig.recordTraces = {
#     'V_soma_pyr_pre': {'sec': 'soma_0','var': 'v', 'conds': {'pop':'pre', 'cellList': [0]}},
#     'V_soma_sst': {'sec': 'soma_0', 'var': 'v', 'conds': {'pop':'sst', 'cellList': [0]}},
#     'V_soma_pyr_post': {'sec': 'soma_0', 'var': 'v', 'conds': {'pop':'post', 'cellList': [0]}}
# }

simConfig.recordTraces = {'V_soma':{'sec':'soma_0','loc':0.5,'var':'v'}}  # Dict with traces to record

# Result control
simConfig.filename = 'Disynaptic_Inhibition_LOOP_MODEL'  # Set data output location
# simConfig.saveFig = '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData'  # Set plot output location (NOT WORKING)
simConfig.analysis['plotRaster'] = {'saveFig': True}                  # Plot a raster"F""
simConfig.analysis['plot2Dnet'] = {'saveFig': True}                   # plot 2D cell positions and connections
# simConfig.analysis['plotTraces'] = {'include': [('pre', 0), ('sst', 0), ('post', 0)], 'saveFig': True}  # Plot recorded traces for this list of cells
simConfig.analysis['plotTraces'] = {'include': [('pre', 0), ('sst', 0), ('post', 0)], 'oneFigPer':'trace','saveFig': True}  # Plot recorded traces for this list of cells
# simConfig.analysis['plotTraces'] = {'include': [('pre', 0), ('sst', 0), ('post', 0)],'saveFig': True}  # Plot recorded traces for this list of cells

# Create network and run simulation
sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig)
#========================================================================


# #==========================TEST CODE of LOADING BIOPHYSICAL PARAMETERS==============================================
# # # Now you can access the biophysical parameters of your cells
# # cell = sim.net.cells[0]  # Get the first cell
# # cellProps = cell.__dict__  # Get the properties of the cell

# # # Now print the properties
# # for prop, value in cellProps.items():
# #     print(f"{prop}: {value}")

# # Check if the section 'somatic' exists for each cell in your desired population.
# for cell in sim.net.cells:
#     if cell.tags['pop'] in ['pre', 'post', 'sst']:  # Replace with your population names
#         if 'soma' in cell.secs:
#             print(f"Cell {cell.gid} in population {cell.tags['pop']} has a 'soma' section.")
#         else:
#             print(f"Cell {cell.gid} in population {cell.tags['pop']} does not have a 'soma' section.")
