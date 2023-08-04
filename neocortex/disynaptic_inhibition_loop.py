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

# SUCCESSFUL import using the HOC template directly
# cellRule = netParams.importCellParams(
#         label='PYR',
#         conds={'cellType': 'PYR'},
#         fileName='PYRtemplate.hoc',
#         cellName='PYRtemplate',
#         importSynMechs=False)

netParams.sizeX = 500.0 # x-dimension (horizontal length) size in um
netParams.sizeY = 500.0 # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = 950.0 # z-dimension (horizontal depth) size in um
netParams.propVelocity = 100.0 # propagation velocity (um/ms)
netParams.probLengthConst = 100.0 # length constant for conn probability (um)

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

'''
20230711
Based on the source code of the paper, changing the synapstic mechanisms to the ones that paper used: excitatory: ProbAMPANMDA, inhibitroy: ProbUFDsyn
'''

# netParams.synMechParams['excite_pyr_pre_AMPA'] = {'mod': 'Exp2Syn', 'tau1': 0.3, 'tau2': 3, 'e': 0}  # excitatory synaptic mechanism
# netParams.synMechParams['inhibit_sst_pre'] = {'mod': 'Exp2Syn', 'tau1': 1, 'tau2': 10, 'e': -80}  # inbitary synaptic mechanism
# netParams.synMechParams['excite_pyr_pre_NMDA'] = {'mod': 'NMDA', 'tau_r_NMDA': 2, 'tau_d_NMDA': 65, 'e': 0}  # excitatory synaptic mechanism

netParams.synMechParams['Excitatory_PYR_SST'] = {'mod': 'ProbAMPANMDA', 'tau_r_AMPA': 0.3, 'tau_d_AMPA': 3, 'tau_r_NMDA': 2, 'tau_d_NMDA': 65, 'e': 0, 'u0':0, 'Dep': 140, 'Fac': 670, 'Use': 0.09, 'gmax': 0.38}
netParams.synMechParams['Inhibitory_SST_PYR'] = {'mod': 'ProbUDFsyn', 'tau_r': 1, 'tau_d': 10, 'e': -80, 'u0':0, 'Dep': 1300, 'Fac': 2, 'Use': 0.3, 'gmax': 1.24}
# netParams.synMechParams['excite_pyr_pre_NMDA'] = {'mod': 'Exp2Syn', 'tau1': 2, 'tau2': 65, 'e': 0}  # excitatory synaptic mechanism
#========================================================================


#==========================STIMULATION PARAMETERS (TEST COMPLETED)==============================================
# Insert a current clamp to generate initial spike in the PYR cell.
netParams.stimSourceParams['CurrentClamp'] = {
    'type': 'IClamp',
    'amp': 1.5,  # amplitude of current, in nA
    'dur': 3000,  # duration of current, in ms
    'delay': 200  # delay before current onset, in ms
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

# Setting up excitatory synapses connection independently with different delay parameters
# Number of connections created
# n_c = 1

# for i in range(n_c):

#     delay = 50 + i*10

#     netParams.connParams[f'PYR->SST_{i}'] = {
#         'preConds': {'pop': 'pre'},
#         'postConds': {'pop': 'sst'},
#         'weight': 0.0001,  # Not really sure about the weight parameter values used for the parameters, waiting for adjustment. It seems like increasing the weight value will decrease the amplitude of the EPSPs
#         'delay': delay, # Millisecond delay between when the pre-synaptic neuron fires and when that signal affects the post-synaptic neuron  
#         'synMech': 'excite_pyr_pre',
#         'sec': 'dend', # section of the postsynaptic cell to connect to 
#         'synsPerConn': 1
#     }
netParams.connParams['PYR->SST'] = {
    'preConds': {'pop': 'pre'},
    'postConds': {'pop': 'sst'},
    'weight':  0.0006,  # Not really sure about the weight parameter values used for the parameters, waiting for adjustment. It seems like increasing the weight value will decrease the amplitude of the EPSPs
    'delay': 0.5, # Millisecond delay between when the pre-synaptic neuron fires and when that signal affects the post-synaptic neuron. Original: 'dist_3D/propVelocity'.
    'synMech': 'Excitatory_PYR_SST', # ['excite_pyr_pre_AMPA', 'excite_pyr_pre_NMDA']
    'sec': 'dend', # section of the postsynaptic cell to connect to 
    'synsPerConn': 8  
}

# Inhibitory synapses (TEST COMPLETED)
netParams.connParams['SST->PYR'] = {
    'preConds': {'pop': 'sst'},
    'postConds': {'pop': 'post'},
    'weight': 0.006,  # Not really sure about the weight parameter values used for the parameters, waiting for adjustment
    'delay': 'dist_3D/propVelocity', # Millisecond delay between when the pre-synaptic neuron fires and when that signal affects the post-synaptic neuron  
    'synMech': 'Inhibitory_SST_PYR',
    'sec': 'apical', # sections of the postsynaptic cell to connect to 
    'synsPerConn': 12 
}
#========================================================================

#==========================STIMULATION PARAMETERS==============================================
# Set up simulation configuration
simConfig = specs.SimConfig()  
simConfig.duration = 4500  # Duration of the simulation, in ms
simConfig.dt = 0.025  # Internal integration timestep to use
simConfig.verbose = False  # Show detailed messages

# Setting up the hParams
simConfig.hParams = {'v_init': -80, 'celcius': 34}  # Set celsius temp and starting v
# simConfig.hParams = {'v_init': -80}  # Set celsius temp and starting v

# simConfig.recordTraces = {
#     'V_soma_pyr_pre': {'sec': 'soma_0','var': 'v', 'conds': {'pop':'pre', 'cellList': [0]}},
#     'V_soma_sst': {'sec': 'soma_0', 'var': 'v', 'conds': {'pop':'sst', 'cellList': [0]}},
#     'V_soma_pyr_post': {'sec': 'soma_0', 'var': 'v', 'conds': {'pop':'post', 'cellList': [0]}}
# }

simConfig.recordTraces = {'V_soma':{'sec':'soma_0','loc':0.5,'var':'v'}}  # Dict with traces to record

# Changing the files and data saving location to the current ones. 
# Result control
simConfig.filename = 'Disynaptic_Inhibition_LOOP_MODEL'  # Set data output location
# simConfig.saveFig = '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData'  # Set plot output location (NOT WORKING)
simConfig.saveFolder = '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/morphonetwork/neocortex/FigureAndData'  # Set plot output location (
simConfig.analysis['plotRaster'] = {'saveFig': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/morphonetwork/neocortex/FigureAndData/raster_plot'}                  # Plot a raster"F""
simConfig.analysis['plot2Dnet'] = {'saveFig': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/morphonetwork/neocortex/FigureAndData/2D_net_plot'}                   # plot 2D cell positions and connections
simConfig.analysis['plotTraces'] = {'include': [('pre', 0), ('sst', 0), ('post', 0)], 'saveFig': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/morphonetwork/neocortex/FigureAndData/Disynaptic_Inhibition_LOOP_MODEL', 'saveData': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/Voltage_trace.pkl'}  

# Trace files
# simConfig.analysis['plotTraces'] = {'include': [('pre',0)], 'saveFig': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_0', 'saveData': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_0.pkl'}  # Plot recorded traces for this list of cells'}
# sim.analysis.plotTraces()

# simConfig.analysis['plotTraces'] = {'include': [('sst',0)], 'saveFig': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_1', 'saveData': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_1.pkl'}  # Plot recorded traces for this list of cells'}

# simConfig.analysis['plotTraces'] = {'include': [('post',0)], 'saveFig': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_2', 'saveData': '/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_2.pkl'}  # Plot recorded traces for this list of cells'}

simConfig.savePickle = True
# Create network and run simulationå
sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig)
#========================================================================

# # Plot and save traces
# for i, cell in enumerate([('pre', 0), ('sst', 0), ('post', 0)]):
#     sim.analysis.plotTraces(include=[cell], saveFig=f'/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_{i}.png', saveData=f'/Users/nelsonwu/Library/CloudStorage/OneDrive-Personal/FridmanLab/Computational_Project/neocortex_ABAN/FigureAndData/V_trace_{i}.pkl')


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

# Then print section names for each cell of desired population
# for cell in sim.net.cells:
#     if cell.tags['pop'] in ['pre', 'post', 'sst']:  # Replace with your population names
#         print(f"Cell {cell.gid} in population {cell.tags['pop']} has the following sections:")
#         for sec in cell.secs:
#             print(f"    {sec}")

# print(sim.net.cells)

# for cell in sim.net.cells:
#     print('Cell GID:', cell.gid)  # print cell global identifier
#     print('Cell Tags:', cell.tags)  # print cell properties
#     print('Cell sections:', cell.secs.keys())  # print section names of the cell
#     print('-------------------')

# #========================================================================