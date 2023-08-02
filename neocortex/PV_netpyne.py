from neuron import h

class PV_netpyne:
    def __init__(self):

        self._setup_morphology()
        self._setup_biophysics()
        
    
    def _setup_morphology(self):
        h.load_file('stdrun.hoc')
        h.nrn_load_dll('nrnmech.dll')
        h.load_file("import3d.hoc")

        h.load_file('NeuronTemplate.hoc')

        self.cell = h.NeuronTemplate("HL23PV.swc")
    
    def _setup_biophysics(self):
        h.load_file('biophys_HL23PV.hoc')
        h.biophys_HL23PV(self.cell)