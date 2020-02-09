import os
from neuron import h, gui
os.chdir('../function_scripts')
from synapse_functions import PutSyns, PutSynsOnSpines, PutSpines
import synapse_functions

import copy
import pdb
import numpy as np

class Simulators():

	def __init__(self, InitDict, GUI_obj):
		self.exc_locs 		= []
		self.inh_locs 		= []
		self.n_exc 			= 0
		self.n_inh 			= 0
		self.GUI 			= GUI_obj
		self.vectors		= {'t': [], 'shaft_v': [], 'spine_v': [], 'soma_v': []}
		self.soma 			= None
		self.spine_heads	= None
		self.spine_necks	= None
		self.spines_exist	= False

		for item in InitDict:
			temp_val = InitDict[item][1]
			
			# Set simulator properties
			if '.' not in item:
				if type(temp_val) == str:
					exec('self.' + item + " = \'%s\'"%temp_val)
				else:
					exec('self.' + item + " = %s"%float(temp_val))
			
			# Set simulation parameters (inside hoc interpreter(=h))
			elif 'h.' in item:				
				exec(item + " = %s"%float(temp_val))

		exc_tstart = 100
		inh_tstart = exc_tstart + InitDict['dEI'][1]
		for att in ['exc_tstart', 'inh_tstart']:
			setattr(synapse_functions, att, eval(att))


	def CreateCompartment(self, sec_name, **kwargs):

		h('create ' + sec_name)
		var = eval('h.'+sec_name)

		var.insert('pas')

		if 'soma' in sec_name:
			var.insert('kv')
			var.insert('na')

		for attribute in kwargs:
			setattr(var, attribute, kwargs[attribute])

		h.define_shape()

		exec('self.'+sec_name+' = var')

		return var

	def Calculate(self, func_name, sec):

		def CalcInputR(sec):
			im = h.Impedance()
			im.loc(0.5, sec=sec)
			im.compute(0, 1) # 1st argument: impedance frequency to compute; 2nd argument: when ==1, the calculation considers the effect of differential gating states
			Rinput = im.input(0.5, sec=sec)
			
			return Rinput

		def CalcAxialR(sec):
			"""
			Calculates compartment resistance given resistivity, diameter and length of compartment.

			Inputs:
				- sec: Section object (should include Ra, diam and L properties for volume calculation)
					* Ra = resistivity (specific resistance for unit area) of compartment [ohm*cm] (=> property of cytoplasm)

			Outputs:
				- Ra_T: Total axial resistance [M-ohm] (=> property of cable)
			"""

			L = sec.L * (1e-4) 			# Convert from um to [cm]
			diam = sec.diam * (1e-4)/2   	# Convert from um to [cm]
			A = np.pi * (diam ** 2)	# in units [cm^2]

			Ra_T = (sec.Ra * L) / A 		# units: [ohm-cm] * [cm] / [cm^2] = [ohm]
			Ra_T = Ra_T * (1e-6)		# Convert to M-ohm

			return Ra_T

		def CalcLambda(sec):
			Rm = 1 / sec.g_pas
			diam = sec.diam
			Ra = sec.Ra

			return np.sqrt((Rm * diam) / (Ra * 4))

		if func_name == 'Rinput':
			return CalcInputR(sec)
		elif func_name == 'Raxial':
			return CalcAxialR(sec)
		elif func_name == 'lambda':
			return CalcLambda(sec)
		else:
			raise Exception('Invalid function name!')

	def CreateSynLocs(self, syn_type, start, interval, density, L, current_locs, dist_ = 'Uniform'):
		'''
		Return synapse locations, given:
			- syn_type: (string) indicating syhapse type ('exc' or 'inh')
			- start: Starting position
			- interval: Interval length which includes synapses 
			- density: synapse density
			- L: branch length
			- dist_: synapse distribution method. Value options are: 'Uniform', 'Random', 'Clustered' (STILL under construction)
			*all given in units of micro-meter
		
		* NOTE: Returned value is synapse locations in arbitrary units of fraction of branch length, 
		and not physical units of micro-meter as given in input (this is what synapse-placing function uses)
		'''
		
		if density == 0:
			locs = []
		else:
			if dist_=='Uniform':
				locs = sorted(np.arange(start, start + interval, 1/density))
				locs = [i / L for i in locs]
			
			elif dist_ == 'Random':
				
				n_syns = int(np.ceil(density * interval))
				locs = start + (np.random.rand(n_syns)) * interval # Put only on relevant section
				locs = sorted([i/L for i in locs])			   # Return to 0-1 range

			elif dist_ == 'Clustered':
				error_msg = 'Error in synapse locations: \'Clustered\' is under Construction. Defaulting to Uniform placement...'
				print(error_msg)
				self.GUI.PopMessageBox(error_msg)

				self.CreateSynLocs(syn_type, start, interval, density, L, current_locs, dist_ = 'Uniform') # Recursive

				return # Exit the function, because self.locs has updated inside recursion

			elif dist_ == 'Freeze':
				if current_locs:
					locs = current_locs

				else:
					# In case someone manages to try and freeze non-existing locations
					error_meg = 'Error in synapse locations: Nothing to freeze. Defaulting to Uniform placement...'
					print(error_msg)
					self.GUI.PopMessageBox(error_msg)

					self.CreateSynLocs(syn_type, start, interval, density, L, current_locs, dist_ = 'Uniform') # Recursive				
					
					return # Exit the function, because self.locs has updated inside recursion

			# Sanity check
			if len([i for i in locs if i>1]) > 0: 
				raise Exception('Synapse location goes beyond cable!')


		if syn_type == 'exc':
			self.exc_locs = locs
			self.GUI.suppress_entry_callback = True
			self.GUI.Entries['synapse_exc']['exc_dist']['Var'].set(dist_)
			self.GUI.suppress_entry_callback = False
		elif syn_type == 'inh':
			self.inh_locs = locs
			self.GUI.suppress_entry_callback = True
			self.GUI.Entries['synapse_inh']['inh_dist']['Var'].set(dist_)
			self.GUI.suppress_entry_callback = False

	def PlaceSynapses(self, syn_type):

		"""
		Places synapses on spines or shaft, according to where_syns argument
		"""
		if syn_type == 'exc':
			self.exc_synapses = None # Delete any existing synapses if there are any
			where_syns = self.GUI.RadioButtons['syn_loc']['Var'].get()

			if where_syns == 1 and self.spines_exist: # Spine
				self.exc_synapses, self.exc_netstim, self.exc_netcon = PutSynsOnSpines(self.spine_heads, 'exc', weight=self.exc_g_max)
				
				n_spines = len(self.exc_synapses)
				spines_v = [h.Vector()] * n_spines
				[spines_v[i].record(self.spine_heads[i](1)._ref_v) for i in range(n_spines)]

				self.spines_v = spines_v
			else:
				if where_syns == 1: # If user wanted to put synapses on spines, but there are no spines
					print('***** \n No Spines! Defaulting to synapses on shaft \n *****')
					self.GUI.RadioButtons['syn_loc']['Var'].set(2)
				
				self.exc_synapses, self.exc_netstim, self.exc_netcon = PutSyns(self.dend, self.exc_locs, 'exc', weight=self.exc_g_max)

			self.n_exc = len(self.exc_synapses)

		elif syn_type == 'inh':
			self.inh_synapses = None
			where_syns = 2 # For now inhibitory synapses always on shaft
			self.inh_synapses, self.inh_netstim, self.inh_netcon = PutSyns(self.dend, self.inh_locs, 'inh', weight=self.inh_g_max)
			self.n_inh = len(self.inh_synapses)
