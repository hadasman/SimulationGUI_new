from __future__ import division
from collections import OrderedDict
import sys, pdb, matplotlib, os
import numpy as np

os.chdir('../function_scripts')
from neuron import gui,h
matplotlib.use("TkAgg")

os.chdir('../function_scripts')
import synapse_functions

# Keep at end of imports
from CableGUI import *
from Simulators import *

# =========================== Updates Functions (connection between GUI and Simulator) ===========================
def UpdateEntryParams(Simulator, GUI):

	'''
	Update parameters according to values in Entries. This is used once upon initialization of the GUI (when mode==
	'initialized'), and then again each time there is an update to one or more of the GUI features (mode=='update').

	Uses global variables:
			Example: {'simulation': {'attribute': value}, 
					  'synapses': {'exc': {'attribute': value}, 'inh': {'attribute': value}}, 
					   etc.}
		- Simulator: Simulators object that is updated according to user choices
		- GUI: GUI object through which user makes choices
	'''

	all_dict = {}
	all_dict.update(GUI.Entries['simulation'])
	all_dict.update(GUI.Entries['synapse_exc'])
	all_dict.update(GUI.Entries['synapse_inh'])
	if Simulator.soma:
		all_dict.update(GUI.Entries['soma'])
	if Simulator.spines_exist:
		all_dict.update(GUI.Entries['spines'])

	for item in all_dict:
		temp_val = all_dict[item]['Var'].get()
		# Set Simulator properties
		if 'h.' not in item:
			try:
				exec('Simulator.' + item + " = %s"%float(temp_val))
			except:
				exec('Simulator.' + item + " = \'%s\'"%temp_val)

		# Set simulation parameters (inside hoc interpreter(=h))
		elif 'h.' in item:				
			exec(item + " = %s"%float(temp_val))

	# Set spine parameters manually
	# !Add Simulator.UpdateSpines()
	if Simulator.spines_exist:
		Simulator.UpdateSpineParams('neck', diam = Simulator.neck_diam,
											L = Simulator.neck_L,
											Ra = Simulator.neck_Ra)
		Simulator.UpdateSpineParams('head', diam = Simulator.head_radius*2,
											L = Simulator.head_radius*2)

	# Set synapse parameters manually

	exc_tstart = Simulator.t_start
	inh_tstart = exc_tstart + Simulator.dEI
	for att in ['exc_tstart', 'inh_tstart']:
		setattr(synapse_functions, att, eval(att))

def UpdateSynLocs():

	synapse_filled_length = Simulator.dend.L - Simulator.start_syns
	
	exc_start = Simulator.start_syns
	exc_dist = GUI.Entries['synapse_exc']['exc_dist']['Var'].get()
	Simulator.CreateSynLocs('exc', exc_start, synapse_filled_length, Simulator.exc_dense, Simulator.dend.L, Simulator.exc_locs, dist_ = exc_dist)

	inh_start = Simulator.start_syns
	inh_dist = GUI.Entries['synapse_inh']['inh_dist']['Var'].get()
	Simulator.CreateSynLocs('inh', inh_start, synapse_filled_length, Simulator.inh_dense, Simulator.dend.L, Simulator.inh_locs, dist_ = inh_dist)

def UpdateErrorLabel(error_dict, error_msg=None):
	label = error_dict['Label']
	var = error_dict['Var']
	
	if error_msg:
		if 'Warning' in error_msg:
			var.set(error_msg)
			label['fg'] = 'darkorange'
		elif 'Error' in error_msg:
			var.set(error_msg)
			label['fg'] = 'red'
		else:
			var.set(error_msg)
			label['fg'] = 'blue'
	else:
		var.set('No error. Good!')
		label['fg'] = 'green'

def UpdatePresentedValues(labels_dict):
	lambda_ = Simulator.Calculate('lambda', Simulator.dend)
	labels_dict['lambda']['Var'].set(u"\u03BB = %.2f cm = %.2f \u03BCm"%(lambda_, lambda_*1e4))
	labels_dict['Ri']['Var'].set(u"R\u1d62 = %.2f M\u2126"%Simulator.Calculate('Rinput', Simulator.dend))
	labels_dict['n_exc']['Var'].set(u"N\u2091\u2093 = %i"%Simulator.n_exc)
	labels_dict['n_inh']['Var'].set(u"N\u2091\u2093 = %i"%Simulator.n_inh)

# ================================= Callback Functions (everything user can activate) =================================
def AddSpines_callback(mode='button_callback'):

	# !! Check if adding spines when dist_='Freeze' keeps synapse locations and just puts spines there
	UpdateEntryParams(Simulator, GUI)
	UpdateSynLocs()

	if Simulator.spines_exist: # Either in the case of removing spines or just updating location (will be False for adding new spines)
		for spine in h.allsec():
			if 'spine' in spine.name():
				h.delete_section(sec=spine)

	if (mode=='update' and Simulator.spines_exist) or (mode=='button_callback' and GUI.Buttons['AddSpines']['text']=='Add Spines'):
		Simulator.spine_heads, Simulator.spine_necks = PutSpines('cable', Simulator.dend, Simulator.exc_locs, 
													neck_diam = float(GUI.Entries['spines']['neck_diam']['Var'].get()), 
													neck_len = float(GUI.Entries['spines']['neck_L']['Var'].get()), 
													head_radius = float(GUI.Entries['spines']['head_radius']['Var'].get()), 
													Ra = float(GUI.Entries['spines']['neck_Ra']['Var'].get()), 
													cm=Simulator.dend.cm, 
													e_pas = Simulator.dend.e_pas, 
													g_pas=4.6716e-5)
		Simulator.spines_exist = True
		GUI.Buttons['AddSpines']['text'] = 'Remove Spines'
		h.define_shape()

		# Unlock spine-related radio buttons
		GUI.RadioButtons['syn_loc']['Buttons'][0].config(state='normal')
		GUI.RadioButtons['volt_loc']['Buttons'][0].config(state='normal')

	elif mode=='button_callback' and GUI.Buttons['AddSpines']['text']=='Remove Spines':

		Simulator.spine_heads, Simulator.spine_necks = [], []
		Simulator.spines_exist = False

		GUI.Buttons['AddSpines']['text'] = 'Add Spines'

		# Make sure user can't put synapses on non-existing spines
		GUI.RadioButtons['syn_loc']['Var'].set(2)
		GUI.RadioButtons['syn_loc']['Buttons'][0].config(state='disabled')
		GUI.RadioButtons['volt_loc']['Var'].set(2)
		GUI.RadioButtons['volt_loc']['Buttons'][0].config(state='disabled')


	Simulator.PlaceSynapses('exc')
	GUI.DrawSections(colors)
	UpdatePresentedValues(GUI.ChangingLabels)

def AddSoma_callback():

	s = ttk.Style()
	s.configure('Gray.TEntry', background='gray')

	if GUI.Buttons['AddSoma']['text'] == 'Add Soma':

		soma_size = float(GUI.Entries['soma']['soma.diam']['Var'].get())
		soma_cm = float(GUI.Entries['soma']['soma.cm']['Var'].get())
		Simulator.CreateCompartment('soma', 
									L = soma_size, 
									diam = soma_size, 
									Ra = 110, 
									e_pas = h.v_init, 
									g_pas = 1.0 / 1500.0, 
									nseg = int(soma_size) * 5, 
									cm = soma_cm)
		
		h.disconnect(sec=Simulator.soma)
		Simulator.dend.connect(Simulator.soma, 1, 0) # Disconnect from parents if exist (following weird bug in which soma was created as child of last created section (spine head))
		h.define_shape() # !Explain this

		GUI.Buttons['AddSoma']['text'] = 'Remove Soma' #!Toggle soma function inside GUI
		GUI.RadioButtons['volt_loc']['Buttons'][2].config(state='normal')

	elif GUI.Buttons['AddSoma']['text'] == 'Remove Soma':
		h.delete_section(sec=Simulator.soma)
		Simulator.soma = None
		GUI.Buttons['AddSoma']['text'] = 'Add Soma'

		GUI.RadioButtons['volt_loc']['Var'].set(2)
		GUI.RadioButtons['volt_loc']['Buttons'][2].config(state='disabled')

	GUI.DrawSections(colors)
	UpdatePresentedValues(GUI.ChangingLabels)

def RunSim_callback():
	# !Insert most of this to Simulator
	for button in GUI.Buttons:
		if button is not 'Reset':
			GUI.Buttons[button].config(state='disabled')

	Simulator.RunSim()

	VoltRadio_callback([], [], [])

def Reset_callback():
	for button in GUI.Buttons:
		GUI.Buttons[button].config(state='normal')

	freeze_plots = GUI.CheckBoxes['freeze_plots']['Var'].get()
	if not freeze_plots:
		GUI.Figures['volt']['ax'].clear()
		GUI.Figures['volt']['graph'].draw()

	Simulator.vectors['t']			= []
	Simulator.vectors['shaft_v']	= []
	Simulator.vectors['spine_v']	= []
	Simulator.vectors['soma_v']		= []

def VoltRadio_callback(a, b, c):
	'''
	Callback to radio button indicating location of voltage traces to be shown in figure.
	Inputs are eventdata (automatically generated by Tkinter).
	'''
	ax = GUI.Figures['volt']['ax']
	graph = GUI.Figures['volt']['graph']
	var = GUI.RadioButtons['volt_loc']['Var']

	freeze_plots = GUI.CheckBoxes['freeze_plots']['Var'].get()
	if not freeze_plots:
		ax.clear()

	if var.get() == 1:
		ax.plot(Simulator.vectors['t'], Simulator.vectors['spine_v'], color=colors['spine_head'])
		v_where = 'Spine'
	
	elif var.get() == 2:
		ax.plot(Simulator.vectors['t'], Simulator.vectors['shaft_v'], color=colors['dend'])
		v_where = 'Shaft'
	
	elif var.get() == 3:
		ax.plot(Simulator.vectors['t'], Simulator.vectors['soma_v'], color=colors['soma'])
		v_where = 'Soma'
	
	elif var.get() == 4:
		ax.plot(Simulator.vectors['t'], Simulator.vectors['shaft_v'], color=colors['dend'], label='Shaft')
		
		if len(Simulator.vectors['spine_v']) > 0:
			ax.plot(Simulator.vectors['t'], Simulator.vectors['spine_v'], color=colors['spine_head'], label='Spine Head')
		
		if len(Simulator.vectors['soma_v']) > 0:
			ax.plot(Simulator.vectors['t'], Simulator.vectors['soma_v'], color=colors['soma'], label='Soma')

		ax.legend()
		v_where = 'All Recorded Locations'

	ax.set_title('Voltage Trace in %s'%v_where)
	graph.draw()

def UpdateMorph_callback():
	"""
	Run an update call on all GUI paraeters and functions. Including:
		- Updating parameters taken from entries
		- Updating synapse locations
		- Placing synapses according to updated locations
		- Re-drawing morphology
		- Updating calculated values presented to user
	"""

	UpdateEntryParams(Simulator, GUI)

	# Update spine locations
	UpdateSynLocs()	
	AddSpines_callback(mode='update')
	Simulator.PlaceSynapses('exc')
	Simulator.PlaceSynapses('inh')
	GUI.DrawSections(colors)

	UpdatePresentedValues(GUI.ChangingLabels)
	UpdateErrorLabel(GUI.ChangingLabels['errors'])

def EntryTracking_callback(a, b, c):
	'''
	Function called whenever an entry value is changed by user. 
	Inputs are eventdata (automatically generated by Tkinter).
	'''

	if not GUI.suppress_entry_callback:
		UpdateErrorLabel(GUI.ChangingLabels['errors'], \
			'Warning: You changed something! Press \'Update Morphology\' to implement changes before running simulation')

# ================================================ Initialize Parameters =================================================
#! Add logs
colors = {	'dend': 'black', 
			'soma': 'lightblue', 
			'spine_neck': 'royalblue', 
			'spine_head': 'darkblue'}

# !Get this inside Simulator and if arguments not given, default to it

# All elements in this dictionary will appear as entries in GUI (separate keys for separate blocks)
UserParamDict = {
'simulation': OrderedDict([
	('dend.L', [u'Branch Length [\u03BCm]', 80]), 
	('dend.diam', [u'Branch Diameter [\u03BCm]', 0.26]),
	('dend.cm', [u'Branch Capacitance [\u03BCF/cm\u00B2]', 1]),
	('start_syns', [u'Synapse Start [\u03BCm]', 4]), 		
	('dEI', ['dt(E, I) [ms]', 10]),
	('h.v_init', ['Resting Potential [mV]', -75]),
	('h.tstop', ['Simulation Time', 250])
	]),
'synapses':  {'exc': OrderedDict([
	('exc_dense', ['Exc. Density', 0.58]),
	('exc_g_max', ['Exc. g_max', 0.4]),
	('exc_dist', ['Exc. Distribution', 'Uniform']),
	]),
'inh': OrderedDict([
	('inh_dense', ['Inh. Density', 0.2]), 
	('inh_g_max', ['Inh. gmax', 0.5]),
	('inh_dist', ['Inh. Distribution', 'Uniform'])
	])},
'spines': OrderedDict([
	('neck_diam', [u'Neck diam [\u03BCm]', 0.0394]),
	('neck_L', [u'Neck L [\u03BCm]', 1]),
	('neck_Ra', [u'Neck Ra (per unit area) [\u03A9-cm]', 50]),
	('head_radius', [u'Head radius [\u03BCm]', 0.297])	
	]),
'soma': OrderedDict([
	('soma.diam', [u'diam [\u03BCm]', 10]),
	('soma.cm', [u'Cm [\u03BCF/cm\u00B2]', 1]),
	('dend.Ra', ['Ra', 110])
	])
}

SimInitDict = {}
SimInitDict.update(UserParamDict['simulation'])
SimInitDict.update(UserParamDict['synapses']['exc'])
SimInitDict.update(UserParamDict['synapses']['inh'])
 
GUI = CableGUI('Cable Simulation', (1500, 750), 'white') # Main window object, from class CableGUI
Simulator = Simulators(SimInitDict, GUI)				 # Simulator object, from class Simulators

Simulator.CreateCompartment('dend', 
	L = UserParamDict['simulation']['dend.L'][1], 
	diam = UserParamDict['simulation']['dend.diam'][1], 
	cm = UserParamDict['simulation']['dend.cm'][1],
	Ra = 150, 
	e_pas = UserParamDict['simulation']['h.v_init'][1], 
	g_pas = 1.0 / 1500.0, 
	nseg = int(UserParamDict['simulation']['dend.L'][1]) * 5)

# ================================================ Create GUI - Main Window (root) =================================================	
GUI.AddLabel(GUI, text='Cable Simulation GUI', font=('TkDefaultFont', 30), foreground='darkred', background='white', sticky='W')

GUI.AddButton(GUI, 'RunSim', 'Run Simulation', command=RunSim_callback, row=1, column=0)
GUI.AddButton(GUI, 'Reset', 'Reset Simulation', command=Reset_callback, row=2, column=0)
GUI.AddButton(GUI, 'UpdateMorph', 'Update Morphology', command=UpdateMorph_callback, row=1, column=1)

GUI.AddTab('param_tab', 'Parameter Set', rowspan=3, columnspan=10)
GUI.AddTab('fig_tab', 'Figures')

# ================================================ Param Tab Design =================================================
# Initialize parameters
params_col, image_row_span, AddColumn = 0, len(UserParamDict['simulation']) + 1, 5
image_dim, reported_dim = (1, 4), (2, 7)

# Simulation- and synapse-related user entries
GUI.AddEntries(GUI.Tabs['param_tab'], 'simulation', UserParamDict['simulation'], 
						command=EntryTracking_callback, labelColumn=params_col, entryTitle='Set Simulation Parameters')


GUI.AddRadioButton(GUI.Tabs['param_tab'], 'syn_loc', ["Synapses on Spines", "Synapses on Shaft"], 
						command=EntryTracking_callback, default_val=2, sticky='WE', column=0)

GUI.RadioButtons['syn_loc']['Buttons'][0].config(state='disabled')

GUI.AddEntries(GUI.Tabs['param_tab'], 'synapse_exc', UserParamDict['synapses']['exc'], 
						command=EntryTracking_callback, labelColumn=params_col, entryTitle='Exc. Synapses')

GUI.AddEntries(GUI.Tabs['param_tab'], 'synapse_inh', UserParamDict['synapses']['inh'], 
						command=EntryTracking_callback, labelColumn=params_col, entryTitle='Inh. Synapses')

# Morphology figure
GUI.PutFigure(GUI.Tabs['param_tab'], 'morph', xlabel='X', ylabel='Y', facecolor='darkgray', 
						row=image_dim[0], column=image_dim[1], rowspan=image_row_span, columnspan=3, padx=45)

# Add compartments buttons
AddSpines_row = image_dim[0] + image_row_span
GUI.AddButton(GUI.Tabs['param_tab'], 'AddSpines', 'Add Spines', 
						command=AddSpines_callback, row=AddSpines_row, column=AddColumn)
GUI.AddEntries(GUI.Tabs['param_tab'], 'spines', UserParamDict['spines'], 
						command=EntryTracking_callback, labelRow=AddSpines_row, labelColumn=AddColumn+1)
Simulator.spines_exist = False

print('!!! Make actual spine params change with entry!!!')

AddSoma_row = AddSpines_row + 1 + len(GUI.Entries['spines'])
GUI.AddButton(GUI.Tabs['param_tab'], 'AddSoma', 'Add Soma', 
						command=AddSoma_callback, row=AddSoma_row, column=AddColumn)
GUI.AddEntries(GUI.Tabs['param_tab'], 'soma', UserParamDict['soma'], 
						command=EntryTracking_callback, labelRow=AddSoma_row, labelColumn=AddColumn+1)

# Add text presented to user
GUI.AddLabel(GUI.Tabs['param_tab'], text='Calculated Values', font=15, row=reported_dim[0], column=reported_dim[1])
for i, name in enumerate(['lambda', 'Ri', 'n_exc', 'n_inh']):
	GUI.AddLabel(GUI.Tabs['param_tab'], varName=name, font=15, fg='darkgreen',
						row=reported_dim[0]+i+1, column=reported_dim[1], padx=5, sticky='We')
UpdatePresentedValues(GUI.ChangingLabels)

GUI.AddLabel(GUI.Tabs['param_tab'], varName='errors', font=15, column=0, columnspan=10, padx=5, sticky='WE')
UpdateErrorLabel(GUI.ChangingLabels['errors'], None)

# ================================================ Figure Tab Design =================================================
GUI.PutFigure(GUI.Tabs['fig_tab'], 'volt', xlabel='T (ms)', ylabel='Voltage (mV)', facecolor = 'gray', row=1)
GUI.AddCheckbox(GUI.Tabs['fig_tab'], "Freeze Plots", 'freeze_plots', False, sticky='W')  

GUI.AddRadioButton(GUI.Tabs['fig_tab'], 'volt_loc', 
	["Show Spine Voltage", "Show Shaft Voltage", "Show Soma Voltage", "Show Overlay of Voltages"], 
	command=VoltRadio_callback, default_val=2, sticky='WE')
GUI.RadioButtons['volt_loc']['Buttons'][0].config(state='disabled')
GUI.RadioButtons['volt_loc']['Buttons'][2].config(state='disabled')

# ================================================ Make Updates =================================================
UpdateSynLocs()
Simulator.PlaceSynapses('exc')
Simulator.PlaceSynapses('inh')
GUI.DrawSections(colors) 

print('*****\nTO DO:\n \
	- After defaulting to uniform locations change exc_dense to \'Uniform\'\
	- Add option for changing spine neck resistance, length, diam and head size!\n\
	- In freeze plots think of option to go freeze in all compartments (maybe put on different axes and show one each time if possible?)\n\
	- IMPORTANT: After updating entry params make sure real values update!\n*****')

