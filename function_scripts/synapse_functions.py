from __future__ import division
import numpy as np
import os
try:
	import cPickle
except:
	import _pickle as cPickle
import pickle
from scipy import stats
import sys
import matplotlib.pyplot as plt; plt.ion()
import pdb
import os
import sys

# home_path = os.path.expanduser('~') 
# path_ = '%s/Documents/GitHubProjects/MIT_spines' %home_path
# os.chdir(path_)

from neuron import h, gui

# ================================================ Define Functions =================================================
def PutSyns(comp_obj, syn_locs, syn_type_string, weight=None):

	"""
	Arguments (by order):
		- section object on which to put synapses
		- array of locations of synapses along the given section
		- string indicating which synapse (exc / inh)
	
	Outputs:
		- synapses: array of synapse objects
		- netstim: array of netstims
		- netcon: array of netcons
	"""

	synapses, netstim, netcon = [], [], []
	for loc in syn_locs:

		if syn_type_string == 'exc':
			if not weight:
				weight = 0.4 # Default value for exc. synapses

			synapses.append(h.ProbAMPANMDA2_RATIO(loc, sec = comp_obj)) # Go back to _RATIO!
			netstim.append(h.NetStim(loc, sec = comp_obj))
			netcon.append(h.NetCon(netstim[-1], synapses[-1]))

			# synapses[-1].mgVoltageCoeff = 0.08

			netstim[-1].number = 1
			netstim[-1].interval = 1
			netstim[-1].noise = 0
			netstim[-1].start = exc_tstart

			netcon[-1].weight[0] = weight # Literature says ~30pS for 1 NMDAR which is 0.03 here (here it's [nS])
			netcon[-1].delay = 0

		elif syn_type_string == 'inh':
			if not weight:
				weight = 0.5 # Default value for inh. synapses

			synapses.append(h.ProbUDFsyn2_lark(loc, sec = comp_obj)) # Go back to _lark!
			netstim.append(h.NetStim(loc, sec = comp_obj))
			netcon.append(h.NetCon(netstim[-1], synapses[-1]))

			synapses[-1].tau_r = 0.18
			synapses[-1].tau_d = 5
			synapses[-1].e = -80
			synapses[-1].Dep = 0
			synapses[-1].Fac = 0
			synapses[-1].Use = 0.25
			synapses[-1].u0 = 0
			synapses[-1].gmax = 0.001 # don't touch - weight conversion factor to (us) times conductance in nS

			netstim[-1].number = 1
			netstim[-1].interval = 1
			netstim[-1].start = inh_tstart
			netstim[-1].noise = 0

			netcon[-1].weight[0] = weight
			netcon[-1].delay = 0

	return synapses, netstim, netcon

def PutSynsOnSpines(spine_objs, syn_type_string, weight=None):
	# Put synapses on spine_head(1)!

	synapses, netstim, netcon = [], [], []
	# syn_locs = [1] * len(spine_objs) # Put spine synapses on edge of head
	
	for spine in spine_objs:
		SYN, NETSTIM, NETCON = PutSyns(spine, [1], syn_type_string, weight=weight)

		synapses.append(SYN); del SYN
		netstim.append(NETSTIM); del NETSTIM
		netcon.append(NETCON); del NETCON
	
	synapses = [j for i in synapses for j in i]
	netstim = [j for i in netstim for j in i]
	netcon = [j for i in netcon for j in i]

	return synapses, netstim, netcon

def PutSpines(section_type, sec_obj, spine_locs, neck_diam = 0, neck_len = 0, head_radius=0, Ra = 0, cm = None, e_pas=None, g_pas=None):
	
	"""
	section type must be 'whole_cell', 'cable' or 'single_branch' (names of template files are '%s_template.hoc'%section_type). 
		- Currently only adapted for 'cable' type because spine heads and necks are accessed from h.allsec() (doesn't distinguish spines on different branches)
		- Currently only adapted for new cable or data from Neurolucida because new hoc template is uploaded and spines added there
	"""

	dend_ref = h.SectionRef(sec_obj)
	neck_diam = neck_diam
	neck_len = neck_len
	head_radius = head_radius
	
	spine_heads, spine_necks = [], []
	for i in range(len(spine_locs)):
		# temp_neck = h.Section(name='spine_neck[%i]'%i)
		# I don't use h.Section() because then h.delete_section() doesn't work
		
		h('create spine_neck_%i'%i)
		exec('h.spine_neck_' + str(i) + '.diam = neck_diam')
		exec('h.spine_neck_' + str(i) + '.L = neck_len')
		exec('h.spine_neck_' + str(i) + '.insert(\'pas\')')
		exec('h.spine_neck_' + str(i) + '.connect(sec_obj(spine_locs[%i]), 0)'%i)

		
		# temp_head = h.Section(name='spine_head[%i]'%i)
		h('create spine_head_%i'%i)
		exec('h.spine_head_' + str(i) + '.diam = 2 * head_radius')
		exec('h.spine_head_' + str(i) + '.L = 2 * head_radius')
		exec('h.spine_head_' + str(i) + '.insert(\'pas\')')
		exec('h.spine_head_' + str(i) + '.connect(h.spine_neck_%i(1), 0)'%i)

		exec('spine_necks.append(h.spine_neck_' + str(i) + ')')
		exec('spine_heads.append(h.spine_head_' + str(i) + ')')

	for param in ['cm', 'e_pas', 'Ra', 'g_pas']:
		for spine in spine_heads+spine_necks:
			setattr(spine, param, eval(param))

	h.define_shape()
	return spine_heads, spine_necks













