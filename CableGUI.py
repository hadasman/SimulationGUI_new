import os, pdb, matplotlib
os.chdir('../function_scripts')
from neuron import h, gui

try: # For python3
	from tkinter import *
	from tkinter import ttk, messagebox
	from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg
except: # For python2
	from Tkinter import *
	import ttk, messagebox
	from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# In main script: from GUI import * - will import tkinter when class is imported
class CableGUI(Tk):

	def __init__(self, title_, dimensions, background_color):

		super(CableGUI, self).__init__()

		# Design GUI layout
		self.title(title_)
		self.geometry("%sx%s"%(dimensions[0], dimensions[1]))
		self.config(background=background_color)

		# Initialize dictionaries
		self.Buttons 		= {}
		self.Tabs 			= {}
		self.RadioButtons 	= {}
		self.Entries 		= {}
		self.Figures 		= {}
		self.ChangingLabels = {}
		self.CheckBoxes		= {}

		self.suppress_entry_callback = False

	def AddTab(self, tab_var_name, tab_text, rowspan=1, columnspan=1, columnConfigure=None):
		'''
		Add new tab to main window. If parent tab doesn't exist in self.parent_tab, the function creates it.
		Inputs:
			- tab_var_name: (string) name of the tab, as will be written in its key inside the self.Tab dictionary.
			- tab_text: (string) text on the tab.
			- rowspan, columnspan, sticky: see .grid() in tkinter.
			- columnConfigure: Tuples indicating arguments for columnconfigure method of tkinter: [(column1, weight1), (column2, weight2)]
		'''

		if not hasattr(self, 'parent_tab'):
			self.parent_tab = ttk.Notebook(self)
			self.parent_tab.grid(rowspan=rowspan, columnspan=columnspan)

		temp_tab = ttk.Frame(self.parent_tab)
		self.parent_tab.add(temp_tab, text = tab_text)
		
		if columnConfigure:
			for column, weight in columnConfigure:
				self.Tabs[tab_var_name].columnconfigure(column, weight=weight)

		self.Tabs[tab_var_name] = temp_tab

	def AddEntries(self, widget_loc, name, param_dict, command=None, entryTitle=None, labelColumn=0, labelRow=None, columnspan=None):
		'''
		Add an entry widget and save in to self.Entries dictionary. 
		Left of the entry will be a label describing its content.
		Inputs:
			- widget_loc: widget location (main window, tab).
			- name: (string) the key name under which the entry will be saved in self.Entries.
			- param_dict: parameter dictionary that includes entry description (string) and default value (string or int)
						  as values, and the variable to be changed by the entry (as string) as key.
						  Example: param_dict = {
					  							'dend.L': ['Branch Length', 80], 
						  						'Simulator.exc_dist': ['Exc. Distribution', 'Uniform']
						  						}
			- entryTitle: if exists, this will be the title above this entry (used for beginning of new entry block)
			- labelColumn, labelRow: row and column of the label describing entry content.
			- columnspan: column span of the entry widget.
		'''

		values = list(param_dict.values())
		labels, defaults = [i[0] for i in values], [i[1] for i in values]

		if entryTitle:
			self.AddLabel(widget_loc, entryTitle, 
						font='Helvetica 15 bold', foreground='darkred', column=labelColumn, padx=5, pady=5, sticky='W')

		self.Entries[name] = {key: {} for key in param_dict}

		for label, default, key in zip(labels, defaults, param_dict.keys()):
			L = Label(widget_loc, text=label, fg='blue')
			L.grid(row=labelRow, column=labelColumn, sticky=EW)
			if labelRow:
				labelRow += 1
			textVar = StringVar()

			temp_entry = ttk.Entry(widget_loc, textvariable=textVar)
			temp_entry.grid(row=L.grid_info()['row'], column=labelColumn+1, columnspan=columnspan)
			textVar.set(default)
			textVar.trace_variable("w", command)

			self.Entries[name][key]['Var'] = textVar
			self.Entries[name][key]['Entry'] = temp_entry

	def AddButton(self, widget_loc, var_name, text, command=None, row=None, column=None):
		'''
		Add a button widget and save in to self.Buttons dictionary.
		Inputs:
			- widget_loc: widget location (main window, tab).
			- var_name: (string) the key name under which the button will be saved in self.Buttons.
			- text: text to be presented on the button.
			- command: the function instance to be called upon pressing the button.
			- row, column: see .grid() in tkinter.
		'''
		
		buttonObj = Button(widget_loc, text=text, command=command, activebackground='gray', padx=20)
		buttonObj.grid(sticky=W, column=column, row=row)

		self.Buttons[var_name] = buttonObj

	def AddLabel(self, widget_loc, text=None, varName=None, font=None, fg=None, foreground=None, background=None, sticky=None, row=None, column=None, columnspan=None, rowspan=None, padx=None, pady=None):
		'''
		Add a label widget. If label text is not constant, i.e. controlled by a string variable, 
		save label pounter in self.ChangingLabels dictionary.
		One of the arguments 'text' or 'varName' must be given, or else an exception will be raised. 
		If both are given, the default will be to present 'text' as a constant text in label.
		Inputs:
			- widget_loc: widget location (main window, tab).
			- text: (string) constant text to (optional).
			- varName: (string) the key name under which the label and associated StringVar will 
								be saved in self.ChangingLabels. If given, and no text is given, 
								the function will create a string variable, associate it with the 
								label and save them both in self.ChangingLabels dictionary.
								Example: self.ChangingLabels[varName] = {'Label': label_widget, 'Variable': StringVar}
			- command: the function instance to be called upon pressing the button.
			- row, column, columnspan, rowspan, sticky, padx, pady: see .grid() in tkinter.
			- font, fg, foreground, background: see Label() properties in tkinter.
		'''
	
		if text or varName:
			if text:
				L = Label(widget_loc, text=text, font=font, foreground=foreground, background=background, fg=fg)
			else:
				var = StringVar()
				L = Label(widget_loc, textvariable=var, font=font, foreground=foreground, background=background, fg=fg)
				self.ChangingLabels[varName] = {'Label': L, 'Var': var}
		else:
			raise Exception('No text argument for label!')
		
		L.grid(row=row, column=column, rowspan=rowspan, columnspan=columnspan, sticky=sticky, padx=padx, pady=pady)
			
	def AddRadioButton(self, widget_loc, name, all_values, command=None, default_val=1, row=None, column=None, sticky=None):
		'''
		Add a radio button widget and save in self.RadioButtons dictionary, along with its associated variable.
		Inputs:
			- widget_loc: widget location (main window, tab).
			- name: (string) the key name under which the radio button will be saved in self.RadioButtons.
			- all_values: list of texts to be presented on each radio button, with respect to order.
			- command: the function instance to be called upon selecting one of the buttons.
			- default_val: default value at GUI initiation.
			- row, column, sticky: see .grid() in tkinter.
		'''

		RadioVar = IntVar()
		RadioVar.set(default_val)		

		self.RadioButtons[name] = {'Var': RadioVar, 'Buttons': []}
		for i, text in enumerate(all_values):
			temp_button = Radiobutton(widget_loc, text=text, variable=RadioVar, value=i+1, indicatoron=1)
			temp_button.grid(row=row, column=column, sticky=sticky)
			self.RadioButtons[name]['Buttons'].append(temp_button)		

		_ = RadioVar.trace_variable("w", command) # Keep at end of function

	def AddCheckbox(self, widget_loc, text, varName, default, row=None, column=None, sticky=None):
		
		var = IntVar()
		Checkbutton(widget_loc, text=text, variable=var).grid(row=row, column=column, sticky=sticky)
		var.set(default)
		self.CheckBoxes[varName] = {'Var': var}

	def DrawSections(self, colors):
		
		if 'morph' not in self.Figures:
			raise Exception('Create morphology axis first!')
		else:
			morph_ax = self.Figures['morph']['ax']
			morph_graph = self.Figures['morph']['graph']

		morph_ax.clear()# morph_graph.draw()

		for sec in h.allsec():
			# Get section name
			key = [k for k in colors.keys() if k in sec.name()][0]
			# Get section color
			color = colors[key]

			# To generalize soma and spine_head, use: x * sum(np.diff([sec.x3d(i) for i in range(sec.n3d())])) +... to figure out which dimension doesn't change
			# statement = 'sec.%s3d(0) - sec.%s3d(sec.n3d()-1)'%x
			# dim_vec = [(sec.x3d(0)- sec.x3d(sec.n3d()-1))!=0, (sec.y3d(0)- sec.y3d(sec.n3d()-1))!=0, (sec.z3d(0)- sec.z3d(sec.n3d()-1))!=0]

			# Plot lines as line, and compartments with shape using fill_between
			if key == 'soma':
				x = [sec.x3d(i) for i in range(sec.n3d())]
				center_ = x[-1] - (x[-1] - x[0])/2
				morph_ax.fill_between(x, y1=x[0]-center_, y2=x[-1]-center_, color=color)
			
			elif key == 'spine_head':
				y = [sec.y3d(i) for i in range(sec.n3d())]
				center_ = y[-1] - (y[-1] - y[0])/2
				x = [sec.x3d(i)+sec.y3d(i)-center_ for i in range(sec.n3d())]

				morph_ax.fill_between(x, y1=y[0], y2=y[-1], color=color)

			else:
				morph_ax.plot([sec.x3d(i) for i in range(sec.n3d())], [sec.y3d(i) for i in range(sec.n3d())], color=color)

		xlim = morph_ax.get_xlim()
		center = xlim[1] - (xlim[1] - xlim[0])/2
		cut = 20
		morph_ax.set_ylim(xlim[0] - center + cut, xlim[1] - center - cut)
		morph_graph.draw()

	def PutFigure(self, widget_loc, fig_name, xlabel='', ylabel='', facecolor = 'white', row=None, column=None, rowspan=None, columnspan=None, padx=0, pady=0):

		fig, ax = plt.subplots()
		fig.set_figheight(2.4)
		fig.set_figwidth(3.2)
		plt.close(fig)

		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.set_facecolor(facecolor)

		graph = FigureCanvasTkAgg(fig, master=widget_loc)
		graph.get_tk_widget().grid(row=row, column=column, rowspan=rowspan, columnspan=columnspan, padx=padx, pady=pady)

		self.Figures[fig_name] = {'ax': ax, 'graph': graph}

	def PopMessageBox(self, message):
		messagebox.showinfo("Error!", message)


