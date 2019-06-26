"""Created by nasim zolaktaf, 2019
This file plots the MSE vs iteration,  (such as fig 4 and fig 5 of the paper).
 To run this file, first run map.py to do parameter estimation once with FPEI, once with SSAI and to generate neccessary files.
 Then in this file set  dc1 and dc2  to correct flies.
 Run 'plot_ssavsFPEI.py'
"""
from __future__ import division, absolute_import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import seaborn as sns
import csv
__all__ = ['minimize', 'minimize_scalar']
import os
import sys
sys.path.insert(0,os.path.realpath('../learndnakinetics'))
import myenums
sns.set()
sns.set_context(rc={"figure.figsize": (8, 4)})

"""Setting options for different plots based on the plots in the literature. """
class PlotOptions(object):

	def __init__(self,  figuresavename, xlabel, ylabel, names, documents  , yax_min , yax_max, title, time):
		self.figuresavename= figuresavename
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.names= names
		self.documents=  documents
		self.yax_min = yax_min
		self.yax_max =  yax_max
		self.title= title
		self.time = time

"""open a csv file"""
def open_document(document) :

	my_CSV = list(csv.reader(open(document, 'rb')))
	return my_CSV

"""Drawing the MSE of FPEI VS SSA"""
def draw_plot(po):

	title_fontsize = 36
	label_fontsize = 28
	tick_fontsize=26
	propsize = 26

	colors=['red','blue', 'green']
	linewidth = 1
	loc =1
	smarkersize= [20,20,20,20]
	FPEIstarcolor= "g"
	FPEIstarstyle= "*"
	linestyle = [ '-' , '--'  , '--']


	fig, ax = plt.subplots()
	ymin = np.inf
	ymax= -np.inf
	lengths= []
	for document, name  in zip( po.documents , po.names) :
		if name =="":
			lengths.append(0)
		else:
			row =  open_document(document)
			row= row[0]
			lengths.append(len(row))

	for document, name , i in zip( po.documents , po.names , range(len(po.names))) :
		if name =="" :
			continue
		row =  open_document(document)
		row = row[0]
		x = range(len(row))
		y = []
		for cc in  x :

			y.append(row[cc])
		if y[-1] =="" :

			x = x[:-1]
			y = y [:-1]
		y = map(float, y)
		ymint= min(y)

		if ymint  < ymin :
			ymin = ymint
		ymaxt=  max(y)
		if ymaxt > ymax and ymaxt!= np.inf:
			ymax = ymaxt

		if myenums.NamesinPlot.SSAI.value in name:
			ind =  0
		elif myenums.NamesinPlot.FPEI.value in name:
			ind= 1
		plt.plot(x, y,  color= colors[i ], linewidth = linewidth, linestyle=  linestyle[i], label = name +" (" +po.time[ind] + " s / iter")

		#########reading only values when new paths are generated
		if  usexnewandynew  == True and   myenums.NamesinPlot.FPEI.value  in name:

			dl = document.split("/")
			documenty  = ""
			print dl
			for i  in range (len(dl)-1):
				print dl[i]
				documenty +=  dl[i]
				documenty += "/"
			documenty =documenty +   myenums.ImportantFiles.MSEYNEW.value

			row =  open_document(documenty )
			row = row[0]
			ynew = []
			lengths = range( len(row))
			for cc in  lengths :
				ynew.append(row[cc])
			dl = document.split("/")
			documentx  = ""
			for i  in range (len(dl)-1):
				documentx += dl[i]
				documentx += "/"
			documentx =documentx +  myenums.ImportantFiles.MSEXNEW.value

			row =  open_document(documentx )
			row = row[0]
			xnew = []
			lengths = range( len(row))
			for cc in  lengths :
				xnew.append(row[cc])
			if ynew[-1] =="" :

				xnew = xnew[:-1]
				ynew = ynew [:-1]

			ynew = map(float, ynew)
			xnew = map(int, xnew)

			xnew= xnew[:8]
			ynew= ynew[:8]
			plt.plot(xnew, ynew,  color= FPEIstarcolor, marker= FPEIstarstyle, markersize=smarkersize[i], linewidth = linewidth, linestyle=  linestyle[1])


		plt.ylim(bottom = 0 )

		if "Modified"in po.figuresavename:
			plt.ylim(top = 7 )
		else:
			plt.ylim(top = 1.6 )

	ttl = ax.title
	ttl.set_position([.5, 1.02])
	plt.title(po.title, fontsize=title_fontsize)
	if po.yax_min!= None  :
		plt.ylim(bottom = po.yax_min)
	if po.yax_max!= None  :
		plt.ylim(top = po.yax_max)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()

	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(tick_fontsize)
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(tick_fontsize)
	plt.legend(loc=loc, borderaxespad=0., prop={'size': propsize})
	plt.xlabel(po.xlabel, fontsize=label_fontsize)
	plt.ylabel(po.ylabel, fontsize=label_fontsize)
	plt.savefig(po.figuresavename, bbox_inches='tight')
	plt.show()
	plt.close(fig)


def main():

	global usexnewandynew 
	usexnewandynew = True

	dc1=   ["../learndnakinetics/testcisse-ssa"] #dc1 should contain  path to the SSAI folder results, same as parameter_folder in config_file.txt
	dc2=["../learndnakinetics/testcisse-fpei"] # This should contain path to the FPEI folder results, same as parameter_folder in config_file.txt
	times  = [ [r"$-$",r"$-$"]]
	for a, b , time  in zip(dc1, dc2 ,times  ) :
		if  "Modified" in a:
			title =  "Initialization: "+r"$\theta_0''$"

		else:
			title =  "Initialization: "+r"$\theta_0'$"

		documents1 = [ a, b]
		names = [ myenums.NamesinPlot.SSAI.value , myenums.NamesinPlot.FPEI.value] # do not change these names

		path = "plotmse/"
		if not os.path.exists(path):
			os.makedirs(path)
		documents2= documents1

		# These files are created from running map.py  The files are created and written to in learndnakinetics.py
		documents1 = [doc+"/" +  myenums.ImportantFiles.MSE.value for doc in documents1]
		documents2 = [doc+"/"+   myenums.ImportantFiles.OVERALLTIME.value for doc in documents2]


		xlabel = "Iterations"
		ylabel1 = "MSE"
		ylabel2 = "Overall time (s)"
		figuresavename1 = path+ "msevsiteration.pdf"
		figuresavename2 = path+ "totaltimevsiteration.pdf"

		plotoptions = PlotOptions( figuresavename= figuresavename1, xlabel = xlabel, ylabel = ylabel1 , names= names, documents = documents1,  yax_min = None  , yax_max= None, title= title, time= time )

		draw_plot(plotoptions)

		plotoptions = PlotOptions( figuresavename= figuresavename2, xlabel = xlabel, ylabel = ylabel2, names= names, documents = documents2,yax_min = None  , yax_max= None , title = title , time= time )
		draw_plot(plotoptions)

if __name__ == "__main__":
	main()
