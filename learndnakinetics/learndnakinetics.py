"""
Created by Nasim Zolaktaf, 2019
This file is used in map.py, it return the GMM to map.py and writes the MSE to file"""
from __future__ import division
import ConfigParser
import numpy as np
import csv
import cPickle as pickle
import timeit
import os
import multiprocessing
import sys
import math
import myenums
import copy
import string
sys.path.insert(0,os.path.realpath('../reactions'))
import parent
from parent import *
from datasets import *

DATASET_PATH = '../dataset'
iter = 0
totaliter= 0
thetalist = []
bestmin = np.inf

NUCLEOTIDES = "ACTG"
TRANSLATION_TABLE = string.maketrans(NUCLEOTIDES, "TGAC")
def Complement(str):
	return ''.join(list(reversed(str))).translate(TRANSLATION_TABLE)

"""A reaction object, general conditions are set here, reaction specific conditions are set in function read_dataset in this file and in datasets.py """
class Reaction(object):

	def __init__(self, **kwargs):
		self.simulation_time  =  100000. # Multistrand simulation time should always be a FLOAT !

		self.use_initialfinalfrompathway =  True # do not change
		for k, v in kwargs.items():
			setattr(self, k, v)

	"""setting other attributes from datasets.py"""
	def set_specific ( self,  **kwargs) :
		for k, v in kwargs.items():
			setattr(self, k, v)

	def __str__(self):
		output =""
		output  += " num_simulations: " + str(self.num_simulations)
		output  += " simulation_time: " + str(self.simulation_time)
		output  += " use_initialfinalfrompathway: " + str(self.use_initialfinalfrompathway)
		output  += " temperature: " + str(self.temperature) + " 1000/T:  " + str(1000/self.temperature)
		output  += " join_concentration: " + str(self.join_concentration)
		output  += " temperature_change: " + str(self.temperature_change)
		output  += " join_concentration_change: " + str(self.join_concentration_change)
		output  += " sodium: " + str(self.sodium)
		output  += " magnesium " + str(self.magnesium)
		output  += " real_rate  " + str(self. real_rate ) + " log10 k:  "  + str( self.real_rate )
		return output

"""reading reactions from database and setting their conditions """
def read_dataset( done_queue ,  dataset ) :
	temperature_column= 1
	rate_column =[ 4, 18]
	flur_position  =  ""
	toehold_length = ""
	sodium_column = 2
	magnesium_column = 3
	timescale_column = 5
	join_concentration_column = [6,7]
	strand_column = [8,9, 10, 11, 12, 13, 14, 15]
	try :
		temperature = float(dataset.row[dataset.counter_cell][temperature_column])
	except :
		raise ValueError("no temperature set in excel file" )
	if dataset.temperature_type ==myenums.DatasetSpecifications.TEMPKELVININV.value  :
		temperature = 1000/ temperature - 273.15
	elif dataset.temperature_type ==myenums.DatasetSpecifications.TEMPKELVIN.value  :
		temperature =  temperature- 273.15
	try :
		sodium = float(dataset.row[dataset.counter_cell][sodium_column])
	except:
		sodium = 0.01
	# sodium can't be less than 0.01 or error will be raised
	if sodium < 0.01 :
		sodium = 0.01
	try :
		magnesium = float(dataset.row[dataset.counter_cell][magnesium_column])
	except:
		magnesium  = 0.
	try :
		join_concentrationsize = 2
		join_concentration = np.max ( (float(dataset.row[dataset.counter_cell][join_concentration_column[0]]) , float(dataset.row[dataset.counter_cell][join_concentration_column[1]]) ))
	except :
		join_concentrationsize = 1
	if join_concentrationsize == 1:
		try :
			join_concentration =float(dataset.row[dataset.counter_cell][join_concentration_column[0]])
		except:
			join_concentration = 0.00000001
	try:
		real_rate = float(dataset.row [dataset.counter_cell][rate_column[0]])
		if dataset.rate_type ==myenums.DatasetSpecifications.LOG10RATECONSTANT.value :
			real_rate =  math.pow(10, real_rate)
		if dataset.rate_type == myenums.DatasetSpecifications.RATECONSTANT10POW5.value:
			real_rate = real_rate *  ( 10 **5 )
	except :
		timescale = float (dataset.row [dataset.counter_cell ][timescale_column])
		if dataset.rate_type == myenums.DatasetSpecifications.TIMESCALE.value  :
			if dataset.bimolecular_reaction == True :
				real_rate  = 1 /  ( timescale  * join_concentration)
			elif dataset.bimolecular_reaction == False  :
				real_rate  = 1 / timescale
	strands_list = []
	for i in strand_column:
		try :
			if dataset.row[dataset.counter_cell][i] != '' :
				strands_list.append(dataset.row[dataset.counter_cell][i].rstrip() )
		except   :
			"no more strings to read"
	if dataset.dataset_name == myenums.DatasetName.HATA.value  :
		strand1 = strands_list[0 ]
		strand2 = Complement(strand1)
		strands_list= [strand1 ,strand2]
	#strands_list  = ["GAAT", "ATTC"]
	dataset.set_specific(cutoff =1 )
	dataset.set_specific(  temperature =  temperature , join_concentration =join_concentration, sodium = sodium, magnesium =  magnesium,  real_rate = real_rate, strands_list = strands_list, toehold_length = toehold_length, flur_position = flur_position)
	complex = ParentComplex(dataset)
	output = parent.main( complex )
	done_queue.put( ( dataset.dataset_name , output.error  , dataset.counter_cell, dataset.document,    output.predicted_log_10_rate, output.real_log_10_rate)  )

"""These values are unacceptable for the Arrhenius and Metropoli smodels"""
def filter_undefined_parameter_set(theta_simulation , alpha, sigma):
	if (parent.rate_method ==Literals.arrhenius and alpha  <= 0)   or sigma <= 0  or (parent.rate_method ==Literals.metropolis and  ( theta_simulation[0] <= 0 or theta_simulation[1]  <= 0 )  ) :
		return True

""" We bound the rates to be between 10** 4 and 10 ** 9 as follows) """
def filter_parameter_set(theta_simulation , alpha ):
	global simulationTimeAVG
	simulationTimeAVG =  0
	mybool = False
	unacceptable_high_unimolecularrate = math.pow (10, 9)
	unacceptable_low_unimolecularrate = math.pow (10, 4)
	if parent.rate_method ==Literals.metropolis :
		simulationTimeAVG =   (1./theta_simulation  [0] )
		if theta_simulation [0 ] > unacceptable_high_unimolecularrate   or theta_simulation [0 ] < unacceptable_low_unimolecularrate or theta_simulation [1 ] > unacceptable_high_unimolecularrate   or theta_simulation [1] < unacceptable_low_unimolecularrate :
			mybool = True
			print "rate is unacceptable!!!!!! " , alpha

	R= 0.0019872036  # kcal / K mol
	RT = R * (23 + 273.15)
	if parent.rate_method == Literals.arrhenius :
		for birate in [1, theta_simulation[-1] ]  :
			n_local_contexts = 7
			count = 0
			for i in range(n_local_contexts):
				for j in range(i, n_local_contexts):
					lnA = theta_simulation[2* i   ] + theta_simulation [2* j  ]
					E = theta_simulation[2* i  +1 ] + theta_simulation [2* j  +1]
					rate = np.e ** (lnA - E / RT) * birate

					simulationTimeAVG +=   ( 1./rate)
					if rate > unacceptable_high_unimolecularrate   or rate < unacceptable_low_unimolecularrate  :

						print "  rate is unacceptable!!!!!" , rate
						mybool = True
					count+=1
			simulationTimeAVG = simulationTimeAVG / count
	if filter_smallandlarge_rates  ==  False  :
		mybool = False
	return mybool

def objective_function_auxilary(   done_queue, dataset_list, n ,dataset):
	for counter_cell in dataset.counter_celllist  :
		print "using " , dataset.document  ,  counter_cell
		datasetx = copy.deepcopy(dataset)
		docID = datasetx.docID + str(counter_cell )
		datasetx.set_specific (  docID=  docID, counter_cell = counter_cell )
		if use_multiprocess == True:
			dataset_list.append( ForMultiProcess( read_dataset, (done_queue, datasetx)))
		else :
			read_dataset(done_queue,  datasetx)
		n +=1
	return  (dataset_list, done_queue, n )

"""This function returns the error of the predictions on four datasets"""
def objective_function(theta_simulationpp):

	global iter , thetalist, totaliter, bestmin
	theta_simulationp = copy.deepcopy (theta_simulationpp )
	if  parent.rate_method == Literals.metropolis:
		for i in range(len(theta_simulationp)):
			theta_simulationp[i]= theta_simulationpp [i]
	elif parent.rate_method == Literals.arrhenius:
		for i in range(len(theta_simulationp)):
			if i %2 == 0 and i <= 12:
				theta_simulationp[i]= math.log(theta_simulationp[i])
	dumptheta  = True

	if iter ==0 :
		#in iter == 0 we  build and save the fixed  paths for FPEI
		load = False
		save = True

	else :
		#in iter >= 1  we load and reuse the fixed paths in FPEI
		load = True
		save = False



	start_time = timeit.default_timer()
	theta_simulation =[]

	#if (parent.rate_method == Literals.arrhenius and len(theta_simulationp)!= 16 ) or (parent.rate_method == Literals.metropolis and len(theta_simulationp)!=3):
	if (parent.rate_method == Literals.arrhenius and len(theta_simulationp)!= 15 ) or (parent.rate_method == Literals.metropolis and len(theta_simulationp)!=2):
		raise ValueError('[paramet set not initialized correctly!, probably forgot to set sigma to learn ')
	for x in theta_simulationp :
		theta_simulation.append(x)
	if parent.rate_method == Literals.arrhenius:
		theta_simulation = [theta_simulationp[0] , theta_simulationp[1] , theta_simulationp[2], theta_simulationp[3] , theta_simulationp[4] , theta_simulationp[5], theta_simulationp[6] ,theta_simulationp[7], theta_simulationp[8], theta_simulationp[9], theta_simulationp[10] , theta_simulationp[11],  theta_simulationp[12] , theta_simulationp[13], theta_simulationp[14] ]
		alpha = theta_simulation [14]
	elif parent.rate_method == Literals.metropolis :
		theta_simulation = [theta_simulationp[0] , theta_simulationp[1]]
		alpha =1
	else:
		raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis!')

	#sigma = theta_simulationp[len(theta_simulationp)-1]
	sigma = 1

	if filter_undefined_parameter_set(theta_simulation , alpha, sigma) == True or filter_parameter_set(theta_simulation, alpha) == True :
		return np.inf
	parameter_file = open(parameter_file_name, 'a')
	MSE_file = open(MSE_file_name, 'a')
	overalltime_file = open(overalltime_file_name, 'a')
	parameter_file.write("Iteration " + str(iter) +" "+str(theta_simulation) +  " " + str(sigma) + '\n')
	error = 0
	n = 0
	done_queue =  multiprocessing.Manager().Queue()
	dataset_list = []
	directories =[]

	#comment or uncomment the dataset you want to consider in parameter estimation.

	#Hata
	#(dataset_list, done_queue, n ) =  read_Hata2017 ( done_queue  =done_queue, dataset_list = dataset_list, n = n  ,  directories = directories, theta_simulation = theta_simulation ,    load  =load , save =save )

	#Cisse2012
	#(dataset_list, done_queue, n ) =  read_Cisse2012( done_queue  =done_queue, dataset_list = dataset_list, n = n  ,  directories = directories, theta_simulation = theta_simulation  ,  load  =load , save =save, )

	#Bonnet
	(dataset_list, done_queue, n ) =  read_Bonnet98 (done_queue  =done_queue,dataset_list =dataset_list, n = n,  directories = directories, theta_simulation = theta_simulation , load  =load , save =save)

	#FPEI currently does not work for symmetric reactions  with the Arrhenius model, eventhough it works for SSA  and RVSSA, because it does not produce correct neighbors for a state in the Builder file. For example for the state where  "AAA" and "TTT" do not have any base pairs it produces four neighbors, but it actually has 9 neighbors

	check_directories (directories)
	dataset_error = dict()
	dataset_count = dict()
	terror  = multi_process(done_queue , dataset_list   , iter , dataset_error , dataset_count  )
	error += terror


	elapsed = timeit.default_timer() - start_time
	thetalist.append(theta_simulation )

	if dumptheta == True :
		with open (thetalistpicklefile, "wb" ) as p:
			pickle.dump( thetalist, p)

	#writing information to file
	parameter_file.write( "Iteration:" + str(iter) +  "		 SE:" + str( error) +   "		iteration time:" + str(format(elapsed, "10.4E")) +  "     Overalltime: " + str( format(timeit.default_timer() - OVERALL_STARTTIME ,"10.4E"))+ "      Size of dataset: "+str(n)+ '\n')
	parameter_file.write( "Iteration:" + str(iter) +  "		 MSE:" + str(error/n) + "    Average Overalltime: " + str(format(  (timeit.default_timer() - OVERALL_STARTTIME)/ (iter + 1) , "10.4E" ))+ "        Size of dataset: "+str(n)+ '\n')

	#This file is used to generate Fig 4 and 5 in the paper
	MSE_file.write(str(error/n) + ",")
	overalltime_file.write(str(timeit.default_timer() - OVERALL_STARTTIME) + ",")
	if iter  ==0 :
		MSExnew_file = open(MSExnew_file_name, 'a')
		MSEynew_file = open(MSEynew_file_name, 'a')
		MSExnew_file.write(str(totaliter ) + ",")
		MSEynew_file.write(str(error/n) + ",")
		MSExnew_file.close()
		MSEynew_file.close()
	parameter_file.write( "SE:  ")
	for cs in dataset_error:
		parameter_file.write(cs + ": " + str(dataset_error[cs])  + ",    ")
	parameter_file.write( "\n")
	parameter_file.write( "MSE:  ")
	for cs in dataset_error:
		parameter_file.write(cs + ": " + str(dataset_error[cs]/ dataset_count[cs])  + ",    ")
	parameter_file.write( "\n\n")
	parameter_file.close()

	print  ( "Iteration:" + str(iter) +  "		,error:" + str(error) + "		iteration time:" + str(elapsed) +  '\n' )
	iter += 1   # Do not move this line or you'll get errors later
	totaliter +=1

	if np.isnan(error)  or error == np.inf:
		error  = np.inf

	if  error  < bestmin :
		bestmin =  error
		bestlnawise = theta_simulation
		parameter_file = open(parameter_file_name, 'a')
		parameter_file.write ("\n best parameter set so far\n\n")
		parameter_file.write(str(bestlnawise) + "\n")
		parameter_file.close()
	return error

"""This class used for multiprocessing"""
class ForMultiProcess(object):
	def __init__(self, function_name, arguments) :
		self.function_name = function_name
		self.arguments = arguments

"""multi processing function """
def multi_process(done_queue ,  dataset_list ,  iter , dataset_error ,  dataset_count) :

	error = 0
	if use_multiprocess== True  :
		pool = multiprocessing.Pool( processes = n_processors)
		for ds in dataset_list:
			compile_error = pool.apply_async( ds.function_name ,  ds.arguments )
		print  ( "Errors: " + str(compile_error.get()) )
		pool.close( )
		pool.join ()
	while not done_queue.empty():
		(dataset_name, s  ,  counter_cell, document  , predicted_log_10_rate ,real_log_10_rate   ) = done_queue.get()

		error += s
		if dataset_name in dataset_error :
			dataset_error [dataset_name ] += s
			dataset_count[dataset_name] += 1
		else :
			dataset_error[dataset_name ] = s
			dataset_count[dataset_name] = 1
	return error

def open_document(document) :
	"""open a csv file"""
	my_CSV = list(csv.reader(open(document, 'rb')))
	return my_CSV

def check_directories (directories) :
	for dir in directories:
		if not os.path.exists(dir):
			os.makedirs(dir)


"""creating required directories"""
def initconf(my_name , directories ) :
	dataset_path =  parameter_folderround +my_name
	document = DATASET_PATH + my_name + '.csv'
	directories +=  [ dataset_path + "/" +  myenums.Permanent_Folder.MYBUILDER.value , dataset_path + "/"+  myenums.Permanent_Folder.TRAJECTORY.value  ]
	if not os.path.exists(dataset_path):
		os.makedirs(dataset_path)
	row =  open_document(document)
	return dataset_path, document , row

"""creating directory for round mew round of of inference"""
def set_folderfiles() :
	global parameter_folderround
	parameter_folderround = parameter_folder+ "/round" +str(map_i)#map_i is initialized in map.py
	if not os.path.exists(parameter_folderround):
		os.makedirs(parameter_folderround)

"""reading config_file.txt and setting filenames"""
def set_configuration():
	#setting some configurations!
	global MSEynew_file_name, MSExnew_file_name ,  thetalistpicklefile, MSE_file_name, overalltime_file_name,  parameter_file_name,  parameter_folder ,   n_processors, use_multiprocess, filter_smallandlarge_rates,  OVERALL_STARTTIME
	OVERALL_STARTTIME = timeit.default_timer()
	configParser = ConfigParser.ConfigParser()
	configParser.readfp(open(r'config_file.txt'))
	CONFIG_NAME = 'learndnakinetics'
	parent.rate_method= configParser.getint(CONFIG_NAME, 'rate_method')
	print "Literals.rate_method is set to ", parent.rate_method

	n_processors =    configParser.getint(CONFIG_NAME, 'n_processors')
	use_multiprocess =  bool(configParser.getint(CONFIG_NAME, 'use_multiprocess') )
	filter_smallandlarge_rates =  bool(configParser.getint(CONFIG_NAME, 'filter_smallandlarge_rates') )
	parameter_folder = configParser.get(CONFIG_NAME, 'parameter_folder')
	if not os.path.exists(parameter_folder):
		os.makedirs(parameter_folder)
	MSE_file_name= parameter_folder +"/" + myenums.ImportantFiles.MSE.value
	MSExnew_file_name= parameter_folder + "/"  + myenums.ImportantFiles.MSEXNEW.value
	MSEynew_file_name= parameter_folder + "/" + myenums.ImportantFiles.MSEYNEW.value
	overalltime_file_name= parameter_folder + "/"+ myenums.ImportantFiles.OVERALLTIME.value
	thetalistpicklefile= parameter_folder + "/thetalistpicklefile"
	parameter_file_name  = parameter_folder + "/parameterlog"
	parameter_file_name+= ".txt"