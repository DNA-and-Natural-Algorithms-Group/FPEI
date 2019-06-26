"""
created by Nasim Zolaktaf, 2019
This file estimates the reaction rate constant of a reactions and returns the squared error between the predicted log reaction rate cosntant and experimental log reaction rate constant
"""
from __future__ import division
import warnings
import gc
import numpy as np
import ConfigParser
import math
import cPickle as pickle
from  scipy.sparse.linalg import *
from subprocess import Popen, PIPE , call
import copy
import timeit
import myenums
import learndnakinetics
from multistrand.options import Options, Literals
from multistrand.builder import Builder, BuilderRate ,  transitiontype,  localtype
from multistrand.concurrent import  FirstStepRate, FirstPassageRate, Bootstrap, MergeSim
from multistrand.system import SimSystem
from builderChild import BuilderRate_FPEI
from paths import  PathHelix, PathHairpin
from multistrand.objects import  StopCondition
import os
configParser = ConfigParser.ConfigParser()
configParser.readfp(open(r'../learndnakinetics/config_file.txt'))
CONFIG_NAME = 'parent'
use_FPEI_MFPT= bool(configParser.getint(CONFIG_NAME, 'use_FPEI_MFPT'))
use_Gillespie_MFPT = bool(configParser.getint(CONFIG_NAME, 'use_Gillespie_MFPT'))

RETURN_MINUS_INF = None
LNA_INDEX= 0
E_INDEX = 1
bimolecular_scaling = "bimolecular_scaling"
unimolecular_scaling = "unimolecular_scaling"

class Output(object):
	#output object returned to learndnakinetics.py
	def __init__(self, **kwargs):
		for k, v in kwargs.items():
			setattr(self, k, v)

	def set_specific ( self,  **kwargs) :
		for k, v in kwargs.items():
			setattr(self, k, v)

class ParentComplex(object):
	#
	def __init__(self , dataset) :
		if rate_method == Literals.arrhenius :
			theta = dataset.theta_simulation
			self.kinetic_parameters_simulation = {localtype.stack: (theta[0] ,theta[1]) ,
			                                  localtype.loop: (theta[2] ,theta[3]),
			                                  localtype.end: (theta[4] ,theta[5]),
			                                  localtype.stackloop: (theta[6] ,theta[7]),
			                                  localtype.stackend: (theta[8] ,theta[9]),
			                                  localtype.loopend: (theta[10] ,theta[11]),
			                                  localtype.stackstack:  (theta[12] ,theta[13]),
			                                  bimolecular_scaling : (theta[14]) }
		elif rate_method == Literals.metropolis :
			theta = dataset.theta_simulation
			self.kinetic_parameters_simulation ={ unimolecular_scaling :theta[0] , bimolecular_scaling :theta[1] }   #used in matrix computations

		else:
			raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis in the configuration file!')
		self.strands_list = dataset.strands_list
		self.use_initialfinalfrompathway= dataset.use_initialfinalfrompathway
		self.dataset_type = dataset.dataset_type
		self.reaction_type =dataset.reaction_type
		self.dataset_name = dataset.dataset_name
		self.startStates = None
		self.load = dataset.load
		self.save =  dataset.save
		self.real_rate = dataset.real_rate
		self.reaction_type  =dataset.reaction_type
		self.bimolecular_reaction = dataset.bimolecular_reaction
		self.dataset_path = dataset.dataset_path
		self.docID = dataset.docID
		self.temperature = dataset.temperature
		self.join_concentration = float(dataset.join_concentration ) # concentration has to be a float  (join_concentration in  Multistrand has to be set to a float)
		self.sodium = dataset.sodium
		self.magnesium = dataset.magnesium
		self.num_simulations= dataset.num_simulations
		self.simulation_time = dataset.simulation_time
		self.join_concentration_change =1 # Do not change this
		self.temperature_change = 0 #Do not change this

		self.loglikelihoods= dict()
		self.MFPTs= dict()
		self.original_loglikelihoods=  dict( )
		if  self.reaction_type == myenums.ReactionType.HELIXDISSOCIATION.value :
			pathway= PathHelix( False , self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type,cutoff = dataset.cutoff )
		elif  self.reaction_type == myenums.ReactionType.HELIXASSOCIATION.value :
			pathway = PathHelix( True, self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type,cutoff =dataset.cutoff )
		elif  self.reaction_type == myenums.ReactionType.HAIRPINCLOSING.value :
			pathway = PathHairpin( True , self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type ,cutoff =dataset.cutoff)
		elif self.reaction_type == myenums.ReactionType.HAIRPINOPENING.value :
			pathway = PathHairpin( False , self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type , cutoff =dataset.cutoff)

		if  self.use_initialfinalfrompathway == True :
			self.startStates =  pathway.get_intial_final()

	def find_meanfirstpassagetime(self):
		attributes_file  = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt"
		attributes_file = open(attributes_file, 'a')
		attributes_file.write("num_simulations: " +str(self.num_simulations) +  "    simulations_time: " + str(self.simulation_time)   +"       concentration_change: "  + str(self.join_concentration_change) + "      concentration: " + str(self.join_concentration )+   "      temperature: " + str(self.temperature )+  "      temprature_change: " + str(self.temperature_change)+ "   sodium: " +str( self.sodium)+  "    magnesium: "+ str(self.magnesium) +"\n")
		attributes_file.close()
		if 	 use_FPEI_MFPT == True :
			return self.find_meanfirstpassagetime_FPEI(  )
		elif use_Gillespie_MFPT == True :
			return self.find_meanfirstpassagetime_Gillespie()

	def print_path(self,o):
		print o.full_path[0][0][3]   # the strand sequence  #if you get  list index out of range,  set this options.output_interval = 1
		print o.start_state[0].structure   # the starting structure
		for i in range(len(o.full_path)):
			time = o.full_path_times[i]
			state = o.full_path[i][0]
			struct = state[4]
			sequence = state[3]
			dG = state[5]
			print struct + ' t=%11.9f seconds, dG=%6.2f kcal/mol' % (time, dG)


	"""Estimating the Mean First Passage Time with FPEI """
	def find_meanfirstpassagetime_Gillespie(self) :
		startime = timeit.default_timer()
		options = doReaction([self.num_simulations,self.simulation_time, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction, self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, rate_method , self.use_initialfinalfrompathway, self.startStates])
		#options.output_interval = 1 #to store all the transitions!!!!!!!!!!!!!!!!!!
		s = SimSystem(options )
		s.start()
		myRates = FirstPassageRate( options.interface.results)
		del s
		finishtime = timeit.default_timer()
		#print "average sampling time is " , str ( (finishtime -startime ) / self.num_simulations )
		mfpt = myRates.k1()
		del myRates
		gc.collect()
		return  mfpt


	def call_Builder(self, num_simulations ) :
		myBuilder = Builder(doReaction, [ num_simulations, self.simulation_time, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction,self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, rate_method , self.use_initialfinalfrompathway, self.startStates])
		start_time = timeit.default_timer()
		if self.num_simulations > 0 :
			myBuilder.genAndSavePathsFile(supplyInitialState= self.startStates[0])
			builderpath = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+"/"+myenums.Permanent_Folder.MYBUILDER.value+str(self.docID)
			lenp = len( myBuilder.protoSpace)
			start  = 0
			myBuilder.protoSpacebackup= copy.deepcopy(myBuilder.protoSpace)
			pathuniqueid = builderpath + myenums.EDITIONAL.UNIQUEID.value
			with open(pathuniqueid   , "wb" ) as p :
				pickle.dump(myBuilder.uniqueID_number,  p )
			pathspace=   builderpath + myenums.EDITIONAL.PROTOSPACEBACKUP.value
			with open(pathspace   , "wb" ) as p :
				pickle.dump(myBuilder.protoSpacebackup,  p )
			pathsequences= builderpath + myenums.EDITIONAL.PROTOSEQUENCES.value
			with open(pathsequences , "wb" ) as p :
				pickle.dump(myBuilder.protoSequences,  p )
			pathoptions = builderpath + myenums.EDITIONAL.PATHOPTIONS.value
			with open(pathoptions , "wb" ) as p :
				pickle.dump(myBuilder.optionsArgs,  p )
			batchsize = 2000
			while start < lenp :
				st =  timeit.default_timer ( )
				end = min(lenp ,  start+batchsize)
				print "progress " ,  str(end)  , " / ", str(lenp) ,  self.docID

				# There was some memory leak issues when I used the fattenStateSpace function in builder.py, so  I added fatten helper to avoid the memory issues by saving intermediate results, restarting Multistrand, and restoring the intermediate results
				command = ["python", "fattenhelper.py"  , str(start), str(end) , builderpath, pathspace, pathsequences, pathoptions, pathuniqueid]
				shell = call(command )
				ft  = timeit.default_timer()
				#print "making fatten state space time" , ft-st
				del shell
				with open( builderpath + myenums.EDITIONAL.TEMPSTATESPACE.value + str(start)+"-"  +str(end) , "rb" ) as p:
					tempstatespace  = pickle.load (  p)
				with open( builderpath + myenums.EDITIONAL.TEMPTRANSITIONS.value + str(start)+"-"  +str(end) , "rb" ) as p:
					temptransitions = pickle.load(  p)
				os.remove(builderpath + myenums.EDITIONAL.TEMPSTATESPACE.value + str(start)+"-"  +str(end))
				os.remove(builderpath + myenums.EDITIONAL.TEMPTRANSITIONS.value + str(start)+"-"  +str(end))
				myBuilder.mergeSet(myBuilder.protoSpace,  tempstatespace)
				myBuilder.mergeSet(myBuilder.protoTransitions_FPEI, temptransitions )
				start  = end
				with open(pathuniqueid   , "rb" ) as p :
					myBuilder.uniqueID_number  = pickle.load(  p )
			os.remove(pathuniqueid)
			os.remove(pathspace)
			os.remove(pathsequences)
			os.remove(pathoptions)
			del myBuilder.protoSpacebackup
			del myBuilder.uniqueID_number
			print "Statistics: " ,  "statespace: ", len(myBuilder.protoSpace), "finalstates: ",   len(myBuilder.protoFinalStates), "initialstates: ", len(myBuilder.protoInitialStates)

		return myBuilder

	def get_builder(self , ignoreSavedBuilder =  False , ignoreNumSimulations = False  , picklepath = ""  ):

		if ignoreNumSimulations == True :
			num_simulations = 1
		else:
			num_simulations  = self.num_simulations
		attributes_file  = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt"
		attributes_file = open(attributes_file, 'a')
		if self.load == False or ignoreSavedBuilder == True  :
			myBuilder = self.call_Builder( num_simulations )
			if  len(myBuilder.protoFinalStates) == 0 or len(myBuilder.protoInitialStates) == 0 :
				raise Exception('"no final states found! multiply simulation_time by 1000. in except block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1"')
		else  :
			start_time = timeit.default_timer()
			attributes_file.write( "Starting to load Builder   ")
			#print "starting to load builder"
			with open(picklepath, "rb") as p:
				myBuilder= pickle.load(p)
			attributes_file.write("finished loading Builder took :  " +str(timeit.default_timer()  -start_time) + "\n")
			if rate_method == Literals.metropolis:
				set_Metropolis_params(myBuilder.options, self.kinetic_parameters_simulation)
			elif  rate_method == Literals.arrhenius:
				set_Arrhenius_params(myBuilder.options,self.kinetic_parameters_simulation)
			else:
				raise ValueError(' Parameter method not set correctly ' )

		return myBuilder

	"""Estimating the Mean First Passage Time with FPEI """
	def find_meanfirstpassagetime_FPEI(self  ) :

		attributes_file  = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt"
		attributes_file = open(attributes_file, 'a')
		attributes_file.write("use_string: " + str(self.use_initialfinalfrompathway)+   "     num_simulations: " +str(self.num_simulations) +  "    simulations_time: " + str(self.simulation_time)   + "      concentration: " + str(self.join_concentration )+   "      temperature: " + str(self.temperature )+ "       concentration_change: "  + str(self.join_concentration_change) + "      temprature_change: " + str(self.temperature_change)+ "   sodium: " +str( self.sodium)+  "    magnesium: "+ str(self.magnesium))
		successful_simulation = False
		while successful_simulation == False :
			random_time =  False
			trial_count =  self.num_simulations
			switch =  1
			if learndnakinetics.iter >= switch:
				original = False
			else :
				original = True
			self.original = original
			path_file =  self.dataset_path+"/"+myenums.Permanent_Folder.TRAJECTORY.value+"/"+myenums.Permanent_Folder.TRAJECTORY.value+str(self.docID)

			def building_paths() :
				start= 1
				num_simulations_temp =start
				builderpath = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+"/"+myenums.Permanent_Folder.MYBUILDER.value+str(self.docID)
				#solutions_builderRate = 0
				while num_simulations_temp <= self.num_simulations :
					myBuilder  =self.get_builder(ignoreSavedBuilder = True  , ignoreNumSimulations = True)
					builderRate_FPEI= BuilderRate_FPEI(myBuilder, trial_count  =trial_count, random_time = random_time, stopStep = 100000 ,path_file =   path_file + "-" +  str( num_simulations_temp   -1 )  + "-" , original= original,  num_simulations_temp = num_simulations_temp, parameter_folder  = learndnakinetics.parameter_folder  )
					num_simulations_temp +=1
					loglikelihood1 , MFPT1, original_loglikelihood1= builderRate_FPEI.get_essentials ( )
					self.loglikelihoods  [ num_simulations_temp   -2]   = loglikelihood1
					self.MFPTs [num_simulations_temp -2] = MFPT1
					self.original_loglikelihoods[ num_simulations_temp -2] = original_loglikelihood1
					dir = builderpath
					if not os.path.exists(dir):
						os.makedirs(dir)
					with open(dir +"/" +str(num_simulations_temp - 2 ) , "wb") as f:
						pickle.dump( myBuilder  , f)
					del myBuilder
					del builderRate_FPEI

				self.lengthPaths = self.num_simulations
				solutions_builderRate = self.averageTime()
				return solutions_builderRate

			if  original == True :
				solutions_builderRate = building_paths()

			else:
				if learndnakinetics.iter  <= 1  :
					builderpath = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+"/"+myenums.Permanent_Folder.MYBUILDER.value+str(self.docID)
					builderratepath =  builderpath + myenums.EDITIONAL.BUILDERRATE.value
					num_simulations_temp = 0
					while num_simulations_temp < self.num_simulations  :
						num_simulations_temp +=1
						myBuilder  = self.get_builder(picklepath=builderpath +"/" + str(num_simulations_temp - 1 )  )

						builderRate_FPEI= BuilderRate_FPEI(myBuilder, trial_count  =trial_count, random_time = random_time, stopStep = 100000 ,path_file =   path_file  + "-" +  str( num_simulations_temp   -1 )  + "-"  , original= original, num_simulations_temp = num_simulations_temp , parameter_folder= learndnakinetics.parameter_folder )
						loglikelihood1 , MFPT1, original_loglikelihood1= builderRate_FPEI.get_essentials ( )
						self.loglikelihoods  [ num_simulations_temp  -1]   = loglikelihood1
						self.MFPTs [num_simulations_temp - 1 ] = MFPT1
						self.original_loglikelihoods[ num_simulations_temp -1 ] = original_loglikelihood1
						dir =builderratepath
						if not os.path.exists(dir ):
							os.makedirs(dir)
						del  builderRate_FPEI.build.neighbors_FPEI
						del builderRate_FPEI.build.protoTransitionsCount
						del builderRate_FPEI.build.protoSequences

						del builderRate_FPEI.build.protoTransitions_FPEI
						builderRate_FPEI.build.protoTransitions_FPEI = copy.deepcopy( builderRate_FPEI.build.protoTransitions )
						del builderRate_FPEI.build.protoTransitions
						newprotoSpace=dict()
						for state1, state2 in builderRate_FPEI.paths[0].edges :
							newprotoSpace[state1] =  builderRate_FPEI.build.protoSpace[state1]
							newprotoSpace[state2] =  builderRate_FPEI.build.protoSpace[state2]
						builderRate_FPEI.build.protoSpace = copy.deepcopy(newprotoSpace)
						del newprotoSpace
						with open (dir  + "/" +  str(num_simulations_temp -1 ) , "wb") as picklerbuilderrate :
							pickle.dump ( builderRate_FPEI, picklerbuilderrate)
						del myBuilder
						del builderRate_FPEI
					self.lengthPaths = num_simulations_temp
					solutions_builderRate= self.averageTime()

				else:
					builderpath = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+"/"+myenums.Permanent_Folder.MYBUILDER.value+str(self.docID)
					builderratepath =  builderpath + myenums.EDITIONAL.BUILDERRATE.value
					countbuilderrates= 0
					num_simulations_temp = 0
					while num_simulations_temp < self.num_simulations  :
						with open(builderratepath +"/"+ str(num_simulations_temp), "rb") as picklerbuilderrate:
							builderRate_FPEI  = pickle.load(picklerbuilderrate )
						countbuilderrates +=1 #do not move this line or the counter which is used to make avg mfpt would be wrong!
						if rate_method == Literals.metropolis:
							set_Metropolis_params(builderRate_FPEI.build.options, self.kinetic_parameters_simulation)
						elif  rate_method == Literals.arrhenius:
							set_Arrhenius_params(builderRate_FPEI.build.options,self.kinetic_parameters_simulation)
						else:
							raise ValueError(' Parameter method not set correctly ' )
						num_simulations_temp +=1
						builderRate_FPEI.reset( trial_count  =trial_count, random_time = random_time, stopStep = 100000 ,path_file =   path_file   + "-" +  str( num_simulations_temp   -1 )  + "-" , original= original, num_simulations_temp = num_simulations_temp , parameter_folder= learndnakinetics.parameter_folder )
						loglikelihood1 , MFPT1, original_loglikelihood1= builderRate_FPEI.get_essentials ( )
						self.loglikelihoods  [ countbuilderrates -1  ]   = loglikelihood1
						self.MFPTs [countbuilderrates-1] = MFPT1
						self.original_loglikelihoods[ countbuilderrates-1] = original_loglikelihood1
						del builderRate_FPEI
					self.lengthPaths = countbuilderrates
					solutions_builderRate = self.averageTime()
			successful_simulation =  True

		if self.save == True :
			attributes_file.write("started saving builder")
			start_time = timeit.default_timer()
			attributes_file.write("finished saving builder" +str(timeit.default_timer() - start_time)+ "\n")

		gc.collect()
		return solutions_builderRate

	def averageTime(self):
		return self.averageTime_normalweighing( )

	"""the Mean First Passage Time is calculated where every path has the same weight"""
	def averageTime_normalweighing(self):
		avg_MFPT = 0
		length =self.lengthPaths
		#print "Individual MFPTs of paths for reaction is ", self.MFPTs
		for i  in  range ( length  )    :
			mfpt =self.MFPTs[i]
			avg_MFPT +=  mfpt
		avg_MFPT = avg_MFPT / length
		#print "MFPT of reaction is " , avg_MFPT
		return avg_MFPT

	""" Calculates Mean First Passage Time, from either FPEI or SSA, and then returns the log reaction rate constant"""
	def find_answers(self):
		concentration  = self.join_concentration
		real_rate = self.real_rate
		bimolecular_reaction = self.bimolecular_reaction
		meanfirstpassagetime = self.find_meanfirstpassagetime()
		attributes_file  = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt"
		attributes_file = open(attributes_file, 'a')
		if meanfirstpassagetime == RETURN_MINUS_INF  or   meanfirstpassagetime<= 0:
			#meanfirstpassagetime should be greater then 0!
			attributes_file.write( "None or Negative MFPT is: "+self.docID  +str(meanfirstpassagetime )+ "\n")
			return   Output( error = np.inf ,  predicted_log_10_rate = np.inf , real_log_10_rate= np.inf  )

		#calculating reaction rate constant from mean first passage time.
		if use_FPEI_MFPT == True :
			if bimolecular_reaction == True :
				predicted_rate= 1.0 / (meanfirstpassagetime * concentration)
			else :
				predicted_rate= 1.0 / meanfirstpassagetime
		elif use_Gillespie_MFPT  == True:
			if bimolecular_reaction == True :
				predicted_rate = meanfirstpassagetime * (1./concentration)
			else:
				predicted_rate = meanfirstpassagetime

		warnings.filterwarnings('error')
		try :
			predicted_log_10_rate =np.log10(predicted_rate)
			real_log_10_rate = np.log10(real_rate)
			error  = math.pow( real_log_10_rate - predicted_log_10_rate, 2)

		except  Exception as e :
			print "exception " , e,   e.args
			return   Output( error = np.inf ,  predicted_log_10_rate = np.inf , real_log_10_rate= np.inf  )
		gc.collect()
		attributes_file.write(  " error: " + str(error) +  "         real_log_10_rate: "  +str (real_log_10_rate)  +  "         predicted_log_10_rate: "  +str (predicted_log_10_rate)  + "\n" )
		attributes_file.write("Next iteration *************************************************************************************************** "+"\n")

		return   Output( error = error ,  predicted_log_10_rate =  predicted_log_10_rate , real_log_10_rate= real_log_10_rate  )

def doReaction(arguments ) :
	# the first argument is always the number of paths
	options =Options(trials = arguments[0])
	#options.output_interval = 1 # do not uncomment ---> might get memory issues
	options.num_simulations= arguments[0]
	options.simulation_time = arguments[1]
	options.sodium = arguments[5]
	options.magnesium = arguments[6]
	options.temperature = arguments[9] +arguments[10]
	options.join_concentration =    arguments[12] * arguments[11]
	if arguments[13] == Literals.metropolis :
		set_Metropolis_params(options, arguments[7] )
	elif arguments[13] == Literals.arrhenius :
		set_Arrhenius_params(options, arguments[7]  )
	options.simulation_mode = Literals.trajectory
	if arguments[14] == True :
		endComplex1 = arguments[15][-1][0]
		stopSuccess = StopCondition(Literals.success, [(endComplex1, Literals.exact_macrostate, 0)])
		options.stop_conditions = [stopSuccess]
		if use_Gillespie_MFPT ==True:
			options.start_state = arguments[15][0]
	return options

def doReaction2(n_trials, arguments ) :
	arguments[0] = n_trials
	return doReaction(arguments)

"""Setting parameters for the Arrhenius kinetic model"""
def set_Arrhenius_params (options, params):
	options.rate_method = Literals.arrhenius
	options.lnAStack = float( params[localtype.stack][LNA_INDEX] )
	options.EStack = float(  params[localtype.stack][E_INDEX] )
	options.lnALoop = float( params[localtype.loop][LNA_INDEX] )
	options.ELoop = float( params[localtype.loop][E_INDEX] )
	options.lnAEnd = float( params[localtype.end][LNA_INDEX] )
	options.EEnd = float(params[localtype.end][E_INDEX] )
	options.lnAStackLoop =  float(params[localtype.stackloop][LNA_INDEX])
	options.EStackLoop =   float( params[localtype.stackloop][E_INDEX] )
	options.lnAStackEnd = float(params[localtype.stackend][LNA_INDEX] )
	options.EStackEnd = float(params[localtype.stackend][E_INDEX] )
	options.lnALoopEnd = float(params[localtype.loopend][LNA_INDEX] )
	options.ELoopEnd =float( params[localtype.loopend][E_INDEX] )
	options.lnAStackStack = float(params[localtype.stackstack][LNA_INDEX] )
	options.EStackStack =float( params[localtype.stackstack][E_INDEX]  )
	options.bimolecular_scaling = float(params[bimolecular_scaling] )


"""Setting parameters for the Arrhenius kinetic model"""
def set_Metropolis_params(options, params):
	options.rate_method = Literals.metropolis
	options.unimolecular_scaling =  float(params [unimolecular_scaling])
	options.bimolecular_scaling  =float( params[bimolecular_scaling])

def main(complex ):
	return complex.find_answers( )
