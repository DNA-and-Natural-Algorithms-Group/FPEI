"""created by Nasim Zolaktaf, 2019
This file calculates the MFPT of a pathway"""
from __future__ import division
from multistrand.builder import transitiontype, Energy, BuilderRate, Builder
import random
import cPickle as pickle
import numpy as  np
import learndnakinetics
import parent
from multistrand.options import Options, Literals
import math

class Path(object):

	def __init__(self):
		self.states= []
		self.edges= []
		self.stateCount = dict()
		self.edgeCount = dict()
		self.MFPT = None
		self.loglikelihood= None

class BuilderRate_FPEI(BuilderRate):

	def __init__(self, builderIn , trial_count, random_time, stopStep ,  path_file , original, num_simulations_temp, parameter_folder )  :
		BuilderRate.__init__(self, builderIn)
		self.trial_count = trial_count
		self.random_time = random_time
		self.stopStep = stopStep
		self.parameter_folder= parameter_folder
		self.paths =  []
		self.holding_time =dict()
		self.rates= dict( )
		self.total_rate = dict(  )
		self.original = original
		self.states_total_rate= dict( )
		self.states_rates = dict()
		self.original = original
		self.path_file = path_file
		self.holding_file = self.path_file+"holdingtime"
		self.rates_file = self.path_file + "rates"
		self.num_simulations_temp= num_simulations_temp
		self.classifyneighbors=  dict( )
		self.lenneighbors= dict()
		self.initialize()

	def reset(self, trial_count, random_time, stopStep ,  path_file , original,  num_simulations_temp, parameter_folder)  :
		self.trial_count = trial_count
		self.random_time = random_time
		self.stopStep = stopStep
		self.parameter_folder= parameter_folder
		self.paths =  []
		self.holding_time =dict()
		self.rates= dict( )
		self.total_rate = dict(  )
		self.original = original
		self.states_total_rate= dict( )
		self.states_rates = dict()
		if learndnakinetics.iter <=1 :
			self.classifyneighbors = dict()
			self.lenneighbors= dict()
		self.original = original
		self.path_file = path_file
		self.holding_file = self.path_file+"holdingtime"
		self.rates_file = self.path_file + "rates"
		self.num_simulations_temp= num_simulations_temp
		self.initialize()

	def initialize(self):

		if  self.original== True:
			self.BuilderPaths()
			#setting the loglikelihoods of paths before saving, because likelihoods will be used in importance sampling
			for path in self.paths:
				self.get_loglikelihood(path)
				self.get_MFPT(path)
			filepickle =  open( self.path_file  , "wb" )
			pickle.dump(self.paths, filepickle)
		else:
			with open(self.path_file , "rb" ) as filepickle:
				self.paths = pickle.load(filepickle)
			for path in self.paths:
				path.MFPT = None
				path.loglikelihood= None
			with open(self.path_file , "rb" ) as filepickle:
				self.paths_original = pickle.load(filepickle)

	"""storing states and edges of path. If you are not implementing importance sampling then there is no need to store edges. In the current FPEI implementation edges is not used"""
	def BuilderPaths(self):
		path = Path ()
		totaltransition = 0

		for edge in self.build.protoTransitions :
			state1, state2= edge
			transitioncount= self.build.protoTransitionsCount[state1, state2]
			path.edgeCount[state1, state2] =transitioncount
			totaltransition +=transitioncount
			path.edges.append((state1,state2) )
			try:
				path.stateCount[state1] +=transitioncount
			except:
				path.states.append(state1)
				path.stateCount[state1] = transitioncount

		self.paths.append(path)


	"""get MFPT of path """
	def get_MFPT(self, path):

		if path.MFPT == None :
			path.MFPT=0
			states= path.states
			for state in states:
				path.MFPT += ( self.get_holding_time(state)  *path.stateCount[state])
		return path.MFPT

	"""Log Likelhoods would be be useful for importance sampling,
	 Currently these values are not used in FPEI!"""
	def get_loglikelihood(self, path):

		if path.loglikelihood == None:
			#path.loglikelihood = 1
			path.loglikelihood = 0
			for edge in path.edges:
				state1, state2= edge
				rate1, rate2 = self.get_rate_stable(state1, state2)
				#path.loglikelihood = path.loglikelihood *  ( ( rate1 / self.get_total_rate(state1) )  ** path.edgeCount[state1,state2] )
				path.loglikelihood +=    ( path.edgeCount[state1,state2]  *    ( math.log10( rate1 )  - math.log10(self.get_total_rate(state1) )  ) )

		return path.loglikelihood

	def get_rate_stable(self , state1, state2 ):
		if not (state1, state2) in  self.rates:
			try :
				rate1 , rate2 = self.get_rate(state1, state2)
			except:
				rate2, rate1=  self.get_rate(state2,state1)
			self.rates[state1, state2]  = rate1
			self.rates[state2, state1] = rate2
		return self.rates[state1, state2],  self.rates[state2, state1]


	"""calculate holding time of state"""
	def get_holding_time(self ,state) :
		if state not in self.holding_time:
			total_rate = self.get_total_rate(state)
			if self.random_time == True :
				self.holding_time[state] = random.expovariate(total_rate)
			else :
				self.holding_time[state] =   1./total_rate

		return self.holding_time[state]

	"""classify neighbors based on type of transition to reduce memory usage"""
	def get_classifyneighbors(self , state):
		myT = self.build.options._temperature_kelvin
		RT = Energy.GAS_CONSTANT * myT
		if state not in self.classifyneighbors:
			self.classifyneighbors[state] =dict( )
			self.lenneighbors[state] =dict( )
			neighbors= set( self.build.neighbors_FPEI[state] )

			for state2 in neighbors:
				structureType = self.build.protoTransitions_FPEI[(state, state2)]
				dG1 = self.build.protoSpace[state].dG(myT)
				dG2 = self.build.protoSpace[state2].dG(myT)
				DeltaG = dG2 - dG1
				DeltaG2= -DeltaG
				dgsum = np.e **  ( -DeltaG/RT)
				dg2sum =  np.e **  ( -DeltaG2/RT)
				if DeltaG > 0 :
					greater = +1
				else:
					greater = -1
				if parent.rate_method ==  Literals.arrhenius:
					try :
						self.classifyneighbors[state][structureType[0], structureType[1], structureType[2], greater][0] += dgsum
						self.lenneighbors[state][ structureType[0], structureType[1], structureType[2], greater ] +=1
					except:
						self.classifyneighbors[state][ structureType[0], structureType[1], structureType[2], greater ]  =  [ dgsum , 0 ]
						self.lenneighbors[state][ structureType[0], structureType[1], structureType[2], greater ]  = 1
					self.classifyneighbors[state][  structureType[0], structureType[1], structureType[2], greater][1] +=  dg2sum

				elif parent.rate_method ==   Literals.metropolis:
					#-4 is not used , dont bother about it
					try :
						self.classifyneighbors[state][structureType[0],  -4,  -4 , greater][0] += dgsum
						self.lenneighbors[state][ structureType[0],  -4,  -4,  greater ] +=1
					except:
						self.classifyneighbors[state][ structureType[0], -4, -4, greater ]  =  [ dgsum , 0 ]
						self.lenneighbors[state][ structureType[0], -4, -4,  greater ]  = 1
					self.classifyneighbors[state][structureType[0], -4, -4, greater][1] +=  dg2sum
		return self.classifyneighbors[state] , self.lenneighbors[state]

	"""calculate total outgoing rate of state"""
	def get_total_rate_classified(self, state):
		if state not in self.total_rate:
			total_rate   = 0
			dgsums , lenneighbors   = self.get_classifyneighbors(state)

			for s0, s1, s2, greater in dgsums :
				[dgsum , dg2sum ] = dgsums[s0, s1, s2 , greater  ]
				n = lenneighbors[s0, s1, s2, greater]
				if parent.rate_method== Literals.arrhenius:
					rate1, rate2 = self.arrhenius_rate (state, "", [s0, s1, s2], greater , dgsum, dg2sum , n  )
				elif parent.rate_method == Literals.metropolis:
					rate1, rate2 = self.metropolis_rate (state, "", [s0, s1, s2], greater , dgsum, dg2sum , n  )
				total_rate += rate1
			self.total_rate [state ] = total_rate
		return self.total_rate[state]

	def get_total_rate(self ,state) :
		return self.get_total_rate_classified(state)

	def get_essentials(self):
		path = self.paths[0]
		try :
			to  = self.paths_original[0].loglikelihood
		except:
			to = None
		return self.get_loglikelihood(path) , self.get_MFPT(path), to

