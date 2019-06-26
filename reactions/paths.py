"""created by Nasim Zolaktaf, 2019
This files creates the start and stop states for a path for hairpin closing and opening, helix association and dissociation reactions"""
from __future__ import division
from multistrand.experiment import standardOptions, makeComplex
import string
NUCLEOTIDES = "ACTG"
TRANSLATION_TABLE = string.maketrans(NUCLEOTIDES, "TGAC")
def Complement(str):
	return ''.join(list(reversed(str))).translate(TRANSLATION_TABLE)

class Path (object) :

	def __init__(self ,  strands_list , reaction_type,  dataset_type, dataset_name   )  :
		self.dataset_type = dataset_type
		self.reaction_type = reaction_type
		self.dataset_name =  dataset_name
		self.strands_list  = strands_list

	def generate_statespace( self   ):
		#Generates the state space
		state = self.initial_final_state_config()[0]
		self.statespace = []
		self.fast_access= dict()
		self.statespace.append(state)
		self.fast_access[state] = len(self.statespace)- 1
		color= dict()
		color [state] = True
		head =  0
		tail = 1
		while head < tail    :
			state = self.statespace[head ]
			pstates = self.possible_states( state )
			for s in pstates :
				if s not in color  :
					color[s] = True
					self.statespace.append(s )
					self.fast_access[s] = len(self.statespace)- 1
					tail +=1
			head += 1
		return self.statespace

	def myFilter(self, state):
		dotparens= dict()
		dotparen_list = self.checkformismatchbasepair(state)
		if str(dotparen_list) in dotparens:
			return False
		else :
			dotparens[str(dotparen_list)]  = True
		if self.checkfordisconnectedstrand(dotparen_list) == True  :
			return False
		sequence_list= self.sequence(state)
		return sequence_list, dotparen_list

	def checkfordisconnectedstrand(self, dotparenlist):

		for dotparen in dotparenlist:
			dotparens = dotparen.split('+')
			if len(dotparens) ==1 :
				continue
			else:
				for dp in dotparens:
					counter = 0
					for i in range(len(dp)):
						if dp[i]!= '.' :
							counter +=1
					if counter == 0 :
						return True
		return False

	def checkformismatchbasepair(self, state):
		def match(sq, sq2):
			if sq == 'A' and sq2 == 'T' or sq == 'T' and sq2 == 'A':
				return True
			elif  sq == 'C' and sq2 == 'G' or sq == 'G' and sq2 == 'C':
				return True
			else:
				return False
		dotparens = self.dot_paren(state)
		sequences=  self.sequence(state , formismatch = True  )
		mystacksq =  []
		mystacki =  []
		newdotparens= []

		for dp , sq in zip ( dotparens, sequences ) :
			dplist= list(dp)
			for  i in range (len(dp)):

				if dp[i] == "(" :
					mystacksq.append (sq[i])
					mystacki.append(i)
				elif dp[i] == ")" :
					sq2 =  mystacksq.pop( )
					ii = mystacki.pop()
					if match(sq[i], sq2)  == False :

						dplist[ii]  = '.'
						dplist[i]  = '.'
			dplist= "".join(dplist)
			newdotparens.append(dplist)

		return newdotparens

	def placeInitialFinal(self, pathway_list):
		#initial state must be in first position of pathway_list  and final state must be in last position of pathway_list
		initial, final= self.initial_final_state_config ()
		i  = self.statespace.index(initial)
		f = self.statespace.index(final)
		if i !=  0 :
			tempi = pathway_list[0 ]
			pathway_list[0 ] = pathway_list[i]
			pathway_list[ i ] = tempi
		size = len(self.statespace) -1
		if f !=  size :

			tempf = pathway_list[size ]
			pathway_list[size] =  pathway_list[f]
			pathway_list[ f ] = tempf

	def get_pathway(self) :
		pathway_list  = self.get_startStates()
		self.placeInitialFinal(pathway_list)
		return pathway_list

	def get_intial_final(self):
		pathway_list = self.get_pathway()
		initial_final = []
		initial_final.append(pathway_list[0])
		initial_final.append(pathway_list[len(pathway_list) - 1])
		return initial_final

class PathHelix(Path):

	def __init__(self ,association ,  strands_list , reaction_type , dataset_name ,  dataset_type ,  cutoff = 1. )  :
		Path.__init__(self,   strands_list , reaction_type ,  dataset_type, dataset_name  )
		self.cutoff = cutoff
		self.association = association
		self.strand1 = self.strands_list [0]
		self.strand2 = self.strands_list[1]
		self.L = len(self.strand1)

		self.botleftdangle = ""
		self.botrightdangle = ""
		self.topleftdangle= ""
		self.toprightdangle= ""

	def get_startStates(self):

		self.initial_final_state_config()
		self.generate_statespace()
		pathway_list = []
		strand_ids = [100,200]
		for state in self.statespace:
			output = self.myFilter(state)
			if output == False :
				continue
			else :
				sequence_list, dotparen_list  = output
			if len(dotparen_list)   == 1 :
				new_complex = makeComplex (sequence_list, dotparen_list[0], strand_ids  )
				pathway_list.append([new_complex] )
			elif len(dotparen_list ) == 2:
				new_complex0 = makeComplex ( [ sequence_list[0] ] , dotparen_list[0], [strand_ids[0]]  )
				new_complex1 = makeComplex ( [ sequence_list[1] ], dotparen_list[1] , [strand_ids[1]]  )
				pathway_list.append( [new_complex1, new_complex0] )


		return pathway_list

	def possible_states(self, state):
		"""Returns the neighbors of state"""
		i, j= state
		if (i == j):
			states = [(n, n + 1) for n in range(0, self.L)]
		else:
			states = [(i - 1, j),
			          (i + 1, j),
			          (i, j - 1),
			          (i, j + 1)]
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and 0 < s[0] <= self.L  :
				removeList.append((s[0],s[1]))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((0,0))
		return filter(self.allowed_state, states)

	def initial_final_state_config(self ):
		"""sets the initial and final state for helix association (association == True) and helix dissociation (association == False ) """
		if self.association == True :
			initialStateConfig = (0, 0 )
			finalStateConfig  = (0, self.L)
		if self.association == False:
			initialStateConfig = (0, self.L )
			finalStateConfig  = (0, 0)
		return [initialStateConfig, finalStateConfig]

	def allowed_state(self, state):
		"""Check that a state is allowed."""
		i, j= state
		allow = 0 <= i <= j <=  self.L
		if  ( j - i ) / self.L > self.cutoff :
			allow =  False
		return allow

	def dot_paren(self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		L = self.L
		i, j = state
		dotparen2 = '.' * len(self.botleftdangle)+'.' * (L - j) + '(' * (j - i) + '.' * i+ '.' * len(self.botrightdangle)
		dotparen1 =   '.' * len(self.topleftdangle) +   '.' * i       + ')' * (j - i) + '.' * (L-j) + '.' * len(self.toprightdangle)

		if i == j:
			return  [  dotparen2  , dotparen1]
		else:
			return  [ dotparen2 + '+' +  dotparen1 ]

	def sequence(self, state , formismatch= False):
		i, j  =state
		strand1 = self.topleftdangle +self.strand1 + self.toprightdangle
		strand2= self.botleftdangle +self.strand2 + self.botrightdangle
		if i ==j :

			return [strand2 , strand1 ]
		else:
			if formismatch == True :
				return [strand2 +'+'+  strand1 ]
			else:
				return [strand2 ,  strand1  ]

class PathHairpin(Path):

	def __init__(self ,association ,  strands_list , reaction_type , dataset_name ,  dataset_type, cutoff =1. )  :
		Path.__init__(self,   strands_list , reaction_type ,  dataset_type, dataset_name  )
		self.cutoff = cutoff
		self.association = association
		self.hairpin = strands_list [0] +strands_list [1]  +Complement(strands_list [0])
		self.m = len(strands_list[0])
		self.L = len(self.hairpin)

	def get_startStates(self):
		self.initial_final_state_config()
		self.generate_statespace()
		pathway_list = []
		strand_ids = [1]
		for state in self.statespace:
			output = self.myFilter(state)
			if output == False :
				continue
			else :
				sequence_list, dotparen_list  = output

			new_complex = makeComplex (sequence_list, dotparen_list[0],  strand_ids  )
			pathway_list.append([new_complex] )

		return pathway_list

	def possible_states(self, state ):
		"""Returns the neighbors of state"""
		i, j = state
		if (i == j):
			states = [(n, n + 1) for n in range(0, self.m)]
		else:
			states = [(i - 1, j),
			          (i + 1, j),
			          (i, j - 1),
			          (i, j + 1)]
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and 0  < s[0]  <=   self.m  :
				removeList.append((s[0],s[1]))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((0,0))
		return filter(self.allowed_state, states)

	def initial_final_state_config(self ):
		"""sets the initial and final state for hairpin closing (association == True) and hairpin opening (association == False ) """
		if self.association == True :
			initialStateConfig = (0, 0 )
			finalStateConfig  = (0, self.m)
		if self.association == False:
			initialStateConfig = (0, self.m )
			finalStateConfig  = (0, 0)
		return [initialStateConfig, finalStateConfig]


	def allowed_state(self, state):
		"""Check that a state is allowed."""
		i, j = state
		allow =0 <= i <= j <= self.m
		if  ( j - i ) / self.m > self.cutoff :
			allow = False
		return allow

	def dot_paren(self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		m, L = self.m, self.L
		i, j = state
		dotPar = '.' * i + '('* (j-i) + '.' * (m - j ) + '.' *  (L - 2*m)  + '.' *( m - j ) + ')'* (j-i) + '.' * i
		return  [dotPar]

	def sequence(self  ,state , formismatch= False):
		return [self.hairpin ]

