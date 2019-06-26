"""Created by nasim zolaktaf, 2019
This files runs parameter estimation using  either FPEI or Gillespie's SSA. In config_file.txt set which one to use.
To do parameter estimation  run  'python map.py' """
from __future__ import division
from scipy.optimize import minimize
import learndnakinetics
import timeit
import sys
import os
sys.path.insert(0,os.path.realpath('../reactions'))
import parent
import ConfigParser
from multistrand.options import Literals

""" This function runs the MAP approach with  the Nelder-Mead optimization technique!"""
def map ( ):

	learndnakinetics.set_configuration() #setting some configuations

	#Initializing the parameters for Nelder-Mead
	if parent.rate_method == Literals.arrhenius :
		theta = [ 468832.10581490654, 3.,   468832.10581490654,3.,  468832.10581490654, 3.,  468832.10581490654 , 3.,   468832.10581490654, 3.,   468832.10581490654,  3.,  468832.10581490654 , 3.,    0.0402  ]
		#theta = [ 468832.10581490654, 2.,   468832.10581490654,2.,  468832.10581490654, 2.,  468832.10581490654 , 2.,   468832.10581490654, 2.,   468832.10581490654,  2.,  468832.10581490654 , 2.,    0.0402  ]
	elif parent.rate_method == Literals.metropolis:
		theta = [8.2 *  (10 **6), 3.3  * (10**5) ]
	else:
		raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis!')

	learndnakinetics.OVERALLTIME= timeit.default_timer()
	options = dict()

	options['disp'] = True

	configParser = ConfigParser.ConfigParser()
	configParser.readfp(open(r'../learndnakinetics/config_file.txt'))
	CONFIG_NAME = 'parent'
	use_FPEI_MFPT= bool(configParser.getint(CONFIG_NAME, 'use_FPEI_MFPT'))
	use_Gillespie_MFPT = bool(configParser.getint(CONFIG_NAME, 'use_Gillespie_MFPT'))

	if use_FPEI_MFPT == True and use_Gillespie_MFPT == False  :

		options['maxfev'] = 200
		options['maxiter']  =200
		rounds= 10
	elif use_Gillespie_MFPT == True and use_FPEI_MFPT == False :
		options['maxfev'] = 2000
		options['maxiter']  =2000
		rounds= 1
	else :
		raise ValueError('Error: set either use_Gillespie_MFPT or use_FPEI_MFPT to 1 and the other one to 0 in in the configuration file')

	for learndnakinetics.map_i in range (rounds):

		#if learndnakinetics.map_i  == 0 :
		#	learndnakinetics.iter= 1
		learndnakinetics.set_folderfiles()
		if  use_FPEI_MFPT == True :
			#in FPEI in iter == 0  we  build and save the fixed  paths
			learndnakinetics.objective_function(theta)
		#optimizing the paramters
		thetaOptimized = minimize(learndnakinetics.objective_function, theta, method='nelder-mead',options=options )

		learndnakinetics.iter = 0 #After each round we need to reset iter so,  so load, save is overwritted
		theta= thetaOptimized.x
		parameter_file = open(learndnakinetics.parameter_file_name, 'a')
		parameter_file.write("thetaoptimized" + str(thetaOptimized.x) + "\n")
		parameter_file.close()
		print(thetaOptimized.x)

if __name__ == "__main__":
	map()
