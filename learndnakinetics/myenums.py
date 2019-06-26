from enum import Enum

"""used for saving computations to prevent repeatable computations"""
class Permanent_Folder(Enum):
	MYBUILDER = "MyBuilder"
	TRAJECTORY = "trajectory"

class EDITIONAL(Enum):
	TEMPTRANSITIONS= "temptransitions"
	TEMPSTATESPACE = "tempstatespace"
	UNIQUEID= "uniqueid"
	PROTOSPACEBACKUP ="protoSpacebackup"
	PROTOSEQUENCES=  "protoSequences"
	PATHOPTIONS= "pathoptions"
	BUILDERRATE= "builderrate"

""" Dataset name"""
class DatasetName(Enum):
	BONNET= "Bonnet98"
	HATA= "Hata2017"
	WETMUR = "Wetmur76"
	CISSE2012  = "Cisse2012"

"""Type of reaction"""
class ReactionType(Enum):
	HELIXASSOCIATION = "HelixAssociation"
	HELIXDISSOCIATION ="HelixDissociation"
	HAIRPINCLOSING = "HairpinClosing"
	HAIRPINOPENING = "HairpinOpening"
	THREEWAYDISPLACEMENT = "ThreeWayDisplacement"
	FOURWAYBRANCHMIGRATION = "FourWayBranchMigration"
	BUBBLECLOSING = "BubbleClosing"

class DatasetType(Enum):
	NODANGLE  = "NoDangle"
	DANGLETWOSUBSTRATE = "DangleTwoSubstrate"
	INCUMBENTDANGLE= "IncumbentDangle"
	MISMATCH= "Mismatch"

"""Type of how temperature, reaction rate constant, timescale  is presented in database"""
class 	DatasetSpecifications(Enum) :
	TEMPKELVININV = "TempKevlinInverse"
	TEMPCELCIUS = "TempCelcius"
	TEMPKELVIN  = "TempKelvin"
	RATECONSTANT= "RateConstant"
	LOG10RATECONSTANT= "Log10RateConstant"
	RATECONSTANT10POW5 = "RateConstant10Pow5"
	TIMESCALE = "TimeScale"

"""used in learndnakinetics.py and plot_ssavsftei.py"""
class ImportantFiles(Enum) :
	MSEYNEW= "MSEynew.csv"
	MSEXNEW= "MSExnew.csv"
	MSE= "MSE.csv"
	OVERALLTIME = "OVERALLTIME.csv"

"""used in plot_ssavsftei.py"""
class NamesinPlot(Enum):
	SSAI =  "SSAI"
	FPEI = "FPEI"