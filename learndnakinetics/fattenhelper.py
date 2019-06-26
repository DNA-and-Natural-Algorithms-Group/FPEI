import sys
import cPickle as pickle
import os
import time
sys.path.insert(0,os.path.realpath('../reactions'))
from multistrand.builder import Builder
import parent
from parent import *
import myenums

start = sys.argv[1]
end = sys.argv[2]
pathbuilder = sys.argv[3]
pathspace= sys.argv[4 ]
pathsequences = sys.argv[5]
pathoptions= sys.argv[6]
pathuniqueid= sys.argv[7]


mytime = open ("times_log_remove.txt", "a")
mytime.write( pathbuilder + "   start " + str( start) + "end "  + str(end) + "\n" )
st = time.time()

with open(pathoptions  , "rb" ) as p:
	optionsArg  = pickle.load( p)

myBuilder=   Builder(parent.doReaction,   optionsArg )

with open(pathspace  , "rb" ) as p:
	myBuilder.protoSpacebackup  = pickle.load( p)


with open(pathuniqueid  , "rb" ) as p:
	myBuilder.uniqueID_number  = pickle.load( p)

with open( pathsequences  , "rb" ) as p :
	myBuilder.protoSequences  = pickle.load( p)


mytime.write( "load time " + str( time.time()  - st )+"\n")

st  = time.time ( )
myBuilder.fattenStateSpace(start = int(start)  , end= int( end))

mytime.write( "fatten time time " + str( time.time()   - st ) +"\n")

st = time.time()


with open( pathuniqueid , "wb" ) as p:
	pickle.dump(myBuilder.uniqueID_number,  p)



with open( pathbuilder + myenums.EDITIONAL.TEMPSTATESPACE.value + str(start)+"-"  +str(end) , "wb" ) as p:
	pickle.dump(myBuilder.tempstatespace,  p)


with open( pathbuilder + myenums.EDITIONAL.TEMPTRANSITIONS.value + str(start)+"-"  +str(end) , "wb" ) as p:
	pickle.dump(myBuilder.temptransitions,  p)

mytime.write( "save time " + str( time.time()  - st ) +"\n")
mytime.close()

