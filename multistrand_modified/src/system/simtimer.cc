//This file is modified by NZ, using RVSSA instead  of SSA

#include "simoptions.h"
#include "simtimer.h"

SimTimer::SimTimer(SimOptions& myOptions) {

	maxsimtime = myOptions.getMaxSimTime();
	stopcount = myOptions.getStopCount();
	stopoptions = myOptions.getStopOptions();

	// saving the pointer to enable access to cotranscriptional timing values
	simOptions = &myOptions;

}

// advances the simulation time according to the set rate
void SimTimer::advanceTime(void) {

	rchoice = rate * drand48();
	stime += 1./rate;  //NZ: stime uses expected  holding time as in RVSSA.  In the original Multistrand, stime = (log(1. / (1.0 - drand48())) / rate)
 	//cout << 1./rate << "\n";
	stime2 += (log(1. / (1.0 - drand48())) / rate); 
	// NZ: stime2 is sampling a  random  holding time. stime2 is not actually used in the sampling. If you want to use SSA  instead of RVSSA, swap the values of stime and stime2 

}

// returns TRUE if this transition needs to be executed.
bool SimTimer::wouldBeHit(const double rate) {
	return rchoice < rate;
}

// returns which Nth collision needs to be used
int SimTimer::checkHitBi(const double collisionRate) {

	return (int) floor(rchoice / collisionRate);

}

// returns TRUE if this transition needs to be executed
// and subtracts the rate
bool SimTimer::checkHit(const double rate) {

	rchoice = rchoice - rate;
	return (rchoice < 0);

}

// returns TRUE if a new nucleotide is to be added to the chain
bool SimTimer::checkForNewNucleotide(void) {

	if (simOptions->cotranscriptional && stime > (nuclAdded + simOptions->initialActiveNT) * simOptions->cotranscriptional_rate) {

		nuclAdded++;
		return true;
	}

	return false;
}

std::ostream& operator<<(std::ostream& ss, SimTimer& timer) {

	ss << "rchoice ";
	ss << timer.rchoice << "  rate  ";
	ss << timer.rate << "  simTime  ";
	ss << timer.stime << "\n";

	return ss;
}
