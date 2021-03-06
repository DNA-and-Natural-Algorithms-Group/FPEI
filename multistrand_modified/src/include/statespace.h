/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

/*
 *  Created on: Feb 23, 2018
 *      Author: Frits Dannenberg
 *
 *      This class collects every visited state, and transitions between visited states.
 *      After the simulation is done, the set is exported to a text file.
 *
 */
#ifndef __STATESPACE_H__
#define __STATESPACE_H__

#include <scomplexlist.h>
#include <unordered_map>
#include <unordered_set>

using std::unordered_map;
using std::unordered_set;

class SimOptions;

/*
 * 		protoSpace : map with states as entries
 *
 *		protoTransitions: map with transitions as entries
 *
 *		key: states. Value: a initCountFlux object that tells how many times the state has been the initial state and the join flux (rate)
 *		self.protoInitialStates = dict()
 *
 *		key: states: Value: the result of this final state is typically SUCCES or FAILURE (tag)
 *		self.protoFinalStates = dict()
 */
class Builder {
public:

	Builder();
	Builder(SimOptions* options);
	friend std::ostream& operator<<(std::ostream&, Builder&);

	void addState(ExportData&, const double arrType);
	void stopResultNormal(double, string);
	void writeToFile(void);
	string filename(string);

	ExportData lastState;

	static const string the_dir;

private:

	SimOptions* simOptions = NULL;

	unordered_set<ExportData> protoSpace;
	unordered_set<ExportTransition> protoTransitions;
	
	unordered_map<ExportTransition, int> protoTransitionsCount; // added by nasim to count the number of time a transition appears
	
	unordered_map<ExportData, ExportFinal> protoFinalStates;
	unordered_map<ExportData, ExportInitial> protoInitialStates;
	//int totalcounter = 0 ; 
};

#endif
