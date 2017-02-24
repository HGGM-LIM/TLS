/*
 * bifurcations.h
 *
 *  Created on: Mar 25, 2011
 *      Author: mario
 */

#ifndef BIFURCATIONS_H_
#define BIFURCATIONS_H_

#include "globals_itk.h"

class BifurcationChecker {

public:

//	bool check(const PropagationInfo *);
	std::vector< PointType >  walk_connected(const PointType, const PropagationInfo *);
	std::vector< std::vector<PointType> > split_connected(PropagationInfo);

};

#endif /* BIFURCATIONS_H_ */
