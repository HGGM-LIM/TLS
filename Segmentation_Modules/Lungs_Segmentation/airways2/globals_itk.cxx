/*
 * globals_itk.cxx
 *
 *  Created on: Apr 6, 2011
 *      Author: mario
 */
#define USE_GLOBALS

#include "globals_itk.h"
//#include <assert.h>

std::vector< PointType > getNeighs(const PointType point, const PointType bb_start, const PointType bb_end) {
	/*
	 * Computes the six-connectivity neighbours of the given point.
	 * Points which are outside the bounding box are not added
	 * In order to check boundary we have a start and a stop point which represents the two vertex
	 * of a bounding box
	 * if p - start has some components less then zero than it is outside the lower bound
	 * if p - end has some components more than zero than it is outside upper bound
	 * (in the following example p1 is ok all the others are outside)
	 *
	 *  p2                        p5
	 *       S-----------------
	 *       -                -
	 *  p3   -         p1     -
	 *       -----------------E
	 *                            p4
	 *
	 */

	std::vector< PointType > neighs;

	for (unsigned int i=0; i<OFFSET_LENGHT; i++){



		IntPointType new_point = point;
		new_point[0] += offsets[i][0];
		new_point[1] += offsets[i][1];
		new_point[2] += offsets[i][2];

		VectorType lower_test = new_point - bb_start;
		VectorType upper_test = new_point - bb_end;


		bool inside = true;
		for (int j=0;j<LDimension;j++){
			if (lower_test[j] <= 0) {
				inside = false;
				break;
			}
		}
		for (int j=0;j<LDimension;j++){
			if (upper_test[j] >= 0) {
				inside = false;
				break;
			}
		}

		if (inside)
			neighs.push_back( new_point );

	}




	return neighs;
}

bool contains(const std::vector< PointType > *bag, const PointType point) {
	/*
	 * Checks if a point is already present in a container.
	 */

	bool found = false;

	std::vector< PointType >::const_iterator it;

	it = std::find(bag->begin(), bag->end(), point);
	if (it != bag->end()) return true;

//	for ( it = bag->begin() ; it < bag->end(); it++ ) {
//		PointType curr = (*it);
//		bool all_comp_equals=true;
//		for (int i=0; i<curr.PointDimension;i++){
//			if (curr[i] != point[i]){
//				all_comp_equals=false;
//				break;
//			}
//		}
//
//		if (all_comp_equals){
//			found = true;
//			break;
//		}
//
//	}
	return found;

}

