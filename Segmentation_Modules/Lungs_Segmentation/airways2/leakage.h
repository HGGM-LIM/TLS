/*
 * ChangeProbabilityLeakageCorrector.h
 *
 *  Created on: May 5, 2011
 *      Author: mario
 */

#ifndef LEAKAGE_H_
#define LEAKAGE_H_

#include "globals_itk.h"
#include "log.h"
#include "segments.h"

class ChangeRequest {
public:
	FloatImagePointer updated_speed_image;
	FloatImagePointer old_speed_image;
	IntImagePointer roi_grayscale;
	RegionType roi;
	bool valid;


};

class LeakageResolutionRequest {


public:
	PropagationInfo bad_prop;
	Segment seg;
	IntImagePointer grayscale_image;
	FloatImagePointer speed_image;

	LeakageResolutionRequest(Segment curr_seg, PropagationInfo prop, IntImagePointer orig_image, FloatImagePointer s_image) {
		this->bad_prop = prop;
		this->seg = curr_seg;
		this->grayscale_image = orig_image;
		this->speed_image = s_image;
	}
};


class ILeakageCorrector
{
    public:
        virtual ~ILeakageCorrector() {}
        virtual ChangeRequest correct(LeakageResolutionRequest) = 0;
        //virtual void set_threshold(int);
};


class ChangeProbabilityLeakageCorrector: public ILeakageCorrector {

	pLogType l;


public:
	int THRESHOLD;
	ChangeProbabilityLeakageCorrector();

	ChangeRequest correct(LeakageResolutionRequest);
	RegionType compute_roi(Segment, const PropagationInfo);
	int get_new_threshold(PropagationInfo);
	void set_threshold(int th) {
		this->THRESHOLD = th;
	}

};




#endif /* LEAKAGE_H_ */
