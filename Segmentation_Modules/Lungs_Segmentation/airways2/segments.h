#ifndef SEGMENT_H_
#define SEGMENT_H_

#include <vector>
#include "propagator.h"
#include "log.h"


class Segment {

	pLogType l;
	bool sealed;

	std::vector<PointType> processed_points;
	std::string status;

public:
	std::vector<PropagationInfo> propagations;
	PropagationInfo get_last_propagation(void);

	bool leaking;
	int segmentNumber;
	int parentSegmentNumber;
	double meanRadius;
	double minRadius;
	int leakage_count;
	void mark_leaking() { seal(); status="leaking"; leaking=true;};

	void setSegmentNumber(int);
	void setParentSegmentNumber(int);
	void SetMeanRadius(double);
	void ComputeMeanRadius(void);
	void ComputeMinimumRadius(void);
	double GetMinimumRadius(void);
	double GetMeanRadius(void);

	IntImagePointer as_image();
	Segment();
	Segment(int, int);
	std::vector<PointType> get_wavefront(void);
	void accept_propagation(PropagationInfo);
	void seal();
	void unseal();

	bool is_sealed() { return sealed; }
	bool has_propagations() { return ( propagations.size() > 0 ); }

	int get_propagations_size() { return propagations.size(); }


	std::string get_status() { return status; }

	std::vector<PointType> get_points_copy(void);

	int get_size(void);
	RegionType get_bbox(int pad = 0);
	RegionType get_bbox(const PropagationInfo, int pad = 0);

};

#endif /* SEGMENT_H_ */
