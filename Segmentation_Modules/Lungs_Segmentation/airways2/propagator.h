/*
 * propagator.h
 *
 *  Created on: Mar 24, 2011
 *      Author: mario
 */
#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_
#include <string>
#include "globals_itk.h"

#include "itkFastMarchingRestartableImageFilter.h"
#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include <algorithm>
#include <assert.h>

class Propagator {
	/*
	 *  Handles all the communications with the fast marching code.
	 *  You can set debug_propagations=true to write the internal label image for
	 *  each time step. Beware that it creates a different image for each timestep!
	 */

	bool firstRun;
public:
	bool debug_propagations;

	double stop_time;
	double get_propagation_time() {
		return this->fastMarching->get_last_arrival_time();
	}
	double timestep;

	std::string input_image_fname;

	typedef itk::BinaryThresholdImageFilter< FloatImageType, UCharImageType	>	ThresholdingFilterType;
	ThresholdingFilterType::Pointer thresholder;

	typedef itk::FastMarchingRestartableImageFilter< FloatImageType, FloatImageType >    FastMarchingRestartableFilterType;
	FastMarchingRestartableFilterType::Pointer  fastMarching;
	FloatImageType::SizeType  image_size;

	FloatImagePointer input_image;
	FloatImageType::SpacingType input_spacing;
	FloatImagePointer current_output;
	typedef itk::ImageRegionConstIterator<FloatImageType> ConstIteratorType;

	typedef FastMarchingRestartableFilterType::NodeContainer           NodeContainer;
	typedef FastMarchingRestartableFilterType::NodeType                NodeType;

	Propagator(){
//		debug_propagations=false;
		debug_propagations=true;
	}
	FloatImagePointer read_float_image(const std::string);
	IntImagePointer read_int_image(const std::string);
	void updateSize(FloatImageType::SizeType);
	void write_image(UCharImagePointer, const std::string);
	void write_image(IntImagePointer, const std::string);
	void write_image(FloatImagePointer, const std::string);
	void set_input_image(std::string path) {
		input_image_fname = path;
	}
	bool isWhole(float num, float epsilon) //could use doubles instead, of course
	{
	    return (fabs(num - (int)num) <= epsilon);
	}
	PropagationInfo grow();
	void initialize();
	void set_timestep(double value) {
		timestep = value;
	}
	void set_stop_time(double value) {
		stop_time = value;
	}
	void set_wavefront(const std::vector<PointType>);
	void save_segmentation(const std::string);
	void save_internal_representation(const std::string path);

	void color_segmentation(const std::string, const std::map<UCharPixelType,std::vector<PointType> > );
	void clean_color_segmentation(const std::string, const std::map<UCharPixelType,std::vector<PointType> > );
	std::map<UCharPixelType, std::vector<PointType> > get_cleaned_colormap(const std::map<UCharPixelType, std::vector<PointType> > );


	void clear_processed_points(void);
	void clear_trials(void);
	void clear_propagation_effects(const PropagationInfo);
	std::vector<PointType> get_processed_points(void);
	int get_stabilization_threshold(const float);
	void update_speed_map(FloatImagePointer, RegionType);


	#ifndef M_PI
	#define M_PI       3.14159265358979323846
	#endif
};


#endif /* PROPAGATOR_H_ */









