/*
 * propagator.h
 *
 *  Created on: Mar 24, 2011
 *      Author: mario
 */


#ifndef MANAGER_H_
#define MANAGER_H_

#include "segments.h"
#include "propagator.h"
#include "bifurcations.h"
#include <vector>
#include "log.h"
#include "leakage.h"


class SegmentationParameters {
	/*
	 * This class holds the parameters of the segmentation manager class
	 */
public:
	float TRACHEA_TYPICAL_DIMENSION; //Minimum length of the trachea, in mm
	std::string role;
	std::string input_gray_image_fname;
	std::string input_speed_image_fname;

	std::string output_image_fname;
	std::string output_label_image_fname;
	std::string output_textfile_fname;

	PointType seed;
	float timestep;
	float stop_time;
	bool leakage_recovery_enabled;
	bool leakage_detection_enabled;
    bool debug_bifurcations;

	/*
	 * Leakage check thresholds:
	 */
	float ST_TRIALS_TH;
	float RADIUS_TH;
	float RADIUS_RATIO_TH;
	int LEAKAGE_CORR_TH;

	SegmentationParameters() {
		TRACHEA_TYPICAL_DIMENSION = 50.0; //Minimum length of the trachea, in mm
		seed.Fill(0);
		timestep = 1.0;
		stop_time = 1000.0;
		leakage_recovery_enabled = false;
		leakage_detection_enabled = false;
		debug_bifurcations = false;
		ST_TRIALS_TH = 1.5;
		RADIUS_TH = 1.5;
		RADIUS_RATIO_TH = 1.1;
		LEAKAGE_CORR_TH = -700;

	}

	bool check() {
		/*
		 * Checks if parameters are correct!
		 */
		return true;
	}

	void print() {

		std::cout << "Mandatory arguments: " << std::endl;
		std::cout << " - Role is: " << role << std::endl;
		std::cout << " - Input speed image is: " << input_speed_image_fname<< std::endl;;
		std::cout << " - Seed point: (" << seed[0] << ", " << seed[1] << ", " << seed[2] << ")" << std::endl;

		std::cout << "General options: " << std::endl;
		std::cout << " - Output segmented image: " << output_image_fname << std::endl;
		std::cout << " - Output colored image: " << output_label_image_fname << std::endl;
		std::cout << " - Output csv file with the segmentation results: " << output_textfile_fname<< std::endl;;
		std::cout << " - debug_bifurcations: " << debug_bifurcations << std::endl;


		std::cout << "Options affecting the propagations: " << std::endl;
		std::cout << " - Timestep: " << timestep  << std::endl;
		std::cout << " - Stop time: " << stop_time << std::endl;
		std::cout << " - Leakage detection status: " << leakage_detection_enabled << std::endl;
		std::cout << " - Leakage recovery status: " << leakage_recovery_enabled << std::endl;
		if (leakage_recovery_enabled) {
			std::cout << " - Maximum accepted wavefront increment between two adjacent propagations: " << ST_TRIALS_TH << std::endl;
			std::cout << " - Maximum accepted wavefront radius increment between two adjacent propagations: " << RADIUS_TH << std::endl;
			std::cout << " - Maximum accepted wavefront ratio between two adjacent propagations: " << RADIUS_RATIO_TH << std::endl;
			std::cout << " - Leakage correction threshold: " << LEAKAGE_CORR_TH << std::endl;
			std::cout << " - Input image used for leakage correction: " << input_gray_image_fname<< std::endl;;
		}


	}

};

class SegmentationManager {
/*
    Segmentation manager takes responsibility over initial setup
    iterate over the segments queue and sending events

*/

	SegmentationParameters params;
public:

	ChangeProbabilityLeakageCorrector* leakage_correction_strategy;
	IntImagePointer grayscale_image;

    std::vector<Segment> segments_queue;
    std::vector<Segment> accepted_segments;


    BifurcationChecker bif_checker;


    int airwayNumber;
    int maxSegmentNumber;


    Propagator propagator;
	int times;
	int speed_map_updated;

	Segment curr_seg;



        SegmentationManager(SegmentationParameters params) {
        	/*
        	 * Check if parameters are consistent!
        	 */

        	pLogType l = get_logger("MANAGER");
        	this->params = params;
       		this->propagator = Propagator();

       		this->times = 0;
       		this->bif_checker = BifurcationChecker();
       		this->leakage_correction_strategy = new ChangeProbabilityLeakageCorrector();

       		l->info("Manager initialized");

        }



        void write_debug_images();
        void write_propagations(const std::string);
        void write_segments_tree(const std::string);
        void write_segments_tree_oldformat(const std::string);
        bool grow_curr_segment(bool, const int);



        bool handle_leakage(Segment, PropagationInfo, IntImagePointer, FloatImagePointer);
        bool handle_change_request(ChangeRequest);
        void do_segmentation();


        IntImagePointer get_grayscale_image();

    	int check_propagation(const PropagationInfo*, const PropagationInfo*, Segment *, Segment *, int);

    	Segment * getNumberedSegment(int number);
    	void print_segments();
    	void get_statistics(std::map<UCharPixelType,std::vector<PointType> >);
    	double calc_angleFromParent(std::vector<PointType>, std::vector<PointType> );
    	double calc_length(std::vector<PointType> );
    	std::vector<PointType> lsqfit_3d(std::vector<PointType> );

    	bool check_cangrow(PropagationInfo*);
    	bool check_leakage(PropagationInfo*, const PropagationInfo*, const Segment*, Segment*);
    	void save_color_segments();
    	std::map<UCharPixelType,std::vector<PointType> > get_color_map();

};


#endif /* MANAGER_H_ */
