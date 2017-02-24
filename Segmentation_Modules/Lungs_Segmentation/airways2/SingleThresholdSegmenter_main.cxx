#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "talk2manager.h"
#include "manager.h"
#include "propagator.h"

#include "log.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>

void usage() {

// input image file, output file, seed x, seed y, seed z
	std::cout << "extract_air inputGrayscaleImageFile outputImageFile outputlabelImageFile seed_x seed_y seed_z [timestep] "
			"[stop_time] [debug_on] [inputSpeedImageFile] [output_segment_file] [leakage_recovery] [max_inc_points] [max_inc_trials] "
			"[max_inc_radius] [max_radius ratio] [leakage_th]" << std::endl;
}

int main ( int argc, char ** argv )
{

	init_logging("SingleThresholdSegmenter");
	pLogType l = get_logger("MAIN");

	SegmentationParameters params;

	try {
    	// Declare the supported options.
    	po::options_description desc("Allowed options");
    	desc.add_options()
    	    ("help,h", "produce help message")
    	    ("seed,s", po::value< std::vector<int> >(), "Input seed (call it many times to set x,y,z in sequence)")
    	    ("write-debug-images", "Write debug images")
    	    ("write-colormap", "Forces the visitor to also write an image with the segments and their labels")
    	    ("enable-leakage-recovery,r", "Enable leakage recovery code")
    	    ("enable-leakage-detection,d", "Enable leakage detection code")
    	    ("quiet,q", "Set log level to WARN")
    	    ("verbose,v", "Set log level to DEBUG (A lot of output!)")

    	    ("visitor", po::value<std::string>(&params.role)->implicit_value("visitor"), "Set visitor role")
    	    ("segmenter", po::value<std::string>(&params.role)->implicit_value("segmenter"), "Set segmenter role")
    	    ("input-grayscale,g", po::value< std::string >(&params.input_gray_image_fname), "input grayscale file")
    	    ("input-speed-map,i", po::value< std::string >(&params.input_speed_image_fname), "input speed map file")
    	    ("output-file,o", po::value< std::string >(&params.output_image_fname)->default_value("air.dcm"), "output air segmentation file")
    	    ("output-labels-file", po::value< std::string >(&params.output_label_image_fname)->default_value("labels.tiff"), "output labels file")
    	    ("output-segment-file", po::value< std::string >(&params.output_textfile_fname)->default_value("segs.csv"), "Write segments statistics to a file")
    	    ("timestep,t", po::value<float>(&params.timestep)->default_value(1.0), "Propagation timestep")
    	    ("stop-time", po::value<float>(&params.stop_time)->default_value(1000.0), "Propagation stop time")
			("max-inc-trials", po::value<float>(&params.ST_TRIALS_TH)->default_value(1.5),
								"Maximum increments in trials points bewteen two subsequent propagations")
			("max-inc-radius", po::value<float>(&params.RADIUS_TH)->default_value(1.5),
								"Maximum increments in radius bewteen two subsequent propagations")
			("max-radius-ratio", po::value<float>(&params.RADIUS_RATIO_TH)->default_value(1.1),
								"Maximum ratio bewteen the radius of two subsequent propagations")
			("leakage-th", po::value<int>(&params.LEAKAGE_CORR_TH)->default_value(-700),
								"Threshold used to recover from leakage")
    	;

    	//po::positional_options_description p;
    	//p.add("role", -1);

    	po::variables_map vm;
    	//po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    	po::notify(vm);


    	/*
    	 * General parameters handling
    	 */

    	if (vm.count("quiet"))
    		set_priority(log4cpp::Priority::WARN);

    	if (vm.count("verbose"))
    	    set_priority(log4cpp::Priority::DEBUG);

    	if (vm.count("help") || (argc == 1)) {
    	    std::cout << desc << std::endl;
    	    return 1;
    	}

    	if ((!vm.count("visitor")) && (!vm.count("segmenter"))) {
    		l->warn("You have selected no role, I will default on visitor!");
    		params.role = "visitor";

    	}

    	if (!vm.count("input-speed-map")) {
			l->error("You have to specify the speed image!");
			return 1;
		}


    	if (vm.count("enable-leakage-detection"))
    		params.leakage_detection_enabled = true;

    	if (vm.count("seed")) {

			std::vector<int> seed = vm["seed"].as<std::vector<int> >();
			if (seed.size() != params.seed.GetPointDimension()) {
				l->errorStream() << "I need a " << params.seed.GetPointDimension() << "-D seed, but you passed a "<< seed.size() << "-D one";
				return -1;
			}
			for (unsigned int i=0; i<params.seed.GetPointDimension(); i++)
				params.seed[i] = seed[i];

		} else {
			l->error("You have to specify the seed!");
			return 1;
		}



    	if (params.role == "visitor") {
    		/*
    		 * Visitor role parameters
    		 */
    		if (vm.count("enable-leakage-recovery")) {
    			l->warn("Turning off leakage recovery code in visitor mode...");
    			params.leakage_recovery_enabled = false;
    		}
    		if (vm.count("input-grayscale") || vm.count("max-inc-trials") || vm.count("max-inc-radius")
    				|| vm.count("max-radius-ratio") || vm.count("leakage-th")) {
				l->warnStream() << "With leakage recovery disabled, any of the options: " <<
						"input-grayscale " <<
						"max-inc-trials " <<
						"max-inc-radius " <<
						"max-radius-ratio " <<
						"leakage-th " <<
						" are of no use.";

			}


    	} else if (params.role == "segmenter") {
    		/*
			 * Segmenter role parameters
			 */

    		if (!vm.count("enable-leakage-detection")) {
    			l->info("The leakage detection code is needed byt the segmenter, so I'll enable it!");
    			params.leakage_detection_enabled = true;
    		}

    		if (vm.count("enable-leakage-recovery")) {
				params.leakage_recovery_enabled = true;
				if (!vm.count("input-grayscale")) {
					l->error("You enabled the leakage recovery, so you'll also need to specify the grayscale image we use in the recover process.");
					return -1;
				}
			}



    	} else {
    		l->errorStream() << "You specified an unimplemented role: " << params.role;
    		return 1;
    	}

    	// Print the parameters
    	params.print();
    	// Create segmentation manager
	SegmentationManager manager = SegmentationManager(params);
	manager.do_segmentation();
	manager.print_segments();
	manager.write_segments_tree(params.output_textfile_fname);
	//manager.write_propagations("propagations_all.csv");
	manager.write_segments_tree_oldformat("segments_all.csv");

	if (params.role == "segmenter") {
		l->infoStream() << "Saving segmentation...";
		manager.propagator.save_segmentation( params.output_image_fname );
		manager.save_color_segments();
	}

	if (vm.count("write-colormap")) {
		l->infoStream() << "Writing colormap...";
		manager.save_color_segments();
	}

	if (vm.count("write-debug-images")) {
		l->infoStream() << "Writing debug images...";
		manager.propagator.save_internal_representation( params.output_label_image_fname );
		manager.write_debug_images();
	}


    }
    catch(std::logic_error &err)
    {
         std::cerr << "std::logic_error:" << err.what() << std::endl;
         for (unsigned int i=0; i< argc; i++)
             std::cerr << argv[i] << std::endl;
    }
    catch(std::exception& e) {
    	std::cout << e.what() << std::endl;
    }



	return 0;

} 

