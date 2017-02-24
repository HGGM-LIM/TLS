#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "talk2manager.h"
#include "manager.h"
#include "propagator.h"

#include "log.h"



int main ( int argc, char ** argv )
{
	init_logging("test_restart");
	pLogType l = get_logger("MAIN");

	std::string input_gray_image = "./1.3.12.2.1107.5.1.4.54680.30000009090206262446800014419/5-lungs.dcm";
	std::string input_speed_image = "../1.3.12.2.1107.5.1.4.54680.30000009090206262446800014419/4-lungs_mask.dcm";
	std::string output_image = "test_air.dcm";
	std::string output_label_image = "test_labels.tif";

	std::string outputFileCSV="segs.csv";

	int x=260; int y=292; int z=0;
	double timestep = 0.5, stop_time = 4.0;
	int leakage_th = -700;

	bool debug_prop=true;
	bool debug_bifurcations=true;
	bool leakage_recovery_enabled=true;

	float max_inc_points=1.5; float max_inc_trials=1.5; float max_inc_radius=1.5; float max_radius_ratio=1.1;

	l->infoStream() << "input_gray_image: " << input_gray_image ;
	l->infoStream() << "input_speed_image: " << input_speed_image;
	l->infoStream() << "output_image: " << output_image ;
	l->infoStream() << "output_label_image: " << output_label_image ;
	l->infoStream() << "Seed point: (" << x << ", " << y << ", " << z << ")";
	l->infoStream() << "timestep: " << timestep ;
	l->infoStream() << "stop_time: " << stop_time ;
	l->infoStream() << "debug_prop: " << debug_prop ;
	l->infoStream() << "debug_bifurcations: " << debug_bifurcations ;
	l->infoStream() << "Leakage control enabled: " << leakage_recovery_enabled ;

	l->infoStream() << "Maximum accepted points increment between two adjacent propagations: " << max_inc_points ;
	l->infoStream() << "Maximum accepted wavefront increment between two adjacent propagations: " << max_inc_trials ;
	l->infoStream() << "Maximum accepted wavefront radius increment between two adjacent propagations: " << max_inc_radius ;
	l->infoStream() << "Maximum accepted wavefront ratio between two adjacent propagations: " << max_radius_ratio ;




		// Activate segmentation manager
		SegmentationManager manager = SegmentationManager(timestep, stop_time);

		manager.ST_POINTS_TH = max_inc_points;
		manager.ST_TRIALS_TH = max_inc_trials;
		manager.RADIUS_TH = max_inc_radius;
		manager.RADIUS_RATIO_TH = max_radius_ratio;

		manager.set_input_gray_image(input_gray_image);
		manager.set_input_speed_image(input_speed_image);
		manager.set_output_label_image(output_label_image);
		manager.set_output_image(output_image);
		manager.set_leakage_recovery_enabled(leakage_recovery_enabled);
		manager.set_output_textfile(outputFileCSV);
		manager.set_leakage_threshold(leakage_th);
		// Give first seed
		PointType seed;
		seed[0]=x;seed[1]=y;seed[2]=z;
		manager.set_seed(seed);

		// Start segmentation
		manager.debug_propagation(debug_prop);
		manager.set_debug_bifurcations(debug_bifurcations);

		manager.do_segmentation();

		l->infoStream() << "Saving segmentation...";
		manager.propagator.save_segmentation( manager.output_image_fname );
		manager.propagator.save_internal_representation( manager.output_label_image_fname );
		manager.save_color_segments();
		manager.print_segments();
		l->infoStream() << "Speed map changed: " << manager.speed_map_updated;
		l->infoStream() << "Segmentation try: " << manager.times;



	return -1;

} 

