#include "itkNumericTraits.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "segments.h"
#include "propagator.h"
#include "log.h"
#include <math.h>
#include "itkImageRegionIterator.h"
#include "dicom_utilities.h"

FloatImagePointer Propagator::read_float_image(const std::string path) {
	typedef  itk::ImageFileReader< FloatImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( path.c_str() );

	 try {
	    reader->Update();
	 }catch( itk::ExceptionObject & excep )
	    {
	    std::cerr << "Exception caught !" << std::endl;
	    std::cerr << excep << std::endl;
	    }

	return reader->GetOutput();
}

IntImagePointer Propagator::read_int_image(const std::string path) {
	typedef  itk::ImageFileReader< IntImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( path.c_str() );

	 try {
	    reader->Update();
	 }catch( itk::ExceptionObject & excep )
	    {
	    std::cerr << "Exception caught !" << std::endl;
	    std::cerr << excep << std::endl;
	    }

	return reader->GetOutput();
}

void Propagator::updateSize(FloatImageType::SizeType _size) {
	image_size = _size;
	if (debug_propagations)
		std::cout << "Size set to be: "<< image_size << "\n";
}

void Propagator::write_image(UCharImagePointer image, const std::string path) {

	typedef  itk::ImageFileWriter< UCharImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName( path.c_str() );
	writer->SetInput(image);

	 try {
		 writer->Update();
	 }catch( itk::ExceptionObject & excep )
	    {
	    std::cerr << "Exception caught !" << std::endl;
	    std::cerr << excep << std::endl;
	    }

}

void Propagator::write_image(FloatImagePointer image, const std::string path) {

	typedef  itk::ImageFileWriter< FloatImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName( path.c_str() );
	writer->SetInput(image);

	 try {
		 writer->Update();
	 }catch( itk::ExceptionObject & excep )
	    {
	    std::cerr << "Exception caught !" << std::endl;
	    std::cerr << excep << std::endl;
	    }

}
void Propagator::write_image(IntImagePointer image, const std::string path) {

	typedef  itk::ImageFileWriter< IntImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName( path.c_str() );
	writer->SetInput(image);

	 try {
		 writer->Update();
	 }catch( itk::ExceptionObject & excep )
	    {
	    std::cerr << "Exception caught !" << std::endl;
	    std::cerr << excep << std::endl;
	    }

}

void Propagator::clear_trials() {
	/*
	 * This will set all trial points to Alive, blocking prop
	 */

	  fastMarching->ClearTrials();
}

void Propagator::set_wavefront(const std::vector<PointType> wf) {
	/*
	 *  Call this function each time you want to switch wavefront
	 *  (e.g.) A new segment started
	 */

	pLogType l = get_logger("PROP::set_wavefront");
	l->infoStream() << "Wavefront contains " << wf.size() << " points";
	std::vector<PointType>::const_iterator it;


	if (!firstRun)
		this->clear_trials();


	double seedValue;

	NodeContainer::Pointer seeds = NodeContainer::New();
	seeds->Initialize();

	int el_count = 0;
	for ( it = wf.begin() ; it < wf.end(); it++ ) {
		FloatImageType::IndexType  seedPosition;

		for (unsigned int i=0;i<seedPosition.GetIndexDimension();i++){
			seedPosition[i] = (*it)[i];
		}



				/*
		 * Initialize seed values to their correct time:
		 * - If output is not ready, just set the time to 0.0
		 * - Otherwise read the correct time from the output image itself
		 */

		if (firstRun){
			seedValue = 0.0;
		} else {
			/*
			 * If trials already have a good value, use it, otherwise set their initial time to propagation time
			 */
			double current_time = static_cast< double >(  this->fastMarching->GetOutput()->GetPixel(seedPosition) );
			if (current_time == itk::NumericTraits< FloatPixelType >::max() / 2.0)
				seedValue =  this->get_propagation_time();
			else
				seedValue = current_time;

		}

		l->debugStream() << "SeedPoint: " << seedPosition << " value: " << seedValue;

		NodeType node;

		node.SetValue( seedValue );
		node.SetIndex( seedPosition );

		seeds->InsertElement( el_count++, node );

	}

	fastMarching->SetTrialPoints(  seeds  );
	if (!firstRun)
		fastMarching->RestartPropagation(); //Empty current m_TrialHeap and write back m_TrialPoints to m_TrialHeap


	/*
	 * FastMarchingRestartableFilterType::LabelImagePointer labels = this->fastMarching->GetLabelImage();
		for ( it = wf.begin() ; it < wf.end(); it++ ) {
			FloatImageType::IndexType  seedPosition;

			for (unsigned int i=0;i<seedPosition.GetIndexDimension();i++){
				seedPosition[i] = (*it)[i];
			}

			labels->SetPixel(seedPosition, FastMarchingRestartableFilterType::TrialPoint);



		}
	*/



}

void Propagator::save_segmentation(const std::string path) {
	/*
	 * Save a binary image of the segmented object
	 */
	if(! this->firstRun) {


		thresholder = ThresholdingFilterType::New();
		thresholder->SetLowerThreshold(	0.0);
		thresholder->SetUpperThreshold( this->get_propagation_time() );

		thresholder->SetOutsideValue( 0 );
		thresholder->SetInsideValue( 255 );
		thresholder->SetInput(this->current_output);

		write_image(thresholder->GetOutput(), path);
	}

}

void Propagator::update_speed_map(FloatImagePointer new_speed_image, RegionType roi) {

	/*
	 * Update speed image values.
	 * You can pass an image with the values and the region where to copy them
	 */
	pLogType l = get_logger("PROP::UpdateSpeedMap");
	  l->infoStream() << "start()";

	  typedef itk::ImageRegionIterator< FloatImageType > SpeedIterator;

	  SpeedIterator oldIt( this->input_image, roi );
	  SpeedIterator newIt( new_speed_image , new_speed_image->GetLargestPossibleRegion());

	  oldIt.GoToBegin(); newIt.GoToBegin();
	  while( !newIt.IsAtEnd() )
		{

		 oldIt.Set(newIt.Get());
		++oldIt;
		++newIt;
		}


	//this->fastMarching->SetInput(new_image);

}

void Propagator::color_segmentation(const std::string path, const std::map<UCharPixelType,std::vector<PointType> > colormap) {
	/*
	 * Save a binary image of the segmented object
	 */


	pLogType l = get_logger("PROP::color_segmentation");
	 // First, segment the image
	thresholder = ThresholdingFilterType::New();
	thresholder->SetLowerThreshold(	0.0);
	thresholder->SetUpperThreshold( this->get_propagation_time() );

	thresholder->SetOutsideValue( 0 );
	thresholder->SetInsideValue( 255 );
	thresholder->SetInput(this->current_output);
	thresholder->Update();
	UCharImageType::Pointer segmentation = thresholder->GetOutput();


	// Then change the color of each pixel to the segments it belongs to

	  std::map<UCharPixelType,std::vector<PointType> >::const_iterator it;
	  std::vector<PointType>::iterator pit;

      // show content:
	  for ( it=colormap.begin() ; it != colormap.end(); it++ ) {
		UCharPixelType color = (*it).first;
	    std::vector<PointType> points = (*it).second;
	    l->infoStream() << "Setting " << points.size() << " points to value " << (int)color;
	    for(pit = points.begin(); pit != points.end(); pit++ ){
	    	//Check if the points overlaps
	    	UCharIndexType idx;
	    	for (unsigned int i=0;i<idx.GetIndexDimension();i++)
	    		idx[i] = (*pit)[i];


	    	segmentation->SetPixel(idx,color);

	    }

	  }

	  /*
	   * Set spacing so we can use this image for further anatomical processing too
	   */

	  //segmentation->SetSpacing(this->input_image->GetSpacing());
	  write_image(segmentation, path);

}


void Propagator::clean_color_segmentation(const std::string path, const std::map<UCharPixelType,std::vector<PointType> > colormap) {
	/*
	 * Save a binary image of the segmented object
	 */


	 // First, segment the image
	thresholder = ThresholdingFilterType::New();
	thresholder->SetLowerThreshold(	0.0);
	thresholder->SetUpperThreshold( this->get_propagation_time() );

	thresholder->SetOutsideValue( 0 );
	thresholder->SetInsideValue( 255 );
	thresholder->SetInput(this->current_output);
	thresholder->Update();
	UCharImageType::Pointer segmentation = thresholder->GetOutput();

	// Then change the color of each pixel to the segments it belongs to

	std::map<UCharPixelType,std::vector<PointType> >::const_iterator it;
	std::vector<PointType>::iterator pit;

	// show content:
	for ( it=colormap.begin() ; it != colormap.end(); it++ ) {
		UCharPixelType color = (*it).first;
		std::vector<PointType> points = (*it).second;
		for(pit = points.begin(); pit != points.end(); pit++ ){
			//Check if the points overlaps
			UCharIndexType idx;
			for (unsigned int i=0;i<idx.GetIndexDimension();i++)
				idx[i] = (*pit)[i];

			UCharPixelType previous_color = segmentation->GetPixel(idx);
			if (previous_color != 0) // if not marked as insideValue by fastMarching, need to discard
				segmentation->SetPixel(idx,color);

		}

	}

	/*
	 * Set spacing so we can use this image for further anatomical processing too
	 */

	//segmentation->SetSpacing(this->input_image->GetSpacing());
	write_image(segmentation, path);

}

std::map<UCharPixelType, std::vector<PointType> > Propagator::get_cleaned_colormap(const std::map<UCharPixelType, std::vector<PointType> > colormap) {
	std::map<UCharPixelType, std::vector<PointType> > cleaned;
	std::map<UCharPixelType,std::vector<PointType> >::const_iterator it;
	for ( it=colormap.begin() ; it != colormap.end(); it++ ) {
		UCharPixelType color = (*it).first;
		cleaned[color].clear();
	}

	pLogType l = get_logger("PROP::get_cleaned_colormap");
	l->debug("start()");
	// First, segment the image -- 'clean' / definitive
	thresholder = ThresholdingFilterType::New();
	thresholder->SetLowerThreshold(	0.0);
	thresholder->SetUpperThreshold( this->get_propagation_time() );

	thresholder->SetOutsideValue( 0 );
	thresholder->SetInsideValue( 255 );
	thresholder->SetInput(fastMarching->GetOutput());
	thresholder->Update();
	UCharImageType::Pointer segmentation = thresholder->GetOutput();

	std::vector<PointType>::iterator pit;
	for ( it=colormap.begin() ; it != colormap.end(); it++ ) {
		UCharPixelType color = (*it).first;
//		l->infoStream() << "color: " << color << "\n";

		std::vector<PointType> points = (*it).second;
		for(pit = points.begin(); pit != points.end(); pit++ ){
			//Check if the points overlaps
			UCharIndexType idx;
			for (unsigned int i=0;i<idx.GetIndexDimension();i++)
				idx[i] = (*pit)[i];

			UCharPixelType previous_color = segmentation->GetPixel(idx);
			if (previous_color != 0) {// if not marked as insideValue by fastMarching, need to discard
//				segmentation->SetPixel(idx,color);
				cleaned[color].push_back(*pit);
			}
		}
	}

	return cleaned;
}

void Propagator::save_internal_representation(const std::string path) {
	/*
	 * Save an image that exposes the information on the Trial and Alive points
	 */

	//FIXME::
	std::cout << "PROP::save_internal_representation: ignoring " << path << std::endl;
	std::cout << "PROP::save_internal_representation: images are written to speed-image.mhd and labels.tiff" << std::endl;
	write_image(fastMarching->GetLabelImage(), "labels.tiff");
	write_image(input_image, "speed-image.mhd");
}


PropagationInfo Propagator::grow() {

	pLogType l = get_logger("PROP::grow");

	l->infoStream() << "Current propagation time is " << this->get_propagation_time() << " and stop time is "<<this->stop_time;

	if (this->get_propagation_time() >= stop_time) {
		l->warnStream() << "We've already reached stop time!";
	}

	//clear_processed_points();
	 // Enable processed point for adding then to segments
	fastMarching->SetCollectPoints(true);
	fastMarching->SetStoppingValue( this->get_propagation_time()  + timestep );

	fastMarching->Update();
	this->current_output = fastMarching->GetOutput(); //Assign it so we don't need to recalculate if we call save_segmentation()
	this->current_output->SetSpacing(input_spacing);

	l->infoStream() << "Creating propagation info...";
	PropagationInfo prop = PropagationInfo();

	// Save information for the walkConnected function
	prop.image_size = image_size;

	/*
	 * We populate the propagation object with the new trial points after the propagation.
	 * This is our next wavefront.
	 */

	prop.trials = fastMarching->GetCurrentTrials();
	prop.points = get_processed_points();

	clear_processed_points();


	if (debug_propagations) {
		/*
		 * If you want to debug the wavefront propagation, just set the "debug_propagation flag"
		 * and I'll write you a copy of the internal label image (far, alive and trial points)
		 * for each timestep
		 */


		std::stringstream s_curr_time;
		s_curr_time << this->get_propagation_time();
		if (isWhole(this->get_propagation_time(), 0.0001)) {
			std::string fname = "prop_" + s_curr_time.str() + ".tiff";
			l->infoStream() << "Writing debug images...";
			save_internal_representation( fname );
		}

	}
	//std::cout << "pp 6\n";


	/*
	 *  After first propagation, we need to be sure that the fastmarching is in
	 * "restart mode" fastMarching->SetRestarted(true);
	 *  Or it will loose memory of previous state
	 *
	 */

	if (firstRun){
		firstRun=false;
		fastMarching->SetRestarted(true);
	}


	return prop;

}

void Propagator::initialize() {

	/*
	 *  Instantiate the fast marching code, and convert the input image to a speed image
	 *  Also holds a reference to the fastmarching potential image to control it and restart.
	 */

	  input_image = this->read_float_image(this->input_image_fname);
	  FloatImageType::IndexType index;
	    index[0] = 255;
	    index[1] = 240;
	    index[2] = 0;
	  std::cout <<"PIXEL VALUE " <<input_image->GetPixel(index) << std::endl;
	  input_spacing = this->input_image->GetSpacing();
	  input_spacing[2] = get_slice_thickness_tag(this->input_image);
	  std::cout << input_spacing[0]  << " "<< input_spacing[1] << " " << input_spacing[2] << std::endl;

	  updateSize(input_image->GetLargestPossibleRegion().GetSize());
	  fastMarching = FastMarchingRestartableFilterType::New();

	  /*
	   * Create speed image
	   */
	  fastMarching->SetInput( input_image );
	  fastMarching->SetOutputSize(  input_image->GetBufferedRegion().GetSize() );

	  /*
	   * After we use the fm once, we need to set the restarted property to ensure proper keeping of the data
	   */
	  firstRun=true;

}

std::vector<PointType> Propagator::get_processed_points() {

	std::vector<PointType> bag;
	NodeContainer::Iterator it= fastMarching->GetProcessedPoints()->Begin();
	NodeContainer::Iterator end= fastMarching->GetProcessedPoints()->End();

	while(it != end) {
		PointType point;
		NodeType node = it.Value();
		UCharIndexType idx = node.GetIndex();
		for (unsigned int i = 0; i<idx.GetIndexDimension();i++)
			point[i]=idx[i];

		bag.push_back(point);
		++it;
	}

	return bag;
}

void Propagator::clear_processed_points() {
	fastMarching->GetProcessedPoints()->Initialize();
}

void Propagator::clear_propagation_effects(const PropagationInfo bad_prop) {
	/*
	 * Remove effects of the propagation.
	 * Effects to be reversed are the following:
	 * 1) Points and trials go back to Far status
	 * 2) Their potential value shold be reversed to the one they have at the beginning of the propagation
	 *   This is more tricky because we don't have a stack to save last propagation status.
	 *   We could take two different approaches:
	 *   a) actually the potential value corresponds to the arrival time of the FM wave.
	 *   	So we could just put it back to an invalid value or to the last known value
	 *   	To do this we should save the current time in the propinfo.
	 *   b) We can change the underlining FM implementation to do a copy-on-write of the potential map
	 *   	Thus allowing for a simple stack solution
	 *   c) Do nothing about it. I don't really know if it affects the results.
	 *
	 *
	 * FIXME: For now we don't reverse the potential, if we see it is a problem, we'll dwelve deeper
	 */

	// Get a reference to the label image
	FastMarchingRestartableFilterType::LabelImagePointer labels = fastMarching->GetLabelImage();
	FastMarchingRestartableFilterType::OutputImagePointer output = fastMarching->GetOutput();
	FastMarchingRestartableFilterType::PixelType large_val = static_cast< FastMarchingRestartableFilterType::PixelType >
		( itk::NumericTraits< FastMarchingRestartableFilterType::PixelType >::max() / 2.0 );

	// Set all points and trials to Far
	std::vector< PointType >::const_iterator it;
	FloatImageType::IndexType  idx;

	for ( it=bad_prop.points.begin() ; it < bad_prop.points.end(); it++ ) {
		for (unsigned int i=0;i<idx.GetIndexDimension();i++)
			idx[i] = (*it)[i];

		labels->SetPixel(idx,FastMarchingRestartableFilterType::FarPoint);
		output->SetPixel(idx,large_val);

	}

	for ( it=bad_prop.trials.begin() ; it < bad_prop.trials.end(); it++ ) {
		for (unsigned int i=0;i<idx.GetIndexDimension();i++)
			idx[i] = (*it)[i];

		labels->SetPixel(idx,FastMarchingRestartableFilterType::FarPoint);
		output->SetPixel(idx,large_val);
	}


}

int Propagator::get_stabilization_threshold(const float typical_lenght) {
	/*
	 * The stabilization threshold is the expected number of propagation steps which are needed
	 * to reach the typical dimension of the object to propagate into
	 */
	pLogType l =  get_logger("get_stabilization_threshold");
	l->debugStream() << "Typical lenght is " << typical_lenght;
	l->debugStream() << "X Spacing is " << input_image->GetSpacing()[0];
	l->debugStream() << "Stable wavefront should be reached in " << typical_lenght/input_image->GetSpacing()[0];
	return ceil( typical_lenght / input_image->GetSpacing()[0] );
}
