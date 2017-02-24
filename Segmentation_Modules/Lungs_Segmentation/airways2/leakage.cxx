/*
 * ChangeProbabilityLeakageCorrector.cpp
 *
 *  Created on: May 5, 2011
 *      Author: mario
 */

#include "leakage.h"



ChangeProbabilityLeakageCorrector::ChangeProbabilityLeakageCorrector() {
	l = get_logger("LeakageCorrector");
	l->info("Change Probability Strategy created");

}

int ChangeProbabilityLeakageCorrector::get_new_threshold(PropagationInfo prop) {
	return this->THRESHOLD;
}

ChangeRequest ChangeProbabilityLeakageCorrector::correct(LeakageResolutionRequest lrr) {
	/*
	 * The correction will work in this way:
	 * 1) We create the ROI adding up to the segment bounding box also the leaking propagation
	 * 	  and expanding it a bit to account for border effects (espcially in CurvatureDiffusion filter)
	 * 1b) Write down the ROI to test it
	 * 2) We set the processing region only to the ROI
	 * 3) We start by implementing only the Ostu re-thredholding
	 * 4) The output of the otsu filter is casted to float and the corresponding region of the speed map is overwritten
	 * 5) the resulting change request will ask for removal of the last propagation and restart the entire segment
	 *
	 */

	ChangeRequest cr = ChangeRequest();
	Segment seg = lrr.seg;
	l->infoStream() << "Reading leakage resolution request...";
	PropagationInfo bad_prop = lrr.bad_prop;
	IntImagePointer grayscale_image = lrr.grayscale_image;
	FloatImagePointer speed_image = lrr.speed_image;


	RegionType bb= bad_prop.get_bbox();
	l->infoStream() << "Leaking segment: " << seg.segmentNumber << " (parent:" << seg.parentSegmentNumber << ")";
	l->infoStream() << "Segment already had: " << seg.leakage_count << " previous leakages.";
	l->infoStream() << "Leaking propagation: #trials " << bad_prop.trials.size() << " #points " << bad_prop.points.size() << " bbox: "<<bb.GetIndex()<<" -> "<<bb.GetSize();

	//1st create the ROI
	RegionType roi = compute_roi(seg, bad_prop);


	// Multi th strategy

	ExtractFilterPointer extractor = ExtractFilterType::New();
	FloatExtractFilterPointer fextractor = FloatExtractFilterType::New();
	BinFilterPointer thresholder = BinFilterType::New();
	int lowerThreshold = -1024;

	int upperThreshold = this->get_new_threshold(bad_prop);

	if (upperThreshold != -1) {//no more tries...
		extractor->SetInput(grayscale_image);
		extractor->SetExtractionRegion(roi);
		thresholder->SetLowerThreshold(lowerThreshold);
		thresholder->SetUpperThreshold(upperThreshold);
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->SetInput(extractor->GetOutput());
		l->infoStream() << "Re-thresholding the roi region...";
		thresholder->Update(); //Ensure output is valid before going out of scope
		l->infoStream() << "Packing output...";
		FloatImagePointer new_map = thresholder->GetOutput();

		cr.updated_speed_image = new_map; //This is only a roi-sized image.
		cr.roi = roi; //We need roi info to be able to paste it again in the speed image.
		cr.roi_grayscale = extractor->GetOutput(); //Let's see what the roi actually took
		fextractor->SetInput(speed_image);
		fextractor->SetExtractionRegion(roi);
		fextractor->Update();
		cr.old_speed_image = fextractor->GetOutput(); //Let's see what the old speed imge was
		cr.valid = true;
		return cr;
	}

	return cr;

}

RegionType ChangeProbabilityLeakageCorrector::compute_roi(Segment seg, const PropagationInfo prop) {
	RegionType reg1 = seg.get_bbox();
	RegionType reg2 = seg.get_bbox(prop);

	l->infoStream() << "The bbox for the segment *without* the leaking prop is " << reg1.GetIndex() << " -> " << reg1.GetSize();
	l->infoStream() << "and with the leaking prop it becomes " << reg2.GetIndex() << " -> " << reg2.GetSize();

	return reg2;
}
