#include "manager.h"
#include "globals_itk.h"
#include "leakage.h"

#include <itkPointSetToImageFilter.h>
#include "itkPointSet.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkImageRegionConstIterator.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkSimplexMesh.h"
#include "itkSimplexMeshToTriangleMeshFilter.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"
#include "itkVector.h"
#include "itkQuadEdgeMesh.h"
#include "itkVTKPolyDataReader.h"

#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkVTKPolyDataWriter.h"
#include "itkImageFileReader.h"
#include "itkQuadEdgeMeshExtendedTraits.h"
#if ITK_VERSION_MAJOR >= 4
#include "itkNormalQuadEdgeMeshFilter.h"
#else
#include "itkQuadEdgeMeshNormalFilter.h"
#endif
#include <stdlib.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vcl_iostream.h>
#include <fstream>
#include <boost/lexical_cast.hpp>

IntImagePointer SegmentationManager::get_grayscale_image() {
	return this->propagator.read_int_image(this->params.input_gray_image_fname);
}

bool SegmentationManager::handle_change_request(ChangeRequest cr) {
	/*
	 * The leakage resolution process goes as follows:
	 * - The map is updated
	 * - The current segment is marked as leaking
	 * - Its parent is identified
	 * - All its propagations are pruned, except the first one
	 * - All its children are removed
	 * - The segment is regrown
	 * -
	 *
	 */

	pLogType l = get_logger("MANAGER::handle_leakage");

	std::string roi_gr_path = "cr_gray.mhd";
	std::string new_speed_map_path = "cr_speed.mhd";
	std::string old_speed_map_path = "cr_old_speed.mhd";
	std::string val;
	bool success = false;

	if (cr.valid) {

		l->info("Leakage correction code returned the following change request:");
		l->infoStream() << "ROI: " << cr.roi;
		/*
		this->propagator.write_image(cr.roi_grayscale, roi_gr_path);
		l->infoStream() << "The ROI portion of the original grayscale image is available at: " << roi_gr_path;
		this->propagator.write_image(cr.updated_speed_image, new_speed_map_path);
		l->infoStream() << "The updated values of the speed map are available at: " << new_speed_map_path;
		this->propagator.write_image(cr.old_speed_image, old_speed_map_path);
		l->infoStream() << "The ROI portion of the old speed map is available at: " << old_speed_map_path;
		 */
		l->info("So we'll try to MERGE the new speed map with the previous one.");

		this->propagator.update_speed_map(cr.updated_speed_image, cr.roi);
		this->speed_map_updated = true;

		this->curr_seg.leaking = true;
		//set_priority(log4cpp::Priority::DEBUG);

		this->propagator.clear_trials();
		// Can segment be grown?
		if (this->curr_seg.has_propagations()) {
			// Pass current wavefront
			propagator.set_wavefront(this->curr_seg.get_wavefront());
		} else {
			l->warnStream() << "Sorry segment has no propagations and thus cannot grown more...";
			this->curr_seg.seal();
		}

		success = true;
		/* TODO:
		 * This part implements the restart of the segmentation from the closest ancestor.
		 * To enable back only when restarting from last propagation works :)
		 *
		//Identify the parent segment
		int p = this->curr_seg.parentSegmentNumber;
		std::vector<PropagationInfo>::iterator propIt;
		if (p > 0) {//not the trachea
			Segment* parent = this->getNumberedSegment(p);
			if (parent != NULL) {
				propIt = parent->propagations.begin();
				l->infoStream() << "Removing all but first propagations from the parent segment " << p;
				++propIt; //Do not remove the first propagation or the segment won't grow anymore
				while(propIt != parent->propagations.end()) {
					this->propagator.clear_propagation_effects(*propIt);
					propIt = parent->propagations.erase(propIt);
				}

				//Find all the children and remove them
				//FIXME: we should not remove good segmentation results.
				std::vector<Segment>::iterator segIt = this->accepted_segments.begin();
				while(segIt != this->accepted_segments.end()) {
					if (segIt->parentSegmentNumber == p) {
						l->infoStream() << "Erasing segment " << segIt->segmentNumber;

						for( propIt = segIt->propagations.begin(); propIt != segIt->propagations.end();) {
							this->propagator.clear_propagation_effects(*propIt);
							propIt = segIt->propagations.erase(propIt);
						}
						segIt = this->accepted_segments.erase(segIt); //Update the iterator once the collection is modified

					} else
						++segIt;
				}
				segIt = this->segments_queue.begin();
				while(segIt != this->segments_queue.end()) {
					if (segIt->parentSegmentNumber == p) {
						l->infoStream() << "Erasing segment " << segIt->segmentNumber;
						for( propIt = segIt->propagations.begin(); propIt != segIt->propagations.end();) {
							this->propagator.clear_propagation_effects(*propIt);
							propIt = segIt->propagations.erase(propIt);
						}
						segIt = this->segments_queue.erase(segIt); //Update the iterator once the collection is modified

					} else
						++segIt;
				}

				//Put the parent segment on top of the queue to regrow it
				parent->unseal();
				this->segments_queue.push_back(*parent);

			} else {
				l->warnStream() << "No parent for segment " << this->curr_seg.segmentNumber;
			}

			success = true;

		} else {
			l->error("Don't know how to recover from a leakage in the trachea!");
			success = false;

		}


*/
	} else
		l->warn("Sorry, no idea on how to solve that!");

	//l->info("Press enter to continue after inspecting the change request");
	//std::cin >> val;

	return success;
}


bool SegmentationManager::handle_leakage(Segment seg, PropagationInfo prop, IntImagePointer grayscale_image, FloatImagePointer speed_image) {

	pLogType l = get_logger("MANAGER::handle_leakage");
	bool recoverable = false;
	if (this->params.leakage_recovery_enabled) {
		if (seg.leakage_count > 3) {
			l->infoStream() << "Too much leakages in this segment: I give it up!";
			recoverable = false;
		} else {

			/*
			l->infoStream() << "Saving debug image just before leakage correction...";
			this->propagator.save_internal_representation( this->output_label_image_fname );
			l->infoStream() << "Also save current segment image...";
			std::ostringstream oss;
			oss << "debug_seg_" << this->curr_seg.segmentNumber << ".mhd";
			std::string seg_image_path = oss.str();
			// Consider the segment with the additional propagation (used by the leakage code)
			this->curr_seg.propagations.push_back(prop);
			this->propagator.write_image(this->curr_seg.as_image(), seg_image_path);
			this->curr_seg.propagations.pop_back();
			*/
			LeakageResolutionRequest lrr = LeakageResolutionRequest(seg, prop, grayscale_image, speed_image);
			ChangeRequest cr = this->leakage_correction_strategy->correct(lrr);

			// Reconsider what to do
			this->propagator.clear_propagation_effects(prop);
			recoverable = this->handle_change_request(cr);
		}

	} else {
		l->warn("No correction strategy enabled: leakage is unrecoverable");
	}

	return recoverable;
}



bool SegmentationManager::check_cangrow(PropagationInfo* prop){
	/*
	 * A propagation can always grow if it has trial points
	 */
	if (prop->trials.empty()) {
		prop->status = PropagationInfo::REJ_NO_POINTS;
		return false;
	}

	return true;
}

bool SegmentationManager::check_leakage(PropagationInfo* prop, const PropagationInfo* prev, const Segment *seg, Segment *parent) {

	// FIXME: We should put the parent into the segment class instead of passing it externally

	pLogType l = get_logger("MGR::check_leakage");


	int st_trials = prop->trials.size();
	int st_l_points = prev->points.size();
	int st_l_trials = prev->trials.size();

	float prevRadius = prev->get_radius();
	float curRadius = prop->get_radius();

	// if previous contains no points, that's most likely because it's the first prop
	// in the segment (thus has only trials but no points)
	if (st_l_points == 0) {
		/*
		 * We cannot compare against a propagation which has grown no points!
		 * let's substitute with the father ones, if available
		 */
		if (parent != NULL) {
			st_l_trials = parent->get_last_propagation().trials.size();
		} else {
			// We accept the propagation because we have no info to judge about!

			/*
			 * FIXME: I'm not completely sure: we could always go on with the tests, no?
			 * For instance the trial size is always available once we are here.
			 */
			return false;
		}
	}




	if (st_trials > this->params.ST_TRIALS_TH*st_l_trials) {
		l->infoStream()  << "rejected because current trials (" << st_trials << ") > "<< this->params.ST_TRIALS_TH << "* prev points (" << st_l_trials<<")";
		prop->status = PropagationInfo::REJ_TRIAL_EXPLOSION;
		return true;
	}

	// check radius by getting minimum radius of parent
	/*
	 * FIXME: for now I went back to last segment only. As we enforce that from the beginning,
	 * it is always true that a child segment has a smalled radius and we don't need to go up
	 * in the ancestors.
	 */
	if (parent != NULL) {

	//	Rules: (van Ginneken et al '08)
	//	1. current radius must be smaller than 1.5*minRadius(parentSegment)

		if (curRadius > this->params.RADIUS_TH * parent->GetMinimumRadius()) {
			l->infoStream() << "rejected because current radius (" << curRadius << ") > "<< this->params.RADIUS_TH << "* minRadius of parent (" << parent->GetMinimumRadius()<< ")";
			prop->status = PropagationInfo::REJ_MINRADIUS;
			return true;
		}

	//	2. front should not touch any other segment
		l->warnStream() << "2. front should not touch any other segment - NOT implemented!";

	//	3. length should not be more than 5*radius
		l->warnStream() << "3. length should not be more than 5*radius - NOT implemented!";

	//	4. avg radii of 2 consecutive wf shouldn't differ by more than 10%
		if (curRadius / prevRadius > this->params.RADIUS_RATIO_TH ) {
			l->infoStream() << "rejected because current radius (" << curRadius << ") > " << this->params.RADIUS_RATIO_TH << "* prevRadius (" << prevRadius<<")";
			prop->status = PropagationInfo::REJ_PREVRADIUS;
			return true;

		}
	}

	// If we are here, we didn't detected any leakage!
	return false;
}

void SegmentationManager::do_segmentation(){

 /*
  * Basic idea:
  * 	- Create a new segment to hold the trachea
  * 	- Init the fast marching propagation method with a seed from the user
  * 	- Put this segment to the segment queue
  * 	- For all the segments in the queue
  * 		- curr_seg <- queue.pop()
  * 		- grow(curr_seg)
  * 		- evaluate(curr_seg)
  * 		- decide if continue growing, stopping segment or rejecting it
  */

   pLogType l = get_logger("MGR::do_segmentation");
  
  /*
   * Initialize propagator
   */


  this->speed_map_updated = false; //Has the speed map been modified?
  this->propagator.set_input_image(this->params.input_speed_image_fname);
  this->propagator.set_stop_time(this->params.stop_time);
  this->propagator.set_timestep(this->params.timestep);
  this->propagator.initialize();

  if (this->params.leakage_recovery_enabled)
	  this->grayscale_image = this->get_grayscale_image(); //FIXME: it is a method! let's call it load_g_image and do its work

  const int STABILIZATION_TH = propagator.get_stabilization_threshold(this->params.TRACHEA_TYPICAL_DIMENSION);
  l->infoStream() << "Propagations will be considered stabilized after " << STABILIZATION_TH << " secs of propagation time";
  bool stabilized = false;

  this->airwayNumber = 1;

  Segment trachea = Segment(airwayNumber, -1);
  PropagationInfo first_wf = PropagationInfo();
  first_wf.trials.push_back(this->params.seed);
  trachea.accept_propagation(first_wf);

  this->segments_queue.push_back(trachea);


  this->maxSegmentNumber = airwayNumber;
  while(!this->segments_queue.empty())
     {
		// get segment from queue
		this->curr_seg = this->segments_queue.front();
		this->segments_queue.erase(this->segments_queue.begin());


		l->infoStream() << "Growing segment [" << this->curr_seg.segmentNumber << "," << this->curr_seg.parentSegmentNumber << "] (still in queue: " << segments_queue.size() << ")";

		// Can segment be grown?
		if (this->curr_seg.has_propagations()) {
			// Pass current wavefront
			this->propagator.set_wavefront(this->curr_seg.get_wavefront());
		} else {
			l->warnStream() << "Sorry segment has no propagations and thus cannot grown more...";
			this->curr_seg.seal();
			continue;
		}


		bool valid = this->grow_curr_segment(stabilized, STABILIZATION_TH);

		/*
		 * FIXME: for now let's accept this. Then we will merge the leakage analysis here
		 */
		if (valid) {
			// accept last segment
			this->accepted_segments.push_back(this->curr_seg);
		}

     }


}

bool SegmentationManager::grow_curr_segment(bool stabilized, const int STABILIZATION_TH) {
	pLogType l = get_logger("MGR::grow_current_segment");

	bool valid = true;
	bool leak_rec = false;

	int count = 0;
	while(this->propagator.get_propagation_time() < this->propagator.stop_time)  {

		if (this->curr_seg.is_sealed()) {
			l->infoStream() << "Current segment has been sealed. Choosing next in the queue...";
			// Get out of the propagation loop and choose another segment
			break;
		}


		if (not stabilized) {

			if (this->propagator.get_propagation_time() >= STABILIZATION_TH) {
				l->warnStream() << "Wavefront should be stable by now...";
				stabilized = true;
				this->propagator.set_timestep(this->params.timestep);

				if (this->propagator.get_propagation_time() > this->params.stop_time) {
					l->warnStream() << "You asked to stop at time: " << this->params.stop_time << " but FM is already at " << this->propagator.get_propagation_time();
					l->warnStream() << "I'll stop instead at: " << this->propagator.get_propagation_time() + this->params.timestep;
					this->params.stop_time = this->propagator.get_propagation_time() + this->params.timestep;
				}
				this->propagator.set_stop_time(this->params.stop_time);

			} else {
				l->warnStream() << "Wavefront still stabilizing...";
				this->propagator.set_timestep(STABILIZATION_TH);
				this->propagator.set_stop_time(STABILIZATION_TH + 1.0);

			}

		}


		l->infoStream() << "About to grow: timestep:  " << propagator.timestep << " stop_time: " << propagator.stop_time;

		// Make the actual propagation
		PropagationInfo prop = propagator.grow();
		count++;

		// Now populate the propagation info
		prop.currSegmentNum = this->curr_seg.segmentNumber;
		prop.parentSegmentNum = this->curr_seg.parentSegmentNumber;
		l->infoStream() << "prop's trial size: " << prop.trials.size() << ", points size: "<< prop.points.size();


		/*
		 * Entering propagation check section. Here we do:
		 * 1) Validity of the propagation. If the propagation has no new trials, it cannot
		 * 		grow anymore. The prop is accepted and the segment sealed
		 * 2) Leakage check. In case of leakage, the last propagation is reversed
		 * 		and the current segment sealed.
		 * 3) Bifurcation check. If the last prop is disconnected, for each disc region
		 * 		we create a new segment. We seal the last one
		 *
		 * If a propagation pass all the previous checks it is a regular propagation and we accept it
		 */

		//1st check
		// WARN: It exits early!
		bool cangrow = this->check_cangrow(&prop);
		if (not cangrow) {
			/*
			 * It is usually ok if you are in smaller branch. If you think this is an error,
			 * maybe you should check if the speed image supports growing of the past wavefront.
			 */
			l->warnStream() <<  "Sorry propagation has produced no trial points and thus we cannot grown anymore...";

			this->curr_seg.accept_propagation(prop);
			break;

		} // end 1st check


		/*
		 * We moved here the accepting of the propagation because we do not want to accept propagation which
		 * can't grow even we are stabilizing!
		 */
		if (not stabilized) {
			l->warnStream() << "Unconditionally accepting last propagation because wavefront is not stable";
			this->curr_seg.accept_propagation(prop);
			continue;
		}


		/*
		 * 2nd: check for leakage, if enabled (visitor will not use it by default!)
		 * This is useful because often constraints that enables a prompt leakage detection and
		 * thus facilitates the recovery hamper the detection of branched in the visitor.
		 * As the visitor is not going to do leakage recovery anyway, the checks are quite pointless
		 * and we can make them stricter for the segmenter instead.
		 */
		if (this->params.leakage_detection_enabled) {
			PropagationInfo last = this->curr_seg.get_last_propagation();
			// if not the first, it should have a parent segment
			Segment * parent = (this->curr_seg.segmentNumber > 1)? getNumberedSegment(this->curr_seg.parentSegmentNumber): NULL;
			bool leaking = check_leakage(&prop, &last, &this->curr_seg, parent);

			if (leaking) {
				/*
				 * Once a leakage is found, the leakage recovery strategy is run, if enabled.
				 */

				this->curr_seg.leakage_count += 1;
				bool recoverable = handle_leakage(this->curr_seg, prop, this->grayscale_image, this->propagator.input_image);

				//If no recovery is possible, default to stop here the prop
				if (not recoverable) {
					l->warn("Leakage is unrecoverable...");
					this->propagator.clear_propagation_effects(prop);
					this->curr_seg.accept_propagation(prop);
					break;
				} else {
					leak_rec = true;
					l->info("Trying to recover from leakage...");

					continue;
				}


			}
		}
		 // end leakage check


		// 3rd: is bifurcation detected?

		bool is_connected = true;
		std::vector< std::vector<PointType> > regions;

		regions =  bif_checker.split_connected(prop);
		is_connected = (regions.size() == 1);

		if (!is_connected) {

			/*
			 * - Split wavefront into connected pieces
			 * - for each piece create a new segment with this as a starting wavefront
			 * - Put segments in the queue - if OK
			 */


			// i think we need to add _this_ disconnected wavefront to the current segment, else
			// we'll have a 'missing' zone between the current segment and the next

			std::vector< std::vector<PointType> >::iterator it;
			std::vector< PointType >::iterator it_wf;
			std::vector< std::vector<PointType> > discarded;

			l->infoStream() << "Found " << regions.size() << " disconnected regions";
			std::vector<PointType> last_region;
			int region_count = 0;
			for ( it = regions.begin() ; it != regions.end(); it++ ) {
				int rsize = it->size();
				l->infoStream() << "Region " << region_count << ": (size = " << rsize << ") - ";

				// use the last propagation instead of the current wavefront as explosion also affects this
				if (rsize < (0.15 * prop.trials.size())) {
					l->infoStream() << "discarding too small region with size: " << rsize << " ("<< prop.trials.size() <<")";
					discarded.push_back(*it);
					continue;
				}
				region_count++;
				last_region = *it;
			}

			/*
			 * It might happen that, even if disconnected regions are found, only one is actually fit to create a new segment.
			 * In this case just merge it into the parent one and continue
			 */
			if(region_count < 2) {
				l->infoStream() << "Propagation is disconnected, but has only generated one segment: merging it into the parent...";
				this->curr_seg.accept_propagation(prop); //We accept the current one as if it was a normal prop

				// Let's now create a new connected wavefront for the segment
				std::vector< PointType >::iterator it_wf;
				PropagationInfo new_wf;
				for ( it_wf = last_region.begin() ; it_wf < last_region.end(); it_wf++ )
					new_wf.trials.push_back(*it_wf);

				this->curr_seg.accept_propagation(new_wf);

				//FIXME: discarded regions should have their points included somewhere... Maybe in last prop.points?

				// Continue with next propagation. But curr_seg has been already modified, so don't re-accept twice last prop
				continue;

			} else {

				l->infoStream() << "Propagation is disconnected";
				prop.status = PropagationInfo::DISCONNECTED;
				this->curr_seg.accept_propagation(prop); //WARNING: accepting a prop with a bifurcation will seal the segment!

				int region_count = 0;
				for ( it = regions.begin() ; it != regions.end(); it++ ) {
					int rsize = it->size();

					// use the last propagation instead of the current wavefront as explosion also affects this
					if (rsize < (0.15 * prop.trials.size())) {
						continue;
					}
					region_count++;

					// Now create new segment
					Segment seg = Segment(++maxSegmentNumber, this->curr_seg.segmentNumber);

					// Setup the disconnected region as it's first wavefront
					PropagationInfo new_wf;

					for ( it_wf = it->begin() ; it_wf < it->end(); it_wf++ )
						new_wf.trials.push_back(*it_wf);

					// Give the new segment its wavefront
					seg.accept_propagation(new_wf);


					// Now push the newly created segment to the queue
					segments_queue.push_back(seg);
				}


				// Continue with next propagation. But curr_seg has been already modified, so don't re-accept twice last prop
				// FIXME: maybe we should just change to continue and seal the segment.
				break;
			}


		}

		// If we are here we always accept the prop and let the segment class choose what to do
		this->curr_seg.accept_propagation(prop);


	}

	l->infoStream() << "Exiting at: " << this->propagator.get_propagation_time() << " " << this->propagator.stop_time;

	return valid;

}

void SegmentationManager::write_propagations(const std::string path) {
	pLogType l = get_logger("write_propagations");
	std::ofstream  data(path.c_str());
	std::string sep = ", ";

	 //[SegmentID,ParentSegmentID,PropagationNum,Points]
	data << "SegmentID" << sep << "ParentSegmentID" << sep << "Propagation#" << sep << "Points" << sep << "Centroid" << sep << "Normals" << std::endl;
	std::vector<Segment>::iterator segIt = this->accepted_segments.begin();

	typedef double        CoordType;
	typedef itk::Mesh< float, LDimension >   MeshType;
	typedef itk::PointSet< CoordType, LDimension > PointSetType;
	typedef itk::QuadEdgeMesh< CoordType, LDimension > InputMeshType;
	typedef itk::BinaryMask3DMeshSource< UCharImageType, MeshType >   MeshSourceType;
	typedef itk::Vector< CoordType, LDimension > VectorType;

	typedef itk::QuadEdgeMeshExtendedTraits <
			VectorType,
			LDimension,
			2,
			CoordType,
			CoordType,
			VectorType,
			bool,
			bool > Traits;
	typedef itk::QuadEdgeMesh < VectorType, LDimension, Traits > OutputMeshType;
	#if ITK_VERSION_MAJOR >= 4
	  typedef itk::NormalQuadEdgeMeshFilter< InputMeshType, OutputMeshType > NormalFilterType;
	#else
	  typedef itk::QuadEdgeMeshNormalFilter< InputMeshType, OutputMeshType > NormalFilterType;
	#endif

	// original idea is to supply binary image & in-value

	while(segIt != this->accepted_segments.end()) {
		Segment seg = *segIt;
		int pcount = 0;
		std::vector<PropagationInfo>::iterator pIt;
		if (seg.get_points_copy().empty()) {
			data << seg.segmentNumber << sep << seg.parentSegmentNumber << sep << "0" << sep << "'[]'" << sep << "[]" << sep << "[Nan, Nan, Nan]" << std::endl ;
			++segIt;
			continue;
		}
		std::string allcentroids = "";
		std::string o_norm = "";
		std::vector<float> norm, norm2;
		std::vector<double> lastcentroid;
		for (pIt=seg.propagations.begin(); pIt!=seg.propagations.end(); pIt++){

			data << seg.segmentNumber << sep << seg.parentSegmentNumber << sep;
			allcentroids = "";
			if (pIt->trials.empty()) {
				data << pcount << sep << "'[]'" << sep << "[]" << sep << "[Nan, Nan, Nan]" << std::endl ;
				continue;
			}
			std::vector<double> centroidV = pIt->get_centroid_index();
			std::string centroid = "[" + boost::lexical_cast<std::string>( centroidV[0]  ) + "," + boost::lexical_cast<std::string>( centroidV[1]  ) + "," + boost::lexical_cast<std::string>( centroidV[2] ) + "]";

			std::vector<PointType>::const_iterator point;
			allcentroids += "'(";
			PointType t;
			int pct = 0;
			PointSetType::Pointer  pointsSet = PointSetType::New();
			typedef PointSetType::PointType PPointType;
			PPointType ppt0;

			for (point=pIt->trials.begin(); point!=pIt->trials.end(); point++) {
				if (point != pIt->trials.begin()) allcentroids += ", ";
				t = *point;
				ppt0[0] = t[0];
				ppt0[1] = t[1];
				ppt0[2] = t[2];
				pointsSet->SetPoint( pct, ppt0 ); // FIXME to see if PointSetType::PointType == Point::PointType
				++pct;

//				allcentroids += "(" + boost::lexical_cast<std::string>( (*point)[0]  ) + "," + boost::lexical_cast<std::string>( (*point)[1]  ) + "," + boost::lexical_cast<std::string>( (*point)[2] ) + ")";
				allcentroids += boost::lexical_cast<std::string> (t) ;
			}
			allcentroids += ")'" ;
			std::cout << "# of points: " << pointsSet->GetNumberOfPoints() << "\n";
//			UCharPixelType objectValue = 1;

			if (pct > 3) { // assumes minimum of 3 points needed to create mesh
				const UCharPixelType objectValue = static_cast<UCharPixelType>( 1 );
				UCharImageType::Pointer OutputImage = UCharImageType::New();

				OutputImage->SetRegions(this->propagator.input_image->GetLargestPossibleRegion());
				OutputImage->Allocate();   // allocate the image
				OutputImage->FillBuffer(0);

				typedef typename PointSetType::PointsContainer::ConstIterator  PointIterator;
				PointIterator pointItr = pointsSet->GetPoints()->Begin();
				PointIterator pointEnd = pointsSet->GetPoints()->End();

				typename UCharImageType::IndexType index;

				int count = 0;
				while( pointItr != pointEnd )
				{
					if(OutputImage->TransformPhysicalPointToIndex(pointItr.Value(),index))
					{
						OutputImage->SetPixel(index,1);
						count++;
					}
					else {
						std::cout << "Point " << pointItr.Value() << " --> " << index << " is outside image region" << std::endl;
					}
					pointItr++;
				}
				std::cout << "OutputImage (setting "<< count <<" pts)... done\n";

				MeshSourceType::Pointer meshSource = MeshSourceType::New();
				meshSource->SetObjectValue( objectValue );
				meshSource->SetInput( OutputImage );
				try
				{
					meshSource->Update();
				}
				catch( itk::ExceptionObject & exp )
				{
					std::cerr << exp << std::endl;
				}
//				std::cout << "meshSource... done\n";

				typedef itk::VTKPolyDataWriter<MeshType> VTKWriterType;
				VTKWriterType::Pointer writer = VTKWriterType::New();
				writer->SetInput(meshSource->GetOutput());
				writer->SetFileName("test.vtk");
				writer->Update();

				typedef itk::VTKPolyDataReader< InputMeshType > VTKReaderType;
				VTKReaderType::Pointer reader = VTKReaderType::New( );
				reader->SetFileName( "test.vtk" );

				try
				{
					reader->Update( );
				}
				catch( itk::ExceptionObject & excp )
				{
					std::cerr << "Exception thrown while reading the input file " << std::endl;
					std::cerr << excp << std::endl;
				}
				InputMeshType::Pointer mesh = reader->GetOutput( );
				NormalFilterType::Pointer normals = NormalFilterType::New( );
				NormalFilterType::WeightType weight_type = NormalFilterType::GOURAUD;
				normals->SetInput( mesh );
				normals->SetWeight( weight_type );

				try
				{
					normals->Update( );
				}
				catch( itk::ExceptionObject & excp )
				{
					std::cerr << excp << std::endl;
				}
				OutputMeshType::Pointer output = normals->GetOutput( );

				//  std::cout << normals << std::endl;

				typedef itk::Vector< double, 3 > InputVectorType;
//				InputVectorType norm_values;
//				norm_values[0] = 0.0;
//				norm_values[1] = 0.0;
//				norm_values[2] = 0.0;
				//			    std::cout << normals << "\n";
//				OutputMeshType::PointDataContainerPointer pointdata = output->GetPointData( );
				std::cout <<"Face Normal" <<std::endl;

				// mean of face normal
				InputVectorType meanFaceNorm;
				meanFaceNorm[0] = 0.0;
				meanFaceNorm[1] = 0.0;
				meanFaceNorm[2] = 0.0;

				int faceCount = 0;
				OutputMeshType::CellDataContainerPointer celldata = output->GetCellData( );
				for( OutputMeshType::CellDataContainerIterator n_it = celldata->Begin( );
						n_it != celldata->End( );
						n_it++ )
				{
					//	  std::cout <<n_it->Index( ) <<"  " <<n_it->Value( ) <<std::endl;
					//	  meanFaceNorm += n_it->Value();
					meanFaceNorm[0] += n_it->Value()[0];
					meanFaceNorm[1] += n_it->Value()[1];
					meanFaceNorm[2] += n_it->Value()[2];
					faceCount++;
				}

				meanFaceNorm[0] = meanFaceNorm[0]*1.0 / faceCount;
				meanFaceNorm[1] = meanFaceNorm[1]*1.0 / faceCount;
				meanFaceNorm[2] = meanFaceNorm[2]*1.0 / faceCount;
				data << pcount << sep << allcentroids << sep << centroid << sep << meanFaceNorm << std::endl;
			} else {
				std::cout << "Too little points ("<< pct <<") to create mesh\n";
				data << pcount << sep << allcentroids << sep << centroid << sep << "[Nan, Nan, Nan]" << std::endl;
			}
			++pcount;
//			data << pcount << sep << allcentroids << sep << norm_values << std::endl;
		}
		++segIt;
	}

	data.close();

}

void SegmentationManager::write_segments_tree(const std::string path) {
	/*
	 * Prints extended statistics on segmentation
	 */

	std::cout << this->propagator.current_output->GetSpacing() << std::endl;
	pLogType l = get_logger("write_segments_tree");
	 std::ofstream  data(path.c_str());
	 std::string sep = ", ";
	data << "SegId" << sep << "ParentId" << sep << "PropId" << sep << "c_x" << sep << "c_y" << sep << "c_z" << sep << "Radius" << std::endl;
	l->debugStream() << "SegId" << sep << "ParentId" << sep << "PropId" << sep << "c_x" << sep << "c_y" << sep << "c_z" << sep << "Radius";
	std::vector<Segment>::iterator segIt = this->accepted_segments.begin();
	while(segIt != this->accepted_segments.end()) {
		Segment seg = *segIt;
		int pcount = 0;
		std::vector<PropagationInfo>::iterator pIt;
		for (pIt=seg.propagations.begin(); pIt!=seg.propagations.end(); pIt++){
			if (pIt->trials.empty())
				continue; //Do not write propagation with no trials (centroid would be 0,0,0)

			/* FIXME: Add the conversion information to allow index2point conversion in the get_centroid method.
			 * this is a very bad fix and should be changed soon in the main code.
			 */

			pIt->space_converter = this->propagator.current_output;

			data << seg.segmentNumber << sep << seg.parentSegmentNumber << sep;
			data << pcount << sep;
			//std::cout << "Spacing: " << pIt->space_converter->GetSpacing() << std::endl;
			std::vector<double> centroid = pIt->get_centroid();
			//std::cout << "Spacing: " << pIt->space_converter->GetSpacing() << std::endl;
			data << centroid[0] << sep << centroid[1] << sep << centroid[2] << sep;
			data << pIt->get_physical_radius();
			data << std::endl;

			l->debugStream() << seg.segmentNumber << sep << seg.parentSegmentNumber << sep << pcount << sep << centroid[0] << ":" << centroid[1] << ":" << centroid[2] << sep << pIt->get_physical_radius();
			// Remove the converter after using it because it is reference counted
			pIt->space_converter = NULL;

			++pcount;
		}


		++segIt;
	}

	data.close();

}

void SegmentationManager::write_segments_tree_oldformat(const std::string path) {

	pLogType l = get_logger("write_segments_tree_oldformat");
	 std::ofstream  data(path.c_str());
	 std::string sep = ", ";
	 //[SegmentID,ParentSegmentID,Length,Radius,PointCount,Centroids,OrientationNorm,FinalNorm]
	data << "SegmentID" << sep << "ParentSegmentID" << sep << "Length" << sep << "Radius" << sep << "PointCount" << sep << "Centroids" << sep << "OrientationNorm" << sep << "FinalNorm"   << std::endl;
	l->debugStream() << "SegmentID" << sep << "ParentSegmentID" << sep << "Length" << sep << "Radius" << sep << "PointCount" << sep << "Centroids" << sep << "OrientationNorm" << sep << "FinalNorm" ;
	std::vector<Segment>::iterator segIt = this->accepted_segments.begin();
	std::map< UCharPixelType, std::vector<PointType> > lsq;

	std::cout << "W1\n";
	//	1,-1,93.1003959193,448.7,10,"['[260.002096436, 237.916142558, 4.4821802935]', '[261.850855746, 245.740831296, 15.3740831296]', '[261.197802198, 250.145054945, 25.043956044]', '[259.602459016, 253.944672131, 34.6393442623]', '[257.424307036, 257.620469083, 44.1108742004]', '[253.868020305, 261.672588832, 53.1116751269]', '[250.855555556, 265.580555556, 61.9222222222]', '[250.863636364, 269.235930736, 71.0584415584]', '[251.245009074, 271.923774955, 80.1542649728]', '[249.82464455, 272.018957346, 88.1895734597]']","(-0.10525300089, 0.415061098923, 0.903684729297)","(-0.117134326457, 0.323140082497, 0.939074031506)"
	while(segIt != this->accepted_segments.end()) {
		Segment seg = *segIt;
		int pcount = 0;
		std::vector<PropagationInfo>::iterator pIt;
		data << seg.segmentNumber << sep << seg.parentSegmentNumber << sep;
		if (seg.get_points_copy().empty()) {
			data << "0" << sep << "0" << sep << "0" << sep << "\"[]\"" << sep << "(0,0,0)" << sep << "(0,0,0)" << std::endl ;
			++segIt;
			continue;
		}
		//lsq[pcount] = lsqfit_3d(seg.get_points_copy());
		// calc length of this lsq
		//double length = calc_length(lsq[pcount]);
		double length = 0;
		data << length << sep;
		double radius = 0.0;
		std::string allcentroids = "";
		std::string o_norm = "";
		std::vector<float> norm, norm2;
		std::vector<double> lastcentroid;
		for (pIt=seg.propagations.begin(); pIt!=seg.propagations.end(); pIt++){

			if (pIt->trials.empty())
				continue;

			std::vector<double> centroid = pIt->get_centroid_index();
			allcentroids += "'[" + boost::lexical_cast<std::string>( centroid[0]  ) + "," + boost::lexical_cast<std::string>( centroid[1]  ) + "," + boost::lexical_cast<std::string>( centroid[2] ) + "]'" ;
			radius += pIt->get_radius();
			++pcount;

			if (pcount < seg.propagations.size())
				allcentroids += sep;

			if (pcount==1) {
				norm.push_back(centroid[0]);
				norm.push_back(centroid[1]);
				norm.push_back(centroid[2]);
			}
			// we dont know what's the last point
			lastcentroid.push_back(centroid[0]);
			lastcentroid.push_back(centroid[1]);
			lastcentroid.push_back(centroid[2]);
		}
		int s = lastcentroid.size();
		if (pcount > 1) {
			norm[0] = lastcentroid[s-3]-norm[0];
			norm[1] = lastcentroid[s-2]-norm[1];
			norm[2] = lastcentroid[s-1]-norm[2];
		} else  {
			norm[0] = 0;
			norm[1] = 0;
			norm[2] = 0;
		}
		if (allcentroids.length() > 2) {
			allcentroids.erase(allcentroids.end() - 2, allcentroids.end() );
		}
		radius = 0; // FIXME right now at 0
		data << radius << sep;
		data << pcount << sep;
		data << "\"[" << allcentroids << "]\"" << sep;
		o_norm = "("+boost::lexical_cast<std::string>(norm[0])+","+boost::lexical_cast<std::string>(norm[1])+","+boost::lexical_cast<std::string>(norm[2])+")";
		data << o_norm << sep;
		data << o_norm ;
		data << std::endl;
		++segIt;
	}

	data.close();

}


Segment * SegmentationManager::getNumberedSegment(int number) {
	for (unsigned int i=0; i<accepted_segments.size(); i++)
		if (accepted_segments.at(i).segmentNumber == number)
			return &accepted_segments.at(i);
	return NULL;
}


std::map<UCharPixelType,std::vector<PointType> > SegmentationManager::get_color_map() {
	std::map<UCharPixelType,std::vector<PointType> > colormap;
	std::vector<Segment>::iterator it;
	//FIXME: We does it need  a copy of the points? He can always iterate it from the segment!
	for ( it = accepted_segments.begin() ; it < accepted_segments.end(); it++ )
		colormap[it->segmentNumber] = it->get_points_copy();

	return colormap;


}


void SegmentationManager::write_debug_images() {
	std::map<UCharPixelType,std::vector<PointType> > colormap, clean;

	std::cout << "Generating debug colormap..." << std::endl;
	colormap=this->get_color_map();
	this->propagator.color_segmentation("debug-colormap.tiff",colormap);


}

void SegmentationManager::save_color_segments() {
	std::map<UCharPixelType,std::vector<PointType> > colormap, clean;

	std::cout << "Generating colormap..." << std::endl;
	colormap=this->get_color_map();

	std::cout << "Generating clean colormap..." << std::endl;
	/*
	 * Use mhd format instead of tiff and dcm because it correcly save both the x-y and z spacing
	 * Moreover, gdcm is still not able to save dicom multiframe CT image correctly...
	 */
	this->propagator.clean_color_segmentation("segments_colormap.mhd",colormap);


	//TODO: Reenable statistics
	//clean = this->propagator.get_cleaned_colormap(colormap);
	//this->get_statistics(clean);
}

void SegmentationManager::print_segments() {

	std::cout << "Id \t Par \t #Prop \t #Points \t Status \t Status of last prop" << std::endl;
	std::vector<Segment>::iterator it;

	for ( it = accepted_segments.begin() ; it < accepted_segments.end(); it++ ) {

		Segment seg = *it;

		std::cout << seg.segmentNumber << " \t " << seg.parentSegmentNumber << " \t " << seg.get_propagations_size() << " \t " << seg.get_size() << " \t " << seg.get_status() << " \t " << seg.get_last_propagation().get_status() <<  std::endl;

	}


}


double SegmentationManager::calc_angleFromParent(std::vector<PointType> parentpoints, std::vector<PointType> points) {

	double angle = 0.0;
	if (points.empty() || parentpoints.empty()) return -1.0;
	double mag1, mag2;
	DirVectorType v1, v2;
	v1 = parentpoints.at(0) - parentpoints.at(parentpoints.size()-1);
	v2 = points.at(0) - points.at(points.size()-1);
	mag1 = sqrt(pow(abs(v1[0]), 2) + pow(abs(v1[1]), 2) + pow(abs(v1[2]), 2));
	mag2 = sqrt(pow(abs(v2[0]), 2) + pow(abs(v2[1]), 2) + pow(abs(v2[2]), 2));
	double dot = (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]);
// v1•v2 = |v1||v2| cos(angle)
	angle = vcl_acos(dot / mag1*mag2);

	return angle;
}

double SegmentationManager::calc_length(std::vector<PointType> input) {
	double length = 0.0;
	if (input.size()==0) return length;

//	Index min_idx, max_idx;
//	min_idx.Fill(INT_MAX);
//	max_idx.Fill(0);
//
//	std::vector< PointType >::iterator it2;
//	for (it2 = input.begin(); it2 != input.end(); ++it2) {
//		PointType p = (*it2);
//		for (unsigned int i=0; i<min_idx.GetIndexDimension();i++) {
//			min_idx[i] = std::min((int) min_idx[i],(int) p[i]);
//			max_idx[i] = std::max((int) max_idx[i], (int) p[i]);
//		}
//	}
	// length = sqrt(x² * y² + z²)
//	length = pow(max_idx[0] - min_idx[0], 2);
//	length += pow(max_idx[1] - min_idx[1], 2);
//	length += pow(max_idx[2] - min_idx[2], 2);

	PointType max, min;
	max = input.at(0);
	min = input.at(input.size()-1);
//	std::cout << "max, min : " << max << ", "<< min << "\n";
	length = vnl_math_sqr(max[0] - min[0]) + vnl_math_sqr(max[1] - min[1]) + vnl_math_sqr(max[2] - min[2]);
	length = sqrt(length);

//	length = pow(abs(max[0] - min[0]), 2);
//	length += pow(abs(max[1] - min[1]), 2);
//	length += pow(abs(max[2] - min[2]), 2);
//	length = sqrt(length);

	return length;

}

std::vector<PointType> SegmentationManager::lsqfit_3d(std::vector<PointType> input) {

	std::vector<PointType> result;

	std::vector<int> x, y, z;
	int n = input.size();
	double meanx = 0, meany = 0, meanz = 0;
	//	m = [mean(x),mean(y),mean(z)];
	int i = 0;
	for (i=0 ; i<n ; i++) {
		meanx += input.at(i)[0];
		meany += input.at(i)[1];
		meanz += input.at(i)[2];
		x.push_back(input.at(i)[0]);
		y.push_back(input.at(i)[1]);
		z.push_back(input.at(i)[2]);
	}
	std::cout << "Pushed "<< i << " times\n";
	meanx /= n;
	meany /= n;
	meanz /= n;
	vnl_matrix<double> W(n,3);
//	std::vector<PointType> w;
//	PointType p;
	for (int i=0 ; i<n; i++) {
//		p[0] = x.at(i) - meanx;
//		p[1] = y.at(i) - meany;
//		p[2] = z.at(i) - meanz;
//		w.push_back(p);
		W(i, 0) = x.at(i) - meanx;
		W(i, 1) = y.at(i) - meany;
		W(i, 2) = z.at(i) - meanz;
	}

//	[n,mx] = size(x); [ny,my] = size(y); [nz,mz] = size(z);
//	if (mx~=1)|(my~=1)|(mz~=1)|(ny~=n)|(nz~=n)
//	 error('The arguments must be column vectors of the same length.')
//	end

//	w = [x-m(1),y-m(2),z-m(3)]; % Use "mean" point as base
//	a = (1/n)*w'*w; // 'a' is a positive definite matrix
	vnl_matrix<double> A(3,3);
	double factor = 1.0 /n;
//	double factor = 0.0012;
	std::cout << "factor: "<< factor << "\n";
	A = factor * W.transpose() * W;
//	A = element_product(factor, W.transpose() * W);
	vnl_svd<double> svd(A);
	vnl_vector<double> p = svd.nullvector();
	vcl_cout << "Mean values: "<< meanx << ", " << meany << ", " << meanz << vcl_endl;
	vcl_cout << "Method1 Eigenvector: "<< p << vcl_endl;
	vcl_cout << "vnl_svd<double> residual = " << (A * p).magnitude() << vcl_endl;

	vnl_symmetric_eigensystem<double>  eig(A.transpose() * A);
	vnl_vector<double> p2 = eig.get_eigenvector(0);
	vcl_cout << "Method2 Eigenvector: "<< p << vcl_endl;
	vcl_cout << "Method2 Eig residual = " << (A * p2).magnitude() << vcl_endl;


//	[u,d,v] = svd(a); % 'eig' & 'svd' get same eigenvalues for this matrix
//	p = u(:,1)'; % Get eigenvector for largest eigenvalue
//	s = d(2,2)+d(3,3); % Sum the other two eigenvalues

	// results? (x,y,z) = m + t*p
	// determine range of t using z as reference (min & max)
	std::vector<int>::iterator where = std::max_element (z.begin(), z.end());
//	std::cout << "maximum z " << (*where) << std::endl;
	double min_t = ((*where) - meanz) / p[2];
	where = std::min_element (z.begin(), z.end());
//	std::cout << "minimum z " << (*where) << std::endl;
	double max_t = (*where - meanz) / p[2];
// check that min_t is smaller than max_t. if not, interchange!
	if (max_t < min_t) {
		double tmp = max_t;
		max_t = min_t;
		min_t = tmp;
	}

	std::cout << "min_t " << min_t << ", max_t" << max_t << ", floor " << floor(min_t) << ", ceil" << ceil(max_t) << std::endl;

	PointType rp;
	for (int i=floor(min_t); i<ceil(max_t); i++) {
		rp[0] = meanx + i*p[0];
		rp[1] = meany + i*p[1];
		rp[2] = meanz + i*p[2];
//		std::cout << "i: " << i << "-> "<<  rp << ", ";
		result.push_back(rp);
	}
	std::cout << "\n";
	return result;
}

void SegmentationManager::get_statistics(std::map<UCharPixelType,std::vector<PointType> >clean) {

	// accepted_segments contain points which are already marked as out. clean contains only those marked in by FM

	// global variables
	int g_tl = 0, g_bc;
	// segmental variables
	double s_length, s_count, s_compactness;
	DirVectorType s_dir;
// GLOBAL - branch size
	g_bc = clean.size();
	std::map< UCharPixelType, std::vector<PointType> > lsq;
// GLOBAL - tree length: calc length of each segment
	for (unsigned int i=0; i<accepted_segments.size(); i++) {
		std::cout << i << ":" << clean[i].size() << "\n";

		int proptotal = 0, proptotalaccepted = 0;
		// calc center of each propagation
		for (unsigned int j=0; j<accepted_segments.at(i).propagations.size(); j++) {
			PropagationInfo pif = accepted_segments.at(i).propagations.at(j);
			std::cout << "Prop " << j << ", status: "<< pif.get_status() <<", size: " << pif.points.size() << "\n";
			if (pif.get_status() == "ACCEPTED") proptotalaccepted += pif.points.size();
			proptotal += pif.points.size();
		}

//		std::cout << "prop total: " << accepted_segments.at(i).get_points_copy().size() << "\n";
		std::cout << "prop total: " << proptotal << ", proptotalaccepted: "<< proptotalaccepted <<"\n";
		// do a least_square fit the center to a line in 3d
		lsq[i] = lsqfit_3d(clean[i]);
		// calc length of this lsq
		s_length = calc_length(lsq[i]);
	}

	std::cout << "GLOBAL STATS\n";
	std::cout << "Segment Count: " << g_bc << "\n";
	std::cout << "Tree Length: " << g_tl << "\n" ;

	std::cout << "SEGMENT STATS\n";

//	std::map<UCharPixelType,std::vector<PointType> >::const_iterator it;
//	for ( it=clean.begin() ; it != clean.end(); it++ ) {
//
//		std::cout << "Segment "<< (int) (*it).first <<": ";
//		std::cout << "Point Count: " << (*it).second.size()  << "\n";
//
//	}
	for (unsigned int i=0; i<accepted_segments.size(); i++) {
		std::cout << i << ":" << clean[i].size() << "\n";
	}
}



