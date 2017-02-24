
#include "segments.h"
#include "itkImageFileWriter.h"

std::vector<PointType> Segment::get_points_copy() {
	std::vector<PointType> bag;
	std::vector<PointType>::const_iterator pit;
	std::vector<PropagationInfo>::const_iterator it;

	for ( it = propagations.begin() ; it < propagations.end(); it++ ) {
		for (pit = it->trials.begin(); pit != it->trials.end(); pit++ )
			bag.push_back(*pit);
		for (pit = it->points.begin(); pit != it->points.end(); pit++ )
			bag.push_back(*pit);
	}


	return bag;

}
int Segment::get_size() {
	/*
	 * Get total number of points contained in the segment
	 */

	int acc=0;

	std::vector<PropagationInfo>::iterator it;
	for ( it = propagations.begin() ; it < propagations.end(); it++ )
		acc += it->points.size();

	return acc;
}


RegionType Segment::get_bbox(const PropagationInfo prop, int pad) {
	/*
	 * Return the bbox with the passed prop added.
	 * Mainly used by leakage correction code
	 */

	this->propagations.push_back(prop);
	RegionType reg = this->get_bbox(pad);
	this->propagations.pop_back();
	return reg;

}

RegionType Segment::get_bbox(int pad) {
	/*
	 *
	 * We need to compute the BB because most of the correction methods are too expensive
	 * for the whole or make sense only in a reduced scope (e.g. re-thresholding the original image)
	 * We compute it from the first to the last propagation
	 */

/*	Index idx; idx.Fill(INT_MAX);
	Size size; size.Fill(0);

	std::vector<PropagationInfo>::iterator it;
	for ( it = propagations.begin() ; it < propagations.end(); it++ ) {
		RegionType t_reg = it->get_bbox();
		Index t_idx = t_reg.GetIndex();
		Size t_size = t_reg.GetSize();
		if (t_size[0]==0 && t_size[1]==0 && t_size[2]==0) {
			//propagation invalid... skipping it
			continue;
		}
		for (unsigned int i=0; i<idx.GetIndexDimension();i++)
			idx[i] = std::min(idx[i],t_idx[i]);
		for (unsigned int i=0; i<size.GetSizeDimension();i++)
			size[i] = std::max(size[i],t_size[i]);

	}*/

	Index idx; idx.Fill(INT_MAX);
	Index max_idx; max_idx.Fill(0);
	Size size; size.Fill(0);
	std::vector<PointType> points, trials;
	std::vector<PropagationInfo>::iterator it;
	for ( it = propagations.begin() ; it != propagations.end(); ++it ) {

		points = (*it).points; trials = (*it).trials;
		if (it->points.empty() && it->trials.empty()){
			l->info("The bounding box of an empty propagation is 0-sized.");
		} else {
			l->debugStream() << "Prop points Size: " << points.size();
			l->debugStream() << "Trial points Size: " << trials.size();

			std::vector< PointType >::iterator it2;
			for (it2 = points.begin(); it2 != points.end(); ++it2) {
				PointType p = (*it2);
				for (unsigned int i=0; i<idx.GetIndexDimension();i++) {
					idx[i] = std::min((int) idx[i],(int) p[i]);
					max_idx[i] = std::max((int) max_idx[i], (int) p[i]);
				}
			}
			for (it2 = trials.begin(); it2 != trials.end(); ++it2) {
				PointType p = (*it2);
				for (unsigned int i=0; i<idx.GetIndexDimension();i++) {
					idx[i] = std::min((int) idx[i],(int) p[i]);
					max_idx[i] = std::max((int) max_idx[i], (int) p[i]);
				}
			}


		}
	}

//	for ( it = propagations.begin() ; it < propagations.end(); it++ ) {
//		RegionType t_reg = it->get_bbox();
//		Index t_idx = t_reg.GetIndex();
//		Size t_size = t_reg.GetSize();
//		if (t_size[0]==0 && t_size[1]==0 && t_size[2]==0) {
//			//propagation invalid... skipping it
//			continue;
//		}
//		for (unsigned int i=0; i<idx.GetIndexDimension();i++) {
//			idx[i] = std::min(idx[i],t_idx[i]);
//			max_idx[i] = std::max(max_idx[i],t_idx[i]);
//		}
////		for (unsigned int i=0; i<size.GetSizeDimension();i++) {
////			size[i] = std::max(size[i],t_size[i]);
////		}
//
//	}
	if (propagations.size() > 1) {
		l->infoStream() << "MAX: " << max_idx;
		l->infoStream() << "MIN: " <<  idx;
		size[0] = max_idx[0] - idx[0];
		size[1] = max_idx[1] - idx[1];
		size[2] = max_idx[2] - idx[2];
	}
	//Add additional padding pixels

	for (unsigned int i=0; i<size.GetSizeDimension();i++)
		size[i] += pad;

	std::cout << "final reg: "<< idx << ", size: "<< size << "\n";
	RegionType reg(idx, size);
	return reg;
}

IntImagePointer Segment::as_image() {


	RegionType bbox = this->get_bbox(10); //Pad with 10 voxel for security
	Index bidx = bbox.GetIndex();

	/*
	 * Image creation
	 */

	Index ridx; ridx.Fill(0);
	Size size; size[0]=bbox.GetSize()[0];size[1]=bbox.GetSize()[1];size[2]=bbox.GetSize()[2];
	RegionType reg(ridx, size);
	int x_l = ridx[0]+size[0]; 	int y_l = ridx[1]+size[1]; 	int z_l = ridx[2]+size[2];
	std::cout << "x_l: " << x_l << ", "<< "y_l: " << y_l << ", "<< "z: " << z_l;
	IntImagePointer seg_img = IntImageType::New();
	seg_img->SetRegions(reg);
	seg_img->Allocate();

	std::vector<PropagationInfo>::iterator it = this->propagations.begin();
	std::vector<PointType>::iterator pit;
	PointType p;
	unsigned int count = 0;
	while(it != this->propagations.end()) {
		count++;
		std::cout << "Writing prop "<<count << " with " << it->trials.size() << " trials and " << it->points.size() << " points\n";

		for (pit = it->trials.begin(); pit !=it->trials.end(); pit++) {
			//Move the pixels back to origin
			Index idx; idx[0] =(*pit)[0] - bidx[0]; idx[1] =(*pit)[1] - bidx[1]; idx[2] =(*pit)[2] - bidx[2];
			if ( not reg.IsInside( idx ) ) {
				std::cout << "Trial index: " << idx << std::endl;
				std::cout << "Image region: " << reg << std::endl;
				std::cout.flush();
				//Abort the program!
				assert( reg.IsInside( idx ) );
			}
			if (seg_img->GetPixel(idx) == 0)//Do not overwrite previous data
				seg_img->SetPixel(idx, count*10 + 1);

		}

		for (pit = it->points.begin(); pit !=it->points.end(); pit++) {
			//Move the pixels back to origin
			Index idx; idx[0] =(*pit)[0] - bidx[0]; idx[1] =(*pit)[1] - bidx[1]; idx[2] =(*pit)[2] - bidx[2];
			if ( not reg.IsInside( idx ) ) {
				std::cout << "Point index: " << idx << std::endl;
				std::cout << "Image region: " << reg << std::endl;
				std::cout.flush();
				//Abort the program!
				assert( reg.IsInside( idx ) );
			}

			if (seg_img->GetPixel(idx) == 0) //Do not overwrite previous data
				seg_img->SetPixel(idx, count*10);


		}


		it++;

	}
/*
	typedef  itk::ImageFileWriter< IntImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName( "debug_seg_3.mhd" );
	writer->SetInput(seg_img);

	 try {
		 writer->Update();
	 }catch( itk::ExceptionObject & excep )
	    {
	    std::cerr << "Exception caught !" << std::endl;
	    std::cerr << excep << std::endl;
	    }


	 */
	return seg_img;

}


void Segment::seal() {
	this->ComputeMinimumRadius();
	this->ComputeMeanRadius();
	l->infoStream() << "Segment " << segmentNumber << " sealed!" ;
	sealed = true;
	status = "sealed";
}

void Segment::unseal() {

	l->infoStream() << "Segment " << segmentNumber << " unsealed!" ;
	sealed = false;
	status = "growing";
}


void Segment::accept_propagation(PropagationInfo prop) {

	if (is_sealed()) {
		l->warnStream() << "Segment is sealed and won't accept a new propagation!" ;
		return;
	}

	switch (prop.status) {
	  case PropagationInfo::NEW:
		  prop.status = PropagationInfo::ACCEPTED;
		  break;
	  case PropagationInfo::DISCONNECTED:
		  seal();
		  break;
	  case PropagationInfo::REJ_NO_POINTS:
		  seal(); //it cannot grow anymore
		  break;
	  case PropagationInfo::REJ_EXPLOSION:
		  mark_leaking();
		  break;
	  case PropagationInfo::REJ_MINRADIUS:
		  mark_leaking();
		  break;
	  case PropagationInfo::REJ_PREVRADIUS:
		  mark_leaking();
		  break;
	  case PropagationInfo::REJ_TRIAL_EXPLOSION:
		  mark_leaking();
		  break;
	  default:
		  l->warnStream() << "DEFAULT " << prop.status;
		break;
    	}
	propagations.push_back(prop);

}

std::vector<PointType> Segment::get_wavefront() {
	return propagations.back().trials;
}

PropagationInfo Segment::get_last_propagation() {
	/*
	 * Do not return a pointer because vector can invalidate references when a new object is added
	 */
	return propagations.back();
}
void Segment::setSegmentNumber(int _ownSegmentNumber) {
	segmentNumber = _ownSegmentNumber;
}

void Segment::setParentSegmentNumber(int _parentSegmentNumber) {
	parentSegmentNumber = _parentSegmentNumber;
}




void Segment::ComputeMinimumRadius( void )
  {
	l->debugStream() << "SEGMENT::ComputeMinimumRadius\n";
	if (propagations.size()==0) return;
    std::vector<PropagationInfo>::iterator it, end;
    it = propagations.begin();
    end  = propagations.end();
    minRadius = 1000;
    double threshold = propagations.size() * 0.25; // to take care of segment 1
    int count = 1;
    while(it!=end)
    {
        float W = (*it).trials.size(); // if shape is dome-like? Why is the formula pi*r² and not 2pi*r²
        if (W < 2.0) { ++it; ++count; continue; } // exclude seed my dear!!
        // special case for segment 1, consider only propagation from 25th percentile onwards
        if (segmentNumber==1) {
        	if (count < threshold) { ++it; ++count; continue; }
        }
        float radius = sqrt(W/(2*M_PI));
        l->debugStream() << "Prop has: "<< W << " points, radius: " << radius << "\n";
    	if (radius < minRadius) {
    		minRadius = radius;
    	}
    	++it;
    	++count;
    }
    l->debugStream() << "minRadius: "<< minRadius;
  }
double Segment::GetMinimumRadius(void) {
	return minRadius ;
}
void Segment::SetMeanRadius(double r) {
	meanRadius = r;
}
double Segment::GetMeanRadius(void) {
	return meanRadius ;
}

void Segment::ComputeMeanRadius(void)
  {
	std::cout << "SEGMENT::ComputeMeanRadius\n";
	if (propagations.size()==0) return ;
    float meanRadius = 0;
    int count = 0;
    std::vector<PropagationInfo>::iterator it, end;
    it = propagations.begin();
    end  = propagations.end();
    double threshold = propagations.size() * 0.25; // to take care of segment 1

    while(it!=end)
      {
    	float W = (*it).trials.size();
    	if (W < 2.0) { ++it; continue; }
        // special case for segment 1, consider only propagation from 25th percentile onwards
        if (segmentNumber==1) {
        	if (count < threshold) { ++it; ++count; continue; }
        }
    	meanRadius = meanRadius + sqrt(W/(2*M_PI));
    	count = count + 1;
    	++it;
      }
    SetMeanRadius(meanRadius/float(count));
  }

Segment::Segment() {
	sealed=false;
	l = get_logger("SEGMENT");
	status = "new";
	leaking = false;
	leakage_count = 0;
}

Segment::Segment(int _ownSegmentNumber, int _parentSegmentNumber) {

	parentSegmentNumber = _parentSegmentNumber;
	segmentNumber = _ownSegmentNumber;
	sealed=false;
	status = "new";
	l = get_logger("SEGMENT");
	leaking = false;
	leakage_count = 0;
}




