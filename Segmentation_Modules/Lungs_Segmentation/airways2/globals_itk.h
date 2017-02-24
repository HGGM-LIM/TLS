/*
 * globals.h
 *
 *  Created on: Mar 25, 2011
 *      Author: mario
 */
#ifdef USE_GLOBALS
    #define EXTERN
#else
    #define EXTERN extern
#endif

#ifndef GLOBALS_ITK_H_
#define GLOBALS_ITK_H_

#include "itkImage.h"
#include "itkVector.h"
#include "itkImageRegion.h"
#include "itkIndex.h"
#include "itkSize.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkExtractImageFilter.h"
#include <assert.h>
#include <limits.h>
#include <algorithm>
#include "log.h"

#define LDimension 3

//static const unsigned int Dimension = 2;
typedef float FloatPixelType;
typedef unsigned char UCharPixelType;
typedef int IntPixelType;

typedef itk::Image< FloatPixelType, LDimension > FloatImageType;
typedef itk::Image< UCharPixelType, LDimension > UCharImageType;
typedef itk::Image< IntPixelType, LDimension > IntImageType;

typedef typename FloatImageType::Pointer FloatImagePointer;
typedef typename UCharImageType::Pointer UCharImagePointer;
typedef typename IntImageType::Pointer IntImagePointer;

typedef typename UCharImageType::IndexType UCharIndexType;

typedef itk::Point< unsigned int, LDimension > PointType;
typedef itk::Point< int, LDimension > IntPointType;
typedef itk::Point< double, LDimension > DoublePointType;
typedef typename FloatImageType::SizeType SizeType;
typedef typename FloatImageType::IndexType IndexType;

/*
 * FIXME: These two ARE NOT the IndexType and the SizeType of the image
 * They should be removed and all its references changed to use the type above instead
 */
typedef itk::Size<LDimension> Size;
typedef itk::Index<LDimension> Index;
typedef itk::ImageRegion<LDimension> RegionType;

typedef itk::Vector< int, LDimension > VectorType;
typedef itk::Vector< double, LDimension > DirVectorType;

typedef itk::BinaryThresholdImageFilter<IntImageType, FloatImageType> BinFilterType;
typedef BinFilterType::Pointer BinFilterPointer;

typedef itk::ExtractImageFilter<IntImageType, IntImageType> ExtractFilterType;
typedef ExtractFilterType::Pointer ExtractFilterPointer;
typedef itk::ExtractImageFilter<FloatImageType, FloatImageType> FloatExtractFilterType;
typedef FloatExtractFilterType::Pointer FloatExtractFilterPointer;

#define OFFSET_LENGHT 1
#define itkversion ITK_VERSION_MAJOR
// Offsets
#if LDimension == 2
#undef OFFSET_LENGHT
#define OFFSET_LENGHT 8
const int up[2]      = { 0, 1};
const int right[2]   = { 1, 0};
const int down[2]    = { 0,-1};
const int left[2]    = {-1, 0};
const int up_r[2]    = { 1, 1};
const int up_l[2]    = { -1, 1};
const int down_r[2]  = { 1,-1};
const int down_l[2]  = {-1,-1};

//const unsigned int OFFSET_LENGHT = 8;
const VectorType offsets[OFFSET_LENGHT] = {VectorType(up),VectorType(right),VectorType(down),
		VectorType(left),VectorType(up_r),VectorType(up_l),VectorType(down_r),VectorType(down_l)};

/*
 * FIXME: Why I cannot conditionally declare the 'offset' variable for both 2D and 3D?
 */

#elif LDimension == 3
#undef OFFSET_LENGHT
#define OFFSET_LENGHT 26

const int i1[3] = {-1,-1,-1};
const int i2[3] = {-1,-1,0};
const int i3[3] = {-1,-1,1};
const int i4[3] = {-1,0,-1};
const int i5[3] = {-1,0,0};
const int i6[3] = {-1,0,1};
const int i7[3] = {-1,1,-1};
const int i8[3] = {-1,1,0};
const int i9[3] = {-1,1,1};
const int i10[3] = {0,-1,-1};
const int i11[3] = {0,-1,0};
const int i12[3] = {0,-1,1};
const int i13[3] = {0,0,-1};
const int i14[3] = {0,0,1};
const int i15[3] = {0,1,-1};
const int i16[3] = {0,1,0};
const int i17[3] = {0,1,1};
const int i18[3] = {1,-1,-1};
const int i19[3] = {1,-1,0};
const int i20[3] = {1,-1,1};
const int i21[3] = {1,0,-1};
const int i22[3] = {1,0,0};
const int i23[3] = {1,0,1};
const int i24[3] = {1,1,-1};
const int i25[3] = {1,1,0};
const int i26[3] = {1,1,1};

//const unsigned int OFFSET_LENGHT = 6;
const VectorType offsets[OFFSET_LENGHT] = {VectorType(i1),VectorType(i2),VectorType(i3),VectorType(i4),
		VectorType(i5),VectorType(i6),VectorType(i7),VectorType(i8),VectorType(i9),VectorType(i10),
		VectorType(i11),VectorType(i12),VectorType(i13),VectorType(i14),VectorType(i15),VectorType(i16),
		VectorType(i17),VectorType(i18),VectorType(i19),VectorType(i20),VectorType(i21),VectorType(i22),
		VectorType(i23),VectorType(i24),VectorType(i25),VectorType(i26)
};

#endif // Offsets

std::vector< PointType > getNeighs(const PointType,const PointType,const PointType);
bool contains(const std::vector< PointType > *, const PointType );

class PropagationInfo {
	std::string error_descr;
public:

	//FIXME: This is an ugly workaround to allow proper index2points conversion in the centroid calculation and should be removed
	FloatImagePointer space_converter;

	int id; //To indentify propagation when first met. It is set to zero initially. Whoever wants to use it can set to a value

	enum STATUS_T {
		ACCEPTED = 0,
		REJ_NO_POINTS = 1, // 1: rejected for not having points
		REJ_EXPLOSION = 2, // 2: rejected because (# of current points) > 1.5*(# of prev points) - eg. explosion
		REJ_MINRADIUS = 3, // 3: rejected because current radius > 1.5* minRadius of parent segment
		REJ_PREVRADIUS = 4,  // 4: rejected because current radius > 1.1* prevRadius
		REJ_TRIAL_EXPLOSION = 5,// 5: rejected because (# of trial points) > 1.5*(# of prev points) - eg. trial explosion
		NEW = 6, //Just created
		DISCONNECTED = 7 //This was a bifurcation
	};

    STATUS_T status;

    int currSegmentNum, parentSegmentNum;
    std::vector< PointType > trials;
    std::vector< PointType > points;
    SizeType image_size;

    PropagationInfo() {
    	status = NEW;
    }

    std::vector<double> get_centroid() const {
    	/*
		 * What is the centroid of a propagation?
		 * For now imagine it is the mean of all of its points...
		 */


    	/*
    	 * FIXME: This is a serious interoperability issue.
    	 * We mix Indexes and Points. Almost all images have spacing other than 1, so the difference matters.
    	 * For now just try and convert them in the output, but we should decide if we operate into Physical or grid coordinate
    	 * also in the FM code and do the appropriate transforms.
    	 */
    	std::vector<double> centroid;
    	centroid.push_back(0.0); centroid.push_back(0.0); centroid.push_back(0.0);
    	float tot_trials = (float) this->trials.size();
    	//= [ np.asarray(a).mean() for a in zip(*self.new_wf)]
    	std::vector< PointType >::const_iterator it;
    	for (it=this->trials.begin(); it!=this->trials.end();it++) {
    		/*
    		 * FIXME:  even if trials are declared as Points, the values we store inside are indexes.
    		 * So, to try and give correct informations to the following stages of the pipeline, for now
    		 * just convert them to indexes and call on them the appropriate transform function
    		 * But we really should fix it in the main code.
    		 */
    		PointType p = *it;

    		// Start the index/point workaround: 1st is point->index casting
    		IndexType temp_fix; temp_fix[0]=p[0]; temp_fix[1]=p[1]; temp_fix[2]=p[2];
    		DoublePointType temp_p;
    		//2nd Call the transform function
    		space_converter->TransformIndexToPhysicalPoint(temp_fix, temp_p);
    		//std::cout << "Index: " << temp_fix << std::endl;
    		//std::cout << "Point: " << temp_p << std::endl;

    		for (int i=0; i<LDimension;i++)
    			centroid[i] = centroid[i] + temp_p.GetElement(i)/tot_trials;
    	}

    	//std::cout << "Spacing: " << space_converter->GetSpacing() << std::endl;
    	return centroid;
    }

    std::vector<double> get_centroid_index() const {
        /*
                 * What is the centroid of a propagation?
                 * For now imagine it is the mean of all of its points...
                 */


        /*
         * FIXME: This is a serious interoperability issue.
         * We mix Indexes and Points. Almost all images have spacing other than 1, so the difference matters.
         * For now just try and convert them in the output, but we should decide if we operate into Physical or grid coordinate
         * also in the FM code and do the appropriate transforms.
         */
        std::vector<double> centroid;
        centroid.push_back(0.0); centroid.push_back(0.0); centroid.push_back(0.0);
        float tot_trials = (float) this->trials.size();
        //= [ np.asarray(a).mean() for a in zip(*self.new_wf)]
        std::vector< PointType >::const_iterator it;
        for (it=this->trials.begin(); it!=this->trials.end();it++) {
                /*
                 * FIXME:  even if trials are declared as Points, the values we store inside are indexes.
                 * So, to try and give correct informations to the following stages of the pipeline, for now
                 * just convert them to indexes and call on them the appropriate transform function
                 * But we really should fix it in the main code.
                 */
                PointType p = *it;

                // Start the index/point workaround: 1st is point->index casting
                IndexType temp_fix; temp_fix[0]=p[0]; temp_fix[1]=p[1]; temp_fix[2]=p[2];
//                DoublePointType temp_p;
//                //2nd Call the transform function
//                space_converter->TransformIndexToPhysicalPoint(temp_fix, temp_p);

                for (int i=0; i<LDimension;i++)
                        centroid[i] = centroid[i] + temp_fix.GetElement(i)/tot_trials;
        }


        return centroid;
    }


    float get_physical_radius() const {
		/*
		 * Returns the radius of the trials
		 */
		FloatImageType::SpacingType sp;
		sp.Fill(1.0);
		if (space_converter.IsNotNull()) {
			sp = space_converter->GetSpacing();
		} else {
			std::cout << "No spacing info available for radius calculation: defaulting to 1.0" << std::endl;

		}

		int curSize = trials.size();
		float curRadius = -1;

		if (LDimension==2) {
			// for 2D propagations, area is circle: pi*r²
			curRadius = sqrt(curSize*sp[0]*sp[1]/M_PI);
		} else if (LDimension==3) {
			// for 3D propagations, 'area' of half sphere: 0.5 * (4*pi*r²)
			curRadius = sqrt(curSize*sp[0]*sp[1]*sp[2]/(2 * M_PI));
//    	} else {
//    		l->errorStream() << "Unsupported LDimensions: " << Dimension;
		}
		return curRadius;
	}


    float get_radius() const {
    	/*
    	 * Returns the radius of the trials
    	 */


    	int curSize = trials.size();
    	float curRadius = -1;

    	if (LDimension==2) {
    		// for 2D propagations, area is circle: pi*r²
    		curRadius = sqrt(curSize/M_PI);
    	} else if (LDimension==3) {
    		// for 3D propagations, 'area' of half sphere: 0.5 * (4*pi*r²)
    		curRadius = sqrt(curSize/(2 * M_PI));
//    	} else {
//    		l->errorStream() << "Unsupported LDimensions: " << Dimension;
    	}
    	return curRadius;
    }

    RegionType get_bbox() const {
    	/*
    	 * FIXME: use array for LDimension
    	 */
    	pLogType l = get_logger("PROP::get_bbox");

    	int xmin, ymin, zmin, xMAX, yMAX, zMAX;
    	xmin = ymin = zmin = INT_MAX; //FIXME: should be the NumericTrait<int>::max() value
    	xMAX = yMAX = zMAX = -1;

    	if (points.empty() && trials.empty()){
    		l->debug("The bounding box of an empty propagation is 0-sized.");
    		Index idx; idx.Fill(0);
    		Size size; size.Fill(0);
    		RegionType reg(idx, size);

    		return reg;
    	}

    	std::vector< PointType >::const_iterator it;
    	for (it = trials.begin(); it != trials.end(); it++) {
			PointType p = (*it);
			xmin = std::min(xmin, (int)p.GetElement(0)); ymin = std::min(ymin, (int)p.GetElement(1)); zmin = std::min(zmin, (int)p.GetElement(2));
			xMAX = std::max(xMAX, (int)p.GetElement(0)); yMAX = std::max(yMAX, (int)p.GetElement(1)); zMAX = std::max(zMAX, (int)p.GetElement(2));
	   }
    	for (it = points.begin(); it != points.end(); it++) {
    		PointType p = (*it);
    		xmin = std::min(xmin, (int)p.GetElement(0)); ymin = std::min(ymin, (int)p.GetElement(1)); zmin = std::min(zmin, (int)p.GetElement(2));
    		xMAX = std::max(xMAX, (int)p.GetElement(0)); yMAX = std::max(yMAX, (int)p.GetElement(1)); zMAX = std::max(zMAX, (int)p.GetElement(2));
    	}

    	l->debugStream() << "x: " << xmin << ", "<< "y: " << ymin << ", "<< "z: " << zmin;
    	l->debugStream() << "X: " << xMAX << ", "<< "Y: " << yMAX << ", "<< "Z: " << zMAX;
    	assert(xmin>=0); assert(ymin>=0);assert(zmin>=0);
    	assert(xMAX>0); assert(yMAX>0);assert(zMAX>0);

    	Index idx; idx.SetElement(0, xmin); idx.SetElement(1, ymin);idx.SetElement(2, zmin);
    	Size size; size.SetElement(0, xMAX-xmin); size.SetElement(1, yMAX-ymin);size.SetElement(2, zMAX-zmin);
    	RegionType reg(idx, size);

    	return reg;
    }

    std::string get_status(){

    	switch (status) {
    	  case ACCEPTED:
    	    error_descr = "ACCEPTED";
    	    break;
    	  case REJ_NO_POINTS:
    		  error_descr = "REJECTED (no points)";
    		  break;
    	  case REJ_EXPLOSION:
    		  error_descr = "REJECTED (explosion)";
    		  break;
    	  case REJ_MINRADIUS:
    		  error_descr = "REJECTED (current radius >> minRadius of parent )";
    		  break;
    	  case REJ_PREVRADIUS:
    		  error_descr = "REJECTED (current radius >> prevRadius)";
    		  break;
    	  case REJ_TRIAL_EXPLOSION:
    		  error_descr = "REJECTED (trial explosion)";
    		  break;
    	  case NEW:
			  error_descr = "NEW";
			  break;
    	  case DISCONNECTED:
			  error_descr = "BIFURCATION";
			  break;
    	  default:
    		std::stringstream ss;
    		ss << "UNKNOWN STATE (" << status << ")";
    		error_descr = ss.str();
    	    break;
    	}



    	return error_descr;
    }
};


#endif /* GLOBALS_ITK_H_ */
