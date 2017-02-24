/*
 * FeatureExtraction.cxx
 *
 *  Created on: Apr 6, 2016
 *      Author: pmacias
 */

#include "FeatureExtraction.h"
#include "itkLabelMapMaskImageFilter.h"

#include "itkMaskImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCropImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkStatisticsLabelObject.h"

#include <iostream>
#include <fstream>

#include "itkAutoCropLabelMapFilter.h"


FeatureExtraction::FeatureExtraction() {}
FeatureExtraction::FeatureExtraction(const std::string CT_path, const std::string CT_labelled_path, const unsigned int bins,const std::string out_features_file) {
	this->CT_path = CT_path;
	this->CT_labelled_path = CT_labelled_path;
	this->bins = bins;
	this->out_features_file = out_features_file;
}

FeatureExtraction::~FeatureExtraction() {}

int FeatureExtraction::getAllFeatures(){



	ImageFileReaderType::Pointer readerCTImage = ImageFileReaderType::New();
	ImageFileReaderLabelType::Pointer readerLabelImage = ImageFileReaderLabelType::New();

	typedef itk::LabelObject< PixelLabelType, Dimension >  LabelObjectType;
	typedef itk::LabelMap< LabelObjectType >          LabelMapType;
	// convert the label image into a LabelMap
	typedef itk::LabelImageToLabelMapFilter< InputLabelImageType, LabelMapType > LabelImageToLabelMapType;
	LabelImageToLabelMapType::Pointer labelMap = LabelImageToLabelMapType::New();

	 typedef itk::CropImageFilter <InputImageType, InputImageType> CropImageFilterType;
	 CropImageFilterType::Pointer cropImageFilter = CropImageFilterType::New();


	//MapMask
	typedef itk::LabelMapMaskImageFilter< LabelMapType,  InputLabelImageType > LabelMapMaskFilterType;
	LabelMapMaskFilterType::Pointer labelMapMaskFilter = LabelMapMaskFilterType::New();

	typedef itk::LabelMapMaskImageFilter< LabelMapType,  InputImageType > LabelMapMaskFilterType2;
	LabelMapMaskFilterType2::Pointer labelMapMaskFilter2 = LabelMapMaskFilterType2::New();

	TextureFeaturesFilterType::Pointer textureFeatsFilter = TextureFeaturesFilterType::New();



	readerCTImage -> SetFileName(CT_path);
	readerLabelImage->SetFileName(CT_labelled_path);
	//readerCTImage->Update();
	//readerLabelImage->Update();
	readerLabelImage->GetOutput()->SetSpacing(readerCTImage ->GetOutput() -> GetSpacing());
	labelMap->SetInput( readerLabelImage->GetOutput());
	//readerCTImage->Update();
	//readerLabelImage->Update();



	  typedef itk::Image<unsigned char,3> auxiliarImage;
	  typedef itk::ImageFileWriter<InputLabelImageType> WriterType;
	  typedef itk::ImageFileWriter<InputImageType> WriterType2;
	  typedef itk::CastImageFilter<InputLabelImageType,auxiliarImage > CastFilterType;
	//Iterate from here (Need to get label indexes)

	labelMapMaskFilter->SetInput(labelMap->GetOutput());
	labelMapMaskFilter->SetFeatureImage(readerLabelImage->GetOutput()); //Workaround due to texture filter does not work with labels


	labelMapMaskFilter2->SetInput(labelMap->GetOutput());
	labelMapMaskFilter2->SetFeatureImage(readerCTImage->GetOutput());


	labelMap->Update();
	std::cout << "labelMap updated" << std::endl;

	unsigned int nLabels = labelMap->GetOutput()->GetNumberOfLabelObjects();


	std::ostringstream file;
			file << this->out_features_file;
			file << this->bins;
			file << ".csv";
			std::ofstream myfile;
			myfile.open (file.str().c_str());
			myfile << "Label,Energy,Entropy_2,InverseDifferentMoment,Inertia,ClusterShade,ClusterProminance,"
					"SRE,LRE,GLNU,RLNU,LGLRE,HGLRE,SRLGLE,SRHGLE,LRLGLE,LRHGLE,Mean,Sigma,Sum,Skewness,Kurtosis,Entropy,Fifth_Percen_Value,"
					"Ninety_Fifth_Percen_Value,Fifth_Percen_Mean,Ninety_Fifth_Percen_Mean,Min,Max,Median\n";

	for(unsigned int i = 1; i < nLabels; i++) {
		//std::cout << i << std::endl;
		unsigned int label = labelMap->GetOutput()->GetLabels().at(i);

		labelMapMaskFilter->SetLabel(label);
		labelMapMaskFilter->SetCrop(true);

		labelMapMaskFilter2->SetLabel(label);
		labelMapMaskFilter2->SetCrop(true);


/*
		std::ostringstream fn;
		std::ostringstream fn2;
		//labelMapMaskFilter->GetOutput()->Print(fn);
		//std::cout << fn.str() << std::endl;
		fn << "/media/Users_UC3M/amunoz/LOTE_1/STUDY_2087/E35_HIGHDOSE/3WEEKS_21APR12/LUNGS_SE3_interfile/SegmentationTest/label_v3_" << label << ".mhd";
		fn2 << "/media/Users_UC3M/amunoz/LOTE_1/STUDY_2087/E35_HIGHDOSE/3WEEKS_21APR12/LUNGS_SE3_interfile/SegmentationTest/label_v3_Cropped_" << label << ".mhd";
		CastFilterType::Pointer caster = CastFilterType::New();
		caster->SetInput(labelMapMaskFilter->GetOutput());
		//caster->SetInput(autoCropLabelMapFilter->GetOutput(label));

		WriterType::Pointer writer = WriterType::New();
		WriterType2::Pointer writer2 = WriterType2::New();
		writer->SetFileName(fn.str());
		writer->SetInput(labelMapMaskFilter->GetOutput());

		writer2->SetFileName(fn2.str());
		writer2->SetInput(labelMapMaskFilter2->GetOutput());
		//writer->SetInput(labelMapMaskFilter->GetCropBorder());
		//writer->SetInput(autoCropLabelMapFilter->GetOutput(label));
		writer->Update();
		writer2->Update();*/



		std::string feats = averagedTextureFeatures(labelMapMaskFilter2->GetOutput(),labelMapMaskFilter->GetOutput());
		myfile << label<<","<< feats<<"\n";


	}
	myfile.close();
	return EXIT_SUCCESS;
	//labelMapMaskFilter->SetLabel();

}


std::string FeatureExtraction::averagedTextureFeatures(InputImageType* CT_Vol, InputLabelImageType* CT_label_Vol) {
	FeatureExtraction::TextureFeaturesFilterType::Pointer imageToFeaturesFilter = TextureFeaturesFilterType::New();
	FeatureExtraction::RunLengthFeaturesFilterType::Pointer imageToRunLengthFilter = RunLengthFeaturesFilterType::New();


	//Mask to calculate min and mmax
	typedef itk::MaskImageFilter<InputImageType,InputImageType> MaskFilterType;
	MaskFilterType::Pointer mask = MaskFilterType::New();

	typedef itk::CastImageFilter<InputLabelImageType,InputImageType > CastFilterType;
	CastFilterType::Pointer caster = CastFilterType::New();

	typedef itk::ImageFileWriter<InputImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();

	typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxCalcType;
	MinMaxCalcType::Pointer minMaxcalc = MinMaxCalcType::New();


	/*CT_label_Vol -> SetOrigin(CT_Vol->GetOrigin());
	CT_label_Vol -> SetSpacing(CT_Vol->GetSpacing());
	CT_label_Vol->Update();*/
	caster->SetInput(CT_label_Vol);

	//caster->GetOutput()->SetSpacing(CT_Vol->GetSpacing());
	//caster->GetOutput()->SetOrigin(CT_Vol->GetOrigin());

	caster->Update();
	minMaxcalc->SetImage(caster->GetOutput());
	minMaxcalc->Compute();
	signed short insideValue = minMaxcalc->GetMaximum();

	std::cout << "Min y max calculated" << std::endl;

  //Double masking outside value to get what is really inside mask. Workaround
	//Todo work with mask region instead whole image

	mask ->SetMaskImage(caster->GetOutput());
	mask->SetInput(CT_Vol);
	mask->SetOutsideValue(15000);
	mask->Update();
	std::cout << "First Mask calculated" << std::endl;
	minMaxcalc ->SetImage(mask->GetOutput());
	minMaxcalc->Compute();
	signed short min  = minMaxcalc->GetMinimum();
	mask->SetOutsideValue(-15000);
	mask->Update();
	std::cout << "Second mask calculated" << std::endl;
	minMaxcalc->Compute();
	signed short max  = minMaxcalc->GetMaximum();


	/*writer->SetFileName("/tmp/masked.mhd");
	writer->SetInput(mask->GetOutput());
	writer->Update();*/
	imageToFeaturesFilter->SetInput(CT_Vol);
	imageToFeaturesFilter->SetMaskImage(caster->GetOutput());
	imageToFeaturesFilter->SetNumberOfBinsPerAxis(this->bins);
	imageToFeaturesFilter->SetPixelValueMinMax(min, max);
	imageToFeaturesFilter->SetInsidePixelValue(insideValue);
	//imageToFeaturesFilter->FastCalculationsOn();

	InputImageType::SizeType size = CT_Vol->GetLargestPossibleRegion().GetSize();
	//std::cout << size << std::endl;
	InputImageType::PointType p0,p1;
	p0[0] = 0; p0[1] = 0; p0[2] = 0;
	p1[0] = size[0]; p1[1] = size[1]; p1[2] = size[2];

	imageToRunLengthFilter->SetInput(CT_Vol);
	imageToRunLengthFilter->SetMaskImage(caster->GetOutput());
	imageToRunLengthFilter->SetNumberOfBinsPerAxis(this->bins);
	imageToRunLengthFilter->SetPixelValueMinMax(min, max);
	imageToRunLengthFilter->SetInsidePixelValue(insideValue);
	//imageToRunLengthFilter->FastCalculationsOn();
	imageToRunLengthFilter->SetDistanceValueMinMax(0,p0.EuclideanDistanceTo(p1));

	//labelStats->SetInput(CT_Vol);
	//labelStats->SetLabelInput(CT_label_Vol);





	imageToFeaturesFilter->Update();
	imageToRunLengthFilter->Update();



	//imageToFeaturesFilter->GetFeatureMeans();

    const TextureFeaturesFilterType::FeatureValueVector* output = imageToFeaturesFilter->GetFeatureMeans();
    const RunLengthFeaturesFilterType::FeatureValueVector* output2 = imageToRunLengthFilter->GetFeatureMeans();

    std::ostringstream feats;



     for(unsigned int i = 0; i < output->size(); ++i)
      feats << (*output)[i]<<',';
     for(unsigned int i = 0; i < output2->size(); ++i) {
    	 std::string end = i == output2->size() - 1 ? "" : ",";
           feats << (*output2)[i] << end;
     }



     //std::cout << feats.str() << std::endl;
     std::string firstStats =  firstOrderStats(CT_Vol,CT_label_Vol,min,max,insideValue);
     //std::cout << feats.str() +","+ firstStats << std::endl;


/*
	TextureFeatureVectorType::Pointer featureVector = imageToFeaturesFilter->GetFeatureMeans();
	itk::VectorContainer<unsigned char,double >::Iterator vectorIt = featureVector->Begin();
	    while(vectorIt != featureVector->End()) {
	        	std::cout << vectorIt->Value() << std::endl;
	        	vectorIt++;
	        }*/

	return feats.str()+","+ firstStats;
}


std::string FeatureExtraction::firstOrderStats(InputImageType* CT_Vol, InputLabelImageType* CT_label_Vol, signed short min, signed short max,signed short label ) {

	typedef float RealType;
	unsigned int numberOfBins = this->bins;


	  /**
	   * Zeroth and first order measurements
	   * These include:
	   *   1. mean
	   *   2. variance
	   *   3. kurtosis
	   *   4. skewness
	   *   5. entropy
	   *   6. fifth percentile value
	   *   7. ninety-fifth percentile value
	   *   8. mean of lower fifth percentile
	   *   9. mean of upper fifth percentile
	   */

	  RealType mean;
	  RealType sigma;
	  RealType sum;
	  RealType variance;
	  RealType skewness;
	  RealType kurtosis;
	  RealType median;
	  RealType entropy = 0;

	  typedef itk::LabelStatisticsImageFilter<InputImageType, InputLabelImageType> HistogramGeneratorType;
	  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
	  stats->SetInput( CT_Vol);
	  stats->SetLabelInput( CT_label_Vol );

	  stats->SetHistogramParameters( numberOfBins, min, max );

	  stats->Update();

	  mean = stats->GetMean( label);
	  sum = stats->GetSum( label );
	  sigma = stats->GetSigma( label);
	  variance = sigma * sigma;
	  median = stats->GetMedian( label);

	  kurtosis = 0.0;
	  skewness = 0.0;

	  RealType N = 0.0;



	  itk::ImageRegionIterator<InputImageType> ItI( CT_Vol,CT_Vol->GetLargestPossibleRegion() );
	  //itk::ImageRegionIterator<InputLabelImageType> ItL( CT_label_Vol,CT_label_Vol->GetLargestPossibleRegion() );
	  itk::ImageRegionIterator<InputLabelImageType> ItM( CT_label_Vol,CT_label_Vol->GetLargestPossibleRegion() );
	  for ( ItI.GoToBegin(), ItM.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItM ){
	    if ( ItM.Get() == label ){
	      RealType value = ItI.Get();

	      RealType diff = value - mean;
	      skewness += ( diff * diff * diff );
	      kurtosis += ( diff * diff * diff * diff );

	      N += 1.0;
	      }
	    }
	  skewness /= ( ( N - 1 ) * variance * sigma );
	  kurtosis /= ( ( N - 1 ) * variance * variance );


	  typedef HistogramGeneratorType::HistogramType  HistogramType;
	  const HistogramType *histogram = stats->GetHistogram( label );

	  double fifthPercentileValue = histogram->Quantile( 0, 0.05 );
	  double ninetyFifthPercentileValue = histogram->Quantile( 0, 0.95 );

			double fifthPercentileMean = 0.0;
			double fifthN = 0.0;
			double ninetyFifthPercentileMean = 0.0;
			double ninetyFifthN = 0.0;
			double quantileValue = 0.0;

			  entropy = 0.0;
			  //std::cout<< "bins: "<< this->bins << " label "<< label << " Min: "<< min << " Max: " << max  <<  std::endl;
			  for( unsigned int i = 0; i < histogram->Size(); i++ ){
			    RealType p = static_cast<RealType>( histogram->GetFrequency( i, 0 )  )/ static_cast<RealType>( histogram->GetTotalFrequency() );
			    //TODO the width
			    entropy += p > 0 ? ( -p * vcl_log( p ))/* / vcl_log( 2.0 ) ) */: 0;
			    }
			  //std::cout<< "width: " << (max - min + 0.0) / numberOfBins << std::endl;
			  entropy /= vcl_log( (max - min +0.0) / numberOfBins );





			  for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItI ){
			    RealType value = ItI.Get();
			    if ( value <= fifthPercentileValue ){
			      fifthPercentileMean += value;
			      fifthN++;
			      }
			    else if ( value >= ninetyFifthPercentileValue ){
			      ninetyFifthPercentileMean += value;
			      ninetyFifthN++;
			      }
			    }

			  fifthPercentileMean /= fifthN;
			  ninetyFifthPercentileMean /= ninetyFifthN;

			std::ostringstream feats;

	  feats << mean << ","<< sigma << ","<< sum << ","<< skewness << ","<< kurtosis << ","<< entropy << ","
	            << fifthPercentileValue << ","<< ninetyFifthPercentileValue << ","<< fifthPercentileMean << ","
	            << ninetyFifthPercentileMean << ","<< min << ","<< max << ","<< median ;

/*
	  std::cout << "mean:        " << mean << std::endl;
	  std::cout << "sigma:       " << sigma << std::endl;
	  std::cout << "sum:         " << sum << std::endl;
	  std::cout << "skewness:    " << skewness << std::endl;
	  std::cout << "kurtosis:    " << kurtosis << std::endl;
	  std::cout << "entropy:     " << entropy << std::endl;
	  std::cout << "5th %:       " << fifthPercentileValue << std::endl;
	  std::cout << "95th %:      " << ninetyFifthPercentileValue << std::endl;
	  std::cout << "5th % mean:  " << fifthPercentileMean << std::endl;
	  std::cout << "95th % mean: " << ninetyFifthPercentileMean << std::endl;
	  std::cout << "Min value:   " << min << std::endl;
	  std::cout << "Max value:   " << max << std::endl;
	  std::cout << "Median:   " << median << std::endl;
*/

	  return feats.str();

}

void FeatureExtraction::hello(std::string password) {
	std::cout << password << std::endl;
}



