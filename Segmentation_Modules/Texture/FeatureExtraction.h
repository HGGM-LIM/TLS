/*
 * FeatureExtraction.h
 *
 *  Created on: Apr 6, 2016
 *      Author: pmacias
 */

#ifndef FEATUREEXTRACTION_H_
#define FEATUREEXTRACTION_H_

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include  "itkScalarImageToTextureFeaturesFilter.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"


using namespace itk;

class FeatureExtraction {

private:


	std::string CT_path;
	std::string CT_labelled_path;
	std::string out_features_file;
	unsigned int bins;

protected:

public:
	typedef signed short PixelType;
	typedef unsigned short PixelLabelType;
	const static unsigned int Dimension=3;

	typedef itk::Image<PixelType,Dimension> InputImageType;
	typedef itk::Image<PixelLabelType,Dimension> InputLabelImageType;
	typedef itk::ImageFileReader<InputImageType> ImageFileReaderType;
	typedef itk::ImageFileReader<InputLabelImageType> ImageFileReaderLabelType;
	typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<InputImageType> TextureFeaturesFilterType;
	typedef TextureFeaturesFilterType::FeatureValueVector TextureFeatureVectorType;
	typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter<InputImageType> RunLengthFeaturesFilterType;






	FeatureExtraction();
	FeatureExtraction(const std::string CT_path, const std::string CT_labelled_path, const unsigned int bins, const std::string out_features_file);

	virtual ~FeatureExtraction();

	//TextureFeatureVectorType getAllFeatures();
	int getAllFeatures();
	/**
	 * Mask required to be same type of Ct_Vol
	 */

	std::string averagedTextureFeatures(InputImageType* CT_Vol,InputLabelImageType* CT_label_Vol);
	std::string firstOrderStats(InputImageType* CT_Vol, InputLabelImageType* CT_label_Vol, signed short min, signed short max, signed short label );



	static void hello(std::string password = "Hola Mundo");





};

#ifndef ITK_MANUAL_INSTANTIATION
#include "FeatureExtraction.cxx"
#endif
#endif /* FEATUREEXTRACTION_H_ */
