#ifndef __itkVisualScoreImageFilter_txx
#define __itkVisualScoreImageFilter_txx

#include "itkVisualScoreImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkExtractImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkShapeLabelMapFilter.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
VisualScoreImageFilter< TInputImage, TOutputImage >
::VisualScoreImageFilter()
{
	m_Verbose = 0;
	m_MinArea = 2000;
	m_LowTh = -1024;
	m_EnfTh = -960;
	m_TissueTh = -1;
	m_LowestScoreTh = 1.0;

  	m_VisualScore = 0.0;
  	m_LowVisualScore = 0.0;
  	m_MedVisualScore = 0.0;
  	m_UppVisualScore = 0.0;

}

template< class TInputImage, class TOutputImage >
void
VisualScoreImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
	//this->AllocateOutputs();

	typename InputImageType::Pointer inputPtr = this->GetInput(0);
	//typename OutputImageType::Pointer outputPtr = this->GetOutput(0);
	RegionType requestedRegion = inputPtr->GetRequestedRegion();
	IndexType requestedIndex = requestedRegion.GetIndex();
	SizeType requestedSize = requestedRegion.GetSize();

    typedef itk::Image<PixelType, NDimensions-1> InternalInputImageType;
    typedef itk::Image<PixelType, NDimensions-1> InternalOutputImageType;
    typedef typename InternalInputImageType::Pointer InternalInputImagePointer;
    typedef typename InternalInputImageType::RegionType InternalRegionType;
    typedef typename InternalRegionType::SizeType InternalSizeType;
    typedef typename InternalRegionType::IndexType InternalIndexType;

    typedef itk::ExtractImageFilter<InputImageType, InternalInputImageType  >  ExtractFilterType;
    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();

    typedef itk::BinaryThresholdImageFilter <InternalInputImageType, InternalInputImageType> ThresholdImageFilterType;
    typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
    typedef unsigned long LabelType;
    typedef itk::ShapeLabelObject<LabelType, NDimensions-1>	LabelObjectType;
    typedef itk::LabelMap<LabelObjectType> LabelCollectionType;
    typedef itk::BinaryImageToLabelMapFilter<InternalInputImageType, LabelCollectionType> ConverterType;
    typedef itk::ShapeLabelMapFilter<LabelCollectionType>  ShapeFilterType;
    typedef itk::ImageRegionConstIterator<InternalInputImageType> ConstIteratorType;


    typedef itk::TileImageFilter< InternalInputImageType, OutputImageType > TilerType;

	RegionType extractionRegion;
	SizeType extractionSize=requestedSize;
	IndexType extractionIndex;

	// Collapse last dimension
	extractionSize[NDimensions-1]=0;
	extractionRegion.SetSize( extractionSize );

	typename TilerType::Pointer tiler = TilerType::New();

	// Configure Tile filter layout for 2D->3D conversion
	itk::FixedArray< unsigned int, NDimensions > layout;
	layout[0] = 1;	layout[1] = 1;	layout[2] = 0;
	tiler->SetLayout( layout );


	const unsigned int sliceRange = static_cast<unsigned int >(
			requestedSize[NDimensions-1] ) + requestedIndex[NDimensions-1];

	std::vector<unsigned int> scores;
	unsigned int firstSlice = 999;
	unsigned int lastSlice = 0;

	for( unsigned int slice = requestedIndex[NDimensions-1]; slice < sliceRange; slice++ )
	{

			extractionIndex[NDimensions-1]=slice;
			extractionRegion.SetIndex( extractionIndex );
			extractFilter->SetInput(inputPtr);
			extractFilter->SetExtractionRegion( extractionRegion );
			extractFilter->Update();

			InternalInputImagePointer gray_slice = extractFilter->GetOutput();

		    /**
		     ** Binarize the grayscale slice and compute area
		     **/

		    thresholdFilter->SetInput(gray_slice);
		    thresholdFilter->SetLowerThreshold(m_LowTh);
		    thresholdFilter->SetUpperThreshold(m_TissueTh);
		    thresholdFilter->SetInsideValue(255);
		    thresholdFilter->SetOutsideValue(0);
		    thresholdFilter->Update();
		    InternalInputImagePointer bin_slice = thresholdFilter->GetOutput();
		    bin_slice->DisconnectPipeline();

		    // FIXME: erode remaining airways

		    // Compute slice area
		    ConstIteratorType in(bin_slice, bin_slice->GetRequestedRegion());
		    unsigned long sliceArea = 0;
		    for (in.GoToBegin(); !in.IsAtEnd(); ++in)
		    	  if (in.Get() != 0)
		    		  sliceArea++;

		    if (m_Verbose)
		    	std::cout << "Slice "<< slice << " area is: "<< sliceArea;

			/**
			 ** Compute emphysema map on the grayscale image
			 **/

			thresholdFilter->SetInput(gray_slice);
			thresholdFilter->SetLowerThreshold(m_LowTh);
			thresholdFilter->SetUpperThreshold(m_EnfTh);
			thresholdFilter->SetInsideValue(255);
			thresholdFilter->SetOutsideValue(0);
			thresholdFilter->Update();
			InternalInputImagePointer enf_map = thresholdFilter->GetOutput();
			enf_map->DisconnectPipeline();

		    if (sliceArea > m_MinArea) {
		    	if (slice < firstSlice)
		    		firstSlice = slice;
		    	if (slice > lastSlice)
		    		lastSlice = slice;

				typename ConverterType::Pointer enf_labeler =  ConverterType::New();
				enf_labeler->SetInput(enf_map);
				enf_labeler->SetInputForegroundValue(255);
				//We don't want clusters that are separeted only by a few pixels to be considered one
				enf_labeler->SetFullyConnected(0);

				typename ShapeFilterType::Pointer shape = ShapeFilterType::New();
				shape->SetInput( enf_labeler->GetOutput() );
				shape->Update();

				unsigned long enfArea = 0;
				typename LabelCollectionType::Pointer enf_labels = shape->GetOutput();
				for(unsigned int label=1; label<enf_labels->GetNumberOfLabelObjects(); label++ ) {
					typename LabelObjectType::Pointer labelObject = enf_labels->GetLabelObject( label );
					enfArea += labelObject->GetNumberOfPixels();
				}

				float enfAff = 100.0 * enfArea / sliceArea;
				unsigned int score = 0;
				if ((enfAff > m_LowestScoreTh) and (enfAff < 25.0))
					score = 1;
				if ((enfAff>=25.0) and (enfAff<50.0))
					score = 2;
				if ((enfAff>=50.0) and (enfAff<75.0))
					score = 3;
				if (enfAff>75.0)
					score = 4;

				if (m_Verbose)
					std::cout << " enf affectation is " << enfAff << " score is: " << score << std::endl;

				scores.push_back(score);


		    } else {
		    	if (m_Verbose)
		    		std::cout << std::endl;
		    }

		    tiler->SetInput( slice, enf_map );
	}

	// Computing visual score as the mean of the observations

	unsigned int sum = 0;
	unsigned int size = scores.size();
	int sector_size = size / 3;

	if (m_Verbose) {
		std::cout << "Processed " << sliceRange << " slices." <<std::endl;
		std::cout << "Lungs was from " << firstSlice << " to " << lastSlice << " slices. " << std::endl;
		std::cout << "For a total of " << lastSlice - firstSlice << " slice processed" << std::endl;
		std::cout << "Total " << size << " scores computed" << std::endl;
		std::cout << "Sector size is " << sector_size << std::endl;
	}

	// Divide into three sector: upper, medium, lower

	// Upper
	sum=0;
	unsigned int startS = 0;
	unsigned int stopS = startS + sector_size;
	for (unsigned int i=startS; i<stopS; i++)
	    sum += scores[i];
	m_UppVisualScore = 1.0*sum / sector_size;

	if (m_Verbose)
		std::cout << "Upper visual score from " << startS + firstSlice << " to " << firstSlice + stopS
			<< " is " << m_UppVisualScore << std::endl;

	// Medium
	sum=0;
	startS = stopS;
	stopS = startS + sector_size;
	for (unsigned int i=startS; i<stopS; i++)
	    sum += scores[i];
	m_MedVisualScore = 1.0*sum / sector_size;

	if (m_Verbose)
		std::cout << "Medium visual score from " << startS + firstSlice << " to " << firstSlice + stopS
			<< " is " << m_MedVisualScore<< std::endl;

	// Lower
	sum=0;
	startS = stopS;
	stopS = startS + sector_size;
	for (unsigned int i=startS; i<stopS; i++)
		sum += scores[i];
	m_LowVisualScore = 1.0*sum / sector_size;

	if (m_Verbose)
		std::cout << "Lower visual score from " << startS + firstSlice << " to " << firstSlice + stopS
			<< " is " << m_LowVisualScore<< std::endl;

	sum = 0;
	while (!scores.empty()) {
	    sum+=scores.back();
	    scores.pop_back();
	}
	m_VisualScore = 1.0*sum / size;

	if (m_Verbose)
		std::cout << "Visual score calculated on " << size << " slices is " << m_VisualScore << std::endl;

	tiler->Update();
	this->GraftOutput(tiler->GetOutput());

}

template< class TInputImage, class TOutputImage >
void
VisualScoreImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
    
}

} // end namespace

#endif
