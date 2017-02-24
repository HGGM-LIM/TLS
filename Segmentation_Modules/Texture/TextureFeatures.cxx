/*
 * TextureFeatures.cxx
 *
 *  Created on: Mar 3, 2016
 *      Author: pmacias
 */

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkDenseFrequencyContainer2.h>
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkVectorContainer.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

#include "itkObjectByObjectLabelMapFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkMath.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkHistogramToImageFilter.h"
#include  "itkScalarImageToTextureFeaturesFilter.h"
#include "FeatureExtraction.h"
#include "Test.h"
#include "itkLabelMapMaskImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkCastImageFilter.h"



//definitions of used types
typedef itk::Image<float, 3> InternalImageType;
typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::Neighborhood<float, 3> NeighborhoodType;
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InternalImageType> Image2CoOccuranceType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType;
typedef InternalImageType::OffsetType OffsetType;
typedef itk::AddImageFilter <InternalImageType> AddImageFilterType;
typedef itk::MultiplyImageFilter<InternalImageType> MultiplyImageFilterType;

//calculate features for one offset
void calcTextureFeatureImage (OffsetType offset,InternalImageType::Pointer inputImage, InternalImageType::Pointer outInertia,
    InternalImageType::Pointer outCorrelation, InternalImageType::Pointer outEnergy)
{
    //allocate output images
    outInertia->CopyInformation(inputImage);
    outInertia->SetRegions(inputImage->GetLargestPossibleRegion());
    outInertia->Allocate();
    outInertia->FillBuffer(0);
    outCorrelation->CopyInformation(inputImage);
    outCorrelation->SetRegions(inputImage->GetLargestPossibleRegion());
    outCorrelation->Allocate();
    outCorrelation->FillBuffer(0);
    outEnergy->CopyInformation(inputImage);
    outEnergy->SetRegions(inputImage->GetLargestPossibleRegion());
    outEnergy->Allocate();
    outEnergy->FillBuffer(0);

    Image2CoOccuranceType::Pointer glcmGenerator=Image2CoOccuranceType::New();
    glcmGenerator->SetOffset(offset);
    glcmGenerator->SetNumberOfBinsPerAxis(16); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(0, 255); //for input UCHAR pixel type
    Hist2FeaturesType::Pointer featureCalc=Hist2FeaturesType::New();

    typedef itk::RegionOfInterestImageFilter<InternalImageType,InternalImageType> roiType;
    roiType::Pointer roi=roiType::New();
    roi->SetInput(inputImage);

    InternalImageType::RegionType window;
    InternalImageType::RegionType::SizeType size;
    size.Fill(3); //window size=3x3x3
    window.SetSize(size);
    InternalImageType::IndexType pi; //pixel index

    //slide window over the entire image
    for (unsigned x=1; x<inputImage->GetLargestPossibleRegion().GetSize(0)-1; x++)
    {
        pi.SetElement(0,x);
        window.SetIndex(0,x-1);
        for (unsigned y=1; y<inputImage->GetLargestPossibleRegion().GetSize(1)-1; y++)
        {
            pi.SetElement(1,y);
            window.SetIndex(1,y-1);
            for (unsigned z=1; z<inputImage->GetLargestPossibleRegion().GetSize(2)-1; z++)
            {
                pi.SetElement(2,z);
                window.SetIndex(2,z-1);
                roi->SetRegionOfInterest(window);
                roi->Update();
                glcmGenerator->SetInput(roi->GetOutput());
                glcmGenerator->Update();
                featureCalc->SetInput( glcmGenerator->GetOutput() );
                featureCalc->Update();

                outInertia->SetPixel(pi, featureCalc->GetFeature(Hist2FeaturesType::Inertia));
                outCorrelation->SetPixel(pi, featureCalc->GetFeature(Hist2FeaturesType::Correlation));
                outEnergy->SetPixel(pi, featureCalc->GetFeature(Hist2FeaturesType::Energy));
            }
        }
        std::cout<<'.';
    }
}

int itkScalarImageToCooccurrenceMatrixFilterTest(){
	  //Data definitions
	  const unsigned int  IMGWIDTH         =  5;
	  const unsigned int  IMGHEIGHT        =  5;
	  const unsigned int  NDIMENSION       =  2;
	  typedef signed short PixelTypeB;


	  //------------------------------------------------------
	  //Create a simple test images
	  //------------------------------------------------------
	  typedef itk::Image<unsigned char, NDIMENSION> InputImageType;
	  //typedef itk::Image<PixelTypeB, NDIMENSION> InputImageType;

	  typedef itk::ImageRegionIterator< InputImageType > InputImageIterator;
	  typedef itk::MinimumMaximumImageCalculator <InputImageType> ImageCalculatorFilterType;

	  InputImageType::Pointer image = InputImageType::New();

	  InputImageType::SizeType inputImageSize = {{ IMGWIDTH, IMGHEIGHT }};

	  InputImageType::IndexType index;
	  index.Fill(0);
	  InputImageType::RegionType region;

	  region.SetSize( inputImageSize );
	  region.SetIndex( index );

	  //--------------------------------------------------------------------------
	  // Set up the image first. It looks like:
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //--------------------------------------------------------------------------

	  image->SetRegions( region );
	  image->Allocate();

	  // setup the iterator
	  InputImageIterator imageIt( image, image->GetBufferedRegion() );

	  imageIt.GoToBegin();

	  for(unsigned int i = 0; i < 5; i++)
	    {
	    for(unsigned int j = 0; j < 5; j++, ++imageIt)
	      {
	      imageIt.Set(j % 2 + 1);
	     //if(i == 1 && j == 1) imageIt.Set(-16000);
	      //if(i == 2 && j == 2) imageIt.Set(16000);
	      }
	    }



	  //--------------------------------------------------------------------------
	  // Generate the histogram. The un-normalized histogram should look like this:
	  //
	  //     0 1  2 ...
	  //     ------
	  //  0 |0 0  0
	  //  1 |0 24 20
	  //  2 |0 20 16
	  //  3 |0 0  0
	  //  .
	  //  .
	  // with zeroes elsewhere.
	  //--------------------------------------------------------------------------

	  bool passed = true;

	  ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	  imageCalculatorFilter->SetImage(image);
	  imageCalculatorFilter->Compute();

	  typedef itk::ImageFileWriter<InputImageType> WriterType;
	  WriterType::Pointer writer = WriterType::New();
	  writer->SetInput(image);
	  writer->SetFileName("/home/pmacias/Projects/MonkeysTuberculosis/imageTestWithNegatiVal.mhd");
	  writer->Update();

	  typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter< InputImageType> FilterType;
	  FilterType::Pointer filter = FilterType::New();

	    filter->SetInput(image);

	    InputImageType::OffsetType offset1 = {{0, 1}};
	    InputImageType::OffsetType offset2 = {{1, 0}};
	    FilterType::OffsetVectorPointer offsetV = FilterType::OffsetVector::New();
	    offsetV->push_back(offset1);
	    offsetV->push_back(offset2);
	    filter->SetOffsets(offsetV);
	    //filter->SetPixelValueMinMax(imageCalculatorFilter->GetMinimum(),imageCalculatorFilter->GetMaximum());
	    unsigned int bins = 256;
	    filter->SetNumberOfBinsPerAxis(bins);
	    std::cout << "Max: " << filter->GetMax() << " Min: " << filter -> GetMin() << " bins: "<< filter->GetNumberOfBinsPerAxis() <<std::endl;

	    filter->Update();




	    const FilterType::HistogramType * hist = filter->GetOutput();
	    //--------------------------------------------------------------------------
	    // Test the histogram.
	    //--------------------------------------------------------------------------
	    // First make sure the bins are sized properly:

	    float max = hist->GetBinMax(0,255);
	    float min = hist->GetBinMin(0,255);
/*
	    if(itk::Math::NotAlmostEquals(max, 256)|| itk::Math::NotAlmostEquals(min, 255)){
	          std::cerr << "Error" << std::endl;
	          std::cerr << "The calculated bin sizes are incorrect" << std::endl;
	          std::cerr << "Expected [255, 256), got [" << min << ", " << max << ")" << std::endl << std::endl;
	          passed = false;
	          }*/

	        // Now make sure the contents of the bins are correct:
	        typedef FilterType::HistogramType::IndexType IndexType;
	        IndexType one_one( hist->GetMeasurementVectorSize() );
	        IndexType one_two( hist->GetMeasurementVectorSize() );
	        IndexType two_one( hist->GetMeasurementVectorSize() );
	        IndexType two_two( hist->GetMeasurementVectorSize() );
	        IndexType zero_zero(2);
	        IndexType zero_one(2);
	        IndexType zero_two(2);
	        IndexType one_zero(2);
	        IndexType two_zero(2);

	        zero_zero[0] = 0;zero_zero[1] = 0;
	        zero_one[0] = 0; zero_one[1] =1;
	        zero_two[0] = 0; zero_two[1] = 2;
	        one_zero[0] = 1; one_zero[1] =0;
	        two_zero[0] = 2; two_zero[1] =0;
	        one_one[0] = 1;one_one[1] = 1;
	        one_two[0] = 1;one_two[1] = 2;
	        two_one[0] = 2;two_one[1] = 1;
	        two_two[0] = 2;two_two[1] = 2;

	        float ooF, otF, toF, ttF, totalF,zzF,zoF,ztF,ozF,tzF;
	        zzF = hist->GetFrequency(zero_zero);
	        zoF = hist->GetFrequency(zero_one);
	        ztF = hist->GetFrequency(zero_two);
	        ozF = hist->GetFrequency(one_zero);
	        ooF = hist->GetFrequency(one_one);
	        otF = hist->GetFrequency(one_two);
	        tzF = hist->GetFrequency(two_zero);
	        toF = hist->GetFrequency(two_one);
	        ttF = hist->GetFrequency(two_two);
	        totalF = hist->GetTotalFrequency();


	        std::cout<< zzF << ", " << zoF << ", "<< ztF << ", "<< otF << ", "<< ozF << ", "<< ooF << ", " << ttF  << ", " << otF  << ", " << toF <<", "<<totalF<< std::endl;

	        for(unsigned int i = 0; i < bins; i++){
	        	for(int j = 0; j< bins; j++) {
	        		IndexType i_j(2);
	        		i_j[0] = i; i_j[1] =j;
	        		float freq = hist ->GetFrequency(i_j);
	        		if(freq > 0)
	        		std::cout << "i " << i << " j "<<j << " hist: "<<freq << std::endl;
	        	}
	        }

	        typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType;
	           	Hist2FeaturesType::Pointer featsCalc = Hist2FeaturesType::New();
	           	featsCalc->SetInput(filter->GetOutput());
	           	featsCalc->Update();

	           	std::cout << "Inertia: " << featsCalc->GetFeature(Hist2FeaturesType::Inertia) << " Energy: " << featsCalc->GetFeature(Hist2FeaturesType::Energy)
	           			<< " Entropy "<< featsCalc->GetFeature(Hist2FeaturesType::Entropy) << std::endl;


	        FilterType::Pointer filter0 = FilterType::New();

	        filter0->SetInput(image);
	        filter0->SetOffsets(offsetV);
	        filter0->NormalizeOn();
	        //filter0->SetPixelValueMinMax(0,1);

	        if ( filter0->GetNormalize() != true ) {
	          std::cerr << "Normalize boolean is not set correctly";
	          passed = false;
	          }

	        filter0->Update();
	        const FilterType::HistogramType * hist0 = filter0->GetOutput();
	        max = hist->GetBinMax(0,255);
	        min = hist->GetBinMin(0,255);

	        //"Expected 0.3, 0.2, 0.25, 0.25"
	        ooF = hist0->GetFrequency(one_one);
	        otF = hist0->GetFrequency(one_two);
	        toF = hist0->GetFrequency(two_one);
	        ttF = hist0->GetFrequency(two_two);
	        std::cout << "Expected 0.3, 0.2, 0.25, 0.25 got " << ooF << ", " << ttF  << ", " << otF  << ", " << toF << std::endl;


}

int itkScalarImageToCooccurrenceMatrixFilterTest2() {
	  //Data definitions
	  const unsigned int  IMGWIDTH         =  5;
	  const unsigned int  IMGHEIGHT        =  5;
	  const unsigned int  NDIMENSION       =  2;
	  //------------------------------------------------------
	  //Create a simple test images
	  //------------------------------------------------------
	  typedef unsigned char PixelType;
	  typedef itk::Image<PixelType, NDIMENSION> InputImageType;

	  typedef itk::ImageRegionIterator< InputImageType > InputImageIterator;
	  typedef itk::Neighborhood<PixelType, NDIMENSION> NeighborhoodType;


	  InputImageType::Pointer image = InputImageType::New();
	  InputImageType::Pointer mask = InputImageType::New();


	  InputImageType::SizeType inputImageSize = {{ IMGWIDTH, IMGHEIGHT }};

	  InputImageType::IndexType index;
	  index.Fill(0);
	  InputImageType::RegionType region;

	  region.SetSize( inputImageSize );
	  region.SetIndex( index );

	  //--------------------------------------------------------------------------
	  // Set up the image first. It looks like:
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //  1 2 1 2 1
	  //--------------------------------------------------------------------------

	  image->SetRegions( region );
	  image->Allocate();

	  // setup the iterator
	  InputImageIterator imageIt( image, image->GetBufferedRegion() );

	  imageIt.GoToBegin();

	  for(unsigned int i = 0; i < 5; i++)
	    {
	    for(unsigned int j = 0; j < 5; j++, ++imageIt)
	      {
	      imageIt.Set(j % 2 + 1);
	      }
	    }

	  //--------------------------------------------------------------------------
	  // Set up the mask next. It looks like:
	  //  0 0 0 0 0
	  //  0 0 1 0 0
	  //  0 0 1 0 0
	  //  0 0 1 0 0
	  //  0 0 0 0 0
	  //--------------------------------------------------------------------------

	  mask->SetRegions( region );
	  mask->Allocate();

	  // setup the iterator
	  InputImageIterator maskIt( mask, mask->GetBufferedRegion() );
	  maskIt.GoToBegin();
	  for(int i = 0; i < 5; i++)
	    for(int j = 0; j < 5; j++, ++maskIt)
	      {
	      if (j == 2 && i > 0 && i < 4)
	        {
	        maskIt.Set(1);
	        }
	      else
	        {
	        maskIt.Set(0);
	        }
	      if(i == 0 && j==0) maskIt.Set(1);
	      if(i == 1 && j==1) maskIt.Set(1);
	      }

	  //--------------------------------------------------------------------------
	  // Generate the histogram. It should look like this:
	  //
	  //     0 1 2 ...
	  //     -----
	  //  0 |0 0 0
	  //  1 |0 4 0
	  //  2 |0 0 0
	  //  .
	  //  .
	  //  .
	  // with zeroes elsewhere.
	  //--------------------------------------------------------------------------



	    typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InputImageType> FilterType;

	    FilterType::Pointer filter = FilterType::New();

	    filter->SetInput(image);

	    NeighborhoodType neighborhood;
        neighborhood.SetRadius(1);
        unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();

	    InputImageType::OffsetType offset1 = {{0, 1}};
	    InputImageType::OffsetType offset2 = {{1, 0}};
	    FilterType::OffsetVectorPointer offsetV = FilterType::OffsetVector::New();
	    offsetV->push_back(offset1);
	    offsetV->push_back(offset2);

	    /*for (unsigned int d = 0; d <= centerIndex - 1; d++){
	    	InputImageType::OffsetType offset = neighborhood.GetOffset(d);
	    	offsetV->push_back(neighborhood.GetOffset(d));
	    }*/


	    filter->SetOffsets(offsetV);
	    filter->SetMaskImage(mask);
	    filter->Update();
	    const FilterType::HistogramType * hist = filter->GetOutput();

	    //--------------------------------------------------------------------------
	    // Test the histogram.
	    //--------------------------------------------------------------------------
	    bool passed = true;

	    typedef FilterType::HistogramType::IndexType IndexType;
	    IndexType one_one( hist->GetMeasurementVectorSize() );
	    IndexType one_two( hist->GetMeasurementVectorSize() );
	    IndexType two_one( hist->GetMeasurementVectorSize() );
	    IndexType two_two( hist->GetMeasurementVectorSize() );
        IndexType zero_zero(2);
        IndexType zero_one(2);
        IndexType zero_two(2);
        IndexType one_zero(2);
        IndexType two_zero(2);

        zero_zero[0] = 0;zero_zero[1] = 0;
        zero_one[0] = 0; zero_one[1] = 1;
        zero_two[0] = 0; zero_two[1] = 2;
        one_zero[0] = 1; one_zero[1] = 0;
        two_zero[0] = 2; two_zero[1] = 0;

	    one_one[0] = 1;
	    one_one[1] = 1;

	    one_two[0] = 1;
	    one_two[1] = 2;

	    two_one[0] = 2;
	    two_one[1] = 1;

	    two_two[0] = 2;
	    two_two[1] = 2;

	    float ooF, otF, toF, ttF, totalF,zzF,zoF,ztF,ozF,tzF;
        zzF = hist->GetFrequency(zero_zero);
        zoF = hist->GetFrequency(zero_one);
        ztF = hist->GetFrequency(zero_two);
        ozF = hist->GetFrequency(one_zero);
	    ooF = hist->GetFrequency(one_one);
	    otF = hist->GetFrequency(one_two);
	    toF = hist->GetFrequency(two_one);
	    ttF = hist->GetFrequency(two_two);
	    totalF = hist->GetTotalFrequency();

	    return EXIT_SUCCESS;
}

int itkScalarImageToCooccurrenceMatrixFilterTest3() {
	typedef signed short PixelType;
	const unsigned int Dimension = 2;


	typedef itk::Image<PixelType,Dimension> InputImageType;
	typedef itk::ImageFileReader<InputImageType> ImageFileReaderType;
	typedef itk::ImageFileWriter<InputImageType> ImageFileWriterType;
	typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InputImageType> FilterType;
	//typedef itk::MaskImageFilter<InputImageType,InputImageType> MaskImageFilter;

	ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
	ImageFileReaderType::Pointer reader2 = ImageFileReaderType::New();
	ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
	FilterType::Pointer glcmGenerator = FilterType::New();

	reader->SetFileName("/media/Users_UC3M/amunoz/LOTE_1/STUDY_4863/R8/16WEEKS_08SEP13/Chest_-_Lung_Reformat_0.625mm_interfile/LUNG/lungsSlice.mhd");
	reader->Update();
	InputImageType* imageIn = reader->GetOutput();


	InputImageType::IndexType testPixel;
	InputImageType::IndexType testPixel2;

	testPixel[0] = 260;testPixel[1]= 180;
	testPixel2[0] = 145;testPixel2[1] = 240;
	std::cout << imageIn->GetPixel(testPixel)<<" " << imageIn->GetPixel(testPixel2)  << std::endl;

	glcmGenerator->SetInput(imageIn);
	/*writer->SetFileName("/media/Users_UC3M/amunoz/LOTE_1/STUDY_4863/R8/16WEEKS_08SEP13/Chest_-_Lung_Reformat_0.625mm_interfile/LUNG/lungSliceTest.mhd");
	  writer->SetInput(reader->GetOutput());
	  writer->Update();*/
	//mask
	reader2->SetFileName("/media/Users_UC3M/amunoz/LOTE_1/STUDY_4863/R8/16WEEKS_08SEP13/Chest_-_Lung_Reformat_0.625mm_interfile/LUNG/lesion34Slice.mhd");
	reader2->Update();
	InputImageType* imageMask = reader2->GetOutput();
	/*writer->SetFileName("/media/Users_UC3M/amunoz/LOTE_1/STUDY_4863/R8/16WEEKS_08SEP13/Chest_-_Lung_Reformat_0.625mm_interfile/LUNG/imageInMaybeNo2.mhd");
	writer->SetInput(imageIn);
	writer->Update();*/
	glcmGenerator->SetMaskImage(imageMask);

	glcmGenerator->SetPixelValueMinMax(-499,414);
	glcmGenerator->SetInsidePixelValue(34);
	glcmGenerator->SetNumberOfBinsPerAxis(16);
	//glcmGenerator->NormalizeOn();
    InputImageType::OffsetType offset1 = {{0, 1}};
    InputImageType::OffsetType offset2 = {{1, 0}};
    FilterType::OffsetVectorPointer offsetV = FilterType::OffsetVector::New();
    offsetV->push_back(offset1);
    //offsetV->push_back(offset2);

    glcmGenerator->SetOffsets(offsetV);
    //glcmGenerator->SetNormalize(true);

	glcmGenerator->Update();
	//glcmGenerator->NormalizeOn();



	/*std::cout << imageIn->GetPixel(testPixel)<<" " << imageIn->GetPixel(testPixel2)  << std::endl;
	std::cout << "Max " <<glcmGenerator->GetMax()<< " Mins "<<glcmGenerator->GetMin()  << std::endl;
	std::cout << "Inputs " <<glcmGenerator->GetNumberOfInputs()<< " Outputs "<<glcmGenerator->GetNumberOfOutputs()  << std::endl;*/

	int bins = glcmGenerator->GetNumberOfBinsPerAxis();
	const FilterType::HistogramType * hist = glcmGenerator->GetOutput();




    for( int i = 0; i < bins; i++){
    	for(int j = 0; j< bins; j++) {
    		typedef FilterType::HistogramType::IndexType IndexType;
    		IndexType i_j(2);
    		i_j[0] = i; i_j[1] =j;
    		float freq = hist ->GetFrequency(i_j);
    		if(freq > 0)
    		std::cout << "i " << i  << " j "<<j  << " hist: "<<freq << std::endl;
    	}
    }


    std::cout << hist->GetTotalFrequency()  << std::endl;

    typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType;
    	Hist2FeaturesType::Pointer featsCalc = Hist2FeaturesType::New();
    	featsCalc->SetInput(glcmGenerator->GetOutput());
    	featsCalc->Update();

    	std::cout << "Inertia: " << featsCalc->GetFeature(Hist2FeaturesType::Inertia) << " Energy: " << featsCalc->GetFeature(Hist2FeaturesType::Energy)
    			<< " Entropy "<< featsCalc->GetFeature(Hist2FeaturesType::Entropy) << std::endl;

    	offsetV->Initialize();
    	offsetV->push_back(offset2);
    	glcmGenerator->SetOffsets(offsetV);
        glcmGenerator->Update();
        Hist2FeaturesType::Pointer featsCalc2 = Hist2FeaturesType::New();
        featsCalc2->SetInput(glcmGenerator->GetOutput());
        featsCalc2->Update();

        std::cout << "Inertia: " << featsCalc2->GetFeature(Hist2FeaturesType::Inertia) << " Energy: " << featsCalc2->GetFeature(Hist2FeaturesType::Energy)
            			<< " Entropy "<< featsCalc2->GetFeature(Hist2FeaturesType::Entropy) << std::endl;


    std::cout << "this is the end "  << std::endl;


}

int itkScalarImageToTextureFeaturesFilter() {
	typedef signed short PixelType;
	const int Dimension = 3;

	typedef itk::Image<PixelType,Dimension> InputImageType;
	typedef itk::Image<unsigned char,Dimension> MaskImageType;
	typedef itk::ImageFileReader<InputImageType> ImageReaderType;
	typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
	typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<InputImageType> FilterType;

	typedef itk::Image<unsigned char, 2> LabelImageType;
	typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

	ImageReaderType::Pointer readerImage = ImageReaderType::New();
	MaskReaderType::Pointer readerMask = MaskReaderType::New();
	LabelReaderType::Pointer readerLabel = LabelReaderType::New();

	typedef itk::ImageFileWriter<InputImageType> TheDefinitiveWriterType;
	TheDefinitiveWriterType::Pointer theDefinitiveWriter = TheDefinitiveWriterType::New();


	readerImage->SetFileName("/tmp/E35HIGHDOSE3WEEKSPOSTINF_E35HIGHDOSE_20150325_181836_.mhd");
	//readerImage->Update();
	//theDefinitiveWriter->SetInput(readerImage->GetOutput());
	//theDefinitiveWriter->SetFileName("/tmp/test3.mhd");
	//theDefinitiveWriter->Update();

	readerMask->SetFileName("/media/Users_UC3M/amunoz/LOTE_1/STUDY_2087/E35_HIGHDOSE/3WEEKS_21APR12/LUNGS_SE3_interfile/SegmentationTest/label112.tif");
	//readerMask->Update();
	readerLabel->SetFileName("/media/Users_UC3M/amunoz/LOTE_1/STUDY_2087/E35_HIGHDOSE/3WEEKS_21APR12/LUNGS_SE3_interfile/LESIONS/mask_lesions_STUDY_2087_E35HIGHDOSE_21APRIL2012_3W-1.tif");
    readerLabel->Update();

    typedef itk::LabelShapeKeepNObjectsImageFilter< LabelImageType > LabelShapeKeepNObjectsImageFilterType;


/* --------------Cagarro--------------------------------------------------------------------------------------------------------*/
	 typedef itk::BinaryImageToLabelMapFilter<LabelImageType> BinaryImageToLabelMapFilterType;
	  BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
	  binaryImageToLabelMapFilter->SetInput(readerLabel->GetOutput());
	  binaryImageToLabelMapFilter->Update();
	  std::cout << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << std::endl;
/*---------------------------------------------------------------------------------------------------------------------------------*/


	  typedef itk::LabelImageToLabelMapFilter<LabelImageType> LabelImageToLabelMapFilterType;

	  LabelImageToLabelMapFilterType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapFilterType::New();

/*------------------------- labelImagetolabelmapfilter si más mola en este caso-------------------------------------------------------- */
	  labelImageToLabelMapFilter->SetInput(readerLabel->GetOutput());
	  labelImageToLabelMapFilter->Update();
	  std::cout << labelImageToLabelMapFilter->GetOutput() ->GetNumberOfLabelObjects() << std::endl;
	 itk::LabelMap<LabelImageToLabelMapFilterType::LabelObjectType >::LabelVectorType lv =  labelImageToLabelMapFilter->GetOutput()->GetLabels();
	 itk::LabelMap<LabelImageToLabelMapFilterType::LabelObjectType >::LabelVectorType::iterator it= lv.begin();
	 typedef LabelImageToLabelMapFilterType::LabelObjectType lobjectType;
	 unsigned char loboject = labelImageToLabelMapFilter->GetOutput()->GetLabels().at(0);
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

/* Enmascaro los enmascarado para ver como fucionan las etiquetas y porque soy tan anormal como molón-----------------------------*/
	 typedef itk::LabelObject< unsigned char, 2 >  LabelObjectType;
	 typedef itk::LabelMap< LabelObjectType >          LabelMapType;
	typedef itk::LabelMapMaskImageFilter< LabelMapType,  LabelImageType > LabelMapMaskFilterType;
	LabelMapMaskFilterType::Pointer labelMapMaskFilter = LabelMapMaskFilterType::New();


	labelMapMaskFilter->SetInput(labelImageToLabelMapFilter->GetOutput());
	labelMapMaskFilter->SetFeatureImage(readerLabel->GetOutput());

	  typedef itk::ImageFileWriter<LabelImageType> WriterType;

	 for (int i = 0; i < labelImageToLabelMapFilter->GetOutput() ->GetNumberOfLabelObjects(); i++) {
		 labelMapMaskFilter->SetLabel(labelImageToLabelMapFilter->GetOutput()->GetLabels().at(i));
			  WriterType::Pointer writer = WriterType::New();
			  std::ostringstream fn;
			    fn << "/media/Users_UC3M/amunoz/LOTE_1/STUDY_2087/E35_HIGHDOSE/3WEEKS_21APR12/LUNGS_SE3_interfile/LESIONS/label" << i << ".tif";
			  writer->SetFileName(fn.str());
			  writer->SetInput(labelMapMaskFilter->GetOutput());
			  writer->Update();
	 }
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	 itk::LabelMap<LabelImageToLabelMapFilterType::LabelObjectType>::LabelObjectVectorType lov = labelImageToLabelMapFilter->GetOutput()->GetLabelObjects();
	 std::cout << labelImageToLabelMapFilter->GetOutput()->GetLabelObject(41)->GetNameOfClass() << std::endl;




/*-----------------------Muela texturil experimentil-----------------------------*/
	  FilterType::Pointer filter = FilterType::New();
	  	filter->SetInput(readerImage->GetOutput());
	  	//filter->SetMaskImage(readerMask->GetOutput());
	  	typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxType;
	  	MinMaxType::Pointer minMax = MinMaxType::New();

	  	typedef itk::CastImageFilter<MaskImageType,InputImageType > CastFilterType2;
	  	CastFilterType2::Pointer caster2 = CastFilterType2::New();
	  	caster2->SetInput(readerMask->GetOutput());


	  	typedef itk::ImageFileWriter<InputImageType> Writer2Type;
	  	Writer2Type::Pointer writer2 = Writer2Type::New();
	  	writer2->SetInput(readerImage->GetOutput());
	  	writer2->SetFileName("/tmp/test2.mhd");



	  	InputImageType::IndexType pix;
	  	pix[0] = 335;
	  	pix[1] = 360;
	  	pix[2] = 115;
	  	//caster2->Update();

	  	filter->SetMaskImage(caster2->GetOutput());

	  	filter->SetPixelValueMinMax(-500,500);
	  	filter->SetInsidePixelValue(112);



    InputImageType::OffsetType offset1 = {{0, 1}};
    //InputImageType::OffsetType offset2 = {{1, 0}};
    FilterType::OffsetVectorPointer offsetV = FilterType::OffsetVector::New();
    offsetV->push_back(offset1);
    //offsetV->push_back(offset2);
    filter->SetOffsets(offsetV);
    int bins = 16;
    filter->SetNumberOfBinsPerAxis(bins);
    clock_t begin = clock();
    filter->Update();
 	std::cout << "Pixel Val " << readerImage-> GetOutput( ) ->GetPixel(pix)  << std::endl;
	//writer2->Update();
    clock_t end = clock();
    std::cout << "Time: "<<double( end - begin)/CLOCKS_PER_SEC << std::endl;

    minMax->SetImage(readerImage->GetOutput());
    	  	minMax->Compute();
    	  	std::cout << "Min: "<< minMax->GetMinimum() << " Max: "<<minMax->GetMaximum() << std::endl;

    typedef FilterType::FeatureValueVector FeatureVectorType;
    //typedef itk::VectorContainer<unsigned char,double >::Iterator VectorItType;
    //itk::VectorContainer<unsigned char,double > *feats= filter->GetFeatureMeans();

    const FilterType::FeatureValueVector* output = filter->GetFeatureMeans();

     for(unsigned int i = 0; i < output->size(); ++i)
      {
      std::cout << (*output)[i] << std::endl;
      }

    FeatureVectorType::Pointer featureVector = filter->GetFeatureMeans();
    itk::VectorContainer<unsigned char,double >::Iterator vectorIt = featureVector->Begin();
    while(vectorIt != featureVector->End()) {
        	std::cout << vectorIt->Value() << std::endl;
        	vectorIt++;
        }

    FeatureVectorType::Iterator feat = featureVector->Begin();
    while(feat != featureVector->End()) {
    	std::cout << feat->Value() << std::endl;
    	feat++;
    }

}




void objectbyobject(InternalImageType::Pointer inputImage) {}

int main(int argc, char*argv[]){
	//itkScalarImageToCooccurrenceMatrixFilterTest3();
	//itkScalarImageToTextureFeaturesFilter(); //Gets features from a segmented lession
	//itkScalarImageToCooccurrenceMatrixFilterTest2();
//std::cout << "Hello Man"<<std::endl;
	 char * infilename  = argv[1];
	  char * labels  = argv[2];
	  unsigned int bins= atoi(argv[3]);
	  char * outputfilename  = argv[4];
	  std::cout << "Processing " << infilename << " bins: " << bins << " Result " << outputfilename << std::endl;
     FeatureExtraction*  fe = new FeatureExtraction(infilename,labels,bins,outputfilename);
    // FeatureExtraction*  fe = new FeatureExtraction("/home/pmacias/Projects/MonkeysTuberculosis/QuatificationTest/imageMedian.mhd","/home/pmacias/Projects/MonkeysTuberculosis/QuatificationTest/Labels570.mhd",32,"/tmp/dummy");
     //fe->hello();
     //fe->hello("Shit");
     fe->getAllFeatures();
	 /* Test * t = new Test();
	  t->hello();
	  t->hello("Mothafucker");
	  FeatureExtraction::hello("Tus muertos");
	  std::list<int> l;
	  l.push_back(10);
	  l.push_front(5);
	  std::list<int>::iterator it = l.begin();
	  std::cout << *it << std::endl;*/








	//itkScalarImageToCooccurrenceMatrixFilterTest();
	/*
  if(argc < 2)
    {
    std::cerr << "Usage: " << argv[0] << " Required image.mha" << std::endl;
    return EXIT_FAILURE;
    }*/

  //std::string fileName = argv[1];
   // std::string fileName = "/home/pmacias/Projects/MonkeysTuberculosis/monoToTextureText.mhd";
/*
  typedef itk::ImageFileReader<InternalImageType> ReaderType;
  ReaderType::Pointer reader=ReaderType::New();
  reader->SetFileName(fileName);
  reader->Update();
  InternalImageType::Pointer image=reader->GetOutput();

  NeighborhoodType neighborhood;
  neighborhood.SetRadius(1);
  unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
  OffsetType offset;

  typedef itk::ImageFileWriter<InternalImageType> WriterType;
  WriterType::Pointer writer=WriterType::New();

  for ( unsigned int d = 0; d < centerIndex; d++ )
  {
      offset = neighborhood.GetOffset(d);
      InternalImageType::Pointer inertia=InternalImageType::New();
      InternalImageType::Pointer correlation=InternalImageType::New();
      InternalImageType::Pointer energy=InternalImageType::New();
      calcTextureFeatureImage(offset, image, inertia, correlation, energy);

      writer->SetInput(inertia);
//	snprintf(buf, 100, "Inertia%u.mha", d); // Warning: call to int __builtin___snprintf_chk will always overflow destination buffer
      std::stringstream ssInertia;
      ssInertia << "Inertia" << d << ".mhd";
      writer->SetFileName(ssInertia.str());
      writer->Update();
      writer->SetInput(correlation);
      std::stringstream ssCorrelation;
      ssCorrelation << "Correlation" << d << ".mhd";
      writer->SetFileName(ssCorrelation.str());
      writer->Update();
      writer->SetInput(energy);
      std::stringstream ssEnergy;
      ssEnergy << "Energy" << d << ".mhd";
      writer->SetFileName(ssEnergy.str());
      writer->Update();
      std::cout<<'\n';
  }
*/
  return EXIT_SUCCESS;
}



