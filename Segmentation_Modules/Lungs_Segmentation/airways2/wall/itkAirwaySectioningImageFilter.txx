/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.14 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAirwaySectioningImageFilter_txx
#define _itkAirwaySectioningImageFilter_txx
#include "itkAirwaySectioningImageFilter.h"
#include <string>
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkBSplineInterpolateImageFunction.h"

#include "itkBinaryFillholeImageFilter.h"

namespace itk
{

template <class TInputImage>
AirwaySectioningImageFilter<TInputImage>
::AirwaySectioningImageFilter()
{
  m_OutputSize.Fill( 0 );
  m_OutputSpacing.Fill( 1.0 );
  m_OutputOrigin.Fill( 0.0 );
  m_Center.Fill( 0.0 );
  m_Direction.Fill( 0.0 );
  m_Up.Fill( 0.0 );
}

template < class TImageType> 
void AirwaySectioningImageFilter< TImageType >::GenerateOutputInformation()
{
  // Do not call the superclass' implementation of this method since
  // the output is of different dimension
 
  // Get pointers to the input and output
  typename Superclass::OutputImagePointer     outputPtr = this->GetOutput();
  typename Superclass::InputImageConstPointer inputPtr  = this->GetInput();

  if ( !outputPtr || !inputPtr )
    {
    return;
    }

  // Set the output image size
  OutputImageRegionType outputRegion;
  typedef typename OutputImageRegionType::SizeType OutputRegionSizeType;
  typedef typename OutputImageRegionType::IndexType OutputRegionIndexType;
  OutputRegionSizeType outSize;
  //OutputRegionIndexType outIndex;
  outSize[0] = inputPtr->GetLargestPossibleRegion().GetSize()[0];
  outSize[1] = inputPtr->GetLargestPossibleRegion().GetSize()[1];  
  outSize[2] = 1;
  outputRegion.SetSize( outSize );
  outputPtr->SetRegions( outputRegion );

  typename InputImageType::SpacingType spacing;
  spacing = inputPtr->GetSpacing();

  // 061010 - try to 'fix' spacing issue
  m_OutputSpacing[0] = spacing[0];
  m_OutputSpacing[1] = spacing[1];
  m_OutputSpacing[2] = spacing[2];

  
}




template< class TInputImage >
void
AirwaySectioningImageFilter< TInputImage>
::GenerateData()
{
 // int file_segment=0;
  // TODO: Add progress accumulator
  // Get pointers to the input and out
  typename Superclass::InputImageConstPointer input  = this->GetInput();
  
  typename TransformType::Pointer transform = TransformType::New();
  
  // Setup output parameters
  typename InputImageType::SizeType size;
  
  //typename InputImageType::PointType origin;
  //double origin[3];
  //origin[0] = origin[1] =  origin[2] = 0; 

  size = input->GetLargestPossibleRegion().GetSize(); //sets up the output size both filter-output and final image.
  size[2] = 1;
//  cout << "itkAirwaySectioning ... 1, ";

  typename InputImageType::SpacingType spacing;
  spacing = input->GetSpacing();

  // 061010 - try to 'fix' spacing issue
  m_OutputSpacing[0] = spacing[0];
  m_OutputSpacing[1] = spacing[1];
  m_OutputSpacing[2] = spacing[2];


//  cout << "2, ";
  // Create transform initializer
  typename TransformInitalizerType::PointType center;
  center[0] = m_Center[0]; 
  center[1] = m_Center[1];
  center[2] = m_Center[2];
//  cout << "3, ";
  // Added by RINA 160910
  // moved here - to before TransformInitalizer is called 310111
//  center[0]= center[0] * spacing[0];// 200 * spacing[0];
//  center[1]= center[1] * spacing[1]; // 200 * spacing[1];
//  center[2]= center[2] * spacing[2];

  typename TransformInitalizerType::VectorType direction;
  direction[0] = m_Direction[0]; 
  direction[1] = m_Direction[1]; 
  direction[2] = m_Direction[2];
//  cout << "4, ";
  
  typename TransformInitalizerType::VectorType up;
  up[0] = m_Up[0]; 
  up[1] = m_Up[1]; 
  up[2] = m_Up[2];
//  cout << "5, ";

  typename TransformInitalizerType::Pointer init = TransformInitalizerType::New();
  init->SetImage( input );
//  cout << "6, ";
  init->SetTransform( transform );
//  cout << "7, ";

//  cout << "in AirwaySectioning, ctr of transform: "<< center[0]  << ", "<< center[1] << ", " << center[2] << "...";
//end of addition
  init->SetPlane( center, direction, up, size );
//  cout << "9, ";
  init->InitializeTransform();
//  cout << "10, ";

  // Create resample filter
  typename ResampleFilterType::Pointer filterResample = ResampleFilterType::New();
  filterResample->SetInput( input );
  filterResample->SetTransform( transform );
  filterResample->SetSize( size );

  // 3degree bspline interpolator (default is 3 if not set explicitly)
//  typedef itk::BSplineInterpolateImageFunction<
//                       InputImageType, double >  InterpolatorType;
//  InterpolatorType::Pointer interpolator = InterpolatorType::New();
//  filterResample->SetInterpolator( interpolator );

//  cout << "11, ";
  // Added by RINA 160910
  typename TransformInitalizerType::PointType outputOrigin;
  outputOrigin[0] = 0;
  outputOrigin[1] = 0;
  outputOrigin[2] = 0;
  filterResample->SetOutputOrigin( outputOrigin );
  //end of addition
//  cout << "12, ";

  filterResample->SetOutputSpacing( spacing );
  filterResample->Update( );

//  cout << "13, ";

   typename OutputImageType::Pointer new_image = OutputImageType::New();
   new_image->SetRegions(filterResample->GetOutput()->GetRequestedRegion());
   // newly added back - RINA - 310111
   new_image->SetSpacing(filterResample->GetOutputSpacing());
   new_image->Allocate();

//   cout << "14, ";
  typename GaussianFilterType::Pointer GaussianFilter = GaussianFilterType::New();
  GaussianFilter->SetInput(filterResample->GetOutput());
//  GaussianFilter->SetUseImageSpacingOn();
//  cout << "15, ";
  //disable it for mice images.
  //GaussianFilter->SetVariance( 1 );
  //GaussianFilter->SetMaximumKernelWidth( 3 );

  
  
  GaussianFilter->Update(); 
//  cout << "16 ";

  if(m_SliceType==1)
    {
//	  cout << "17, ";
      
 /*   typename WriterType_3D::Pointer writer4 = WriterType_3D::New();
    writer4->SetInput(new_image);
    writer4->SetFileName("prova2.tiff");
    writer4->Update();
    cout << filterResample;*/
//    int tmp;
//    imagecin >> tmp;
    ConstIteratorType constIterator( GaussianFilter->GetOutput(), GaussianFilter->GetOutput()->GetRequestedRegion() );
    IteratorType iterator( new_image, new_image->GetRequestedRegion() );
    OutputPixelType Pixel;
    typename TInputImage::IndexType Pos;
    for(int i=0; i<(size[0]*size[1]);i++)
      {
      Pos = constIterator.GetIndex();
      
      if((Pos[0]>= (size[0]/2) - m_RadiusToExtract-8)&&(Pos[0] <=((size[0]/2) - m_RadiusToExtract-8)+4.25*m_RadiusToExtract)&&(Pos[1]>= (size[1]/2) - m_RadiusToExtract-8)&&(Pos[1] <=((size[1]/2) - m_RadiusToExtract-8)+4.25*m_RadiusToExtract))
        {
	Pixel = constIterator.Get();
	//cout << Pixel << endl;
        iterator.Set(Pixel);
        }
      else
        {
        iterator.Set(-2000);
        }

      iterator.operator++();
      constIterator.operator++();
      //std::cout << "El lector apunta a: " << Pos << std::endl;
      
      }

    // TODO verify - WHY SMOOTHING?
      //We smooth the new image beacuse of the low dose.

    
/*

  typename HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(new_image);


//calcular maximo y minimo con el nuevo filtro.

  int nob = 256; 
  double ms = 10.0;
  double hmin=1000;
  double hMAX =0;


  ConstIteratorType it( new_image, new_image->GetRequestedRegion() );
  it = it.Begin();
  for (; !it.IsAtEnd(); ++it)
   {
    if(it.Get() > hMAX)
      {
      hMAX=double(it.Get());
      }
    if(it.Get() < hmin)
      {
      hmin=double(it.Get());
      }
   
   }

 //  cout << "HMAX: " << hMAX << " Hmin: " << hmin << endl;
  
  histogramGenerator->SetNumberOfBins( nob );
  histogramGenerator->SetMarginalScale( ms );

  histogramGenerator->SetHistogramMin(  hmin );
  histogramGenerator->SetHistogramMax( hMAX );
  histogramGenerator->Compute();



  const HistogramType * histogram = histogramGenerator->GetOutput();


  unsigned int binNumber = 0;
  unsigned long freq=0;
  unsigned long cumfreq=0;
  unsigned long totfreq=histogram->GetTotalFrequency();
  long th = 0;
  double percent = 0;

  typename HistogramType::ConstIterator itr = histogram->Begin();
  typename HistogramType::ConstIterator end = histogram->End();
  m_Histogram.str("");
  m_Histogram << "," << hmin << "," << hMAX << std::endl;
  while( itr != end )
    {
    freq = (unsigned long)itr.GetFrequency();
    cumfreq = cumfreq + freq;
    percent = (cumfreq * 100.0) /totfreq;
    th = (long)(hmin + binNumber*((hMAX-hmin)/nob));
    m_Histogram << percent << "%, " << binNumber << ", " << freq << ", " << th << std::endl;     
 //   cout << percent << "%, " << binNumber << ", " << freq << ", " << th << std::endl;     
    ++itr;
    ++binNumber;
    }
*/
      this->GraftOutput( new_image );
   }
 
  else if(m_SliceType==0)
    {
//    cout << "17B, ";
      typename ConnectedFilter::Pointer connected = ConnectedFilter::New();
      connected->SetInput(filterResample->GetOutput() );
    

   /* typename WriterType_3D::Pointer writer4 = WriterType_3D::New();
      writer4->SetFileName("corte_prueba.tiff");
      writer4->SetInput(filterResample->GetOutput());
      writer4->Update();*/

 

      connected->Update();
      OutputImageType::Pointer conn_image = connected->GetOutput();

      typedef itk::MinimumMaximumImageCalculator< InputImageType > MinCalcType;
      MinCalcType::Pointer minCalc2 = MinCalcType::New();
      minCalc2->SetImage(connected->GetOutput());
      minCalc2->Compute();
      InputPixelType globalMin, globalMax;
      globalMin = minCalc2->GetMinimum();
      globalMax = minCalc2->GetMaximum();

      typename InputImageType::IndexType pixelIndex;

      pixelIndex[0] = size[0]/2;   // x position 282
      pixelIndex[1] = size[1]/2;   // y position 223
      pixelIndex[2] = 0; // z position
    
      // RINA original was InputImageType::PixelType - which is unsigned int. check for effects in other experiments!!!!!!!
//      typename InputImageType::PixelType pixelValue = conn_image->GetPixel( pixelIndex );
      InputPixelType pixelValue = conn_image->GetPixel( pixelIndex );
//      cout << "DEBUG "<< "Init pixelValue: " << pixelValue << endl;
      int i = 0;
      int direction_p = 0; //0->left x, 1->right x, 2->down y, 3->up 
      int distance=5;

      bool emergency = false;
//      while(pixelValue==0)
      while(pixelValue<=0)
      {

    	  //First set new direction

    	  if((direction_p==0)&&(i==distance))
    	  {
    		  pixelIndex[0] = pixelIndex[0] +distance;
    		  direction_p=1;
    		  //std::cout << "Change direction x+ to x-" << std::endl;
    		  i=0;
    	  }
    	  if((direction_p==1)&&(i==distance))
    	  {
    		  pixelIndex[0] = pixelIndex[0] -distance;
    		  direction_p=2;
    		  //std::cout << "Change direction x- to y+" << std::endl;
    		  i=0;
    	  }
    	  if((direction_p==2)&&(i==distance))
    	  {
    		  pixelIndex[1] = pixelIndex[1] +distance;
    		  direction_p=3;
    		  //std::cout << "Change direction y+ to y-" << std::endl;
    		  i=0;
    	  }
    	  if((direction_p==3)&&(i==distance))
    	  {
    		  pixelIndex[1] = pixelIndex[1] -distance;
    		  direction_p=0;
    		  distance++;
    		  //std::cout << "Change direction y+ to y-" << std::endl;
    		  i=0;
    	  }

    	  //then get the pixel.
    	  if(direction_p==0)
    	  {
    		  pixelIndex[0]--;
    	  }
    	  if(direction_p==1)
    	  {
    		  pixelIndex[0]++;
    	  }
    	  if(direction_p==2)
    	  {
    		  pixelIndex[1]--;
    	  }
    	  if(direction_p==3)
    	  {
    		  pixelIndex[1]++;
    	  }

    	  //std::cout << "i: " << i << " direction: " << direction_p <<  std::endl;
    	  //cout << conn_image;
    	  i++;
    	  //cout  << pixelIndex;
    	  cout.flush();

    	  /* comment by RINA -- added over & under conditions to check for indexOutOfBound error
    	   * Why hasn't it bombed before?!
    	   */

    	  //	if( (pixelIndex[0] <= 0) || (pixelIndex[1] <= 0) || (pixelIndex[2] < 0) ){
    	  if( (pixelIndex[0] < 0) || (pixelIndex[1] < 0) || (pixelIndex[2] < 0) || (pixelIndex[0] >= size[0]) || (pixelIndex[1] >= size[1]) ){
//    		  cout << "Pixel index out of bounds! " << pixelIndex[0] << "," << pixelIndex[1] << endl;

    		  //		typename TInputImage::Pointer new_image = filterResample->GetOutput();
    		  //
    		  //		typename WriterType_3D::Pointer writer4 = WriterType_3D::New();
    		  // 280910 - Rina commented the next 2 lines
    		  //	  cout << "Please input a pixel value to continue" << endl;
    		  //	  cin >> pixelValue;
    		  //		pixelValue = 1;
    		  // if the index is out of bounds, the value shd be 0 so that loop looks for other values!
//    		  pixelValue = 0;
    		  // 290910 - decided to skip this & use max value instead
    		  // Q: what if globalMax is not the first non-zero value encountered?
    		  // Solution is not correct. When the thresholded image is not found in the center of the image, something is wrong
//    		  pixelValue = globalMax;
    		  pixelValue = 1; // why not 1?
    		  emergency = true;
    	  } else
    		  pixelValue = conn_image->GetPixel( pixelIndex );

      }
//      cout << "DEBUG "<< "GlobalMin: " << globalMin << ", GlobalMax: " << globalMax << endl;
//      std::cout << "m_SliceType: "<< m_SliceType <<", from index: "<< pixelIndex[0] << "," << pixelIndex[1] << ", Pixel Value: " << pixelValue << std::endl;
    
      typename FilterType::Pointer filter = FilterType::New();

      filter->SetInput(connected->GetOutput());
      filter->SetOutsideValue( 0 );
      filter->SetInsideValue( 1 ); 
      if (emergency)
    	  pixelValue = 1;

      filter->SetLowerThreshold( pixelValue);
      filter->SetUpperThreshold( pixelValue);
      filter->Update();
      //Now i have the central airway segmented and labeled as 1.

// to fill holes

//      typedef itk::BinaryFillholeImageFilter< InputImageType > FillHoles;
//      FillHoles::Pointer filler =  FillHoles::New();
//      filler->SetForegroundValue(255);
//      filler->SetInput(filter->GetOutput());
//      filler->Update();


      typename ReLabelType::Pointer relabelfilter = ReLabelType::New();
      relabelfilter->SetInput(filter->GetOutput());
      relabelfilter->Update();

    
      m_Area= relabelfilter->GetSizeOfObjectsInPhysicalUnits()[0];
     // std::cout << "Slice area: " << m_Area << std::endl;

      if (emergency) {
//    	  cout << "ALERT!!! TO CHECK MANUALLY IF SlicesB is centered! DISCARDING....." << endl;
    	  this->GraftOutput(connected->GetOutput());
    	  return;
      } else
//    	  this->GraftOutput(filler->GetOutput());
    	  this->GraftOutput( filter->GetOutput());

     }

}

template <class TInputImage> int AirwaySectioningImageFilter<TInputImage>::GetArea()
{
  return m_Area;
}

/*
template <class TInputImage> string AirwaySectioningImageFilter<TInputImage>::GetHistogram()
{
  return (m_Histogram.str());
}
*/

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage>
void
AirwaySectioningImageFilter<TInputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  
}

} // end namespace itk

#endif
