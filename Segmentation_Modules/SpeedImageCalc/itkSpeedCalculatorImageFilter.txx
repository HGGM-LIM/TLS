/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSpeedCalculatorImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/07/17 12:26:43 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkSpeedCalculatorImageFilter_txx
#define _itkSpeedCalculatorImageFilter_txx

#include "itkSpeedCalculatorImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkGradientMagnitudeImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkSigmoidImageFilter.h"

#include "itkOtsuThresholdImageCalculator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkAverageImageFilter.h"

#include "itkHessian3DToAirwaynessMeasureImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"

#include "itkNumericTraits.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
SpeedCalculatorImageFilter<TInputImage, TOutputImage>
::SpeedCalculatorImageFilter()
{
 m_Scale = 0.1; 
}

template <typename TInputImage, typename TOutputImage>
void
SpeedCalculatorImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

template< typename TInputImage, typename TOutputImage >
void
SpeedCalculatorImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  typedef ImageRegionConstIterator< TInputImage > IteratorType;
  typedef ImageRegionIterator< TOutputImage> OutIteratorType;

  // get the input and output images
  InputImagePointer input = this->GetInput();

  //typename TOutputImage::Pointer output = this->GetOutput();
  //this->AllocateOutputs();


  
  // gradient 
  typedef GradientMagnitudeImageFilter< TInputImage, TInputImage > GradientType;
  typename GradientType::Pointer gradient = GradientType::New() ;
  
  gradient->SetInput( input );
  gradient->SetUseImageSpacingOn();
  
  // rescale and sigmoid. The output of the sigmoid can be seen as a probability of being non-edge voxel
  typedef RescaleIntensityImageFilter< TInputImage, TOutputImage > RescaleType;
  typename RescaleType::Pointer rescaler = RescaleType::New() ;
  rescaler->SetInput( gradient->GetOutput() );
  rescaler->SetOutputMaximum( 1 );  
  rescaler->SetOutputMinimum( 0 ); 

  typedef SigmoidImageFilter< TOutputImage, TOutputImage > SigmoidType;
  typename SigmoidType::Pointer sigmoid = SigmoidType::New();
  sigmoid->SetInput(rescaler->GetOutput());
  sigmoid->SetOutputMaximum( 1 );
  sigmoid->SetOutputMinimum( 0 );
  // K1=0.2 (value in airway), K2=0.7 (value in border)
  // alpha=(k2-k1)/6; beta=mean(k1,k2)
  sigmoid->SetAlpha(-0.083);
  sigmoid->SetBeta(0.45);
  sigmoid->Update();
  // use the following line to view result of sigmoid
  //this->GraftOutput( sigmoid->GetOutput() );

  //separately, we can compute the intensity probability
  // we compute the Otsu threshold and use this as reference
  // to compute probability of being air/tissue
  typedef OtsuThresholdImageCalculator< TInputImage > OtsuCalculatorType;
  typename OtsuCalculatorType::Pointer otsu = OtsuCalculatorType::New();
  otsu->SetImage( input );
  otsu->Compute();
  std::cout << "Otsu threshold: " << otsu->GetThreshold() << std::endl;
  InputPixelType threshold = otsu->GetThreshold();
  // compute intensity probability point by point
  // we need iterators for that
  typedef ImageRegionConstIterator< TInputImage > InputIteratorType;
  InputIteratorType it = InputIteratorType( input, input->GetLargestPossibleRegion());
  typename TOutputImage::Pointer intensityProb = TOutputImage::New();
  intensityProb->SetRegions( input->GetLargestPossibleRegion() );
  intensityProb->Allocate();
  intensityProb->FillBuffer(0);
  typedef ImageRegionIterator< TOutputImage> OutputIteratorType;
  OutputIteratorType intIt = OutputIteratorType( intensityProb, input->GetLargestPossibleRegion());
  
  // we compute the parameters for the probability line
  float a = 1.0/(-950.0-threshold);
  float b = -a*threshold;  

  // loop to compute probability of each voxel
  for(it.GoToBegin(),intIt.GoToBegin(); !it.IsAtEnd(); ++it, ++intIt)  
    {
    if(it.Get()<-950)
      {
      intIt.Set(1); // it is air for sure
      }
    else if(it.Get()>threshold)
      {
      intIt.Set(0); // it is body for sure
      }
    else
      {
      intIt.Set(a*it.Get()+b); // line between 1 and 0.5 probability
      }
    }
  // use this line to visualize the output of the intensity probability
  //this->GraftOutput( intensityProb );
  
  // now we compute the airwayness, first computing the hessian and then
  // with the airwayness filter, which considers the eigenvalues
  typedef   itk::HessianRecursiveGaussianImageFilter< 
                            TInputImage >              HessianFilterType;
  typedef Hessian3DToAirwaynessMeasureImageFilter<
              typename TOutputImage::PixelType > AirwaynessMeasureFilterType;
  typedef RescaleIntensityImageFilter< OutputImageType, OutputImageType >  RescaleFilterType;
  
  typename HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
  typename AirwaynessMeasureFilterType::Pointer airwaynessFilter = 
                            AirwaynessMeasureFilterType::New();
  typename RescaleFilterType::Pointer rescaler_airwayness = RescaleFilterType::New();
  hessianFilter->SetSigma( this->GetScale() );
  rescaler_airwayness->SetOutputMinimum( 0 );
  rescaler_airwayness->SetOutputMaximum( 1 ); 
  hessianFilter->SetInput( input );
  airwaynessFilter->SetInput( hessianFilter->GetOutput() );
  rescaler_airwayness->SetInput( airwaynessFilter->GetOutput() );
  rescaler_airwayness->Update();
  
  // to view tubularness
  //this->GraftOutput( rescaler_airwayness->GetOutput() );

  // now we can combine the three by simply averaging
  typedef AverageImageFilter< TOutputImage, TOutputImage > AverageFilterType;
  typename AverageFilterType::Pointer averager = AverageFilterType::New();
  averager->SetInput( 0, sigmoid->GetOutput());
  averager->SetInput( 1, intensityProb);
  averager->SetInput( 2, rescaler_airwayness->GetOutput());
  averager->Update();
  this->GraftOutput( averager->GetOutput() );
  
}

} // end namespace itk

#endif
