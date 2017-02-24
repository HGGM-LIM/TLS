/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanAbsoluteSurfaceDistanceImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:50 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMeanAbsoluteSurfaceDistanceImageFilter_txx
#define _itkMeanAbsoluteSurfaceDistanceImageFilter_txx
#include "itkMeanAbsoluteSurfaceDistanceImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkProgressAccumulator.h"
#include "itkMeanAbsoluteSurfaceDistanceImageFilter.h"

#include "itkNotImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"

#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"


namespace itk {


template<class TInputImage1, class TInputImage2>
MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::MeanAbsoluteSurfaceDistanceImageFilter()
{

  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );

  m_MeanAbsoluteSurfaceDistance = NumericTraits<RealType>::Zero;      
}


template<class TInputImage1, class TInputImage2>
void
MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::SetInput2( const TInputImage2 * image )
{
  this->SetNthInput(1, const_cast<TInputImage2 *>( image ) );      
}


template<class TInputImage1, class TInputImage2>
const typename MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::InputImage2Type *
MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::GetInput2()
{
  return static_cast< const TInputImage2 * >
    (this->ProcessObject::GetInput(1));
}



template<class TInputImage1, class TInputImage2>
void
MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // this filter requires:
  // - the largeset possible region of the first image
  // - the corresponding region of the second image
  if ( this->GetInput1() )
    {
    InputImage1Pointer image =
      const_cast< InputImage1Type * >( this->GetInput1() );
    image->SetRequestedRegionToLargestPossibleRegion();

    if ( this->GetInput2() )
      {
      InputImage2Pointer image =
        const_cast< InputImage2Type * >( this->GetInput2() );
      image->SetRequestedRegion( 
        this->GetInput1()->GetRequestedRegion() );
      }

    }
}


template<class TInputImage1, class TInputImage2>
void
MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::EnlargeOutputRequestedRegion(DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}


template<class TInputImage1, class TInputImage2>
void
MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::GenerateData()
{
  // first get the two images
  InputImage1Pointer image1 =
        const_cast< InputImage2Type * >( this->GetInput1() );
  
  InputImage2Pointer image2 =
        const_cast< InputImage2Type * >( this->GetInput2() );
	
  // define the necessary filters  
  typedef itk::NotImageFilter< TInputImage1,  TInputImage1>  NotFilterType1;
  typedef itk::NotImageFilter< TInputImage2,  TInputImage2>  NotFilterType2;
  
  typedef itk::DanielssonDistanceMapImageFilter<
               TInputImage1,  TInputImage1>  DistanceFilterType1;
  typedef itk::DanielssonDistanceMapImageFilter<
               TInputImage1,  TInputImage1>  DistanceFilterType2;	       
  
  typedef itk::AddImageFilter<
               TInputImage1, TInputImage1, TInputImage1 >  AddFilterType1;
  typedef itk::AddImageFilter<
               TInputImage1, TInputImage1, TInputImage1 >  AddFilterType2;
  
  typename DistanceFilterType1::Pointer distanceFilter1 = DistanceFilterType1::New();
  typename DistanceFilterType1::Pointer distanceFilter2 = DistanceFilterType1::New();
  distanceFilter1->InputIsBinaryOn();
  distanceFilter2->InputIsBinaryOn();
  
  typename NotFilterType1::Pointer notFilter1 = NotFilterType1::New();
  
  typename AddFilterType1::Pointer addFilter1 = AddFilterType1::New();
  
  // calculate distance map for image 1 
  distanceFilter1->SetInput( image1 );
  typename TInputImage1::Pointer distMapExt1 = distanceFilter1->GetOutput();
  distanceFilter1->Update();
    
  notFilter1->SetInput( image1 );
  distanceFilter2->SetInput( notFilter1->GetOutput() );
  typename TInputImage1::Pointer distMapInt1 = distanceFilter2->GetOutput();
  distanceFilter2->Update();
  
  //combine the two distance maps
  addFilter1->SetInput1( distanceFilter1->GetOutput() );
  addFilter1->SetInput2( distanceFilter2->GetOutput() );
  typename TInputImage1::Pointer distMap1 = addFilter1->GetOutput();
  addFilter1->Update();
  
  // calculate distance map for image 2 
  typename DistanceFilterType2::Pointer distanceFilter3 = DistanceFilterType1::New();
  typename DistanceFilterType2::Pointer distanceFilter4 = DistanceFilterType1::New();
  distanceFilter3->InputIsBinaryOn();
  distanceFilter4->InputIsBinaryOn();
  
  typename NotFilterType2::Pointer notFilter2 = NotFilterType2::New();
  
  typename AddFilterType2::Pointer addFilter2 = AddFilterType2::New();
  
  distanceFilter3->SetInput( image2 );
  typename TInputImage2::Pointer distMapExt2 = distanceFilter3->GetOutput();
  distanceFilter3->Update();
    
  notFilter2->SetInput( image2 );
  distanceFilter4->SetInput( notFilter2->GetOutput() );
  typename TInputImage2::Pointer distMapInt2 = distanceFilter4->GetOutput();
  distanceFilter4->Update();
  
  //combine the two distance maps
  addFilter2->SetInput1( distanceFilter3->GetOutput() );
  addFilter2->SetInput2( distanceFilter4->GetOutput() );
  typename TInputImage1::Pointer distMap2 = addFilter2->GetOutput();
  addFilter2->Update();
  
  // define iterators
  typedef itk::ImageRegionIterator<TInputImage1> ImageIteratorType1;
  typedef itk::ImageRegionIterator<TInputImage2> ImageIteratorType2;
  
  ImageIteratorType1 distExt1It(distMapExt1, distMapExt1->GetLargestPossibleRegion());
  ImageIteratorType1 distInt1It(distMapInt1, distMapInt1->GetLargestPossibleRegion());
  ImageIteratorType2 distExt2It(distMapExt2, distMapExt2->GetLargestPossibleRegion());
  ImageIteratorType2 distInt2It(distMapInt2, distMapInt2->GetLargestPossibleRegion());
  ImageIteratorType1 distMap1It(distMap1, distMap1->GetLargestPossibleRegion());
  ImageIteratorType2 distMap2It(distMap2, distMap2->GetLargestPossibleRegion());
  
  // iterate through images, finding the distance from image 2 to image 1 and from image 1 to to
  float border_voxels1 = 0;
  float border_voxels2 = 0;
  float sum1 = 0;
  float sum2 = 0;
  
  for (distExt1It.GoToBegin(),distInt1It.GoToBegin(),distMap1It.GoToBegin(),distExt2It.GoToBegin(),distInt2It.GoToBegin(),distMap2It.GoToBegin();!distExt1It.IsAtEnd()||!distInt1It.IsAtEnd();++distExt1It,++distInt1It,++distExt2It,++distInt2It,++distMap1It,++distMap2It)
    {
    // First if loop to measure distance from 2 to 1
    //check if we are in the border of image 2
    if (distExt2It.Get()==0 && distInt2It.Get() == 1)
      {
      //check if we are in border of image 1. If not, get the distance and add to the sum
      if (distExt1It.Get()!=0 || distInt1It.Get() != 1)
        {
        sum1 = sum1 + int(distMap1It.Get());
	}
      border_voxels1++;
      }
         
    // Second if loop to measure distance from 1 to 2
    if(distExt1It.Get()==0 && distInt1It.Get() == 1)
      {
      //check if we are in border of image 1. If not, get the distance and add to the sum
      if (distExt2It.Get()!=0 || distInt2It.Get() != 1)
        {
        sum2 = sum2 + int(distMap2It.Get());
	}
      
      border_voxels2++;
      } 
      
    }
  
  float result1 = sum1/border_voxels1;
  float result2 = sum2/border_voxels2;
  m_MeanAbsoluteSurfaceDistance= 0.5*( result1 + result2 );
  
}



template<class TInputImage1, class TInputImage2>
void 
MeanAbsoluteSurfaceDistanceImageFilter<TInputImage1, TInputImage2>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Mean Absolute Surface Distance (MASD): "  
     << m_MeanAbsoluteSurfaceDistance << std::endl;
}


}// end namespace itk
#endif
