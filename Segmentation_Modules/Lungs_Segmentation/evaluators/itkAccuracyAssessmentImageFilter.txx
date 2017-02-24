/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAccuracyAssessmentImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:50 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAccuracyAssessmentImageFilter_txx
#define _itkAccuracyAssessmentImageFilter_txx
#include "itkAccuracyAssessmentImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkProgressAccumulator.h"
#include "itkAccuracyAssessmentImageFilter.h"

#include "itkShiftScaleImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"


namespace itk {


template<class TInputImage1, class TInputImage2>
AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::AccuracyAssessmentImageFilter()
{

  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );

  m_TruePositiveVF = NumericTraits<RealType>::Zero;      
  m_TrueNegativeVF = NumericTraits<RealType>::Zero;      
  m_FalsePositiveVF = NumericTraits<RealType>::Zero;      
  m_FalseNegativeVF = NumericTraits<RealType>::Zero;      
  
}


template<class TInputImage1, class TInputImage2>
void
AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::SetTestSegmentation( const TInputImage2 * image )
{
  this->SetNthInput(1, const_cast<TInputImage2 *>( image ) );      
}


template<class TInputImage1, class TInputImage2>
const typename AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::InputImage2Type *
AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::GetTestSegmentation()
{
  return static_cast< const TInputImage2 * >
    (this->ProcessObject::GetInput(1));
}



template<class TInputImage1, class TInputImage2>
void
AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // this filter requires:
  // - the largeset possible region of the first image
  // - the corresponding region of the second image
  if ( this->GetReferenceSegmentation() )
    {
    InputImage1Pointer image =
      const_cast< InputImage1Type * >( this->GetReferenceSegmentation() );
    image->SetRequestedRegionToLargestPossibleRegion();

    if ( this->GetTestSegmentation() )
      {
      InputImage2Pointer image =
        const_cast< InputImage2Type * >( this->GetTestSegmentation() );
      image->SetRequestedRegion( 
        this->GetReferenceSegmentation()->GetRequestedRegion() );
      }

    }
}


template<class TInputImage1, class TInputImage2>
void
AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::EnlargeOutputRequestedRegion(DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}


template<class TInputImage1, class TInputImage2>
void
AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::GenerateData()
{
  // first get the two images
  InputImage1Pointer refSegmentation =
        const_cast< InputImage1Type * >( this->GetReferenceSegmentation() );
  
  InputImage2Pointer testSegmentation =
        const_cast< InputImage2Type * >( this->GetTestSegmentation() );
	
  // define the necessary filters  
  typedef itk::ShiftScaleImageFilter< TInputImage1,  TInputImage1>  ShiftScaleFilterType;
  typename ShiftScaleFilterType::Pointer shiftFilter = ShiftScaleFilterType::New();
  
  shiftFilter->SetScale(2.0);
  shiftFilter->SetInput(refSegmentation);
  
  typedef itk::AddImageFilter< TInputImage1, TInputImage2, TInputImage1 >  AddFilterType;
  typename AddFilterType::Pointer addFilter = AddFilterType::New();
  
  addFilter->SetInput1(shiftFilter->GetOutput());
  addFilter->SetInput2(testSegmentation);
  addFilter->Update();
  InputImage1Pointer composedImage = const_cast< InputImage1Type * >( addFilter->GetOutput() );
  
  typedef itk::ImageRegionIterator<TInputImage1> ImageIteratorType;
  ImageIteratorType imageIt(composedImage, composedImage->GetLargestPossibleRegion());
  
  unsigned int FP = 0;
  unsigned int FN = 0;
  unsigned int TP = 0;
  unsigned int TN = 0;
  unsigned int totalVoxels =0;
  for (imageIt.GoToBegin();!imageIt.IsAtEnd();++imageIt,++totalVoxels)
    {
    switch ((int)imageIt.Get())
      {
      case 0:
    	++TN;
	break;
      case 1:
    	++FP;
	break;
      case 2:
    	++FN;
	break;
      case 3:
    	++TP;
	break;
      default: 
    	std::cout<< "Unexpected value. Please check that two inputs are correct (they must be binary!)";
	break;
      }
    }

  m_TruePositiveVF = (double)TP/(TP+FN);
  m_TrueNegativeVF = (double)TN/(totalVoxels-TP-FN);
  //m_TrueNegativeVF = (double)TN/(TP+FN);
  //m_FalsePositiveVF = 1 - m_TrueNegativeVF;
  m_FalsePositiveVF = (double)FP/(TP+FN);
  m_FalseNegativeVF = 1 -m_TruePositiveVF;
  
  //////////FOR DEBUGGING///////////////////////////////////////
  
  std::cout<< "FP: " << FP << std::endl;
  std::cout<< "FN: " << FN << std::endl;
  std::cout<< "TP:" << TP << std::endl;
  std::cout<< "TN: "<<TN << std::endl;
  std::cout<< m_TrueNegativeVF << std::endl;
  std::cout<< m_TruePositiveVF << std::endl;std::cout<< totalVoxels << std::endl;
  
  /*typedef itk::Image< unsigned short, 3 >  OutputImageType;
  typedef itk::ImageFileWriter < OutputImageType >  WriterType;

  typename WriterType::Pointer writer = WriterType::New();
  

  typedef itk::CastImageFilter< TInputImage1, OutputImageType > CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput(composedImage);
  writer->SetFileName( "/home/xabiarta/ITK/AccuracyTest.tif" );
  writer->SetInput(caster->GetOutput());
  writer->Update();*/	       
}



template<class TInputImage1, class TInputImage2>
void 
AccuracyAssessmentImageFilter<TInputImage1, TInputImage2>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "True Positive Volume Fraction): " 
     << m_TruePositiveVF << std::endl;
  os << indent << "True Negative Volume Fraction): " 
     << m_TrueNegativeVF << std::endl;     
  os << indent << "False Positive Volume Fraction): " 
     << m_FalsePositiveVF << std::endl;     
  os << indent << "False Negative Volume Fraction): " 
     << m_FalseNegativeVF << std::endl;
}


}// end namespace itk
#endif
