/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkObliqueSectionImageFilter.txx,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkObliqueSectionImageFilter_txx
#define _itkObliqueSectionImageFilter_txx

#include "itkObliqueSectionImageFilter.h"
#include "itkNumericTraits.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkLookAtTransformInitializer.h"

namespace itk
{

/**
 *
 */
template < class TInputPixel, class TOutputPixel, class TCoordRep >
ObliqueSectionImageFilter< TInputPixel, TOutputPixel, TCoordRep >
::ObliqueSectionImageFilter()
{
  m_OutputSize.Fill( 0 );
  m_OutputSpacing.Fill( 1.0 );
  m_OutputOrigin.Fill( 0.0 );
  m_Interpolator = DefaultInterpolatorType::New();
  m_DefaultPixelValue = NumericTraits< TOutputPixel >::Zero;
  m_Center.Fill( 0.0 );
  m_Direction.Fill( 0.0 );
  m_Direction[2] = 1.0;
  m_Up.Fill( 0.0 );
  m_Up[1] = 1.0;
}


/**
 *
 */
template < class TInputPixel, class TOutputPixel, class TCoordRep >
void
ObliqueSectionImageFilter< TInputPixel, TOutputPixel, TCoordRep >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "OutputSize: " << m_OutputSize << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "DefaultPixelValue: " << m_DefaultPixelValue << std::endl;
  os << indent << "Interpolator: " << m_Interpolator << std::endl;
  os << indent << "Transform: " << m_Transform << std::endl;
}


/**
 *
 */
template < class TInputPixel, class TOutputPixel, class TCoordRep >
void 
ObliqueSectionImageFilter< TInputPixel, TOutputPixel, TCoordRep >
::GenerateOutputInformation()
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
  typename OutputImageType::RegionType outputRegion( m_OutputSize );
  outputPtr->SetLargestPossibleRegion( outputRegion );

  // Set the spacing
  outputPtr->SetSpacing( m_OutputSpacing );
  
  // Set the origin
  outputPtr->SetOrigin( m_OutputOrigin );
}


/**
 *
 */
template < class TInputPixel, class TOutputPixel, class TCoordRep >
void
ObliqueSectionImageFilter< TInputPixel, TOutputPixel, TCoordRep >
::GenerateData()
{
  // TODO: Add progress accumulator
  
  // Get pointers to the input and output
  typename Superclass::OutputImagePointer     outputPtr = this->GetOutput();
  typename Superclass::InputImageConstPointer inputPtr  = this->GetInput();

  // Setup resample output parameters
  typename InputImageType::SizeType resampleSize;
  resampleSize[0] = m_OutputSize[0];
  resampleSize[1] = m_OutputSize[1];
  resampleSize[2] = 1;
  typename InputImageType::SpacingType resampleSpacing;
  resampleSpacing[0] = m_OutputSpacing[0];
  resampleSpacing[1] = m_OutputSpacing[1];
  resampleSpacing[2] = m_OutputSpacing[1];
  typename InputImageType::PointType resampleOrigin;
  resampleOrigin[0] = m_OutputOrigin[0];
  resampleOrigin[1] = m_OutputOrigin[1];
  resampleOrigin[2] = m_OutputOrigin[1];

  // Create transform
  m_Transform = TransformType::New();

  // Initialize transform
  typedef LookAtTransformInitializer< TransformType, InputImageType >
    TransformInitalizerType;
  typename TransformInitalizerType::Pointer init = TransformInitalizerType::New();
  init->SetTransform( m_Transform );
  init->SetImage( inputPtr );
  init->SetPlane( m_Center, m_Direction, m_Up, resampleSize );
  init->InitializeTransform();
  
  // Apply resample filter
  typedef ResampleImageFilter< InputImageType, InputImageType >
    ResampleFilterType;
  typename ResampleFilterType::Pointer filterResample = ResampleFilterType::New();
  filterResample->SetInput( inputPtr );
  filterResample->SetTransform( m_Transform );
  filterResample->SetInterpolator( m_Interpolator );
  filterResample->SetSize( resampleSize );
  filterResample->SetOutputSpacing( resampleSpacing );
  filterResample->SetOutputOrigin( resampleOrigin );
  filterResample->Update( );

  // Apply extract filter
  typedef ExtractImageFilter< InputImageType, OutputImageType >
    ExtractFilterType;
  typename ExtractFilterType::Pointer filterExtract = ExtractFilterType::New();
  InputSizeType extractSize( resampleSize );
  extractSize[2] = 0;
  typename ExtractFilterType::InputImageRegionType extractRegion( extractSize );
  filterExtract->SetExtractionRegion( extractRegion );
  filterExtract->SetInput( filterResample->GetOutput() );

  // Graft output
  filterExtract->GraftOutput( outputPtr );
  filterExtract->Update( );
  this->GraftOutput( filterExtract->GetOutput() );
}

} // end namespace itk

#endif
