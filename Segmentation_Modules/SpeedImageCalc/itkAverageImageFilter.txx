/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAverageImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/07/17 12:26:43 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAverageImageFilter_txx
#define _itkAverageImageFilter_txx

#include "itkAverageImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkNumericTraits.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
void
AverageImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

template< typename TInputImage, typename TOutputImage >
void
AverageImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread,
                        int itkNotUsed(threadId))
{
  typedef ImageRegionConstIterator< TInputImage > IteratorType;
  typedef ImageRegionIterator< TOutputImage> OutIteratorType;

  typename TOutputImage::Pointer output = this->GetOutput();
  
  // Record the number of input images.
  const unsigned int numberOfInputFiles = this->GetNumberOfInputs();

  //  create and initialize all input image iterators
  IteratorType *it = new IteratorType[numberOfInputFiles];
  for ( unsigned int i = 0; i < numberOfInputFiles; ++i)
    {
    it[i] = IteratorType( this->GetInput( i ), 
                          outputRegionForThread );
    }

  OutIteratorType out = OutIteratorType( output, outputRegionForThread );
  for( out.GoToBegin(); !out.IsAtEnd(); ++out )
    {
    typename NumericTraits<OutputPixelType>::RealType sum = it[0].Get();
    ++it[0];
	
    for( unsigned int i = 1; i < numberOfInputFiles; ++i)
      {
      sum += it[i].Get();
      ++(it[i]);
      }
    
    sum *= (1.0 / numberOfInputFiles);
    // instead of casting to output type, we should support rounding if the
    // output type is an integer type. Unfortunately, this does not seem to
    // be supported by the NumericTraits so far.
    out.Set( (OutputPixelType) sum );
    }
  
  delete[] it;  
}

} // end namespace itk

#endif
