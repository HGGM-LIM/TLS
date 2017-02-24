/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHessian3DToAirwaynessMeasureImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-16 16:45:10 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHessian3DToAirwaynessMeasureImageFilter_txx
#define __itkHessian3DToAirwaynessMeasureImageFilter_txx

#include "itkHessian3DToAirwaynessMeasureImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Constructor
 */
template < typename TPixel >
Hessian3DToAirwaynessMeasureImageFilter< TPixel >
::Hessian3DToAirwaynessMeasureImageFilter()
{
  m_Alpha1 = 2;
  m_Alpha2 = 0.5;

  // Hessian( Image ) = Jacobian( Gradient ( Image ) )  is symmetric
  m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
  m_SymmetricEigenValueFilter->SetDimension( ImageDimension );
  m_SymmetricEigenValueFilter->OrderEigenValuesBy( 
      EigenAnalysisFilterType::FunctorType::OrderByValue );
  
}


template < typename TPixel >
void 
Hessian3DToAirwaynessMeasureImageFilter< TPixel >
::GenerateData()
{
  itkDebugMacro(<< "Hessian3DToAirwaynessMeasureImageFilter generating data ");

  m_SymmetricEigenValueFilter->SetInput( this->GetInput() );
  
  typename OutputImageType::Pointer output = this->GetOutput();

  typedef typename EigenAnalysisFilterType::OutputImageType
    EigenValueOutputImageType;

  m_SymmetricEigenValueFilter->Update();
  
  const typename EigenValueOutputImageType::ConstPointer eigenImage = 
                    m_SymmetricEigenValueFilter->GetOutput();
  
  // walk the region of eigen values and get the airwayness measure
  EigenValueArrayType eigenValue;
  ImageRegionConstIterator<EigenValueOutputImageType> it;
  it = ImageRegionConstIterator<EigenValueOutputImageType>(
      eigenImage, eigenImage->GetRequestedRegion());
  ImageRegionIterator<OutputImageType> oit;
  this->AllocateOutputs();
  oit = ImageRegionIterator<OutputImageType>(output,
                                             output->GetRequestedRegion());
  oit.GoToBegin();
  it.GoToBegin();
  while (!it.IsAtEnd())
    {
    // Get the eigen value
    eigenValue = it.Get();
    
    // normalizeValue >= 0 for dark line structures
    double normalizeValue = vnl_math_min( eigenValue[2], eigenValue[1]);
    //std::cout << "EigenValues:" << eigenValue[2] << ", " << eigenValue[1] << ", "<< eigenValue[0] << ", " << std::endl;
    // Similarity measure to a line structure
    if( normalizeValue > 0 )
      {
      double lineMeasure;
      if( eigenValue[0] <= 0 )
        {
        lineMeasure = 
          vcl_exp(-0.5 * vnl_math_sqr( eigenValue[0] / (m_Alpha1 * normalizeValue)));
        }
      else
        {
        lineMeasure = 
          vcl_exp(-0.5 * vnl_math_sqr( eigenValue[0] / (m_Alpha2 * normalizeValue)));
        }
      
      lineMeasure *= normalizeValue;
      //std::cout << lineMeasure << std::endl;
      oit.Set( static_cast< OutputPixelType >(lineMeasure) );
      }
    else
      {
      oit.Set( NumericTraits< OutputPixelType >::Zero );
      }

    ++it;
    ++oit;
    }
    
}

template < typename TPixel >
void
Hessian3DToAirwaynessMeasureImageFilter< TPixel >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "Alpha1: " << m_Alpha1 << std::endl;
  os << indent << "Alpha2: " << m_Alpha2 << std::endl;
}


} // end namespace itk
  
#endif
