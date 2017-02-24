/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSpeedCalculatorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2005/07/15 11:48:09 $
  Version:   $Revision: 1.2 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSpeedCalculatorImageFilter_h
#define __itkSpeedCalculatorImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{
/** \class SpeedCalculatorImageFilter
 *
 * \brief This filter computes a speed image for the fast marching in chest images for airway segmentation.
 * It uses gradient and intensity information to compute a probability
 * per voxel. It expects images in Hounsfield Units.
 *
 * \par INPUT 
 *  Original image
 
 * \par OUTPUT 
 *  Speed image
 
 * \author Xabier Artaechevarria, CIMA
 *
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_EXPORT SpeedCalculatorImageFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef SpeedCalculatorImageFilter Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(SpeedCalculatorImageFilter, ImageToImageFilter);
  
  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TInputImage::PixelType InputPixelType;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  
  /** Image typedef support */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;
  typedef typename InputImageType::ConstPointer InputImagePointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  
  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  
  /**Set and get parameters*/
  // empty for the time being
  itkSetMacro( Scale, double);
  itkGetMacro( Scale, double);

protected:
  SpeedCalculatorImageFilter();
  virtual ~SpeedCalculatorImageFilter() {}  

  void GenerateData();
  
  void PrintSelf(std::ostream&, Indent) const;

private:
  SpeedCalculatorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  float m_Scale;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpeedCalculatorImageFilter.txx"
#endif

#endif
