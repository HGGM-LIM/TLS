/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAirwaySectioningImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/29 14:53:40 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAirwaySectioningImageFilter_h
#define __itkAirwaySectioningImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include <string>
#include <sstream>
#include "itkLinearInterpolateImageFunction.h"
#include "itkRigid3DTransform.h"
#include "itkObliqueSectionImageFilter.h"
#include "itkLookAtTransformInitializer.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkResampleImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
//#include "itkExtractLabelImageFilter.h"
#include <itkConnectedComponentImageFilter.h>
#include <itkMinimumMaximumImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>

using namespace std;
namespace itk
{
/** \class AirwaySectioningImageFilter
 * \brief Extracts section...
 *
 * Extracts..
 *
 * 
 * 
 */
 
template <class TImageType>
class ITK_EXPORT AirwaySectioningImageFilter :
    public ImageToImageFilter< TImageType,  TImageType >
{
public:

  int GetArea();
  string GetHistogram();
   
  /** Extract dimension from input and output image. */
  //itkStaticConstMacro(InputImageDimension, unsigned int,
  //                    TImageType::ImageDimension);
  //itkStaticConstMacro(OutputImageDimension, unsigned int,
  //                    TImageType::ImageDimension);

  itkStaticConstMacro(InputImageDimension, signed short,
                      TImageType::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, signed short,
                      TImageType::ImageDimension);
  /** Convenient typedefs for simplifying declarations. */
  typedef Image<signed short, 3> InputImageType;
  typedef Image<signed short, 3> OutputImageType;
  //typedef Image<typename TImageType::PixelType, 3> InputImageType;
  //typedef Image<typename TImageType::PixelType, 3> OutputImageType;
  
  /** Standard class typedefs. */
  typedef AirwaySectioningImageFilter Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AirwaySectioningImageFilter, ImageToImageFilter);
  
  /** Image typedef support. */
  //typedef typename InputImageType::PixelType InputPixelType;
  //typedef typename OutputImageType::PixelType OutputPixelType;
  typedef  signed short InputPixelType;
  typedef  signed short OutputPixelType;

  typedef typename NumericTraits<InputPixelType>::RealType InputRealType;
  
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;
  typedef typename OutputImageType::SizeType    OutputSizeType;
  typedef typename OutputImageType::SpacingType    OutputSpacingType; 
  typedef typename OutputImageType::PointType   OutputPointType;
 
  // Typedefs for ObliqueSection and TransformInitializer.
  typedef itk::Rigid3DTransform< double > TransformType;
 // typedef itk::Rigid3DTransform< signed short > TransformType;
  typedef itk::LookAtTransformInitializer< TransformType, InputImageType > 	TransformInitalizerType;
  typedef itk::ResampleImageFilter< InputImageType, OutputImageType > 		ResampleFilterType;
  typedef itk::ObliqueSectionImageFilter< InputPixelType, OutputPixelType > 	ObliqueFilterType;


  // Typedefs for Extracting Label Image.

  typedef itk::ConnectedComponentImageFilter< InputImageType, OutputImageType, InputImageType >ConnectedFilter;   //3er param. TMaskImage = TInputImage
  typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > FilterType;
  typedef itk::RelabelComponentImageFilter< InputImageType, OutputImageType > ReLabelType;


  /** Typedefs for extracting the region of interest.*/
  typedef itk::ExtractImageFilter< TImageType, TImageType > ExtractFilterType;
  typedef itk::ImageRegionConstIterator< TImageType > ConstIteratorType;
  typedef itk::ImageRegionIterator< TImageType > IteratorType;

  /**Typedefs for obtaining the image histogram */
/*
  typedef itk::Statistics::ScalarImageToHistogramGenerator< TImageType > HistogramGeneratorType;
  typedef typename HistogramGeneratorType::HistogramType  HistogramType;
  typedef itk::MinimumMaximumImageFilter< TImageType > MinMaxType;
*/
  typedef itk::ImageFileWriter< TImageType > WriterType_3D;  
	
  /** Typedef for vector type. */
  typedef typename ObliqueFilterType::CoordRepType CoordType;
  typedef Vector< CoordType, InputImageDimension > InputVectorType;
  
  /** Typedef for gaussian filter. */
  typedef itk::DiscreteGaussianImageFilter< TImageType, TImageType > GaussianFilterType;
  
  /** Set the plane using a center, direction to "look", and  "up" vector. */
  virtual void SetPlane( InputVectorType center, InputVectorType direction, InputVectorType up, int SliceType)
  {
    m_Center = center;
    m_Direction = direction;
    m_Up = up;
    m_SliceType = SliceType;
    this->Modified();
  }
  virtual void SetRadiusToExtract( int Radius)
  {
    m_RadiusToExtract = Radius;
    this->Modified();
  }
  /** AirwaySectioningImageFilter needs a larger input requested region than
   * the output requested region.  As such, AirwaySectioningImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
 // virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<InputPixelType>));
  /** End concept checking */
#endif

protected:
  AirwaySectioningImageFilter();
  virtual ~AirwaySectioningImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

private:
  AirwaySectioningImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  virtual void GenerateOutputInformation();

  OutputPixelType     m_DefaultPixelValue;
  OutputSizeType      m_OutputSize;
  OutputSpacingType   m_OutputSpacing;
  OutputPointType     m_OutputOrigin;
  InputVectorType     m_Center;
  InputVectorType     m_Direction;
  InputVectorType     m_Up;
  int		      m_Area;
  int		      m_SliceType; //0 binary, 1 grayScale
  int		      m_RadiusToExtract;
//  stringstream	      m_Histogram;

};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAirwaySectioningImageFilter.txx"
#endif

#endif
