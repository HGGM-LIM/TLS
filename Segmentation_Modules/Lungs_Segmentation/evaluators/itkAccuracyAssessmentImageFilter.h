/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAccuracyAssessmentImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/28 19:59:05 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAccuracyAssessmentImageFilter_h
#define __itkAccuracyAssessmentImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"

namespace itk {

/** \class AccuracyAssessmentImageFilter 
 * \brief Computes...
 * 
 * This filter requires the largest possible region of the first image
 * and the same corresponding region in the second image. 
 * It behaves as filter with
 * two input and one output. Thus it can be inserted in a pipeline with 
 * other filters. The filter passes the first input through unmodified.
 *
 * This filter is templated over the two input image type. It assume
 * both image have the same number of dimensions.
 *
 * \sa AccuracyAssessmentImageFilter
 *
 * \ingroup MultiThreaded
 */
template<class TInputImage1, class TInputImage2>
class ITK_EXPORT AccuracyAssessmentImageFilter : 
    public ImageToImageFilter<TInputImage1, TInputImage1>
{
public:
  /** Standard Self typedef */
  typedef AccuracyAssessmentImageFilter Self;
  typedef ImageToImageFilter<TInputImage1,TInputImage1>  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(AccuracyAssessmentImageFilter, ImageToImageFilter);
  
  /** Image related typedefs. */
  typedef TInputImage1 InputImage1Type;
  typedef TInputImage2 InputImage2Type;
  typedef typename TInputImage1::Pointer InputImage1Pointer;
  typedef typename TInputImage2::Pointer InputImage2Pointer;
  typedef typename TInputImage1::ConstPointer InputImage1ConstPointer;
  typedef typename TInputImage2::ConstPointer InputImage2ConstPointer;

  typedef typename TInputImage1::RegionType RegionType ;
  typedef typename TInputImage1::SizeType SizeType ;
  typedef typename TInputImage1::IndexType IndexType ;

  typedef typename TInputImage1::PixelType InputImage1PixelType;
  typedef typename TInputImage2::PixelType InputImage2PixelType;
  
  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage1::ImageDimension);

  /** Type to use form computations. */
  //typedef typename NumericTraits<InputImage1PixelType>::RealType RealType;
  typedef float RealType;

  /** Set the first input. */
  void SetReferenceSegmentation( const InputImage1Type * image )
  { this->SetInput( image ); }

  /** Set the second input. */
  void SetTestSegmentation( const InputImage2Type * image );

  /** Get the first input. */
  const InputImage1Type * GetReferenceSegmentation(void)
  { return this->GetInput(); }
  
  /** Get the second input. */
  const InputImage2Type * GetTestSegmentation(void);
  
  /** Return the computed values. */
  itkGetMacro(TruePositiveVF,RealType);
  itkGetMacro(FalsePositiveVF,RealType);
  itkGetMacro(TrueNegativeVF,RealType);
  itkGetMacro(FalseNegativeVF,RealType);
  
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(Input1HasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImage1PixelType>));
  /** End concept checking */
#endif

protected:
  AccuracyAssessmentImageFilter();
  ~AccuracyAssessmentImageFilter(){};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** GenerateData. */
  void  GenerateData();

  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

  // Override since the filter produces all of its output
  void EnlargeOutputRequestedRegion(DataObject *data);

private:
  AccuracyAssessmentImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  RealType                            m_TruePositiveVF;
  RealType                            m_TrueNegativeVF;
  RealType                            m_FalsePositiveVF;
  RealType                            m_FalseNegativeVF;
  

} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAccuracyAssessmentImageFilter.txx"
#endif

#endif
