/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelVotingImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/04/07 11:42:31 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFindTracheaImageFilter_h
#define __itkFindTracheaImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

//#include "itkArray.h"
#include <vnl/vnl_matrix.h>

// structure "Label" defition
#include "Label.h"

namespace itk
{
/** \class FindTracheaImageFilter
 *
 * \brief This filter looks for the trachea in a single 2D slice.
 *
 * \author Xabier Artaechevarria, CIMA, University of Navarra,
 * 
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_EXPORT FindTracheaImageFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef FindTracheaImageFilter Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(FindTracheaImageFilter, ImageToImageFilter);
  
  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TInputImage::PixelType InputPixelType;
  
  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  
  /** Image typedef support */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  
  /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
 
  /** Get the trachea center. */
  typename TOutputImage::IndexType GetTracheaCenter(void)
  { return this->m_tracheaCenter(); }

  // check is label has already been checked
  bool IsNotComputed( InputPixelType );
  
  // compute label statistics and add data to vector
  void ComputeAndAddLabelStatistics( InputPixelType, InputImageConstPointer );
  
  /**Set n1 and n2*/
  //itkSetMacro( N1, unsigned int );
  //itkGetMacro( N1, unsigned int );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<InputPixelType, OutputPixelType>));
  itkConceptMacro(IntConvertibleToInputCheck,
    (Concept::Convertible<int, InputPixelType>));
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<InputImageDimension, ImageDimension>));
  itkConceptMacro(InputConvertibleToUnsignedIntCheck,
    (Concept::Convertible<InputPixelType, unsigned int>));
  itkConceptMacro(IntConvertibleToOutputPixelType,
    (Concept::Convertible<int, OutputPixelType>));
  itkConceptMacro(InputPlusIntCheck,
    (Concept::AdditiveOperators<InputPixelType, int>));
  itkConceptMacro(InputIncrementDecrementOperatorsCheck,
    (Concept::IncrementDecrementOperators<InputPixelType>));
  itkConceptMacro(OutputOStreamWritableCheck,
    (Concept::OStreamWritable<OutputPixelType>));
  /** End concept checking */
#endif

protected:
  FindTracheaImageFilter();
  virtual ~FindTracheaImageFilter() {}  

  virtual void GenerateInputRequestedRegion ();
  void PrintSelf(std::ostream&, Indent) const;
  void GenerateData ();
  
private:
  FindTracheaImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  // vector to save all candidate regions to be trachea
  std::vector<Label> m_labelVector;
  
  // label that trachea has in label image
  int m_tracheaLabel;
  
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFindTracheaImageFilter.txx"
#endif

#endif
