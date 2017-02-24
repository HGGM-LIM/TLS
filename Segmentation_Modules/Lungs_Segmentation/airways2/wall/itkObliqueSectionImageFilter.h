/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkObliqueSectionImageFilter.h,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkObliqueSectionImageFilter_h
#define __itkObliqueSectionImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRigid3DTransform.h"

namespace itk
{

/** \class ObliqueSectionImageFilter
 * \brief Resamples a 3-D image producing a 2-D image along an arbitrary plane.
 *
 * ObliqueSectionImageFilter resamples a 3-D image to produce a 2-D oblique
 * section. An oblique section is an arbitrary plane through a 3-D image.
 * The plane is specified as a center point, a direction to "look", and
 * an "up" vector. The "up" vector specifies the direction to be considered
 * as the y-axis after the transformation. A warning is raised if the up
 * vector is the same as the direction vector.
 *
 * \author Dan Mueller, Queensland University of Technology, dan.muel[at]gmail.com
 *
 * \ingroup GeometricTransforms
 */
template < class TInputPixel, class TOutputPixel, class TCoordRep=double >
class ITK_EXPORT ObliqueSectionImageFilter:
    public ImageToImageFilter< Image<TInputPixel,3>, Image<TOutputPixel,2> >
{
public:
  /** Standard class typedefs. */
  typedef ObliqueSectionImageFilter Self;
  typedef ImageToImageFilter<
	  Image<TInputPixel,3>,
	  Image<TOutputPixel,2> > Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(ObliqueSectionImageFilter, ImageToImageFilter);

  /** Interpolator precision typedef. */
  typedef TCoordRep CoordRepType;

  /** Image type information. */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;

  /** Typedef to describe the type of pixel. */
  typedef TInputPixel InputPixelType;
  typedef TOutputPixel OutputPixelType;
  
  /** ImageDimension enumeration */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      InputImageType::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      OutputImageType::ImageDimension);

  /** Typedef to describe index, size, and spacing types. */
  typedef typename InputImageType::IndexType    InputIndexType;
  typedef typename InputImageType::SizeType     InputSizeType;
  typedef typename InputImageType::SpacingType  InputSpacingType;
  typedef typename InputImageType::PointType    InputPointType;
  typedef typename OutputImageType::IndexType   OutputIndexType;
  typedef typename OutputImageType::SizeType    OutputSizeType;
  typedef typename OutputImageType::SpacingType OutputSpacingType;
  typedef typename OutputImageType::PointType   OutputPointType;
  
  /** Typedef for vector type. */
  typedef Vector< CoordRepType, InputImageDimension >
    InputVectorType;

  /** Transform typedef. */
  typedef Rigid3DTransform< CoordRepType > TransformType;
  typedef typename TransformType::Pointer TransformPointer;

  /** Interpolator typedef. */
  typedef InterpolateImageFunction< InputImageType, CoordRepType >
    InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;
  typedef LinearInterpolateImageFunction< InputImageType, CoordRepType >
    DefaultInterpolatorType;

  /** Get a pointer to the resultant transform.
   *  The transform is only valid after calling GenerateData(). */
  itkGetConstObjectMacro( Transform, TransformType );

  /** Set the interpolator function.  The default is
   * itk::LinearInterpolateImageFunction<InputImageType, CoordRepType>. Some
   * other options are itk::NearestNeighborInterpolateImageFunction
   * (useful for binary masks and other images with a small number of
   * possible pixel values), and itk::BSplineInterpolateImageFunction
   * (which provides a higher order of interpolation).  */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );
  
  /** Set the pixel value when a transformed pixel is outside of the
   * image. The default default pixel value is 0. */
  itkSetMacro( DefaultPixelValue, OutputPixelType );

  /** Get the pixel value when a transformed pixel is outside of the image */
  itkGetMacro( DefaultPixelValue, OutputPixelType );

  /** Set the size of the output image. */
  itkSetMacro( OutputSize, OutputSizeType );

  /** Get the size of the output image. */
  itkGetConstReferenceMacro( OutputSize, OutputSizeType );
  
  /** Set the spacing of the output image. */
  itkSetMacro( OutputSpacing, OutputSpacingType );

  /** Get the spacing of the output image. */
  itkGetConstReferenceMacro( OutputSpacing, OutputSpacingType );
  
  /** Set the origin of the output image. */
  itkSetMacro( OutputOrigin, OutputPointType );

  /** Get the spacing of the output image. */
  itkGetConstReferenceMacro( OutputOrigin, OutputPointType );
  
  /** Helper method to set the output parameters based on this image */
  void SetOutputParametersFromImage ( typename OutputImageType::Pointer image )
    {
    this->SetOutputOrigin ( image->GetOrigin() );
    this->SetOutputSpacing ( image->GetSpacing() );
    this->SetOutputSize ( image->GetLargestPossibleRegion().GetSize() );
    }
    
  /** Set the plane using a center and direction to "look".
   *  The "up" vector defaults to the y-axis [0,1,0]. */
  virtual void SetPlane( InputPointType center,
                         InputVectorType direction)
  {
    m_Center = center;
    m_Direction = direction;
    this->Modified( );
  }

  /** Set the plane using a center, direction to "look", and
   *  "up" vector. */
  virtual void SetPlane( InputPointType center,
                         InputVectorType direction,
                         InputVectorType up )
  {
    m_Center = center;
    m_Direction = direction;
    m_Up = up;
    this->Modified( );
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputCovertibleToOutputCheck,
    (Concept::Convertible<InputPixelType, OutputPixelType>));
  /** End concept checking */
#endif

protected:
  ObliqueSectionImageFilter();
  ~ObliqueSectionImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** ObliqueSectionImageFilter produces an image which is a different
   * resolution than its input image.  As such, the implementation for
   * GenerateOutputInformation() must be provided in order to inform the
   * pipeline execution model.
   *
   * \sa ProcessObject::GenerateOutputInformaton()  */
  virtual void GenerateOutputInformation( );

  /** This method implements the actual reampling. */
  void GenerateData( );

private:
  ObliqueSectionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  TransformPointer    m_Transform;
  InterpolatorPointer m_Interpolator;
  OutputPixelType     m_DefaultPixelValue;
  OutputSizeType      m_OutputSize;
  OutputSpacingType   m_OutputSpacing;
  OutputPointType     m_OutputOrigin;
  InputPointType      m_Center;
  InputVectorType     m_Direction;
  InputVectorType     m_Up;
};

  
} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkObliqueSectionImageFilter.txx"
#endif
  
#endif
