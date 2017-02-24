#ifndef __itkVisualScoreImageFilter_h
#define __itkVisualScoreImageFilter_h

#include "itkImageToImageFilter.h"

#include <stdio.h>
#include <string>
#include <time.h>
#include <list>
#include <iostream>
#include <fstream>
#include <vector>

#include <sys/types.h>
#include <sys/timeb.h>

namespace itk
{

template<class TInputImage, class TOutputImage>
class ITK_EXPORT VisualScoreImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef VisualScoreImageFilter Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( VisualScoreImageFilter, ImageToImageFilter );

  /** Number of dimensions */
  itkStaticConstMacro(NDimensions, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(NOutputDimensions, unsigned int, TOutputImage::ImageDimension);

  /** typedef for images */
  
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;
  typedef typename InputImageType::RegionType RegionType;
  typedef typename RegionType::IndexType IndexType;
  typedef typename RegionType::SizeType SizeType;
  typedef typename TInputImage::PixelType PixelType;

  itkSetMacro(Verbose, int);
  itkGetMacro(Verbose, int);
  itkSetMacro(MinArea, unsigned long);
  itkGetMacro(MinArea, unsigned long);
  itkSetMacro(LowTh, int);
  itkGetMacro(LowTh, int);
  itkSetMacro(EnfTh, int);
  itkGetMacro(EnfTh, int);
  itkSetMacro(TissueTh, int);
  itkGetMacro(TissueTh, int);
  itkSetMacro(LowestScoreTh, float);
  itkGetMacro(LowestScoreTh, float);

  // Read only parameters
  itkGetMacro(VisualScore, float);
  itkGetMacro(LowVisualScore, float);
  itkGetMacro(MedVisualScore, float);
  itkGetMacro(UppVisualScore, float);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<itkGetStaticConstMacro(NDimensions),
                            itkGetStaticConstMacro(NOutputDimensions)>));
  itkConceptMacro(InputConvertibleToDoubleCheck,
    (Concept::Convertible<typename TInputImage::PixelType, double>));
  itkConceptMacro(DoubleConvertibleToOutputCheck,
    (Concept::Convertible<double, PixelType>));
  /** End concept checking */
#endif

protected:
  VisualScoreImageFilter();
  virtual ~VisualScoreImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method for evaluating the implicit function over the image. */
  void GenerateData();

private:
  VisualScoreImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  	float m_VisualScore;
  	float m_LowVisualScore;
  	float m_MedVisualScore;
  	float m_UppVisualScore;
  	int m_Verbose;
	unsigned long m_MinArea;
	int m_LowTh;
	int m_EnfTh;
	int m_TissueTh;
	float m_LowestScoreTh;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVisualScoreImageFilter.txx"
#endif

#endif
