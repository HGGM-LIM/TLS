#ifndef __itkEmphysemaCalculationsImageFilter_h
#define __itkEmphysemaCalculationsImageFilter_h

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkScalarImageToHistogramGenerator.h"
//#include "itkListSampleToHistogramGenerator.h"

#include "itkSubtractImageFilter.h"
#include "itkShapeLabelObject.h"
#include "itkExtractImageFilter.h"

#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
//#include "itkLabelMapMaskImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkLabelExtractorFilter.h"
#include "itkMedianImageFilter.h"

#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"

#include <stdio.h>
#include <string>
#include <time.h>
#include <list>
#include <iostream>
#include <fstream>
#include <vector>

//#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

namespace itk
{

template<class TInputImage, class TOutputImage>
class ITK_EXPORT EmphysemaCalculationsImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef EmphysemaCalculationsImageFilter Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( EmphysemaCalculationsImageFilter, ImageToImageFilter );

  /** Number of dimensions */
  itkStaticConstMacro(NDimensions, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(NOutputDimensions, unsigned int, TOutputImage::ImageDimension);

  /** typedef for images */
  
  typedef TInputImage                                                   InputImageType;
  typedef TOutputImage                                                  OutputImageType;
  typedef signed short 													PixelType;
  typedef unsigned short 												BinaryPixelType;
  typedef unsigned long 												LabelType;
  typedef itk::Image< PixelType, 3 >                       		 		ImageType;
  typedef itk::Image< BinaryPixelType, 3 >                       		BinaryImageType;
  typedef itk::Image<PixelType, 2>                                	  	SliceType;
  typedef typename OutputImageType::Pointer                             OutputImagePointer;
  typedef typename InputImageType::Pointer                              InputImagePointer;
  typedef typename InputImageType::ConstPointer                         InputImageConstPointer;
  typedef typename InputImageType::RegionType														RegionType;
  typedef typename InputImageType::IndexType 													  IndexType;
  typedef typename InputImageType::SizeType 													  SizeType;
  
	typedef itk::ImageFileWriter<ImageType>                         	  	WriterType;
	typedef itk::ShapeLabelObject<LabelType, 3>            		 				  	LabelObjectType;
  typedef itk::LabelMap<LabelObjectType>                         		  	LabelCollectionType; 
	typedef itk::BinaryImageToLabelMapFilter<BinaryImageType, LabelCollectionType> ConverterType;
	typedef itk::ShapeLabelMapFilter<LabelCollectionType> 								ShapeFilterType;
	
	typedef itk::ImageRegionIterator<ImageType>                    				IteratorType;
	typedef itk::ImageRegionIterator<BinaryImageType>                    		BinaryIteratorType;
	typedef itk::ImageRegionConstIterator<ImageType>                    	ConstIteratorType;
	

	typedef itk::BinaryCrossStructuringElement<BinaryPixelType, 3> StructuringElementType;
    typedef itk::BinaryErodeImageFilter <BinaryImageType, BinaryImageType, StructuringElementType> BinaryErodeImageFilterType;


  /** Image size typedef */
  /*typedef Size<itkGetStaticConstMacro(NDimensions)> SizeType;*/

  /** Image index typedef */
  /*typedef typename TOutputImage::IndexType IndexType;*/

  /** Image pixel value typedef */
  /*typedef typename TOutputImage::PixelType PixelType;*/

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;
  

  itkSetMacro(LowEmphValue, short);
  itkGetMacro(LowEmphValue, short);
  
  itkSetMacro(Backup, bool);
  itkGetMacro(Backup, bool);  
  itkSetMacro(tEmphysema, double);
  itkGetMacro(tEmphysema, double);
  itkSetMacro(tConverter, double);
  itkGetMacro(tConverter, double);
  itkSetMacro(tShape, double);
  itkGetMacro(tShape, double);
  itkSetMacro(tLavFile, double);
  itkGetMacro(tLavFile, double);
  itkSetMacro(tHistogram, double);
  itkGetMacro(tHistogram, double);
	
	itkSetMacro(NumberOfBins, unsigned int);
  itkGetMacro(NumberOfBins, unsigned int);
  itkSetMacro(MarginScale, double);
  itkGetMacro(MarginScale, double);
  
  itkSetMacro(VolLung, unsigned long);
  itkGetMacro(VolLung, unsigned long);
  itkSetMacro(VolEmph, unsigned long);
  itkGetMacro(VolEmph, unsigned long);
  
  itkSetMacro(Smooth, int);
  itkGetMacro(Smooth, int);

  itkSetMacro(SmoothRadius, int);
  itkGetMacro(SmoothRadius, int);

  itkSetMacro(ClusterFile, std::string);
  itkGetMacro(ClusterFile, std::string);

  //Not aleady implemented but in mind
  itkSetMacro(Size, SizeType);
  itkGetMacro(Size, SizeType);
  itkSetMacro(Index, IndexType);
  itkGetMacro(Index, IndexType);
  itkSetMacro(AutomaticallyResize, bool);
  itkGetMacro(AutomaticallyResize, bool);
  
  


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
  EmphysemaCalculationsImageFilter();
  virtual ~EmphysemaCalculationsImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method for evaluating the implicit function over the image. */
  void GenerateData();

private:
  EmphysemaCalculationsImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  int    				m_LowEmphValue;
  bool      		m_Backup;
  bool      		m_AutomaticallyResize;
  unsigned long m_VolLung;
  unsigned long m_VolEmph;
  long 					m_Mean;
  
  double    		m_tEmphysema;
  double    		m_tConverter;
  double    		m_tShape;
  double    		m_tLavFile;
  double    		m_tHistogram;
  
  unsigned int  m_NumberOfBins;
  double    	  m_MarginScale;
  int           m_Smooth;
  int           m_SmoothRadius;
  std::string	m_ClusterFile;
  
  SizeType  m_Size;
  IndexType m_Index;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEmphysemaCalculationsImageFilter.txx"
#endif

#endif
