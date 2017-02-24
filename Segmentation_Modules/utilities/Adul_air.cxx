
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <itkCastImageFilter.h>

#include "itkMask2ColorBorderImageFilter.h"
#include "itkBackgroundImageFilter.h"
#include "itkMeanStdImageFilter.h"

#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>

#include <limits.h>


using namespace itk;
/** Pixels types */
typedef signed short                                        PixelType;
typedef unsigned char                                       RGBType;
typedef itk::RGBPixel<RGBType>                              RGBPixelType;

/** Images types */
typedef Image<PixelType, 3>                                 ShortImageType;
typedef Image<unsigned short, 3>                            USImageType;
typedef Image<float, 3>                            			FImageType;
typedef Image<bool, 3>                                      BinImageType;
typedef Image<signed long, 3>                               SLImageType;
typedef itk::Image<RGBPixelType, 3>                         RGBImageType;

/** Pointers to images types */
typedef ShortImageType::Pointer                             ShortImagePointer;
typedef BinImageType::Pointer                               BinImagePointer;

/** Filters types */
typedef ImageFileReader<ShortImageType>                     ShortReaderType;
typedef ImageFileReader<USImageType>                        UShortReaderType;
typedef ImageFileWriter<ShortImageType>                     ShortWriterType;
typedef ImageFileWriter<BinImageType>                       BinWriterType;
typedef ImageFileReader<BinImageType>                       BinReaderType;
typedef ImageFileWriter<USImageType>                        USWriterType;
typedef ImageFileWriter<ShortImageType>                     SSWriterType;
typedef CastImageFilter<BinImageType, USImageType>          CasterType;
typedef CastImageFilter<SLImageType, USImageType>           CasterSL2USType;
typedef CastImageFilter<SLImageType, ShortImageType>        CasterSL2SSType;
typedef Mask2ColorBorderImageFilter<ShortImageType,
    RGBImageType>                                           BorderFilterType;
typedef BackgroundImageFilter<BinImageType, BinImageType>   BackgroundFilterType;
typedef MeanStdImageFilter<ShortImageType, SLImageType>     MeanStdFilter;

/** Iterators types */
typedef ImageRegionConstIterator<BinImageType>              InBinIteratorType;
typedef ImageRegionConstIterator<ShortImageType>            InShortIteratorType;
typedef ImageRegionConstIterator<SLImageType>               InSLIteratorType;
typedef ImageRegionConstIterator<USImageType>               InUSIteratorType;
typedef ImageRegionIterator<BinImageType>                   OutBinIteratorType;
typedef ImageRegionIterator<ShortImageType>                 OutShortIteratorType;

/** Other types */
typedef ShortImageType::SizeType                            SizeType;
typedef ShortImageType::RegionType                          RegionType;


int main(int argc, char * argv[])
{
  if (argc<3)
  {
    std::cout << "Usage: " << argv[0] << " input_file binary_image [region]" << std::endl;
    return EXIT_FAILURE;    
  }


 if ( (strcmp(argv[1],"--help") == 0) || ( (strcmp(argv[1],"-h") == 0) ) ){
    std::cout << "Usage: " << argv[0] << " input_file binary_image [region]" << std::endl;
    return EXIT_SUCCESS;
  }


  std::string inDirName = "./";
  std::string inFileName = inDirName + argv[1];
  std::string huFileName = inDirName + argv[2];

  unsigned short squareRegionSize = 10;
  if (argc==4)
    squareRegionSize = atoi(argv[3]);

  USImageType::SizeType regionSize = {squareRegionSize, squareRegionSize, 1};

  //Reader for the cropped input image
  ShortReaderType::Pointer inputReader = ShortReaderType::New();
  inputReader->SetFileName(inFileName.c_str());

  BinReaderType::Pointer huThresholdReader = BinReaderType::New();
  huThresholdReader->SetFileName(huFileName.c_str());

  try {
    inputReader->Update();
    huThresholdReader->Update();
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << "  Problems reading the input image"   << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  BackgroundFilterType::Pointer backgroundFilter = BackgroundFilterType::New();
  backgroundFilter->SetInput(huThresholdReader->GetOutput());
  backgroundFilter->SetGrayScaleImage(inputReader->GetOutput());
  try {
    backgroundFilter->Update();
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << "Error while creating the mask Image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  MeanStdFilter::Pointer meanStdFilter = MeanStdFilter::New();
  meanStdFilter->SetRegionSize(regionSize);
  meanStdFilter->SetInput(inputReader->GetOutput());
  meanStdFilter->SetMaskImage(backgroundFilter->GetOutput());

  try {
    meanStdFilter->Update();
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << "Error while creating the mean and std image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  CasterSL2SSType::Pointer caster = CasterSL2SSType::New();
  caster->SetInput(meanStdFilter->GetOutput());

  SSWriterType::Pointer meanWriter = SSWriterType::New();
  meanWriter->SetInput(caster->GetOutput());
  meanWriter->SetFileName("MeanImage.tif");

  USWriterType::Pointer stdWriter = USWriterType::New();
  stdWriter->SetInput(meanStdFilter->GetStdImage());
  stdWriter->SetFileName("StdImage.tif");

  try {
    meanWriter->Update();
    stdWriter->Update();
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << "  Problems writting the std and mean image"   << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  InSLIteratorType meanIt = InSLIteratorType(meanStdFilter->GetOutput(),
                            meanStdFilter->GetOutput()->GetLargestPossibleRegion());
  InUSIteratorType stdIt  = InUSIteratorType(meanStdFilter->GetStdImage(),
                            meanStdFilter->GetOutput()->GetLargestPossibleRegion());

  PixelType meanValue = 0;
  typedef ShortImageType::IndexType IndexType;
  IndexType minIndex;
  minIndex.Fill(0);
  MeanStdFilter::USImageType::PixelType stdMin = USHRT_MAX;
  for (meanIt.GoToBegin(), stdIt.GoToBegin();!meanIt.IsAtEnd(), !stdIt.IsAtEnd();++meanIt, ++stdIt)
    if (stdIt.Get()&&stdIt.Get()<stdMin) {
      stdMin = stdIt.Get();
      meanValue = meanIt.Get();
      minIndex = stdIt.GetIndex();
    }

  // Writes down the air mean in a format that can be easily parsed by a calling script
  std::cout << "Air mean: " << meanValue << ", " << "Standard dev: " << stdMin << " at index: " << minIndex << std::endl;



  return EXIT_SUCCESS;
}

