
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
typedef ImageFileWriter<USImageType>                        UShortWriterType;
typedef ImageFileWriter<BinImageType>                       BinWriterType;
typedef ImageFileReader<BinImageType>                       BinReaderType;
typedef ImageFileWriter<USImageType>                        USWriterType;
typedef ImageFileWriter<ShortImageType>                     SSWriterType;
typedef CastImageFilter<BinImageType, USImageType>          CasterType;
typedef CastImageFilter<SLImageType, USImageType>           CasterSL2USType;
typedef CastImageFilter<SLImageType, ShortImageType>        CasterSL2SSType;
typedef Mask2ColorBorderImageFilter<ShortImageType,
    RGBImageType>                                           BorderFilterType;
typedef BackgroundImageFilter<BinImageType, USImageType>   BackgroundFilterType;


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


  UShortWriterType::Pointer writer = UShortWriterType::New();
  writer->SetFileName("air_mask.tiff");
  writer->SetInput(backgroundFilter->GetOutput());

  try {
    writer->Update();
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << "Error while creating the mask Image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

