
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>

#include "itkMask2ColorBorderImageFilter.h"

#include <stdio.h>
#include <string>
#include <fstream>

using namespace std;
using namespace itk;
int main(int argc, char * argv[])
{

  const int Dimension = 3;

  typedef signed short                                        PixelType;
  typedef Image<PixelType, Dimension>                         ShortImageType;
  typedef Image<bool, Dimension>                              BinImageType;
  typedef RGBPixel<unsigned char>                             RGBPixelType;
  typedef Image<RGBPixelType, Dimension>                      RGBImageType;
  typedef ImageFileReader<ShortImageType>                     ShortReaderType;
  typedef ImageFileWriter<ShortImageType>                     ShortWriterType;
  typedef ImageFileReader<BinImageType>                       BinReaderType;
  typedef ShortImageType::SizeType                            SizeType;
  typedef ShortImageType::RegionType                          RegionType;
  typedef ImageRegionConstIterator<BinImageType>              InBinIteratorType;
  typedef ImageRegionConstIterator<ShortImageType>            InShortIteratorType;
  typedef ImageRegionIterator<BinImageType>                   OutBinIteratorType;
  typedef ImageRegionIterator<ShortImageType>                 OutShortIteratorType;

  typedef Mask2ColorBorderImageFilter<ShortImageType,
      RGBImageType>                                         BorderFilterType;


 if (argc<3)
  {
    std::cout << "Usage: " << argv[0] << " input_image mask_image [region]" << std::endl;
    return EXIT_FAILURE;    
  }

 if ( (strcmp(argv[1],"--help") == 0) || ( (strcmp(argv[1],"-h") == 0) ) ){
    std::cout << "Usage: " << argv[0] << " input_image mask_image [region]" << std::endl;
    return EXIT_SUCCESS;
  }

  std::string voiString = argv[1];
  std::string hessString = argv[2];
  
  bool user_region = false;
  int  ur_x, ur_X, ur_y, ur_Y, ur_z, ur_Z;
  RegionType uRegion;
  RegionType::IndexType origin;
  RegionType::SizeType size;
  
  origin.Fill(0);
  size.Fill(0);
  
  ur_x = ur_X = ur_y = ur_Y = ur_z = ur_Z = 0;
    
  if (argc > 3) {
     if (argc < 9) {
        std::cout << "Usage: " << argv[0] << " input_image mask_image [region]" << std::endl;
        return EXIT_FAILURE;    
     } 
     
     ur_x = atoi(argv[3]); ur_X = atoi(argv[4]); 
     ur_y = atoi(argv[5]); ur_Y = atoi(argv[6]); 
     ur_z = atoi(argv[7]); ur_Z = atoi(argv[8]);
     
     user_region = true;
     
     size[0] = ur_X - ur_x; size[1] = ur_Y - ur_y; size[2] = ur_Z - ur_z;
     uRegion.SetSize(size);
     uRegion.SetIndex(origin);
     
  }
  

  


  //Reader for the cropped input image
  ShortReaderType::Pointer inReader = ShortReaderType::New();
  inReader->SetFileName(voiString.c_str());

  //Reader for the hessian of the cropped input image
  ShortReaderType::Pointer hessReader = ShortReaderType::New();
  hessReader->SetFileName(hessString.c_str());
  
  try {
    inReader->Update();
    hessReader->Update();
  }
  catch (itk::ExceptionObject & excp) {
    cerr << "  Problems reading an image"   << endl;
    cerr << excp << endl << endl;
    return EXIT_FAILURE;
  }

  RegionType inRegion;
  
  inRegion = inReader->GetOutput()->GetLargestPossibleRegion();

  BinImageType::Pointer outImage = BinImageType::New();
  outImage->SetRegions(inRegion);
  outImage->Allocate();
  outImage->FillBuffer(0);

  if (user_region)
    inRegion = uRegion;

  InShortIteratorType hess = InShortIteratorType(hessReader->GetOutput(), inRegion);
  InShortIteratorType voi = InShortIteratorType(inReader->GetOutput(), inRegion);
  OutBinIteratorType out = OutBinIteratorType(outImage, inRegion);



  try {
    float sum = 0.0;
    float dev = 0.0; 
    int N = 0;
    for (hess.GoToBegin(), voi.GoToBegin(), out.GoToBegin();
        !hess.IsAtEnd(), !voi.IsAtEnd(),!out.IsAtEnd();++hess, ++voi, ++out) {
      if (hess.Get() == 1) {
        sum += voi.Get();

        N++;
      }
      out.Set( hess.Get() );
    }

    sum /= N;



    for (hess.GoToBegin(), voi.GoToBegin();!hess.IsAtEnd(), !voi.IsAtEnd(); ++hess, ++voi) {

      if (hess.Get() == 1) {
        dev += pow(voi.Get() - sum,2);

      }
      
    }

    dev = sqrt( dev / N );

    cout << "Mean: " << sum << ", " << "std: " << dev << endl;



  }
  catch (itk::ExceptionObject & excp) {
    cerr << "  Problem getting the mean of the region" << endl;
    cerr << excp << endl << endl;
    return EXIT_FAILURE;
  }
  BorderFilterType::Pointer borderFilter = BorderFilterType::New();
  borderFilter->SetInput(inReader->GetOutput());
  borderFilter->SetMaskImage(outImage);
  borderFilter->SetOutName("aorta.tif");
  borderFilter->SetRedMask(255);
  try {
    borderFilter->Update();
  }
  catch (itk::ExceptionObject & excp) {
    cerr << "  Problem drawing mask image" << endl;
    cerr << excp << endl << endl;
    return EXIT_FAILURE;
  }
  /*ShortWriterType::Pointer writer = ShortWriterType::New();
  writer->SetInput(outImage);
  writer->SetFileName("aorta.tiff");
  writer->Update();*/
  return EXIT_SUCCESS;
}

