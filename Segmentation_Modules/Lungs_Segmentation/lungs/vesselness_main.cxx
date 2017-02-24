
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


int main( int argc, char * argv[] )
{

  float sigma = 0.009;

  if (argc < 3){
    std::cout << "Usage: " << argv[0] << " InputFileName OutputFileName [sigma]" << std::endl;
    return EXIT_SUCCESS;
  }
  
  if ( (strcmp(argv[1],"--help") == 0) || ( (strcmp(argv[1],"-h") == 0) ) ){
    std::cout << "Usage: " << argv[0] << " InputFileName OutputFileName [sigma]" << std::endl;
    return EXIT_SUCCESS;
  }

  typedef signed short PixelType;
  typedef unsigned short USPixelType;
  const unsigned int Dimension=3;

  typedef itk::Image< PixelType, Dimension >                       	ImageType;
  typedef itk::Image< USPixelType, Dimension >                      USImageType;

  typedef itk::ImageFileReader< ImageType >                        	ReaderType;
  ReaderType::Pointer reader =                                     	ReaderType::New();
  typedef itk::ImageFileWriter< ImageType >                        	WriterType;
  typedef itk::ImageFileWriter< USImageType >                       USWriterType;
  
  typedef itk::ImageRegionConstIterator< ImageType > SourceIteratorType;
  typedef itk::ImageRegionIterator< USImageType> DestIteratorType;


  char * infilename  = argv[1];
  char * outfilename  = argv[2];

  if (argc == 4)
	sigma = atof(argv[3]);

  reader->SetFileName( infilename );

  try {
	reader->Update();
  } catch( itk::ExceptionObject & err ) {
    std::cerr << "ExceptionObject caught: can't read serie" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::HessianRecursiveGaussianImageFilter< ImageType > HessianType;

  HessianType::Pointer hessian = HessianType::New();

  hessian->SetSigma( sigma );
  hessian->SetInput( reader->GetOutput() );

  typedef itk::Hessian3DToVesselnessMeasureImageFilter< PixelType > VesselnessType;

  VesselnessType::Pointer vesselness = VesselnessType::New();

  vesselness->SetInput( hessian->GetOutput() );
  
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( vesselness->GetOutput() );
  writer->SetFileName( outfilename );
  writer->Update();
  
  USImageType::Pointer hessian_mask = USImageType::New();
  hessian_mask->SetRegions(vesselness->GetOutput()->GetLargestPossibleRegion());
  hessian_mask->Allocate();
  
  SourceIteratorType inputIt( vesselness->GetOutput(), vesselness->GetOutput()->GetLargestPossibleRegion() );
  DestIteratorType outputIt(hessian_mask, vesselness->GetOutput()->GetLargestPossibleRegion() );

  for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
    outputIt.Set( inputIt.Get() == 0 );


  USWriterType::Pointer uswriter = USWriterType::New();
  uswriter->SetInput( hessian_mask );
  uswriter->SetFileName( "hess_mask.tiff" );
  uswriter->Update();
  
  

  return EXIT_SUCCESS;

}
