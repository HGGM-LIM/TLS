#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"


int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage sx sy sz" << std::endl;
    return 1;
    }

  
  typedef   float           InternalPixelType;
  typedef   unsigned char  OutputPixelType;
  const     unsigned int    Dimension = 3;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter< InternalImageType, OutputImageType >  CastingFilterType;
  CastingFilterType::Pointer caster = CastingFilterType::New();

  // Software Guide : BeginCodeSnippet
  typedef  itk::ImageFileReader< InternalImageType > ReaderType;
  typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
  // Software Guide : EndCodeSnippet


  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );

  typedef itk::ConfidenceConnectedImageFilter<InternalImageType, InternalImageType> ConfidenceConnType;
  ConfidenceConnType::Pointer confConn = ConfidenceConnType::New();

  confConn->SetInput( reader->GetOutput() );
  
  InternalImageType::IndexType seed;
  seed[0] = atoi(argv[3]);
  seed[1] = atoi(argv[4]);
  seed[2] = atoi(argv[5]);

  confConn->SetSeed( seed );
  confConn->SetReplaceValue( 255 );
  confConn->SetNumberOfIterations( 5 );

  caster->SetInput( confConn->GetOutput() );

  typedef itk::VotingBinaryIterativeHoleFillingImageFilter< OutputImageType > FilterType;
  FilterType::Pointer holeFill = FilterType::New();

  holeFill->SetInput( caster->GetOutput() );

  OutputImageType::SizeType radius;
  radius[0] = 1;
  radius[1] = 1;
  radius[2] = 1;

  holeFill->SetRadius( radius );
  holeFill->SetBackgroundValue( 0 );
  holeFill->SetForegroundValue( 255 );
  holeFill->SetMajorityThreshold( 2 );

//filter->SetMaximumNumberOfIterations( 5 );
  writer->SetInput( holeFill->GetOutput() );
  writer->SetFileName( argv[2] );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
 
 
  return 0;
}
