
#include "itkImage.h"
#include "filters/itkFastMarchingImageFilter.h"
//#include "itkFastMarchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage seedX seedY seedZ [minOut MAXOut]" << std::endl;
    return 1;
    }

  int minOut = 0;
  int maxOut = 100;

  if (argc == 8)
     minOut = atoi(argv[6]);
     maxOut = atoi(argv[7]);

  typedef   float           InternalPixelType;
  const     unsigned int    Dimension = 3;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;


  // Software Guide : BeginCodeSnippet
  typedef  itk::ImageFileReader< InternalImageType > ReaderType;
  typedef  itk::ImageFileWriter<  InternalImageType  > WriterType;
  // Software Guide : EndCodeSnippet


  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );


  typedef  itk::FastMarchingImageFilter< InternalImageType, InternalImageType >    FastMarchingFilterType;
  
  FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();

  
  fastMarching->SetInput( reader->GetOutput() );

  typedef FastMarchingFilterType::NodeContainer           NodeContainer;
  typedef FastMarchingFilterType::NodeType                NodeType;
  NodeContainer::Pointer seeds = NodeContainer::New();
  
  InternalImageType::IndexType  seedPosition;
  
  seedPosition[0] = atoi( argv[3] );
  seedPosition[1] = atoi( argv[4] );
  seedPosition[2] = atoi( argv[5] );


  NodeType node;
  const double seedValue = 0.0;
  
  node.SetValue( seedValue );
  node.SetIndex( seedPosition );

  seeds->Initialize();
  seeds->InsertElement( 0, node );

  fastMarching->SetTrialPoints(  seeds  );


  fastMarching->SetOutputSize( 
           reader->GetOutput()->GetBufferedRegion().GetSize() );


  typedef itk::MinimumMaximumImageCalculator< InternalImageType > MiniMaxType;
  MiniMaxType::Pointer miniMax = MiniMaxType::New();

  miniMax->SetImage( fastMarching->GetOutput() );
  miniMax->Compute();

  typedef itk::IntensityWindowingImageFilter< InternalImageType,InternalImageType > RescalerFilterType;
  RescalerFilterType::Pointer rescaler = RescalerFilterType::New();

  rescaler->SetInput( fastMarching->GetOutput() );
  rescaler->SetOutputMaximum( maxOut );
  rescaler->SetOutputMinimum( minOut );
  rescaler->SetWindowMaximum( maxOut );
  rescaler->SetWindowMinimum( 0 );

//  writer->SetInput( fastMarching->GetOutput() );
  writer->SetInput( rescaler->GetOutput() );
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
