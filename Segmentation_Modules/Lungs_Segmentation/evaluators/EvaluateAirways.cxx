
#include <iostream>
#include <fstream>

#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkSimilarityIndexImageFilter.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkMeanAbsoluteSurfaceDistanceImageFilter.h"
#include "itkAccuracyAssessmentImageFilter.h"

int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile1 inputImageFile2 resultsFile" << std::endl;
    return EXIT_FAILURE;
    }


// We start by defining the PixelType and ImageType

  
  typedef double PixelType;
  
  const unsigned int Dimension = 3;
    

  typedef itk::Image< PixelType, Dimension >  ImageType;
    
// The image type is used as a template parameter to instantiate
// the reader.

  typedef itk::ImageFileReader < ImageType >  ReaderType;
  

  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );
  
// Declare SI and Hausdorff distance Image Filters
  typedef itk::SimilarityIndexImageFilter< ImageType, ImageType > SIFilterType;
  SIFilterType::Pointer SIFilter = SIFilterType::New();
  typedef itk::HausdorffDistanceImageFilter< ImageType, ImageType > HDFilterType;
  HDFilterType::Pointer HDFilter = HDFilterType::New();  
  typedef itk::MeanAbsoluteSurfaceDistanceImageFilter< ImageType,ImageType >  MASDFilterType;
  MASDFilterType::Pointer MASDFilter = MASDFilterType::New();
  typedef itk::AccuracyAssessmentImageFilter< ImageType,
  ImageType >  AccuracyFilterType;
  AccuracyFilterType::Pointer accuracyFilter = AccuracyFilterType::New();
  
  
  SIFilter->SetInput1( reader1->GetOutput());
  SIFilter->SetInput2( reader2->GetOutput());
  HDFilter->SetInput1( reader1->GetOutput());
  HDFilter->SetInput2( reader2->GetOutput());
  MASDFilter->SetInput1( reader1->GetOutput());
  MASDFilter->SetInput2( reader2->GetOutput());
  accuracyFilter->SetReferenceSegmentation( reader1->GetOutput() );
  accuracyFilter->SetTestSegmentation( reader2->GetOutput() );

  
     
//  Finally, execution of the pipeline can be triggered by invoking the
//  Update() method in the SI filter. This call must be placed in a try/catch
//  block since exceptions be potentially be thrown in the process of reading
//  or writing the images

  try 
    { 
    //SIFilter->Update();
    //HDFilter->Update();
    //MASDFilter->Update();
    accuracyFilter->Update();
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught!" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
// Software Guide : EndCodeSnippet
  
  std::cout << "Similarity Index: " << SIFilter->GetSimilarityIndex() << std::endl;
  std::cout << "Haussdorf distance: " << HDFilter->GetHausdorffDistance() << std::endl;
  std::cout << "Mean Absolute Surface Distance: " << MASDFilter->GetMeanAbsoluteSurfaceDistance() << std::endl;
  std::cout << "TPVF: " << accuracyFilter->GetTruePositiveVF() << std::endl;
  std::cout << "TNVF: " << accuracyFilter->GetTrueNegativeVF() << std::endl;
  std::cout << "FPVF: " << accuracyFilter->GetFalsePositiveVF() << std::endl;
  std::cout << "FNVF: " << accuracyFilter->GetFalseNegativeVF() << std::endl;
  
  //write results in a file
  std::ofstream myfile;
  myfile.open(argv[3]);
  myfile << "# Segmentation evaluation\n";
  myfile << "Image 1: " << argv[1] << "\n";
  myfile << "Image 2: " << argv[2] << "\n";
  myfile << "Similarity Index: " << SIFilter->GetSimilarityIndex() << "\n";
  myfile << "Haussdorf distance: " << HDFilter->GetHausdorffDistance() << "\n";
  myfile << "Mean Absolute Surface Distance: " << MASDFilter->GetMeanAbsoluteSurfaceDistance() <<"\n";
  myfile << "TPVF: " << accuracyFilter->GetTruePositiveVF() << "\n";
  myfile << "TNVF: " << accuracyFilter->GetTrueNegativeVF() << "\n";
  myfile << "FPVF: " << accuracyFilter->GetFalsePositiveVF() << "\n";
  myfile << "FNVF: " << accuracyFilter->GetFalseNegativeVF() << "\n";
  myfile.close();
  
  return EXIT_SUCCESS;
}
