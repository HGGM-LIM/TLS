#include <iostream>
#include <fstream>

#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkSimilarityIndexImageFilter.h"
//#include "itkHausdorffDistanceImageFilter.h"
//#include "itkMeanAbsoluteSurfaceDistanceImageFilter.h"
//#include "itkAccuracyAssessmentImageFilter.h"

using namespace itk;

typedef unsigned short PixelType;
const unsigned int Dimension = 3;
typedef itk::Image< PixelType, Dimension >  ImageType;
typedef itk::Point<float,3> PointType;
typedef itk::ImageRegionIterator<ImageType>        IteratorType;
typedef ImageType::IndexType IndexType;
typedef itk::ImageFileReader < ImageType >  ReaderType;
typedef ImageType::Pointer ImagePointer;

double calcOverlap(std::vector<IndexType> v1, std::vector<IndexType> v2) {
	// use v1 as reference
	std::vector<IndexType>::iterator it1, it2;
	int count = v1.size(); // start from full
	for (it1=v1.begin() ; it1!=v1.end(); ++it1) {
//		std::cout << "searching for "<< (*it1) << " in v2";
		it2 = find (v2.begin(), v2.end(), (*it1));
		if (it2 == v2.end()) count--; // not found!
	}
//	std::cout << "overlap count: " << count << "\n";
	return  100 * (count / v1.size());
}

void matchTrees(std::map<int, std::vector<IndexType> > refTree, std::map<int, std::vector<IndexType> > targetTree, std::map<int, int> pcTree1, std::map<int, int> pcTree2) {
/**
 * - traverse both trees,
 * 		- establish parent-child r/s (read from CSV instead?)
 * 		- gather stats about each segment (length, points, radius - CSV provided too?)
 * - calculate best overlaps (which target segment is closest to which ref segment), considering the following
 * - create map linking ref_segmentNumber to target (if corresponding not found - use criteria, highlight)
 *
 */
	// #1 traverse both trees & gather stats - perhaps done by CSV already

	// if overlap is above acceptable threshold, use this as a match
	// if overlap is below certain threshold (as can be the case in mismatched label),
}

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

  
//  typedef double PixelType;

// The image type is used as a template parameter to instantiate
// the reader.

  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  
  reader1->SetFileName( argv[1] );
  reader2->SetFileName( argv[2] );
  
// Declare SI and Hausdorff distance Image Filters
  typedef itk::SimilarityIndexImageFilter< ImageType, ImageType > SIFilterType;
  SIFilterType::Pointer SIFilter = SIFilterType::New();
//  typedef itk::HausdorffDistanceImageFilter< ImageType, ImageType > HDFilterType;
//  HDFilterType::Pointer HDFilter = HDFilterType::New();
//  typedef itk::MeanAbsoluteSurfaceDistanceImageFilter< ImageType,ImageType >  MASDFilterType;
//  MASDFilterType::Pointer MASDFilter = MASDFilterType::New();
//  typedef itk::AccuracyAssessmentImageFilter< ImageType,
//  ImageType >  AccuracyFilterType;
//  AccuracyFilterType::Pointer accuracyFilter = AccuracyFilterType::New();
  
  SIFilter->SetInput1( reader1->GetOutput());
  SIFilter->SetInput2( reader2->GetOutput());
//  HDFilter->SetInput1( reader1->GetOutput());
//  HDFilter->SetInput2( reader2->GetOutput());
//  MASDFilter->SetInput1( reader1->GetOutput());
//  MASDFilter->SetInput2( reader2->GetOutput());
//  accuracyFilter->SetReferenceSegmentation( reader1->GetOutput() );
//  accuracyFilter->SetTestSegmentation( reader2->GetOutput() );

  
     
//  Finally, execution of the pipeline can be triggered by invoking the
//  Update() method in the SI filter. This call must be placed in a try/catch
//  block since exceptions be potentially be thrown in the process of reading
//  or writing the images

  try 
    { 
    SIFilter->Update();
    //HDFilter->Update();
    //MASDFilter->Update();
//    accuracyFilter->Update();
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught!" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
// Software Guide : EndCodeSnippet
  
  unsigned long volume1=0, volume2=0, intersection=0;

  ImagePointer image1 = reader1->GetOutput();
  ImagePointer image2 = reader2->GetOutput();

  std::map<int, std::vector<IndexType> > segmentMap1, segmentMap2;

  IteratorType resIt( image1, image1->GetLargestPossibleRegion() ), valIt( image2, image2->GetLargestPossibleRegion() );

  for ( resIt.GoToBegin(), valIt.GoToBegin(); !resIt.IsAtEnd(); ++resIt, ++valIt ) {
    if (resIt.Get()!=0) {
      volume1++;
      segmentMap1[resIt.Get()].push_back(resIt.GetIndex());
      if (valIt.Get()!=0) {
        volume2++;
        segmentMap2[valIt.Get()].push_back(valIt.GetIndex());
        intersection++;
      }
    }
    else {
      if (valIt.Get()!=0) {
    	  volume2++;
    	  segmentMap2[valIt.Get()].push_back(valIt.GetIndex());
      }
    }
  }
  ImageType::SpacingType spacing1 = image1->GetSpacing();
  ImageType::SpacingType spacing2 = image2->GetSpacing();
  double volumeFactor1 = 0.001*spacing1[0]*spacing1[1]*spacing1[2];
  double volumeFactor2 = 0.001*spacing2[0]*spacing2[1]*spacing2[2];
//  std::cout << "(Img1) Volume/pix: " << volumeFactor1 << ", (Img2) Volume/pix: "<< volumeFactor2 <<" mm3\n";

  double tanimotoVal = 100.0 * (double)(intersection) / ((double)(volume1+volume2-intersection));
  std::cout << "Similarity Index: " << SIFilter->GetSimilarityIndex() << std::endl;
  std::cout << "Volume (Image 1): " << volume1 << " pixels, "<< volumeFactor1*volume1<<" mm3\n";
  std::cout << "Volume (Image 2): " << volume2 << " pixels, "<< volumeFactor2*volume2<<" mm3\n";
  std::cout << "Intersection: " << intersection << " pixels\n";
  std::cout << "Tanimoto (both / [count1 + count2 - both]): " << tanimotoVal << " %\n";

  std::cout << "Similarity by segments: (img1 vs img2)\n";
  std::map< int, std::vector<IndexType> >::iterator it, it2;
  for (it = segmentMap1.begin(), it2 = segmentMap2.begin(); it!=segmentMap1.end(); ++it, ++it2) {
//	  std::cout << "Segment " << (int) (*it).first << ": "<< (*it).second.size() << " vs "<< (*it2).second.size() << "\n";
	  std::cout << "Segment " << (int) (*it).first << ": "<< (*it).second.size() << " vs "<< (*it2).second.size() << ", overlap: "<< calcOverlap(it->second, it2->second)<<"\n";
  }


//  std::cout << "Haussdorf distance: " << HDFilter->GetHausdorffDistance() << std::endl;
//  std::cout << "Mean Absolute Surface Distance: " << MASDFilter->GetMeanAbsoluteSurfaceDistance() << std::endl;
//  std::cout << "TPVF: " << accuracyFilter->GetTruePositiveVF() << std::endl;
//  std::cout << "TNVF: " << accuracyFilter->GetTrueNegativeVF() << std::endl;
//  std::cout << "FPVF: " << accuracyFilter->GetFalsePositiveVF() << std::endl;
//  std::cout << "FNVF: " << accuracyFilter->GetFalseNegativeVF() << std::endl;
  
  //write results in a file
  std::ofstream myfile;
  myfile.open(argv[3]);
  myfile << "# Segmentation evaluation\n";
  myfile << "Image 1: " << argv[1] << "\n";
  myfile << "Image 2: " << argv[2] << "\n";
  myfile << "Similarity Index: " << SIFilter->GetSimilarityIndex() << "\n";
//  myfile << "Haussdorf distance: " << HDFilter->GetHausdorffDistance() << "\n";
//  myfile << "Mean Absolute Surface Distance: " << MASDFilter->GetMeanAbsoluteSurfaceDistance() <<"\n";
//  myfile << "TPVF: " << accuracyFilter->GetTruePositiveVF() << "\n";
//  myfile << "TNVF: " << accuracyFilter->GetTrueNegativeVF() << "\n";
//  myfile << "FPVF: " << accuracyFilter->GetFalsePositiveVF() << "\n";
//  myfile << "FNVF: " << accuracyFilter->GetFalseNegativeVF() << "\n";
  myfile.close();
  
  return EXIT_SUCCESS;
}
