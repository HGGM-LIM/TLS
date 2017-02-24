
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVisualScoreImageFilter.h"

#include <stdio.h>
#include <string>
#include <time.h>
#include <list>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
namespace po = boost::program_options;

int main( int argc, char * argv[] )
{ 

	std::string infilepath;
	std::string outfilepath;
	int smoothMe  = 0; //Whether to smooth or not
	int smoothRadius = 1;

	unsigned int minArea=2000;
	int lowTh=-1024;
	int enfTh=-960;
	int tissueTh=-1;
	float lowestScoreTh=1.0;

	int verbose = 0;

	try {
    	// Declare the supported options.
    	po::options_description desc("Allowed options");
    	desc.add_options()
    	    ("help,h", "produce help message")
    	    ("input,i", po::value< std::string >(&infilepath), "input file")
    	    ("output,o", po::value< std::string >(&outfilepath)->default_value("enf.dcm"), "output emphysema mask")
    	    ("erode,e", po::value< int >(&smoothMe)->default_value(0), "Whether to erode or not the mask before calculating")
    	    ("radius,r", po::value< int >(&smoothRadius)->default_value(1), "Erosion radius")
    	    ("threshold,t", po::value< int >(&enfTh)->default_value(-960), "Threshold for emphysema quantification")
    	    ("low_threshold,l", po::value< int >(&lowTh)->default_value(-1024), "Low threshold for lung binarization")
    	    ("high_threshold,h", po::value< int >(&tissueTh)->default_value(-1), "High threshold for lung binarization")
    	    ("lowest_score,s", po::value< float >(&lowestScoreTh)->default_value(1.0), "If emphysema affectation is below this value, score is 0")
    	    ("verbose,v", "Show verbose output")
    	;

    	po::variables_map vm;
    	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    	po::notify(vm);


    	/*
    	 * General parameters handling
    	 */

    	if (vm.count("help") || (argc == 1)) {
    	    std::cout << desc << std::endl;
    	    return 1;
    	}

    	if (vm.count("verbose"))
    		verbose = 1;



    }
    catch(std::logic_error &err)
    {
         std::cerr << "std::logic_error:" << err.what() << std::endl;
         for (int i=0; i< argc; i++)
             std::cerr << argv[i] << std::endl;
    }
    catch(std::exception& e) {
    	std::cout << e.what() << std::endl;
    }
  
  typedef signed short PixelType;
  const unsigned int Dimension=3;
  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  
  typedef itk::VisualScoreImageFilter< ImageType, ImageType > VisualScoreType;
  VisualScoreType::Pointer radiologist = VisualScoreType::New();
  
  //Reading the lungs data
  reader->SetFileName(infilepath.c_str());

  try {
	  reader->Update();
  } catch (itk::ExceptionObject & excp) {
    std::cerr << "  Problems reading the input image"   << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  radiologist->SetVerbose(verbose);
  radiologist->SetInput(reader->GetOutput());
  radiologist->SetMinArea(minArea);
  radiologist->SetLowTh(lowTh);
  radiologist->SetEnfTh(enfTh);
  radiologist->SetTissueTh(tissueTh);
  radiologist->SetLowestScoreTh(lowestScoreTh);

  radiologist->Update();

  std::cout << "Visual Scores:" << std::endl;
  std::cout << "\tUpper " << radiologist->GetUppVisualScore() << std::endl;
  std::cout << "\tMedium " << radiologist->GetMedVisualScore() << std::endl;
  std::cout << "\tLower " << radiologist->GetLowVisualScore() << std::endl;
  std::cout << "\tTotal " << radiologist->GetVisualScore() << std::endl;

  writer->SetInput(radiologist->GetOutput());
  writer->SetFileName(outfilepath.c_str());
  writer->Update();

    
  return EXIT_SUCCESS;
  
}

