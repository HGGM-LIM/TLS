
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkEmphysemaCalculationsImageFilter.h"
#include "dicom_utilities.h"

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
	std::string clusterFile;
	int smoothMe  = 0; //Whether to smooth or not
	int smoothRadius = 1;
	int enf_th = -960;

	try {
    	// Declare the supported options.
    	po::options_description desc("Allowed options");
    	desc.add_options()
    	    ("help,h", "produce help message")
    	    ("input,i", po::value< std::string >(&infilepath), "input file")
    	    ("output,o", po::value< std::string >(&outfilepath)->default_value("enf.dcm"), "output emphysema mask")
    	    ("cluster,c", po::value< std::string >(&clusterFile)->default_value("lavdata.txt"), "Where to save cluster data")
    	    ("erode,e", po::value< int >(&smoothMe)->default_value(0), "Whether to erode or not the mask before calculating")
    	    ("radius,r", po::value< int >(&smoothRadius)->default_value(1), "Erosion radius")
    	    ("threshold,t", po::value< int >(&enf_th)->default_value(-960), "Threshold for emphysema quantification")
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


    }
    catch(std::logic_error &err)
    {
         std::cerr << "std::logic_error:" << err.what() << std::endl;
         for (unsigned int i=0; i< argc; i++)
             std::cerr << argv[i] << std::endl;
    }
    catch(std::exception& e) {
    	std::cout << e.what() << std::endl;
    }

  std::cout << "Reading input serie  " << infilepath << std::endl;
  
  struct timeb initialTime;
  //struct timeb actualTime;
  struct timeb time0;
	ftime(&initialTime);
	ftime(&time0);

  typedef signed short PixelType;
  typedef unsigned short BinaryPixelType;
  const unsigned int Dimension=3;

  typedef unsigned long LabelType;

  typedef itk::Image<BinaryPixelType, Dimension>                       	BinaryImageType;
  typedef itk::Image<PixelType, Dimension>                       		ImageType;
  typedef itk::ShapeLabelObject<LabelType, Dimension>            		LabelObjectType;
  typedef itk::LabelMap<LabelObjectType>                         		LabelCollectionType;
	
	typedef ImageType::IndexType                                     	IndexType;
  typedef ImageType::SizeType                                      	SizeType;
  typedef ImageType::RegionType                                    	OutputImageRegionType;
  
  	typedef itk::ImageFileReader< ImageType >                      	 	ReaderType;
	ReaderType::Pointer reader =                                     	ReaderType::New();
 
	typedef itk::ImageFileWriter< BinaryImageType >                        	WriterType;
	WriterType::Pointer writer =                                     	WriterType::New();
	
	typedef itk::EmphysemaCalculationsImageFilter<ImageType, BinaryImageType> CalculationsType;
	CalculationsType::Pointer calculus =															CalculationsType::New( );
  
  
  //Reading the lungs data
  reader->SetFileName(infilepath.c_str());

  try {
	  reader->Update();
  } catch (itk::ExceptionObject & excp) {
    std::cerr << "  Problems reading the input image"   << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  calculus->SetInput(reader->GetOutput());
  calculus->SetBackup(1);
  calculus->SetAutomaticallyResize(1);
  calculus->SetNumberOfBins(255);
  calculus->SetMarginScale(10.0);
  calculus->SetLowEmphValue(enf_th);
  calculus->SetSmooth(smoothMe);
  calculus->SetSmoothRadius(smoothRadius);
  calculus->SetClusterFile(clusterFile);
  calculus->Update();
  BinaryImageType::Pointer outputImagePointer = calculus->GetOutput();



  /*
   * FIXME:
  copy_dicom_data(reader->GetOutput(),outputImagePointer);
  gives a segfault!!
  */

  writer->SetInput(outputImagePointer);
  writer->SetFileName(outfilepath.c_str());
  writer->Update();

    
  return EXIT_SUCCESS;
  
}

