#include "itkImage.h"
#include "itkImageFileReader.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
namespace po = boost::program_options;

int main( int argc, char ** argv )
{

	std::vector<float> seeds;
	std::vector<float>::const_iterator cit;
	bool isIndex=false;
	bool isPoint=false;

	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("seed,s", po::value< std::vector<float> >(), "Coordinates to probe (call it many times to set x,y,z in sequence)")
		("index,i", "Interpret specified coordinates as image grid indexes")
		("point,p", "Interpret specified coordinates as physical points")
		("in", po::value< std::string >(), "Image to probe")
	;


	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);


	/*
	 * General parameters handling
	 */

	if (vm.count("index") && vm.count("point")) {
		std::cout << "Please specify only --index or --point option" << std::endl;
		return 1;
	}

	if (vm.count("index") || vm.count("point")) {
		if (vm.count("index")) {
			isIndex = true;
			isPoint = false;
		} else {
			isIndex = false;
			isPoint = true;
		}

	}

	if (vm.count("help") || (argc == 1)) {
		std::cout << desc << std::endl;
		return 1;
	}


	if (vm.count("seed")) {

		seeds = vm["seed"].as<std::vector<float> >();

	} else {
		std::cout << "Please specify the coordinates to probe" << std::endl ;
		return 1;
	}


	if ( ! vm.count("in") ) {
		std::cout << "Please specify the image to probe" << std::endl;
		return 1;
	}



// We start by defining the PixelType and ImageType

  
  typedef double PixelType;
  
  const unsigned int Dimension = 3;
    

  typedef itk::Image< PixelType, Dimension >  ImageType;
    
// The image type is used as a template parameter to instantiate
// the reader.

  typedef itk::ImageFileReader < ImageType >  ReaderType;
  

  ReaderType::Pointer reader = ReaderType::New();
  
  reader->SetFileName( vm["in"].as<std::string >() );

  try 
    { 
    reader->Update();
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught!" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
  
  if (isIndex) {

	  ImageType::IndexType idx;
	  idx.Fill(0);

	  int i=0;
	  for(cit=seeds.begin(); cit!=seeds.end(); cit++)
		  idx[i++]=(int)*cit; // Convert to int because it is an index

	  std::cout << reader->GetOutput()->GetPixel(idx) << std::endl;
  }

  if (isPoint) {
	  typedef itk::Point< double, ImageType::ImageDimension > PointType;
	  PointType point;
	  ImageType::IndexType idx;

	  int i=0;
	  for(cit=seeds.begin(); cit!=seeds.end(); cit++)
		  point[i++]=*cit; // Do not convert: it is a real point

	  bool isInside = reader->GetOutput()->TransformPhysicalPointToIndex( point, idx );
	  if ( isInside )
		  std::cout << reader->GetOutput()->GetPixel(idx) << std::endl;
	  else {
		  std::cout << "Point " << point << " gave an index of " << idx << " which is outside image region:" << std::endl;
		  std::cout << reader->GetOutput()->GetLargestPossibleRegion() << std::endl;
	  }

  }

  return EXIT_SUCCESS;
}

