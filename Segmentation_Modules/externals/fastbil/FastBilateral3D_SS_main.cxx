#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkFastBilateralImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
namespace po = boost::program_options;

/**
 * This test was originally taken from the tests for the itkBilateralImageFilter
 * and modified for the itkFastBilateralImageFilter.
 */
int main(int argc, char* argv[] )
{
    std::string infile;
    std::string outfile;
    // Initialization values from standard bilateral filter
    float sigma_domain = 4.0;
    float sigma_range = 50.0;
    
      try {
    	// Declare the supported options.
    	po::options_description desc("Allowed options");
    	desc.add_options()
    	    ("help,h", "produce help message")
    	    ("domain-kernel,d", po::value<float>(&sigma_domain)->default_value(4.0), "Sigma of gaussian filter used for smoothing on spatial domain")
            ("range-kernel,r", po::value<float>(&sigma_range)->default_value(50.0), "Sigma of gaussian filter used for smoothing on range domain")
            ("infile,i", po::value<std::string>(&infile), "Input image file")
            ("outfile,o", po::value<std::string>(&outfile), "Output image file")
    	;

    	po::variables_map vm;
    	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    	po::notify(vm);
    	
    	if (vm.count("help") || (argc == 1)) {
    	    std::cout << desc << std::endl;
    	    return 1;
    	}
    	
    	
    } catch(std::logic_error &err) {
         std::cerr << "std::logic_error:" << err.what() << std::endl;
         for (int i=0; i< argc; i++)
             std::cerr << argv[i] << std::endl;
    } catch(std::exception& e) {
    	std::cout << e.what() << std::endl;
    }

  typedef signed short PixelType;
  const unsigned int dimension = 3;
  typedef itk::Image<PixelType, dimension> myImage;
  itk::ImageFileReader<myImage>::Pointer input 
    = itk::ImageFileReader<myImage>::New();
  input->SetFileName(infile.c_str());
  
  // Create a filter
  typedef itk::FastBilateralImageFilter<myImage,myImage> FilterType;
  
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput(input->GetOutput());
  
  // these settings reduce the amount of noise by a factor of 10
  // when the original signal to noise level is 5
  filter->SetDomainSigma( sigma_domain );
  filter->SetRangeSigma( sigma_range );
  
  
  try
    {
    input->Update();
    filter->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception detected: "  << e.GetDescription();
    return -1;
    }
  catch (...)
    {
    std::cerr << "Some other exception occurred" << std::endl;
    return -2;
    }

  // Generate test image
  itk::ImageFileWriter<myImage>::Pointer writer;
    writer = itk::ImageFileWriter<myImage>::New();
    writer->SetInput( filter->GetOutput() );
    writer->SetFileName( outfile.c_str() );
    writer->Update();

  return EXIT_SUCCESS;
}
