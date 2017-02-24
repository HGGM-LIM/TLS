/*
 * NOVAFilter.cxx
 *
 *  Created on: May 25, 2012
 *      Author: mceresa
 */


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

#include <math.h>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
namespace po = boost::program_options;


int main( int argc, char ** argv )
{

	std::string infilepath;
	std::string outfilepath;
	int smoothRadius = 1;
	float sigma = 1.0;
	float p = 0.0;
	float d = 3.0;

	try {
		// Declare the supported options.
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message")
			("input,i", po::value< std::string >(&infilepath), "input file")
			("output,o", po::value< std::string >(&outfilepath),"output file")
			("radius,r", po::value< int >(&smoothRadius)->default_value(1), "Erosion radius")
			("sigma,s", po::value< float >(&sigma), "Noise level")
			("plateau,p", po::value< float >(&p)->default_value(0.0), "Plateau level")
			("cutoff,d", po::value< float >(&d)->default_value(3.0), "Cutoff level")
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

		/*
		 * Test values of p and d
		 */
		if ((p>d) or (p==d)) {
			std::cout << "Plateau parameter p must be smaller than the cutoff d" << std::endl;
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


  std::cout << "Parameters:" << std::endl;
  std::cout << "Radius: " << smoothRadius << std::endl;
  std::cout << "Sigma: " << sigma << std::endl;

  const unsigned int DIM = 3;
  typedef float						PixelType;
  typedef itk::Image< PixelType, DIM >        ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>        IteratorType;
  typedef itk::ImageRegionIterator<ImageType>       ImageIterator;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( infilepath.c_str() );
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  ImageType::Pointer output = ImageType::New();
  output->SetRegions(reader->GetOutput()->GetRequestedRegion());
  output->Allocate();
/*
  ImageType::Pointer noise_image = ImageType::New();
  noise_image->SetRegions(reader->GetOutput()->GetRequestedRegion());
  noise_image->Allocate();
*/
    // A radius of 1 in all axial directions gives a 3x3x3x3x... neighborhood.
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(smoothRadius);

 // Initializes the iterators on the input & output image regions
 NeighborhoodIteratorType it(radius, reader->GetOutput(),
                                 output->GetRequestedRegion());
 ImageIterator out(output, output->GetRequestedRegion());
 //ImageIterator noise(noise_image, noise_image->GetRequestedRegion());

 // Iterates over the input and output
 for (it.Begin(), out = out.Begin(); ! it.IsAtEnd(); ++it, ++out )
   {
     float accum = 0.0;
     float acc_w = 0.0;
     //float mean = 0.0;
     for (unsigned int i = 0; i < it.Size(); ++i)
       {
    	 // Calculate weight
    	 float I = vnl_math_abs( ( it.GetCenterPixel() - it.GetPixel(i) ) / sigma );
    	 float NOVAweight = 0.0;
    	 if (I <= p)
    		 NOVAweight = 1.0;
    	 else if (I >= d)
    		 NOVAweight = 0.0;
    	 else
    		 NOVAweight = -1.0/(p-d) * I + d/(p-d);

    	 acc_w += NOVAweight;
         accum += it.GetPixel(i)*NOVAweight;
         //mean += it.GetPixel(i);
       }
     out.Set(accum/(float)(acc_w)); //Normalize to preserve the mean
     /*
     mean += it.GetCenterPixel();
     mean = mean/(float)(it.Size());
     float sigma_est = 0.0;
     for (unsigned int i = 0; i < it.Size(); ++i)
        {
    	 sigma_est += pow(( it.GetPixel(i) - mean),2);
        }
      sigma_est = sqrt(sigma_est/(float)(it.Size()));
      noise.Set(sigma_est); //Write the estimated noise

      */
   }



  typedef signed short                          WritePixelType;
  typedef itk::Image< WritePixelType, DIM >        WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType > WriterType;

  typedef itk::CastImageFilter< ImageType, WriteImageType > CasterType;

  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(output);

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outfilepath.c_str() );
  writer->SetInput(caster->GetOutput());
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

/*
  caster->SetInput(noise_image);
  writer->SetFileName( "noise.tiff" );
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

*/
  return 0;
}




