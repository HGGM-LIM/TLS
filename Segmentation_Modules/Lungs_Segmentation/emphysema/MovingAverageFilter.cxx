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

	try {
		// Declare the supported options.
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message")
			("input,i", po::value< std::string >(&infilepath), "input file")
			("output,o", po::value< std::string >(&outfilepath),"output file")
			("radius,r", po::value< int >(&smoothRadius)->default_value(1), "Erosion radius")

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


  std::cout << "Parameters:" << std::endl;
  std::cout << "Radius: " << smoothRadius << std::endl;

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

    // A radius of 1 in all axial directions gives a 3x3x3x3x... neighborhood.
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(smoothRadius);

 // Initializes the iterators on the input & output image regions
 NeighborhoodIteratorType it(radius, reader->GetOutput(),
                                 output->GetRequestedRegion());
 ImageIterator out(output, output->GetRequestedRegion());

  // Iterates over the input and output
 for (it.Begin(), out = out.Begin(); ! it.IsAtEnd(); ++it, ++out )
   {
     float accum = 0.0;
     for (unsigned int i = 0; i < it.Size(); ++i)
       {
         accum += it.GetPixel(i);
       }
     out.Set(accum/(float)(it.Size())); //Normalize to preserve the mean

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




