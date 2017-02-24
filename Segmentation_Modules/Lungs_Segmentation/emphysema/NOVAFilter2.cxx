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

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodOperatorImageFilter.h"
#include "itkNOVAOperator.h"

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
			("output,o", po::value< std::string >(&outfilepath)->default_value("enf.dcm"), "output emphysema mask")
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


  const unsigned int DIM = 2;
  typedef float						PixelType;
  typedef itk::Image< PixelType, DIM >        ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>        IteratorType;

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

  typedef itk::NOVAOperator<float, DIM> NOVAOperatorType;
  NOVAOperatorType novaOperator;
  itk::Size<DIM> radius;
  radius.Fill(1); // a radius of 1x1 creates a 3x3 operator
  novaOperator.SetDirection(0); // Create the operator for the X axis derivative
  novaOperator.CreateToRadius(radius);

  typedef itk::NeighborhoodOperatorImageFilter<ImageType, ImageType> NeighborhoodOperatorImageFilterType;
  NeighborhoodOperatorImageFilterType::Pointer filter = NeighborhoodOperatorImageFilterType::New();
  filter->SetOperator(novaOperator);
  filter->SetInput(reader->GetOutput());
  filter->Update();

  typedef unsigned char                          WritePixelType;
  typedef itk::Image< WritePixelType, DIM >        WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType > WriterType;

  typedef itk::RescaleIntensityImageFilter<
               ImageType, WriteImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput(filter->GetOutput());

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outfilepath.c_str() );
  writer->SetInput(rescaler->GetOutput());
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



  return 0;
}

