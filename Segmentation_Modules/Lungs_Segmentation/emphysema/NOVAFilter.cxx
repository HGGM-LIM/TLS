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


  typedef float						PixelType;
  typedef itk::Image< PixelType, 2 >        ImageType;
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

  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(smoothRadius);
  NeighborhoodIteratorType it( radius, reader->GetOutput(),
                               reader->GetOutput()->GetRequestedRegion() );

  ImageType::Pointer output = ImageType::New();
  output->SetRegions(reader->GetOutput()->GetRequestedRegion());
  output->Allocate();

  IteratorType out(output, reader->GetOutput()->GetRequestedRegion());

  NeighborhoodIteratorType::OffsetType offset1 = {{-1,-1}};
  NeighborhoodIteratorType::OffsetType offset2 = {{0,-1}};
  NeighborhoodIteratorType::OffsetType offset3 = {{1,-1}};
  NeighborhoodIteratorType::OffsetType offset4 = {{-1,0 }};
  NeighborhoodIteratorType::OffsetType offset5 = {{1,0}};
  NeighborhoodIteratorType::OffsetType offset6 = {{-1,1}};
  NeighborhoodIteratorType::OffsetType offset7 = {{0,1}};
  NeighborhoodIteratorType::OffsetType offset8 = {{1,1}};


  for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out)
    {
    float sum;
    it.GetCenterPixel();
    sum = it.GetPixel(offset1) + it.GetPixel(offset2) + it.GetPixel(offset3) +
    		it.GetPixel(offset4) + it.GetPixel(offset5) + it.GetPixel(offset6) +
    		it.GetPixel(offset7) + it.GetPixel(offset8);
    out.Set(sum/8.0);
    }


  typedef unsigned char                          WritePixelType;
  typedef itk::Image< WritePixelType, 2 >        WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType > WriterType;

  typedef itk::RescaleIntensityImageFilter<
               ImageType, WriteImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput(output);

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



