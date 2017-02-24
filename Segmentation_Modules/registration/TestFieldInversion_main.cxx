#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkVector.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPoint.h"


#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
namespace po = boost::program_options;

// We start by defining the PixelType and ImageType
  const unsigned int Dimension = 3;

  typedef itk::Vector< float, Dimension > VPixelType;
  typedef itk::Image< VPixelType, Dimension > VImageType;
  typedef int IPixelType;
  typedef itk::Image< IPixelType, Dimension > IImageType;
  typedef itk::Point<float, Dimension> PointType;

// The image type is used as a template parameter to instantiate
// the reader.

  typedef itk::ImageFileReader < VImageType >  FieldReaderType;
  typedef itk::ImageFileReader < IImageType >  ImageReaderType;

  typedef struct {
  	PointType f_p;						// The point in the fixed image
  	IImageType::IndexType f_p_idx;		// Index of the point in the fixed image
  	PointType m_p;						// The point in the moving image
  	IImageType::IndexType m_p_idx;		// Index of the point in the moving image
  	VPixelType t_vector;				// Transform vector

  } TransformType;



TransformType transform(IImageType::IndexType, IImageType::Pointer, VImageType::Pointer, IImageType::Pointer);
IImageType::Pointer load_image(std::string, std::string);
VImageType::Pointer load_field(std::string, std::string);





int main( int argc, char ** argv )
{

	std::vector<float> seeds;
	std::vector<float>::const_iterator cit;

	bool testPoint = false;
	bool testImage = false;
	bool testTransform = false;

	bool withMask = false;
	// Declare the supported options.
	po::options_description desc("Test inversion of a deformation field");
	desc.add_options()
		("help,h", "produce help message")
		("point,p", "Test a single point")
		("isPoint", "Interpret coordinates as physical point")
		("isIndex", "Interpret coordinates as grid index")
		("image,i", "Test the full image")
		("values,v", "Test values")
		("transform,t", "Test transformation")
		("convert,c", "Apply the specified deformation field to a point/index and print the result")
		("def,d", po::value< std::string >(), "Deformation field")
		("definv", po::value< std::string >(), "Inverse deformation field")
		("fixed,f", po::value< std::string >(), "Fixed image")
		("moving,m", po::value< std::string >(), "Moving image")
		("mask,x", po::value< std::string >(), "Mask to exclude points")
		("seed,s", po::value< std::vector<float> >(), "Only compute at this point")
	;


	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);


	/*
	 * General parameters handling
	 */

	if (vm.count("help") || (argc == 1)) {
		std::cout << desc << std::endl;
		return EXIT_SUCCESS;
	}

	if (vm.count("point") && vm.count("image")) {
		std::cerr << "You might want to test either a single point or the whole image but not both." << std::endl;
		return EXIT_FAILURE;
	}

	if (vm.count("isPoint") && vm.count("isIndex")) {
			std::cerr << "I can consider coordinates either as a physical point or as an index but not both." << std::endl;
			return EXIT_FAILURE;
		}

	if (vm.count("convert") && ( (!vm.count("isPoint")) && (!vm.count("isIndex")) )) {
		std::cerr << "Please specify if to interpret the coordinates as a physical point or as an index" << std::endl;
		return EXIT_FAILURE;
	}


	// We always need the deformation field, the fixed image and the moving one
	if ( ! vm.count("def") ) {
			std::cerr << "Please specify the deformation field" << std::endl;
			return EXIT_FAILURE;
		}

	if ( ! vm.count("fixed") ) {
		std::cerr << "Please specify the fixed image" << std::endl;
		return EXIT_FAILURE;
	}

	if ( ! vm.count("moving") ) {
		std::cerr << "Please specify the moving image" << std::endl;
		return EXIT_FAILURE;
	}

	// If we operate on the entire image, let's assure that the user specifies one
	if (vm.count("image")) {

		testImage = true;

		if ( ! vm.count("fixed") ) {
				std::cout << "Please specify the fixed image" << std::endl;
				return EXIT_FAILURE;
			}

	}

	// If we operate on a single point, let's assure that the user specifies one
	if (vm.count("point") || vm.count("convert")) {

		if (! vm.count("seed") ){
			std::cerr << "Please specify the coordinates of the point with the --seed option" << std::endl;
			return EXIT_FAILURE;
		}

		seeds = vm["seed"].as<std::vector<float> >();
		testPoint = true;
	}


	VImageType::Pointer p_field = load_field( vm["def"].as<std::string >(), "deformation field" );
	IImageType::Pointer p_original_im = load_image( vm["fixed"].as<std::string >(), "fixed image (A)" );
	IImageType::Pointer p_warped_im = load_image( vm["moving"].as<std::string >(), "moving image (B)" );
	IImageType::Pointer p_mask_im;
	VImageType::Pointer p_inv_field;


	/*
	 * EXECUTE THE ANALYSIS
	 *
	 *
	 * We can perform three types of analysis:
	 * 	1) We can check that the values of the warped moving image and the
	 * 		original moving image are the same.
	 * 		In order to do this specify the warped moving as the fixed one. We will iterate
	 * 		over all it's pixel and check
	 * 	2) We can transform each point with the forward and the inverse deformation field
	 * 		and check that the final coordinates are equal to the original one
	 *  3) We simply apply the deformation field to a user-given point and print back the result
	 */


	  if (vm.count("values")) {



		  std::cout << "Index_A; " <<
				  "PhysicalPoint_A; "<<
				  "A(Index_A); "<<
				  "DF(PhysicalPoint_A); " <<
				  "PhysicalPoint_B = A(PhysicalPoint_A) + DF(PhysicalPoint_A); " <<
				  "Index_B; " <<
				  "B(Index_B)" << std::endl;



		  if (testPoint ) {


			  /*
			   * Let's get the index of the fixed image from the command line
			   */
			  IImageType::IndexType f_p_idx;

			  int i=0;
			  for(cit=seeds.begin(); cit!=seeds.end(); cit++)
				  f_p_idx[i++]=(int)*cit; // Convert to int because it is an index

			  TransformType t = transform(f_p_idx, p_original_im, p_field, p_warped_im);

			  std::cout << t.f_p_idx <<
					  "; " << t.f_p <<
					  "; " << p_original_im->GetPixel(t.f_p_idx) <<
					  "; " << t.t_vector <<
					  "; " << t.m_p <<
					  "; " << t.m_p_idx <<
					  "; " << p_warped_im->GetPixel(t.m_p_idx)
					  << std::endl;



		  } else if (testImage) {

				if (withMask)
				  p_mask_im = load_image( vm["mask"].as<std::string >(), "mask" );

			  //Make the entire image
				typedef itk::ImageRegionConstIteratorWithIndex< IImageType > ImageConstIteratorType;
				typedef itk::ImageRegionConstIteratorWithIndex< VImageType > FieldConstIteratorType;

				ImageConstIteratorType oIt( p_original_im, p_original_im->GetLargestPossibleRegion() );


			  int limit=0;
			  for(oIt.Begin(); !oIt.IsAtEnd();++oIt) {

				  /*
				   * Let's get the index of the fixed image
				   */
				  IImageType::IndexType f_p_idx = oIt.GetIndex();

				  //Avoid non interesting points
				  if (withMask)
					  if (p_mask_im->GetPixel(f_p_idx) == 0)
						  continue;

				  if (limit++ > 10)
					  break;

				  TransformType t = transform(f_p_idx, p_original_im, p_field, p_warped_im);

				  std::cout << t.f_p_idx <<
						  "; " << t.f_p <<
						  "; " << p_original_im->GetPixel(t.f_p_idx) <<
						  "; " << t.t_vector <<
						  "; " << t.m_p <<
						  "; " << t.m_p_idx <<
						  "; " << p_warped_im->GetPixel(t.m_p_idx)
						  << std::endl;

			  }
		  }
	  }


	if (vm.count("transform")) {



		  if (! vm.count("definv")) {
			  std::cerr << "Please specify the inverse deformation field with the --definv option" << std::endl;
			  return EXIT_FAILURE;
		  } else {
			  //Load the inverse deformation field
			  p_inv_field = load_field( vm["definv"].as<std::string >(), "inverse deformation field" );

		  }

		  if (testPoint) {
			  /*
			   * Let's get the index of the fixed image from the command line
			   */
			  IImageType::IndexType f_p_idx;

			  int i=0;
			  for(cit=seeds.begin(); cit!=seeds.end(); cit++)
				  f_p_idx[i++]=(int)*cit; // Convert to int because it is an index

			  std::cout << "Starting point: " << f_p_idx << std::endl;

			  /*
			   * We go from the fixed image index to the corresponding point into the
			   * moving image.
			   */
			  TransformType direct = transform(f_p_idx, p_original_im, p_field, p_warped_im);

			  std::cout << "Direct transform: " << std::endl;
			  std::cout << "Index_A; " <<
			 			   "PhysicalPoint_A; "<<
			 			   "DF(PhysicalPoint_A); " <<
			 			   "PhysicalPoint_B = A(PhysicalPoint_A) + DF(PhysicalPoint_A); " <<
			 			   "Index_B; " << std::endl;

			  std::cout << direct.f_p_idx <<
					  "; " << direct.f_p <<
					  "; " << direct.t_vector <<
					  "; " << direct.m_p <<
					  "; " << direct.m_p_idx  << std::endl;

			  /*
			   * Now we want to go back from the index of the moving image to the corresponding
			   * point into the fixed image and see if the coordinates are the same.
			   */

			  // First, test if the point is inside or outside the domain of the transform
			  if (! p_inv_field->GetLargestPossibleRegion().IsInside(direct.m_p_idx)) {

				  std::cerr << "Sorry, index " << direct.m_p_idx <<
						  " is outside the largest possible region of the inverse deformation field." <<  std::endl;
				  return EXIT_FAILURE;
			  }

			  TransformType inverse = transform(direct.m_p_idx, p_warped_im, p_inv_field, p_original_im);

			  std::cout << inverse.f_p_idx <<
					  "; " << inverse.f_p <<
					  "; " << inverse.t_vector <<
					  "; " << inverse.m_p <<
					  "; " << inverse.m_p_idx  << std::endl;



		  } else if (testImage) {

			  std::cerr << "Full image transformation test not implemented!" << std::endl;
		  }

	  }

	if (vm.count("convert")) {


		IImageType::IndexType f_p_idx;

		if (vm.count("isIndex")) {

			  /*
			   * Let's get the index of the fixed image from the command line
			   */

			  int i=0;
			  for(cit=seeds.begin(); cit!=seeds.end(); cit++)
				  f_p_idx[i++]=(int)*cit; // Convert to int because it is an index


		} else if (vm.count("isPoint")) {

			  PointType point;

			  int i=0;
			  for(cit=seeds.begin(); cit!=seeds.end(); cit++)
				  point[i++]=*cit; // Do not convert: it is a real point

			  bool isInside = p_original_im->TransformPhysicalPointToIndex( point, f_p_idx );
			  if (! isInside ) {
				  std::cout << "Point " << point << " corresponds to index " << f_p_idx << " which is outside image region:" << std::endl;
				  return EXIT_FAILURE;

			  }
		}


		TransformType direct = transform(f_p_idx, p_original_im, p_field, p_warped_im);

		std::cout << direct.f_p << " --> " << direct.m_p << std::endl;

	}

  return EXIT_SUCCESS;
}

IImageType::Pointer load_image(std::string fname, std::string descr) {
	/*
	 * Load a scalar image, given it's name and description
	 */

	ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->SetFileName( fname );

	try
	  {
	  std::cout << "Reading " << descr << "...";
	  reader->Update();
	  std::cout << "done!" << std::endl;
	  }
	catch( itk::ExceptionObject & err )
	  {
	  std::cerr << "Cannot read " << descr << "!" << std::endl;
	  std::cerr << err << std::endl;

	  }

	std::cout << descr << " region: " << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

	return reader->GetOutput();
}

VImageType::Pointer load_field(std::string fname, std::string descr) {
	/*
	 * Load a scalar image, given it's name and description
	 */

	FieldReaderType::Pointer reader = FieldReaderType::New();
	reader->SetFileName( fname );

	try
	  {
	  std::cout << "Reading " << descr << "...";
	  reader->Update();
	  std::cout << "done!" << std::endl;
	  }
	catch( itk::ExceptionObject & err )
	  {
	  std::cerr << "Cannot read " << descr << "!" << std::endl;
	  std::cerr << err << std::endl;

	  }

	std::cout << descr << " region: " << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

	return reader->GetOutput();
}

TransformType transform(IImageType::IndexType f_p_idx, IImageType::Pointer p_fixed_im, VImageType::Pointer p_field, IImageType::Pointer p_moving_im) {
	  /*
	   1) Get the index of the fixed image
	   2) Read the corresponding index of the def field
	   3) Add the displacement vector to the point and get the transformed point in the fixed image space
	   4) Transform this point into the moving image space
	   */
	TransformType t;
	t.f_p_idx = f_p_idx;

	/*
	* Now we transform the index in a point in the space
	* of the fixed image
	*/

	p_fixed_im->TransformIndexToPhysicalPoint(f_p_idx, t.f_p);

	/*
	* Let's take the transformation from the def field image
	* and generate a new point in the physical space of the
	* moving image
	*/
	t.t_vector = p_field->GetPixel(f_p_idx);
	t.m_p = t.f_p + t.t_vector;

	/*
	* Convert back from the point to the index coordinates
	* of the moving image grid
	*/
	p_moving_im->TransformPhysicalPointToIndex( t.m_p, t.m_p_idx );

	return t;

}



