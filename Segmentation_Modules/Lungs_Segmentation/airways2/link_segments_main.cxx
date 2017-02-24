#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <locale>
#include <string>

#include "log.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>

namespace po = boost::program_options;

typedef unsigned short PixelType;
const unsigned int Dimension = 3;
typedef itk::Image< PixelType, Dimension >  ImageType;
typedef itk::Point<double, ImageType::ImageDimension> PointType;

PointType linearNSS(const std::vector<PointType> points, PointType p) {
	/*
	 * Implements linear Nearest neighbor search
	 */

	pLogType l = get_logger("MAIN::linearNSS");
	double min_dist=10000;
	PointType closest;

	std::vector<PointType>::const_iterator it;
	for(it = points.begin(); it!=points.end(); it++) {
		double distance = p.EuclideanDistanceTo(*it);
		if (distance < min_dist) {
			min_dist = distance;
			closest = *it;
			l->debugStream() << "New closest is: " << closest << " with dist=" << min_dist;
		}
	}

	return closest;
}


/*
 * Basic function to trim a string.
 * Copied from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
 * FIXME: we could use boost
	#include <boost/algorithm/string.hpp>
	using namespace std;
	using namespace boost;
	string str1(" hello world! ");
	trim(str1);

 */
// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}


// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}




/*
 * FIXME: We don;t use it: remove.
 */
class Segment {
public:
	int id;
	int parent_id;
	Segment(int id, int pid) {
		this->id = id;
		this->parent_id = pid;
	}

};

class PropagationInfo {
	/*
	 * Here we store the info about the prop we read in the csv file. It is simpler than the corresponding class
	 * in global_itk.h but we could think to put it there and convert the actual Propinfo into a subclass
	 * so we have a single library which reads and writes the csv file of the segmentation results
	 */
public:
	int id;
	PointType centroid;
	float radius;
	int from_label;
	int to_label;
	PropagationInfo(int id, double cx, double cy, double cz, float radius, int from_label) {

		this->id = id;
		this->radius = radius;
		this->centroid[0] = cx;
		this->centroid[1] = cy;
		this->centroid[2] = cz;


	}




};


std::map<int, std::vector<PropagationInfo> > read_prop_info(const std::string fname ) {
	/*
	 * Parse the csv with the propagation info generated from the visitor code
	 */
	//pLogType l = get_logger("MAIN::read_prop_info");
	std::map<int, std::vector<PropagationInfo> > map;
	std::map<int, std::vector<PropagationInfo> >::iterator it;

    std::ifstream data(fname.c_str());

	std::string line;
	std::getline(data,line); //Remove header!
	while(std::getline(data,line))
	{
		std::stringstream  lineStream(line);
		std::string        segid;
		std::string        pid;
		std::string        prop_id;
		std::string        cx;
		std::string        cy;
		std::string        cz;
		std::string        radius;

		std::getline(lineStream,segid,',');
		std::getline(lineStream,pid,',');
		std::getline(lineStream,prop_id,',');
		std::getline(lineStream,cx,',');
		std::getline(lineStream,cy,',');
		std::getline(lineStream,cz,',');
		std::getline(lineStream,radius,',');

		int sid = atoi(segid.c_str());
		if (map.find(sid) == map.end()) {
			std::vector<PropagationInfo> props;
			map[sid] = props;
		}

		map[sid].push_back(PropagationInfo(atoi(prop_id.c_str()), atof(cx.c_str()), atof(cy.c_str()), atof(cz.c_str()), atof(radius.c_str()),atoi(segid.c_str())));


	}

	return map;

}

std::vector<PointType> as_pointset(std::map<int, std::vector<PropagationInfo> > data) {

	std::vector<PointType> pointset;
	std::map<int, std::vector<PropagationInfo> >::const_iterator cit;
	std::vector<PropagationInfo>::const_iterator cit2;
	for (cit=data.begin(); cit != data.end(); cit++) {
	  for (cit2=cit->second.begin(); cit2 != cit->second.end(); cit2++)
		  pointset.push_back(cit2->centroid);

	}

	return pointset;
}

int main( int argc, char ** argv )
{

	/* NOTICE: image should be the colormap processed by the visitor, because it usually relabel the segments.
	* Also it is important to use a timestep which does not create too much bifurcation or
	* use a second pass code that links togheter the segments with only a bifurcation
	*/

  init_logging("Segment_linker");
  pLogType l = get_logger("MAIN");

  try {
      	// Declare the supported options.
      	po::options_description desc("Allowed options");
      	desc.add_options()
      	    ("help,h", "produce help message")
      	    ("simple,s", "disable matching algorithms")
      	    ("reference-image,r", po::value<std::string>(), "specify refrence image")
      	    ("test-image,t", po::value<std::string>(), "specify test image")
      	    ("ref-csv,1", po::value<std::string>(), "read reference segments data from file")
      	    ("test-csv,2", po::value<std::string>(), "read test segments data from file")
      	    ("first-label,l", po::value<int>()->default_value(1), "set the label of the first segment into reference image")
      	    ("save-output,o", po::value<std::string>()->default_value("link_data.txt"), "save output to a csv file")


    	;

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
		po::notify(vm);

		if (vm.count("help") || (argc == 1)) {
		  std::cout << desc << std::endl;
	  	  return EXIT_SUCCESS;
        }

		if (!vm.count("test-image")) {
		  std::cout << "You must specify both the test image" << std::endl;
		  return EXIT_FAILURE;
		}

		if (!vm.count("ref-csv")) {
		  std::cout << "You must specify the reference csv file" << std::endl;
		  return EXIT_FAILURE;
		}

		if (!vm.count("test-csv")) {
		  std::cout << "You must specify the test csv file" << std::endl;
		  return EXIT_FAILURE;
		}


		  typedef itk::ImageFileReader < ImageType >  ReaderType;
		  typedef ImageType::Pointer ImagePointer;
		  std::string sep = ", ";


		  l->info("Loading propagation info..");

		  std::string bs = vm["ref-csv"].as< std::string >();
		  std::map<int, std::vector<PropagationInfo> > ref_data = read_prop_info(bs);

		  bs = vm["test-csv"].as< std::string >(); 
		  std::map<int, std::vector<PropagationInfo> > test_data = read_prop_info(bs);

		  std::vector<PointType> test_points = as_pointset(test_data);

		  l->infoStream() << "Loaded " << ref_data.size() << " segments with a total of " << test_points.size() << " points";

		  l->info("Loading test image...");
		  ReaderType::Pointer reader = ReaderType::New();
		  reader->SetFileName( vm["test-image"].as< std::string >() );
		  try {
			  reader->Update();
		  } catch( itk::ExceptionObject & err ) {
			  l->fatalStream() << "ExceptionObject caught!";
			  std::cerr << err;
			  return EXIT_FAILURE;
		  }

		  ImagePointer test_img = reader->GetOutput();

		  /*
		   * Here we iterate over all the propagation and look at which is the corresponding segment label into the
		   * test image.
		   */



		  std::ofstream outfile (vm["save-output"].as<std::string>().c_str());

		  outfile <<  "from_seg" << sep << "pid" << sep << "c_x" << sep << "c_y" << sep << "c_z" <<  sep << "to_seg" << std::endl;

		  std::map<int, std::vector<PropagationInfo> >::iterator cit;
		  std::vector<PropagationInfo>::iterator cit2;

		  for (cit=ref_data.begin(); cit != ref_data.end(); cit++) {
			  l->infoStream() << "Segment " << cit->first;
			  for (cit2=cit->second.begin(); cit2 != cit->second.end(); cit2++) {

				  ImageType::IndexType idx;idx.Fill(0);
				  bool isInside = test_img->TransformPhysicalPointToIndex(cit2->centroid, idx);
				  if (isInside) {

					  int label = test_img->GetPixel(idx);
					  cit2->to_label=label;

					  outfile <<  cit->first << sep << cit2->id << sep << cit2->centroid[0] << sep << cit2->centroid[1] << sep << cit2->centroid[2] << sep;
					  outfile <<  label << std::endl;
					  if (label == 0) {
						  if ( (!vm.count("simple")) ) {
							  /* No immediate match in the other image and matching algorithms enabled.
							   * Search for the closest propagation
							   */

							  PointType closest = linearNSS(test_points, cit2->centroid);
							  bool isInside = test_img->TransformPhysicalPointToIndex(closest, idx);
							  if (isInside) {
								  label = test_img->GetPixel(idx);
								  l->infoStream() << "Propagation " << cit2->id << " centroid " << cit2->centroid << " is close to " << closest << " -> " << idx << " which belongs to segment: " <<label;
							  } else
								  l->infoStream() << "Propagation " << cit2->id << " centroid " << cit2->centroid << " -> " << idx << " is outside ";

						  }
					  }
					  l->infoStream() << "Propagation " << cit2->id << " centroid " << cit2->centroid << " -> " << idx << " belongs to segment: " <<label;




				  }
			  }
		  }

		  outfile.close();


		  /*
		   * Gather statistics
		   */


		  l->infoStream() << "segid " << sep << "tot_prop" << sep << "link" << sep << "#props" << sep << "prob";
		  for (cit=ref_data.begin(); cit != ref_data.end(); cit++) {
			  std::ostringstream oss;
			  //Writes "segid, tot_prop"
			  oss << cit->first << sep << cit->second.size() << sep;
			  std::map<int, int> counts;
			  std::map<int, int>::const_iterator counts_it;
			  for (cit2=cit->second.begin(); cit2 != cit->second.end(); cit2++) {
				  if (counts.find(cit2->to_label) == counts.end())
					  counts[cit2->to_label] = 1;
				  else
					  ++counts[cit2->to_label];
			  }
			  //Writes "link, #prop, prob" for each possibility
			  for (counts_it=counts.begin(); counts_it!=counts.end(); counts_it++) {
				  oss << counts_it->first << sep << counts_it->second << sep << counts_it->second*1.0/cit->second.size() << sep;
			  }

			  l->infoStream() << oss.str();

		  }


  } catch(std::exception& e) {
  	std::cout << e.what() << std::endl;
  }


  return EXIT_SUCCESS;
}
