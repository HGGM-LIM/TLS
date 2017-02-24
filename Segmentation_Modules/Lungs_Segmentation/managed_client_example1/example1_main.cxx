#include <string>
#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>

using namespace std;

string getmodel() {

    string res;
            
    // '*#' is the 'ring' bell to grap manager attention
    std::cout << "*#model" << std::endl;
    std::getline(std::cin,res);
    std::cin.sync();

	return res;
}

string query_manager(string xpath) {

    string res;

    // '*#' is the 'ring' bell to grap manager attention
    std::cout << "*#query" << std::endl;
    std::cout << xpath.c_str() << std::endl;
    
    std::getline(std::cin,res);
    std::cin.sync();
	return res;
}


string exec_manager(string command) {

    string res;

    /* '*#' is the 'ring' bell to grap manager attention
        Notice that you only have access to the model local of the manager.
        Commands that you send will be passed to the os.eval function.
        Please be gentle and don't try to hack the parser :)
    */
    
    std::cout << "*#exec" << std::endl;
    std::cout << command.c_str() << std::endl;
    
    std::getline(std::cin,res);
    std::cin.sync();
	return res;

}

int main ( int argc, char ** argv )
{
    std::cout << "This is an example of a segmentation client"  << std::endl;
    std::cout << "you should run this example from a Segmenter Solidify class..."  << std::endl;
    
    string model1 = getmodel();
    std::cout << "Starting model = " << model1 << std::endl;
    
    string query_image="/SegmentationModel/source_image/@fullPath";
    string fullPath = query_manager(query_image);
    std::cout << "Source image path = " << fullPath << std::endl;        
    
    string query_seedX="/SegmentationModel/Trachea/seed/@x";
    string query_seedY="/SegmentationModel/Trachea/seed/@y";
    string query_seedZ="/SegmentationModel/Trachea/seed/@z";
    int x = atoi( query_manager(query_seedX).c_str() );
    int y = atoi( query_manager(query_seedY).c_str() );    
    int z = atoi( query_manager(query_seedZ).c_str() );    
    std::cout << "Trachea seed is in (x,y,z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
        
    string store="model.add('Segment1')";
    string res1 = exec_manager(store);
    
    string model2 = getmodel() ;
    std::cout << "End model after segment creation is = " << model2 << std::endl;    
    
    
            
	return -1;

	}
