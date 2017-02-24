#include "talk2manager.h"

std::string getmodel() {

    std::string res;
            
    // '*#' is the 'ring' bell to grap manager attention
    std::cout << "*#model" << std::endl;
    std::getline(std::cin,res);
    std::cin.sync();

	return res;
}

std::string query_manager(std::string xpath) {

    std::string res;

    // '*#' is the 'ring' bell to grap manager attention
    std::cout << "*#query" << std::endl;
    std::cout << xpath.c_str() << std::endl;
    
    std::getline(std::cin,res);
    std::cin.sync();
	return res;
}


std::string exec_manager(std::string command) {

    std::string res;

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


