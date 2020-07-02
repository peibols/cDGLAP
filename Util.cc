#include "Util.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace Util {

int IsFile(std::string file_name)
{
 FILE *temp;

 if( (temp = fopen(file_name.c_str(),"r")) == NULL) return 0;
 else
  {
   fclose(temp);
   return 1;
  }
}/* IsFile */

// support comments in the parameters file
// comments need to start with #
std::string StringFind4(std::string file_name, std::string str_in) {
    std::string inputname = file_name;
    std::string str = str_in;

    std::string tmpfilename;
    tmpfilename = "input.default";
    
    // check whether the input parameter file is exist or not
    if (!IsFile(file_name)) {
        if (file_name == "") {
            fprintf(stderr, "No input file name specified.\n");
            fprintf(stderr, "Creating a default file named input.default\n");
        } else {
            std::cerr << "The file named " << file_name << " is absent." << std::endl;
            std::cout << "Creating " << file_name << "..." << std::endl;
            tmpfilename = file_name;
        }
        std::ofstream tmp_file(tmpfilename.c_str());
        tmp_file << "EndOfData" << std::endl;
        tmp_file.close();
        exit(1);
    }/* if isfile */
  
    // pass checking, now read in the parameter file
    std::string temp_string;
    std::ifstream input(inputname.c_str());
    getline(input, temp_string);  // read in the first entry

    int ind = 0;
    std::string para_name;
    std::string para_val;
    while (temp_string.compare("EndOfData") != 0) {
        // check whether it is the end of the file
        std::string para_string;
        std::stringstream temp_ss(temp_string);
        getline(temp_ss, para_string, '#');  // remove the comments
        if (para_string.compare("") != 0
                && para_string.find_first_not_of(' ') != std::string::npos) {
            // check the read in string is not empty
            std::stringstream para_stream(para_string);
            para_stream >> para_name >> para_val;
            if (para_name.compare(str) == 0) {
                // find the desired parameter
                ind++;
                input.close();
                return(para_val);
            }  /* if right, return */
        }
        getline(input, temp_string);  // read in the next entry
    }/* while */
    input.close(); // finish read in and close the file
    
    // the desired parameter is not in the parameter file, then return "empty"
    if (ind == 0) {
        return("empty");
    }
    // should not cross here !!!
    std::cout << "Error in StringFind4 !!!\n";
    return("empty");
}/* StringFind4 */

} // end namespace Util
