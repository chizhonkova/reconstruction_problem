#include "EvolutionTree.h"

#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
    std::string input_filename, output_filename;
    int i = 1;
	while(i < argc)
	{
		std::string a(argv[i]);
		if(a == "--input"){
			input_filename = argv[++i];
        } else if (a == "--output") {
            output_filename = argv[++i];
        }
		i++;
	}
    
    if (input_filename.empty() || output_filename.empty()) {
        std::cout << "usage: ./reconstruction --input input-filename --output output-filename" << std::endl;
        return -1;
    }

    std::ifstream in(input_filename);
    std::string bracket_representation;
    if (!std::getline(in, bracket_representation)) {
        std::cout << "no bracket sequence is given" << std::endl;
        return -1;
    }

    std::cout << "Building raw tree ..." << std::endl;

    auto tree = BuildRawTree(bracket_representation);

    std::cout << "Saving data in " << output_filename << " ..." << std::endl;

    std::ofstream out(output_filename);
    PrintBracketRepresentation(out, tree);

    std::cout << "Done!" << std::endl;

    return 0;
}
