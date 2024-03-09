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
        throw "no bracket sequence is given";
    }

    int edge_count = 0;
    in >> edge_count;
    std::cout << "Edge count: " << edge_count << std::endl;

    std::cout << "Building raw tree ..." << std::endl;

    auto tree = BuildRawTree(bracket_representation, edge_count);

    Reconstruction reconstruction(tree);

    std::cout << "Building structures ..." << std::endl;

    int leaf_id, structe_count;
    std::string structure;
    while(in >> leaf_id >> structe_count) {
        std::getline(in, structure);
        for (int i = 0; i < structe_count; ++i) {
            if (!std::getline(in, structure)) {
                std::cout << "no structure is given" << std::endl;
                throw "no structure is given";
            }
            FillStructure(structure, reconstruction.id_to_subree[leaf_id]);
        }
    }

    InitInnerStructes(tree, edge_count);

    std::cout << "Saving data in " << output_filename << " ..." << std::endl;

    std::ofstream out(output_filename);
    PrintBracketRepresentation(out, tree);
    out << std::endl;

    for (const auto& [id, subtree] : reconstruction.id_to_subree) {
        PrintStructure(out, subtree);
    }

    std::cout << "Done!" << std::endl;

    return 0;
}
