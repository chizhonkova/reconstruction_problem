#include "EvolutionTree.h"

#include <algorithm>
#include <filesystem>
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

    double cut_cost, join_cost, insertion_cost, deletion_cost;
    in >> cut_cost >> join_cost >> insertion_cost >> deletion_cost;

    std::cout << "Cut cost: " << cut_cost << std::endl;
    std::cout << "Join cost: " << join_cost << std::endl;
    std::cout << "Insertion cost: " << insertion_cost << std::endl;
    std::cout << "Deletion cost: " << deletion_cost << std::endl;

    std::string bracket_representation;
    // Skip line.
    std::getline(in, bracket_representation);
    if (!std::getline(in, bracket_representation)) {
        throw "no bracket sequence is given";
    }

    int edge_count = 0;
    in >> edge_count;
    std::cout << "Edge count: " << edge_count << std::endl;

    std::cout << "Building raw tree ..." << std::endl;

    auto tree = BuildRawTree(bracket_representation, edge_count);

    Reconstruction reconstruction(tree);
    reconstruction.cut_cost = cut_cost;
    reconstruction.join_cost = join_cost;
    reconstruction.insertion_cost = insertion_cost;
    reconstruction.deletion_cost = deletion_cost;

    std::cout << "Building structures ..." << std::endl;

    std::vector<int> leaves;

    int leaf_id, structe_count;
    std::string structure;
    while(in >> leaf_id >> structe_count) {
        leaves.emplace_back(leaf_id);
        // Skip line.
        std::getline(in, structure);
        for (int i = 0; i < structe_count; ++i) {
            if (!std::getline(in, structure)) {
                std::cout << "no structure is given" << std::endl;
                throw "no structure is given";
            }
            reconstruction.FillStructure(structure, reconstruction.id_to_subtree[leaf_id]);
        }
    }

    for (int leaf_id : leaves) {
        reconstruction.FillLoopStructures(reconstruction.id_to_subtree[leaf_id]);
    }

    std::cout << "Solving problem ..." << std::endl;

    reconstruction.Solve();
    double cost = reconstruction.CalculateFinalCost(tree);

    std::cout << "Saving data in " << output_filename << " ..." << std::endl;

    std::ofstream out(output_filename);
    PrintBracketRepresentation(out, tree);
    out << std::endl;

    std::vector<int> ids;
    ids.reserve(reconstruction.id_to_subtree.size());
    for (const auto& [id, _] : reconstruction.id_to_subtree) {
        ids.push_back(id);
    }
    std::sort(ids.begin(), ids.end());

    for (int id : ids) {
        reconstruction.PrintStructure(out, reconstruction.id_to_subtree[id]);
    }

    out << "Final cost: " << cost << std::endl;
    std::cout << "Final cost: " << cost << std::endl;
    std::cout << "Done!" << std::endl;

    out.flush();
    out.close();
    
    std::cout << "Output file content: " << std::endl;
    std::ifstream f(output_filename);
    if (f.is_open()) {
        std::cout << f.rdbuf() << std::endl;
    }

    std::cout << "Output absolute path: " << std::filesystem::absolute(output_filename) << std::endl;

    return 0;
}
