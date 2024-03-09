#pragma once

#include "Structure.h"

#include <ostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

class EvolutionTree
{
public:
    EvolutionTree(int root_id)
        : root_id(root_id)
    { }

    int root_id;
    std::shared_ptr<EvolutionTree> parent;
    std::unordered_map<int, std::shared_ptr<EvolutionTree>> children;
    Structure structure;
};

struct Reconstruction
{
    Reconstruction(std::shared_ptr<EvolutionTree> evolution_tree);

    std::shared_ptr<EvolutionTree> evolution_tree;
    std::unordered_map<int, std::shared_ptr<EvolutionTree>> id_to_subree;

private:
    void SaveSubtree(std::shared_ptr<EvolutionTree> evolution_tree);
};

std::shared_ptr<EvolutionTree> BuildRawTree(
    const std::string& bracket_representation,
    int edge_count,
    std::shared_ptr<int> pos = nullptr,
    std::shared_ptr<Graph> complete_graph = nullptr);

void FillStructure(
    const std::string& structure,
    std::shared_ptr<EvolutionTree> evolution_subtree);

void InitInnerStructes(
    std::shared_ptr<EvolutionTree> evolution_tree,
    int edge_count);

void PrintBracketRepresentation(
    std::ostream& stream,
    std::shared_ptr<EvolutionTree> evolution_tree);

void PrintStructure(
    std::ostream& stream,
    std::shared_ptr<EvolutionTree> evolution_subtree);
