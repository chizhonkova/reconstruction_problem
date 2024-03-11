#pragma once

#include "Structure.h"
#include "matching/BinaryHeap.h"

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

    void CalculateInitialCost(std::shared_ptr<EvolutionTree> evolution_tree);
    void CalculatePotentialMatchings();
    void Solve();
    double CalculateFinalCost(std::shared_ptr<EvolutionTree> tree);

    std::shared_ptr<EvolutionTree> evolution_tree;
    std::unordered_map<int, std::shared_ptr<EvolutionTree>> id_to_subtree;
    BinaryHeap heap;
    double cut_cost;
    double join_cost;
    double insertion_cost;
    double deletion_cost;

private:
    void SaveSubtree(std::shared_ptr<EvolutionTree> evolution_tree);
    std::vector<double> CalculateWeights(std::shared_ptr<EvolutionTree> evolution_tree);

    double GetPYesFromChild(int edge_index, std::shared_ptr<EvolutionTree> child);
    double GetPNoFromChild(int edge_index, std::shared_ptr<EvolutionTree> child);
    double GetPYesFromParent(int edge_index, std::shared_ptr<EvolutionTree> parent);
    double GetPNoFromParent(int edge_index, std::shared_ptr<EvolutionTree> parent);

    void CalculatePotentialMatching(std::shared_ptr<EvolutionTree> subtree);
    void CalculateCost(std::shared_ptr<EvolutionTree> subtree);

    double CalculateCostFromChild(
        const std::unordered_set<int>& node_matching,
        const std::unordered_set<int>& child_matching);
};

std::shared_ptr<EvolutionTree> BuildRawTree(
    const std::string& bracket_representation,
    int edge_count,
    std::shared_ptr<int> pos = nullptr,
    std::shared_ptr<Graph> complete_graph = nullptr);

void FillStructure(
    const std::string& structure,
    std::shared_ptr<EvolutionTree> evolution_subtree);

void PrintBracketRepresentation(
    std::ostream& stream,
    std::shared_ptr<EvolutionTree> evolution_tree);

void PrintStructure(
    std::ostream& stream,
    std::shared_ptr<EvolutionTree> evolution_subtree);
