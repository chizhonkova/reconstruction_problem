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

    void FillStructure(
        const std::string& structure,
        std::shared_ptr<EvolutionTree> evolution_subtree, 
        bool without_erase = false);
    void FillLoopStructures(std::shared_ptr<EvolutionTree> evolution_subtree);

    void CalculateInitialCost(std::shared_ptr<EvolutionTree> evolution_tree);
    void CalculatePotentialMatchings();
    void Solve();
    double CalculateFinalCost(std::shared_ptr<EvolutionTree> tree);

    void PrintStructure(
        std::ostream& stream,
        std::shared_ptr<EvolutionTree> evolution_subtree);

    std::shared_ptr<EvolutionTree> evolution_tree;
    std::unordered_map<int, std::shared_ptr<EvolutionTree>> id_to_subtree;
    BinaryHeap heap;
    double cut_cost;
    double join_cost;
    double insertion_cost;
    double deletion_cost;

    int next_id = 1;
    std::unordered_map<std::string, int> edge_name_to_id;
    std::unordered_map<int, std::string> edge_id_to_name;
    std::unordered_set<std::string> all_edges;

private:
    std::vector<std::string> PrepareTokens(
        std::shared_ptr<EvolutionTree> evolution_subtree,
        std::vector<std::string> tokens);

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

    void AddPrevEdge(
        int end,
        std::string& repr);
    void AddNextEdge(
        int start,
        std::string& repr);
    std::string NextStructure(
        const Graph& graph,
        std::unordered_set<int>& not_visited);
};

std::shared_ptr<EvolutionTree> BuildRawTree(
    const std::string& bracket_representation,
    int edge_count,
    std::shared_ptr<int> pos = nullptr,
    std::shared_ptr<Graph> complete_graph = nullptr);

void PrintBracketRepresentation(
    std::ostream& stream,
    std::shared_ptr<EvolutionTree> evolution_tree);
