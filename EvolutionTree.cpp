#include "EvolutionTree.h"

#include <iostream>

Reconstruction::Reconstruction(std::shared_ptr<EvolutionTree> evolution_tree)
    : evolution_tree(evolution_tree)
{
    SaveSubtree(evolution_tree);
}

void Reconstruction::SaveSubtree(std::shared_ptr<EvolutionTree> evolution_tree)
{
    if (evolution_tree == nullptr) {
        return;
    }

    id_to_subree[evolution_tree->root_id] = evolution_tree;

    for (const auto& [_, child] : evolution_tree->children) {
        SaveSubtree(child);
    }
}

std::shared_ptr<Graph> BuildCompleteGraph(int edgeCount)
{
    auto graph = std::make_shared<Graph>();
    for (int i = 0; i < 2 * edgeCount; ++i) {
        graph->AddVertex();
    }
    for (int i = 0; i < 2 * edgeCount; ++i) {
        for (int j = i + 1; j < 2 * edgeCount; ++j) {
            graph->AddEdge(i, j);
        }
    }
    return graph;
}

std::shared_ptr<EvolutionTree> BuildRawTree(
    const std::string& bracket_representation,
    int edgeCount,
    std::shared_ptr<int> pos,
    std::shared_ptr<Graph> completeGraph)
{
    if (pos == nullptr) {
        pos = std::make_unique<int>(0);
    }

    if (completeGraph == nullptr) {
        completeGraph = BuildCompleteGraph(edgeCount);
    }

    if (bracket_representation.empty() ||
        *pos >= bracket_representation.size())
    {
        return nullptr;
    }

    // Extract root ID.
    int root_id = 0;
    while (*pos < bracket_representation.size() &&
        bracket_representation[*pos] != '(' &&
        bracket_representation[*pos] != ')'
    ) {
        root_id = root_id * 10 + (int)(bracket_representation[*pos] - '0');
        *pos = *pos + 1;
    }

    auto root = std::make_shared<EvolutionTree>(root_id);
    root->structure.graph = completeGraph;
    root->structure.costs = std::vector<double>(edgeCount * (2 * edgeCount - 1), 0);

    // Build children.
    while (*pos < bracket_representation.size() &&
        bracket_representation[*pos] == '(')
    {
        *pos = *pos + 1;
        auto child = BuildRawTree(bracket_representation, edgeCount, pos, completeGraph);
        if (child != nullptr) {
            root->children[child->root_id] = child;
            child->parent = root;
        }
    }

    if (*pos < bracket_representation.size() &&
        bracket_representation[*pos] == ')')
    {
        *pos = *pos + 1;
        return root;
    }

    return root;
}

void PrintBracketRepresentation(
    std::ostream& os,
    std::shared_ptr<EvolutionTree> evolution_tree)
{
    if (evolution_tree == nullptr) {
        return;
    }

    os << evolution_tree->root_id;

    for (const auto& [_, child] : evolution_tree->children) {
        os << '(';
        PrintBracketRepresentation(os, child);
        os << ')';
    }
}

void FillStructure(
    const std::string& structure,
    std::shared_ptr<EvolutionTree> evolution_subtree)
{
    std::vector<std::string> tokens;
    std::string delimiter = "|";
    size_t start_pos = 0, end_pos = 0;

    while ((end_pos = structure.find(delimiter, start_pos)) != std::string::npos) {
        tokens.emplace_back(structure.substr(start_pos, end_pos - start_pos));
        start_pos = end_pos + delimiter.length();
    }
    tokens.emplace_back(structure.substr(start_pos, end_pos - start_pos));

    if (tokens.size() <= 1) {
        std::cout << "tokens size is too small" << std::endl;
        throw "tokens size is too small";
    }

    // For cycle structure first_start and last_end points should be united.
    if (tokens.front() == "C" || tokens.size() == 2) {
        int first = std::atoi(tokens[1].c_str());
        bool firstReversed = (tokens[1].back() == '\'');
        int last = std::atoi(tokens.back().c_str());
        bool lastReversed = (tokens.back().back() == '\'');

        auto first_start = 2 * first - 2;
        auto first_end = 2 * first - 1;
        if (firstReversed) {
            std::swap(first_start, first_end);
        }

        auto last_start = 2 * last - 2;
        auto last_end = 2 * last - 1;
        if (lastReversed) {
            std::swap(last_start, last_end);
        }

        auto edge_index = evolution_subtree->structure.graph->GetEdgeIndex(first_start, last_end);
        evolution_subtree->structure.matching.insert(edge_index);
    }

    // Build other edges.
    for (int i = 2; i < tokens.size(); ++i) {
        int first = std::atoi(tokens[i - 1].c_str());
        bool firstReversed = (tokens[i - 1].back() == '\'');
        int last = std::atoi(tokens[i].c_str());
        bool lastReversed = (tokens[i].back() == '\'');

        auto first_start = 2 * first - 2;
        auto first_end = 2 * first - 1;
        if (firstReversed) {
            std::swap(first_start, first_end);
        }

        auto last_start = 2 * last - 2;
        auto last_end = 2 * last - 1;
        if (lastReversed) {
            std::swap(last_start, last_end);
        }

        auto edge_index = evolution_subtree->structure.graph->GetEdgeIndex(first_end, last_start);
        evolution_subtree->structure.matching.insert(edge_index);
    }
}

void PrintStructure(
    std::ostream& stream,
    std::shared_ptr<EvolutionTree> evolution_subtree)
{
    stream << evolution_subtree->root_id << std::endl;
    for (int edge_index : evolution_subtree->structure.matching) {
        auto [u, v] = evolution_subtree->structure.graph->GetEdge(edge_index);
        stream << u << " " << v << std::endl;
    }
}
