#include "EvolutionTree.h"

#include "matching/Matching.h"

#include <iostream>
#include <queue>

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
    id_to_subtree[evolution_tree->root_id] = evolution_tree;

    for (const auto& [_, child] : evolution_tree->children) {
        SaveSubtree(child);
    }
}

std::shared_ptr<Graph> BuildCompleteGraph(int edge_count)
{
    auto graph = std::make_shared<Graph>();
    for (int i = 0; i < 2 * edge_count; ++i) {
        graph->AddVertex();
    }
    for (int i = 0; i < 2 * edge_count; ++i) {
        for (int j = i + 1; j < 2 * edge_count; ++j) {
            graph->AddEdge(i, j);
        }
    }
    return graph;
}

std::shared_ptr<EvolutionTree> BuildRawTree(
    const std::string& bracket_representation,
    int edge_count,
    std::shared_ptr<int> pos,
    std::shared_ptr<Graph> complete_graph)
{
    if (pos == nullptr) {
        pos = std::make_shared<int>(0);
    }

    if (complete_graph == nullptr) {
        complete_graph = BuildCompleteGraph(edge_count);
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
    root->structure.graph = complete_graph;

    // Build children.
    while (*pos < bracket_representation.size() &&
        bracket_representation[*pos] == '(')
    {
        *pos = *pos + 1;
        auto child = BuildRawTree(bracket_representation, edge_count, pos, complete_graph);
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

std::vector<std::string> Reconstruction::PrepareTokens(
    std::shared_ptr<EvolutionTree> evolution_subtree, 
    std::vector<std::string> tokens)
{
    for (int i = 1; i < tokens.size(); ++i) {
        std::string token = tokens[i];
        bool reversed = false;
        if (token.back() == '\'') {
            reversed = true;
            token = {token.begin(), token.end() - 1};
        }
        if (!edge_name_to_id.contains(token)) {
            auto id = next_id++;
            edge_name_to_id[token] = id;
            edge_id_to_name[id] = token;
        }
        if (evolution_subtree == nullptr) {
            std::cout << "wtf" << std::endl;
        }
        evolution_subtree->structure.used_edges.insert(token);
        all_edges.insert(token);
        tokens[i] = std::to_string(edge_name_to_id[token]);
        if (reversed) {
            tokens[i] += '\'';
        }
    }
    return tokens;
}

void Reconstruction::FillStructure(
    const std::string& structure,
    std::shared_ptr<EvolutionTree> evolution_subtree, 
    bool without_erase)
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

    tokens = PrepareTokens(evolution_subtree, std::move(tokens));

    // For cycle structure first_start and last_end points should be united.
    if (tokens.front() == "C") {
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

void Reconstruction::FillLoopStructures(std::shared_ptr<EvolutionTree> evolution_subtree)
{
    auto dif = all_edges;
    for (const auto& edge : evolution_subtree->structure.used_edges) {
        dif.erase(edge);
    }
    for (const auto& edge : dif) {
        std::string structure = "C|" + edge;
        FillStructure(structure, evolution_subtree, true);
    }
}

int GetEnd(int start) {
    int edge = start / 2 + 1;
    if (start % 2 == 0) {
        return 2 * edge - 1;
    }
    return 2 * edge - 2;
}

void Reconstruction::AddNextEdge(
    int start,
    std::string& repr)
{
    int end = GetEnd(start);
    int edge = start / 2 + 1;

    repr += "|" + edge_id_to_name[edge];
    if (start > end) {
        repr += "'";
    }
}

void Reconstruction::AddPrevEdge(
    int end,
    std::string& repr)
{
    int start = GetEnd(end);
    int edge = start / 2 + 1;
    std::string token = "|" + edge_id_to_name[edge];
    if (start > end) {
        token += "'";
    }
    repr = token + repr;
}

void AddType(const Graph& graph, int start, int end, std::string& repr)
{
    if (graph.AdjMat()[start][end]) {
        repr = "C" + repr;
    } else {
        repr = "L" + repr;
    }
}

std::string Reconstruction::NextStructure(
    const Graph& graph,
    std::unordered_set<int>& not_visited)
{
    std::string repr;

    std::queue<int> queue;
    int u = *not_visited.begin();
    int v = GetEnd(u);
    if (u > v) {
        std::swap(u, v);
    }
    queue.push(u);

    // Move forward.
    while (!queue.empty()) {
        int start = queue.front();
        queue.pop();
        AddNextEdge(start, repr);
        not_visited.erase(start);

        int end = GetEnd(start);
        v = end;
        not_visited.erase(end);

        const auto& adjList = graph.AdjList(end);
        if (adjList.size() > 1) {
            std::cout << "bad adjList size" << std::endl;
            throw "bad adjList size";
        }
        for (int next : adjList) {
            if (not_visited.contains(next)) {
                queue.push(next);
            }
        }
    }

    // Move backward.
    const auto& adjList = graph.AdjList(u);
    if (adjList.size() == 0) {
        AddType(graph, u, v, repr);
        return repr;
    }

    int prev = adjList.front();
    if (!not_visited.contains(prev)) {
        AddType(graph, u, v, repr);
        return repr;
    }

    queue.push(prev);

    while(!queue.empty()) {
        int start = queue.front();
        queue.pop();
        AddPrevEdge(start, repr);
        not_visited.erase(start);

        int end = GetEnd(start);
        u = end;
        not_visited.erase(end);

        const auto& adjList = graph.AdjList(end);
        if (adjList.size() > 1) {
            std::cout << "bad adjList size" << std::endl;
            throw "bad adjList size";
        }
        for (int next : adjList) {
            if (not_visited.contains(next)) {
                queue.push(next);
            }
        }
    }

    AddType(graph, u, v, repr);
    return repr;
}

void Reconstruction::PrintStructure(
    std::ostream& stream,
    std::shared_ptr<EvolutionTree> evolution_subtree)
{
    stream << evolution_subtree->root_id << std::endl;

    Graph graph;
    for (int i = 0; i < evolution_subtree->structure.graph->GetNumVertices(); ++i) {
        graph.AddVertex();
    }
    for (int edge_index : evolution_subtree->structure.matching) {
        auto [u, v] = evolution_subtree->structure.graph->GetEdge(edge_index);
        graph.AddEdge(u, v);
    }

    std::unordered_set<int> not_visited;
    for (int i = 0; i < graph.GetNumVertices(); ++i) {
        not_visited.insert(i);
    }

    stream << "Structure:" << std::endl;
    while (!not_visited.empty()) {
        auto repr = NextStructure(graph, not_visited);
        stream << repr << std::endl;
    }

    stream << std::endl;
}

double Reconstruction::CalculateCostFromChild(
    const std::unordered_set<int>& node_matching,
    const std::unordered_set<int>& child_matching)
{
    std::unordered_set<int> matching_union = node_matching;
    matching_union.insert(child_matching.begin(), child_matching.end());

    double cost = 0;

    // Here we don't need all edges, just union of both matchings.
    for (int edge_index : matching_union) {
        if (node_matching.contains(edge_index) &&
            !child_matching.contains(edge_index))
        {
            if (evolution_tree->structure.graph->IsLoop(edge_index)) {
                cost += insertion_cost;
            } else {
                cost += cut_cost;
            }
        } else if (!node_matching.contains(edge_index) &&
            child_matching.contains(edge_index))
        {
            if (evolution_tree->structure.graph->IsLoop(edge_index)) {
                cost += deletion_cost;
            } else {
                cost += join_cost;
            }
        }
    }

    return cost;
}

void Reconstruction::CalculateInitialCost(std::shared_ptr<EvolutionTree> evolution_tree)
{
    if (evolution_tree == nullptr || evolution_tree->children.empty()) {
        return;
    }

    for (const auto& [_, child] : evolution_tree->children) {
        CalculateInitialCost(child);
    }

    CalculateCost(evolution_tree);
}

void Reconstruction::CalculateCost(std::shared_ptr<EvolutionTree> subtree)
{
    if (subtree->children.empty()) {
        return;
    }

    subtree->structure.current_cost = 0;
    for (const auto& [_, child] : subtree->children) {
        subtree->structure.current_cost += CalculateCostFromChild(
            subtree->structure.matching,
            child->structure.matching);
    }

    if (subtree->parent != nullptr) {
        subtree->structure.current_cost += CalculateCostFromChild(
            subtree->parent->structure.matching,
            subtree->structure.matching);
    }
}

void Reconstruction::CalculatePotentialMatching(std::shared_ptr<EvolutionTree> subtree)
{
    if (subtree->children.empty()) {
        return;
    }

    bool noNegative = true;
    auto weights = CalculateWeights(subtree);
    for (double weight : weights) {
        if (weight < 0) {
            noNegative = false;
        }
    }

    subtree->structure.potential_matching = {};
    if (!noNegative) {
        Matching m(*subtree->structure.graph);
        std::pair<std::list<int>, double> solution = m.SolveMinimumCostPerfectMatching(weights);

        for (int i : solution.first) {
            // subtree->structure.potential_matching.insert(i);
            if (weights[i] < 0) {
                subtree->structure.potential_matching.insert(i);
            }
        }
    }

    subtree->structure.potential_cost = 0;
    for (const auto& [_, child] : subtree->children) {
        subtree->structure.potential_cost += CalculateCostFromChild(
            subtree->structure.potential_matching,
            child->structure.matching);
    }

    if (subtree->parent != nullptr) {
        subtree->structure.potential_cost += CalculateCostFromChild(
            subtree->parent->structure.matching,
            subtree->structure.potential_matching);
    }
}

void Reconstruction::CalculatePotentialMatchings()
{
    for (const auto& [_, subtree] : id_to_subtree) {
        CalculatePotentialMatching(subtree);

        heap.Insert(
            subtree->structure.potential_cost - subtree->structure.current_cost,
            subtree->root_id);
    }
}

double Reconstruction::GetPYesFromChild(int edge_index, std::shared_ptr<EvolutionTree> child)
{
    if (child->structure.matching.contains(edge_index)) {
        return 0;
    } else {
        if (child->structure.graph->IsLoop(edge_index)) {
            return insertion_cost;
        }
        return cut_cost;
    }
}

double Reconstruction::GetPNoFromChild(int edge_index, std::shared_ptr<EvolutionTree> child)
{
    if (child->structure.matching.contains(edge_index)) {
        if (child->structure.graph->IsLoop(edge_index)) {
            return deletion_cost;
        }
        return join_cost;
    } else {
        return 0;
    }
}

double Reconstruction::GetPYesFromParent(int edge_index, std::shared_ptr<EvolutionTree> parent)
{
    if (parent->structure.matching.contains(edge_index)) {
        return 0;
    } else {
        if (parent->structure.graph->IsLoop(edge_index)) {
            return deletion_cost;
        }
        return join_cost;
    }
}

double Reconstruction::GetPNoFromParent(int edge_index, std::shared_ptr<EvolutionTree> parent)
{
    if (parent->structure.matching.contains(edge_index)) {
        if (parent->structure.graph->IsLoop(edge_index)) {
            return insertion_cost;
        }
        return cut_cost;
    } else {
        return 0;
    }
}

std::vector<double> Reconstruction::CalculateWeights(std::shared_ptr<EvolutionTree> evolution_tree)
{
    auto edge_count = evolution_tree->structure.graph->GetNumEdges();
    std::vector<double> weights(edge_count, 0);

    for (int i = 0; i < edge_count; ++i) {
        double p_yes = 0;
        double p_no = 0;

        for (const auto& [_, child] : evolution_tree->children) {
            p_yes += GetPYesFromChild(i, child);
            p_no += GetPNoFromChild(i, child);
        }

        if (evolution_tree->parent != nullptr) {
            p_yes += GetPYesFromParent(i, evolution_tree->parent);
            p_no += GetPNoFromParent(i, evolution_tree->parent);
        }

        weights[i] = p_yes - p_no;
        if (weights[i] >= 0) {
            weights[i] = 0;
        }
    }

    return weights;
}

void Reconstruction::Solve()
{
    CalculateInitialCost(evolution_tree);
    CalculatePotentialMatchings();

    for (int i = 0; i < 1000; ++i) {
        if (heap.Size() == 0) {
            break;
        }

        auto [dif, id] = heap.GetMin();
        if (dif == 0) {
            break;
        }

        auto subtree = id_to_subtree[id];

        // Change matching.
        subtree->structure.matching = subtree->structure.potential_matching;
        subtree->structure.current_cost = subtree->structure.potential_cost;
        heap.ChangeKey(0, id);

        // Recalculate structure for neighbours.
        for (const auto& [_, child] : subtree->children) {
            if (child->children.empty()) {
                continue;
            }
            CalculateCost(child);
            CalculatePotentialMatching(child);
            heap.ChangeKey(
                child->structure.potential_cost - child->structure.current_cost,
                child->root_id);
        }

        if (subtree->parent != nullptr) {
            CalculateCost(subtree->parent);
            CalculatePotentialMatching(subtree->parent);
            heap.ChangeKey(
                subtree->parent->structure.potential_cost - subtree->parent->structure.current_cost,
                subtree->parent->root_id);
        }
    }

    auto [dif, id] = heap.GetMin();
    if (dif != 0) {
        std::cout << "There are some nodes with non-null potential cost" << std::endl;
    }
}

double Reconstruction::CalculateFinalCost(std::shared_ptr<EvolutionTree> tree)
{
    if (evolution_tree == nullptr || evolution_tree->children.empty()) {
        return 0;
    }

    double cost = 0;

    for (const auto& [_, child] : tree->children) {
        cost += CalculateFinalCost(child);

        cost += CalculateCostFromChild(
            tree->structure.matching,
            child->structure.matching);
    }

    return cost;
}
