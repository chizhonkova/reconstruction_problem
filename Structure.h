#pragma once

#include "matching/Graph.h"

#include <unordered_set>
#include <utility>

struct Structure
{
    // Graph M is the same for each structer.
    std::shared_ptr<Graph> graph;

    std::vector<double> costs;
    std::unordered_set<int> matching;
    std::unordered_set<int> potential_matching;
    std::unordered_set<int> loops;
    double current_cost = 0;
    double potential_cost = 0;
};
