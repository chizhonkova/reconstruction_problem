#pragma once

#include "matching/Graph.h"

#include <unordered_set>
#include <utility>

class Structure
{
public:
    // Graph M is the same for each structer.
    std::shared_ptr<Graph> graph;
    std::vector<double> costs;
    std::unordered_set<int> matching;
private:
};
