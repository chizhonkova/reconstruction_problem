#include "EvolutionTree.h"

std::shared_ptr<EvolutionTree> BuildRawTree(
    const std::string& bracket_representation,
    std::shared_ptr<int> pos)
{
    if (pos == nullptr) {
        pos = std::make_unique<int>(0);
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

    auto root = std::make_unique<EvolutionTree>(root_id);

    // Build children.
    while (*pos < bracket_representation.size() &&
        bracket_representation[*pos] == '(')
    {
        *pos = *pos + 1;
        auto child = BuildRawTree(bracket_representation, pos);
        if (child != nullptr) {
            root->children[child->root_id] = child;
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

void PrintBracketRepresentation(std::ostream& os, std::shared_ptr<EvolutionTree> evolution_tree)
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

void FillStructures(std::shared_ptr<EvolutionTree> evolution_tree)
{
    
}
