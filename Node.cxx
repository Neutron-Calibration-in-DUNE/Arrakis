/**
 * @file Node.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-20
 */
#include "Node.h"

namespace arrakis
{
    Node::Node()
    {
    }
    Node::~Node()
    {
    }

    Node::Node(size_t rank, NodeType node_type)
    : mRank(rank), mNodeType(node_type)
    {
    }

    Node* Node::GetChild(size_t child)
    {
        return mChildren[child];
    }
    std::vector<Node*> Node::GetChildren()
    {
        return mChildren;
    }
}