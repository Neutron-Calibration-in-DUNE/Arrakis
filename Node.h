/**
 * @file Node.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-24
 */
#pragma once
#include <memory>

#include "NodeData.h"

namespace arrakis
{
    namespace mctree
    {
        class Node : public std::enable_shared_from_this<Node>
        {
        public: 
            explicit Node(NodeData node_data)
            : mNodeData(node_data)
            {
            }
            ~Node()
            {
            }
            explicit Node(NodeData node_data, std::weak_ptr<Node> parent)
            : mNodeData(node_data), mParent(parent)
            {
            }

            void AddChild(std::shared_ptr<Node> child)
            {
                child->mParent = std::weak_ptr<Node>(shared_from_this());
                mChild = child;
            }

            inline NodeData& GetNodeData()  { return mNodeData; }
            std::weak_ptr<Node> Parent()    { return mParent; }
            std::shared_ptr<Node> Child()   { return mChild; }
            std::shared_ptr<Node> Sibling() { return mSibling; }

        private:
            NodeData mNodeData;
            std::weak_ptr<Node> mParent;
            std::shared_ptr<Node> mChild;
            std::shared_ptr<Node> mSibling;
        };
    }
}