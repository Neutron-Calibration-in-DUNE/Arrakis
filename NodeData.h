/**
 * @file NodeData.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-24
 */
#pragma once
#include <memory>

namespace arrakis
{
    namespace mctree
    {
        enum class NodeType
        {
            Empty = -1,
            Primary = 0,
            Daughter = 1,
            TrajectoryPoint = 2,
            EnergyDeposition = 3
        };
        using NodeTypeInt = std::underlying_type<NodeType>::type;

        class NodeData
        {
        public:
            NodeData(){}
            virtual ~NodeData(){}
            const NodeType& Type()  { return mNodeType; }
            Int_t TypeInt()         { return static_cast<NodeTypeInt>(Type()); }

            NodeType mNodeType;
        };

        class Empty : public NodeData
        {
        public:
            const NodeType& Type()  { return mNodeType; }
        private:
            const NodeType mType = NodeType::Empty;
        };
    }
}