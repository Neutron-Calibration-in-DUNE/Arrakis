/**
 * @file NodeData.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-24
 */
#pragma once
#include <memory>

#include "nusimdata/SimulationBase/MCParticle.h"

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

            virtual const Double_t T() = 0;

            NodeType mNodeType;
        };

        class Empty : public NodeData
        {
        public:
            const NodeType& Type()  { return mNodeType; }

            const Double_t T()  { return 0.0; }
        private:
            const NodeType mType = NodeType::Empty;
        };

        class Primary : public NodeData
        {
        public:
            Primary(simb::MCParticle& particle)
            : mParticle(particle)
            {
            }

            const NodeType& Type()  { return mNodeType; }

            const Double_t T()  { return mParticle.T(); }

        private:
            const simb::MCParticle& mParticle;
            const NodeType mType = NodeType::Primary;
        };

        class TrajectoryPoint : public NodeData
        {
        public:
            TrajectoryPoint(simb::MCParticle& particle, Int_t point)
            : mParticle(particle), mPoint(point)
            {
            }

            const NodeType& Type()  { return mNodeType; }

            const Double_t T()  { return mParticle.T(mPoint); }

        private:
            const simb::MCParticle& mParticle;
            const Int_t mPoint;
            const NodeType mType = NodeType::TrajectoryPoint;
        };
    }
}