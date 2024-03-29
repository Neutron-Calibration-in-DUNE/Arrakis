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

        /**
         * Class for containing data for various types of nodes.
         * Each node has a type dictated by the enum "NodeType".
         * Other functions must be overriden by the specific 
         * child class.
        */
        class NodeData
        {
        public:
            NodeData(){}
            virtual ~NodeData(){}
            const NodeType& Type()  { return mNodeType; }
            Int_t TypeInt()         { return static_cast<NodeTypeInt>(Type()); }

            virtual const Double_t T() { return 0.0; }

            NodeType mNodeType = NodeType::Empty;
        };

        class PrimaryData : public NodeData
        {
        public:
            PrimaryData(Int_t index)
            : mParticle(index)
            {
            }

            const NodeType& Type()  { return mNodeType; }

            const Double_t T()  { 
                return mcdata::MCData::GetInstance()->GetMCParticle(mParticle).T(); 
            }

        private:
            const Int_t mParticle;
            const NodeType mType = NodeType::Primary;
        };

        class TrajectoryPointData : public NodeData
        {
        public:
            TrajectoryPointData(simb::MCParticle& particle, Int_t point)
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