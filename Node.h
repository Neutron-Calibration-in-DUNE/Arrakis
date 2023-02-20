/**
 * @file Node.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-20
 */
#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <map>

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

// necessary ROOT libraries
#include <TTree.h>

#include "Generators.h"
#include "Logger.h"

namespace arrakis
{
    enum NodeType
    {
        Empty = -1,
        Primary = 0,
        TrajectoryStep = 1,
        EnergyDeposition = 2,
        DetectorSimulation = 3
    };

    class Node
    {
    public:
        Node();
        ~Node();

        Node(size_t rank, NodeType node_type);

        NodeType Type() { return mType; }
        size_t Rank()   { return mRank; }

        Node* GetParent()   { return mParent; }

        size_t NumberOfChildren()   { return mChildren.size(); }
        Node* GetChild(size_t child);
        std::vector<Node*> GetChildren();

    private:
        NodeType mNodeType = {NodeType::Empty};
        size_t mRank = {0};
        Node* mParent = {0};
        std::vector<Node*> mChildren = {};
    };
}