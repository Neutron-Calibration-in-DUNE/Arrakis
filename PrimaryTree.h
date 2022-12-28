/**
 * @file PrimaryTreemc_part.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12/21
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

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "ParticleMaps.h"

namespace arrakis
{
    enum class NodeType =
    {
        "none" = 0,
        "mc_particle" = 1,
        "edep" = 2,
        "det_sim" = 3,
        "reco" = 4
    };

    class Node
    {
    public:
        Node();
        ~Node();

        enum NodeType GetNodeType() { return mNodeType; }
        
    private:
        enum NodeType mNodeType;
    };

    

    class PrimaryTree
    {

    };
}