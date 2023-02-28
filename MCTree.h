/**
 * @file MCTree.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-20
 */
#pragma once
#include <vector>
#include <memory>
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
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/sim.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "Configuration.h"
#include "Core.h"
#include "Logger.h"
#include "MCData.h"
#include "Particle.h"

namespace arrakis
{
    namespace mcdata
    {
        class MCTree
        {
        public:
            MCTree();
            ~MCTree();

            void ResetEvent();

            /**
             * @brief Main function for processing event data.
             * The individual primary trees are built first,
             * with all the trajectory points and daughter
             * trees.
             * 
             * @param config 
             * @param event 
             */
            void ProcessEvent(const Parameters& config, art::Event const& event);
            void ProcessMCParticles(const Parameters& config, art::Event const& event);

            Particle CreatePrimary(const simb::MCParticle& particle, Int_t index);
        
        private:
            std::vector<Particle> sPrimaries;
        };
    }
}