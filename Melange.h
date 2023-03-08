/**
 * @file Melange.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#pragma once
#include <mutex>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "Configuration.h"
#include "Core.h"
#include "DetectorSimulation.h"
#include "Logger.h"

namespace arrakis
{
    namespace melange
    {
        class Melange
        {
        public:
            Melange();
            ~Melange();

            // methods for processing event data
            void ResetEvent();
            void ProcessEvent(const Parameters& config, art::Event const& event);

        private:
            // Output TTree
            art::ServiceHandle<art::TFileService> sTFileService;
            TTree *sMCDataTree;
        };
    }
}