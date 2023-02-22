/**
 * @file MCData.h
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
#include "lardataobj/lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "Configuration.h"
#include "Logger.h"

namespace arrakis
{
    namespace mcdata
    {
        class MCData
        {
        public:
            // this singleton cannot be cloned
            MCData(MCData &other) = delete;
            // singleton should also not be assignable
            void operator=(const MCData &) = delete;

            // static method that controls access to 
            // the singleton instance
            static MCData* GetInstance();

            // methods for processing event data
            void ProcessEvent(const Parameters& config, art::Event const& event);
            void ProcessMCTruth(art::Event const& event, art::InputTag input_tag);
            void ProcessMCParticle(art::Event const& event, art::InputTag input_tag);
            void ProcessSimEnergyDeposit(art::Event const& event, art::InputTag input_tag);

        protected:
            MCData() { mLogger = Logger::GetInstance("mcdata"); }
            ~MCData() {}

        private:
            static MCData * sInstance;
            static std::mutex sMutex;

            Logger* mLogger;

            art::Handle<std::vector<simb::MCTruth>>         mMCTruthHandle;
            art::Handle<std::vector<simb::MCParticle>>      mMCParticleHandle;
            art::Handle<std::vector<sim::SimEnergyDeposit>> mMCSimEnergyDepositHandle;
            art::Handle<std::vector<sim::SimChannel>>       mMCSimChannelHandle;
            art::Handle<std::vector<raw::RawDigit>>         mMCRawDigitHandle;

        };
    }
}