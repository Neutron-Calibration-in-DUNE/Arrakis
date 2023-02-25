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
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "Configuration.h"
#include "Core.h"
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
            void ProcessMCTruth(art::Event const& event, fhicl::ParameterSet const& generator_labels);
            void ProcessMCParticles(art::Event const& event, art::InputTag input_tag);
            void ProcessSimEnergyDeposits(art::Event const& event, art::InputTag input_tag);
            void ProcessSimChannels(art::Event const& event,
                art::InputTag producer_label, art::InputTag instance_label
            );
            void ProcessRawDigits(art::Event const& event,
                art::InputTag producer_label, art::InputTag instance_label
            );

            art::Handle<std::vector<simb::MCTruth>> GetMCTruth()        { return sMCTruthHandle; }
            art::Handle<std::vector<simb::MCParticle>> GetMCParticles() { return sMCParticleHandle; }
            art::Handle<std::vector<sim::SimEnergyDeposit>> GetSimEnergyDeposits() { return sMCSimEnergyDepositHandle; }
            art::Handle<std::vector<sim::SimChannel>> GetSimChannels()  { return sMCSimChannelHandle; }
            art::Handle<std::vector<raw::RawDigit>> GetRawDigits()      { return sMCRawDigitHandle; }

            const simb::MCParticle& GetMCParticle(Int_t index)  { return (*sMCParticleHandle)[index]; }
            const sim::SimEnergyDeposit& GetMCSimEnergyDeposit(Int_t index) { return (*sMCSimEnergyDepositHandle)[index]; }
            const sim::SimChannel& GetMCSimChannel(Int_t index) { return (*sMCSimChannelHandle)[index]; }
            const raw::RawDigit& GetMCRawDigit(Int_t index)     { return (*sMCRawDigitHandle)[index]; }

            // particle maps from track id
            inline GeneratorLabel GetGeneratorLabel(Int_t trackID) { return sGeneratorLabelMap[trackID]; }
            inline Int_t GetPDGCode(Int_t trackID)          { return sPDGMap[trackID]; }
            inline Int_t GetParentPDG(Int_t trackID)        { return sParentPDGMap[trackID]; }
            inline Int_t GetParentTrackID(Int_t trackID)    { return sParentTrackIDMap[trackID]; }
            inline Int_t GetParticleEnergy(Int_t trackID)   { return sParticleEnergyMap[trackID];}
            inline Int_t GetAncestorPDG(Int_t trackID)      { return sAncestorPDGMap[trackID]; }
            inline Int_t GetAncestorTrackID(Int_t trackID)  { return sAncestorTrackIDMap[trackID]; }
            inline Int_t GetAncestorLevel(Int_t trackID)    { return sAncestorLevelMap[trackID]; }
            inline Int_t GetAncestorEnergy(Int_t trackID)   { return sAncestorEnergyMap[trackID]; }

        protected:
            MCData() {}
            ~MCData() {}

        private:
            static MCData * sInstance;
            static std::mutex sMutex;

            // handles
            std::vector<art::Handle<std::vector<simb::MCTruth>>> sMCTruthHandles;
            art::Handle<std::vector<simb::MCTruth>>         sMCTruthHandle;
            art::Handle<std::vector<simb::MCParticle>>      sMCParticleHandle;
            art::Handle<std::vector<sim::SimEnergyDeposit>> sMCSimEnergyDepositHandle;
            art::Handle<std::vector<sim::SimChannel>>       sMCSimChannelHandle;
            art::Handle<std::vector<raw::RawDigit>>         sMCRawDigitHandle;

            std::map<Int_t, GeneratorLabel> sGeneratorLabelMap;
            std::map<art::InputTag, GeneratorLabel> sGeneratorMap;
            std::map<Int_t, Int_t>      sPDGMap;
            std::map<Int_t, Int_t>      sParentPDGMap;
            std::map<Int_t, Int_t>      sParentTrackIDMap;
            std::map<Int_t, Double_t>   sParticleEnergyMap;
            std::map<Int_t, Int_t>      sAncestorPDGMap;
            std::map<Int_t, Int_t>      sAncestorTrackIDMap;
            std::map<Int_t, Int_t>      sAncestorLevelMap;
            std::map<Int_t, Double_t>   sAncestorEnergyMap;
        };
    }
}