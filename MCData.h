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
#include "DetectorSimulation.h"
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
            void ResetEvent();
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

            const simb::MCParticle& GetMCParticleTrackID(Int_t TrackID)  { return (*sMCParticleHandle)[sParticleMap[TrackID]]; }

            // particle maps from track id
            inline GeneratorLabel GetGeneratorLabel(Int_t trackID) { return sGeneratorLabelMap[trackID]; }
            inline Int_t GetParticleIndex(Int_t trackID)    { return sParticleMap[trackID]; }
            inline Int_t GetPDGCode(Int_t trackID)          { return sPDGMap[trackID]; }
            inline Int_t GetParentPDG(Int_t trackID)        { return sParentPDGMap[trackID]; }
            inline Int_t GetParentTrackID(Int_t trackID)    { return sParentTrackIDMap[trackID]; }
            inline Int_t GetParticleEnergy(Int_t trackID)   { return sParticleEnergyMap[trackID];}
            inline Int_t GetAncestorPDG(Int_t trackID)      { return sAncestorPDGMap[trackID]; }
            inline Int_t GetAncestorTrackID(Int_t trackID)  { return sAncestorTrackIDMap[trackID]; }
            inline Int_t GetAncestorLevel(Int_t trackID)    { return sAncestorLevelMap[trackID]; }
            inline Int_t GetAncestorEnergy(Int_t trackID)   { return sAncestorEnergyMap[trackID]; }
            inline std::vector<Int_t> GetDaughters(Int_t trackID)   { return sDaughterMap[trackID]; }
            inline std::vector<Int_t> GetProgeny(Int_t trackID)     { return sProgenyMap[trackID]; }
            inline std::vector<Int_t> GetAncestry(Int_t trackID)    { return sAncestryMap[trackID]; }
            inline std::vector<Int_t> GetParticleEdep(Int_t trackID){ return sParticleEdepMap[trackID]; }
            inline std::vector<ProcessType> GetParticleEdepProcess(Int_t trackID)  { return sParticleEdepProcessMap[trackID]; }
            inline std::vector<Int_t> GetParticleDetectorSimulation(Int_t trackID) { return sParticleDetectorSimulationMap[trackID]; }
            inline std::vector<Int_t> GetRandomDetectorSimulation(Int_t trackID)   { return sRandomDetectorSimulationMap[trackID]; }

            // maps from edep to process
            inline ProcessType GetEdepProcess(Int_t edepID)   { return sEdepProcessMap[edepID]; }

            // helper functions for organizing data
            ProcessType DetermineEdepProcess(const sim::SimEnergyDeposit& edep);
            std::vector<Int_t> DetermineDetectorSimulationEdeps(const std::vector<sim::IDE>& det_ide, Int_t detsim_id);



            // fill TTree
            void FillTTree();

        protected:
            MCData();
            ~MCData() {}

        private:
            static MCData * sInstance;
            static std::mutex sMutex;

            // Output TTree
            art::ServiceHandle<art::TFileService> sTFileService;
            TTree *sMCDataTree;

            // handles
            std::vector<art::Handle<std::vector<simb::MCTruth>>> sMCTruthHandles;
            art::Handle<std::vector<simb::MCTruth>>         sMCTruthHandle;
            art::Handle<std::vector<simb::MCParticle>>      sMCParticleHandle;
            art::Handle<std::vector<sim::SimEnergyDeposit>> sMCSimEnergyDepositHandle;
            art::Handle<std::vector<sim::SimChannel>>       sMCSimChannelHandle;
            art::Handle<std::vector<raw::RawDigit>>         sMCRawDigitHandle;

            std::map<Int_t, GeneratorLabel> sGeneratorLabelMap;
            std::map<art::InputTag, GeneratorLabel> sGeneratorMap;

            // List of primary track IDs
            std::vector<Int_t> sPrimaries;

            // List of detector simulation structs
            std::vector<DetectorSimulation> sDetectorSimulation;
            DetectorSimulationNoise sDetectorSimulationNoise;

            // MCParticle TrackID maps
            std::map<Int_t, Int_t>      sParticleMap;
            std::map<Int_t, Int_t>      sPDGMap;
            std::map<Int_t, Int_t>      sParentPDGMap;
            std::map<Int_t, Int_t>      sParentTrackIDMap;
            std::map<Int_t, Double_t>   sParticleEnergyMap;
            std::map<Int_t, Int_t>      sAncestorPDGMap;
            std::map<Int_t, Int_t>      sAncestorTrackIDMap;
            std::map<Int_t, Int_t>      sAncestorLevelMap;
            std::map<Int_t, Double_t>   sAncestorEnergyMap;
            std::map<Int_t, std::vector<Int_t>> sDaughterMap;
            std::map<Int_t, std::vector<Int_t>> sProgenyMap;
            std::map<Int_t, std::vector<Int_t>> sAncestryMap;
            std::map<Int_t, std::vector<Int_t>> sParticleEdepMap;
            std::map<Int_t, std::vector<ProcessType>> sParticleEdepProcessMap;
            std::map<Int_t, std::vector<Int_t>> sParticleDetectorSimulationMap;
            std::map<Int_t, std::vector<Int_t>> sRandomDetectorSimulationMap;

            // maps from edepID
            std::map<Int_t, ProcessType> sEdepProcessMap;
            std::map<Int_t, std::vector<Int_t>> sEdepDetectorSimulationMap;

            // maps from detsimID
            std::map<Int_t, std::vector<Int_t>> sDetectorSimulationEdepMap;
            
        };
    }
}