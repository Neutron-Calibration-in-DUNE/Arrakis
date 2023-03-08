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
        using TrackID_t = Int_t;
        using EdepID_t = Int_t;
        using ParticleID_t = Int_t;
        using DetSimID_t = Int_t;

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

            const simb::MCParticle& GetMCParticle(ParticleID_t index)  { return (*sMCParticleHandle)[index]; }
            const sim::SimEnergyDeposit& GetMCSimEnergyDeposit(EdepID_t index) { return (*sMCSimEnergyDepositHandle)[index]; }
            const sim::SimChannel& GetMCSimChannel(Int_t index) { return (*sMCSimChannelHandle)[index]; }
            const raw::RawDigit& GetMCRawDigit(Int_t index)     { return (*sMCRawDigitHandle)[index]; }

            const simb::MCParticle& GetMCParticleTrackID(TrackID_t TrackID)  { return (*sMCParticleHandle)[sParticleMap[TrackID]]; }

            // particle maps from track id
            inline GeneratorLabel GetGeneratorLabel(TrackID_t trackID) { return sGeneratorLabelMap[trackID]; }
            inline ParticleID_t GetParticleIndex(TrackID_t trackID)    { return sParticleMap[trackID]; }
            inline Int_t GetPDGCode(TrackID_t trackID)          { return sPDGMap[trackID]; }
            inline Int_t GetParentPDG(TrackID_t trackID)        { return sParentPDGMap[trackID]; }
            inline Int_t GetParentTrackID(TrackID_t trackID)    { return sParentTrackIDMap[trackID]; }
            inline Int_t GetParticleEnergy(TrackID_t trackID)   { return sParticleEnergyMap[trackID];}
            inline Int_t GetAncestorPDG(TrackID_t trackID)      { return sAncestorPDGMap[trackID]; }
            inline Int_t GetAncestorTrackID(TrackID_t trackID)  { return sAncestorTrackIDMap[trackID]; }
            inline Int_t GetAncestorLevel(TrackID_t trackID)    { return sAncestorLevelMap[trackID]; }
            inline Int_t GetAncestorEnergy(TrackID_t trackID)   { return sAncestorEnergyMap[trackID]; }
            inline std::vector<TrackID_t> GetDaughters(TrackID_t trackID)   { return sDaughterMap[trackID]; }
            inline std::vector<TrackID_t> GetProgeny(TrackID_t trackID)     { return sProgenyMap[trackID]; }
            inline std::vector<TrackID_t> GetAncestry(TrackID_t trackID)    { return sAncestryMap[trackID]; }
            inline std::vector<EdepID_t> GetParticleEdep(TrackID_t trackID){ return sParticleEdepMap[trackID]; }
            inline std::vector<ProcessType> GetParticleEdepProcess(TrackID_t trackID)  { return sParticleEdepProcessMap[trackID]; }
            inline std::vector<DetSimID_t> GetParticleDetectorSimulation(TrackID_t trackID) { return sParticleDetectorSimulationMap[trackID]; }
            inline std::vector<DetSimID_t> GetRandomDetectorSimulation(TrackID_t trackID)   { return sRandomDetectorSimulationMap[trackID]; }

            // maps from edep to process
            inline ProcessType GetEdepProcess(EdepID_t edepID)   { return sEdepProcessMap[edepID]; }

            // helper functions for organizing data
            ProcessType DetermineEdepProcess(const sim::SimEnergyDeposit& edep);
            std::vector<EdepID_t> DetermineDetectorSimulationEdeps(const std::vector<sim::IDE>& det_ide, DetSimID_t detsim_id);

            // functions for collecting track ids
            //std::vector<Int_t> GetPrimariesByProcess(ProcessType process_type);
            std::vector<TrackID_t> GetPrimariesByPDG(Int_t pdg);
            std::vector<TrackID_t> GetDaughtersByPDG(TrackID_t track_id, Int_t pdg);

            std::vector<TrackID_t> FilterParticlesByProcess(std::vector<TrackID_t> track_ids, ProcessType process_type);

            std::vector<EdepID_t> FilterEdepsByVolume(std::vector<EdepID_t> edep_ids, geometry::VolumeType volume_type);

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

            std::map<TrackID_t, GeneratorLabel> sGeneratorLabelMap;
            std::map<art::InputTag, GeneratorLabel> sGeneratorMap;

            // List of primary track IDs
            std::vector<TrackID_t> sPrimaries;

            // List of detector simulation structs
            std::vector<DetectorSimulation> sDetectorSimulation;
            DetectorSimulationNoise sDetectorSimulationNoise;

            // MCParticle TrackID maps
            std::map<TrackID_t, ParticleID_t> sParticleMap;
            std::map<TrackID_t, Int_t>      sPDGMap;
            std::map<TrackID_t, Int_t>      sParentPDGMap;
            std::map<TrackID_t, TrackID_t>  sParentTrackIDMap;
            std::map<TrackID_t, Double_t>   sParticleEnergyMap;
            std::map<TrackID_t, Int_t>      sAncestorPDGMap;
            std::map<TrackID_t, Int_t>      sAncestorTrackIDMap;
            std::map<TrackID_t, Int_t>      sAncestorLevelMap;
            std::map<TrackID_t, Double_t>   sAncestorEnergyMap;
            std::map<TrackID_t, std::vector<TrackID_t>> sDaughterMap;
            std::map<TrackID_t, std::vector<TrackID_t>> sProgenyMap;
            std::map<TrackID_t, std::vector<TrackID_t>> sAncestryMap;
            std::map<TrackID_t, std::vector<EdepID_t>>  sParticleEdepMap;
            std::map<EdepID_t, std::vector<ProcessType>> sParticleEdepProcessMap;
            std::map<TrackID_t, std::vector<DetSimID_t>> sParticleDetectorSimulationMap;
            std::map<TrackID_t, std::vector<DetSimID_t>> sRandomDetectorSimulationMap;

            // maps from edepID
            std::map<EdepID_t, ProcessType> sEdepProcessMap;
            std::map<EdepID_t, std::vector<DetSimID_t>> sEdepDetectorSimulationMap;

            // maps from detsimID
            std::map<DetSimID_t, std::vector<EdepID_t>> sDetectorSimulationEdepMap;
            
        };
    }
}