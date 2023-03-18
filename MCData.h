/**
 * @file MCData.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#pragma once
#include <mutex>
#include <iostream>
#include <iomanip>

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
        /**
         * To reduce confusion throughout MCData code, we will use different names
         * for integers that mean different things.  TrackID_t refers to the actual 
         * track id of a particle, while ParticleID_t refers to the index of a particle
         * in the simb::MCParticle vector.  EdepID_t and DetSimID_t are also the indices
         * of the sim::SimEnergyDeposit and arrakis::DetectorSimulation vectors respectively.
         */
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
            void ProcessSimEnergyDeposits(art::Event const& event, 
                art::InputTag producer_label, art::InputTag instance_label
            );
            void ProcessSimChannels(art::Event const& event,
                art::InputTag producer_label, art::InputTag instance_label
            );
            void ProcessRawDigits(art::Event const& event,
                art::InputTag producer_label, art::InputTag instance_label
            );

            void SetADCThreshold(Double_t ADCThreshold);
            Double_t GetADCThreshold()  { return sADCThreshold; }

            art::Handle<std::vector<simb::MCTruth>> GetMCTruth()        { return sMCTruthHandle; }
            art::Handle<std::vector<simb::MCParticle>> GetMCParticles() { return sMCParticleHandle; }
            art::Handle<std::vector<sim::SimEnergyDeposit>> GetSimEnergyDeposits() { return sMCSimEnergyDepositHandle; }
            art::Handle<std::vector<sim::SimChannel>> GetSimChannels()  { return sMCSimChannelHandle; }
            art::Handle<std::vector<raw::RawDigit>> GetRawDigits()      { return sMCRawDigitHandle; }

            const simb::MCParticle& GetMCParticle(ParticleID_t index)  { return sMCParticleHandle->at(index); }
            const sim::SimEnergyDeposit& GetMCSimEnergyDeposit(EdepID_t index) { return sMCSimEnergyDepositHandle->at(index); }
            const sim::SimChannel& GetMCSimChannel(Int_t index) { return sMCSimChannelHandle->at(index); }
            const raw::RawDigit& GetMCRawDigit(Int_t index)     { return sMCRawDigitHandle->at(index); }

            const simb::MCParticle& GetMCParticleTrackID(TrackID_t TrackID)  { return sMCParticleHandle->at(sTrackIDParticleIDMap[TrackID]); }

            std::vector<DetectorSimulation> GetDetectorSimulation() { return sDetectorSimulation; }
            DetectorSimulationNoise GetDetectorSimulationNoise()    { return sDetectorSimulationNoise; }

            inline Double_t GetEnergyTrackID(TrackID_t trackID, Int_t precision)
            {
                Double_t factor = pow(10, precision);
                Double_t particle_energy = round(sTrackID_EnergyMap[trackID]*factor)/factor;
                return particle_energy;
            }

            /**
             * Various accessors from TrackID.  The convention for the function names are
             * "Get<Value>_<Key>", so for example "GetPDGCode_TrackID" will return the PDGCode
             * for a given TrackID.
             */
            inline ParticleID_t     GetParticleID_TrackID(TrackID_t trackID)     { return sTrackID_ParticleIDMap[trackID]; }
            inline GeneratorLabel   GetGeneratorLabel_TrackID(TrackID_t trackID) { return sTrackID_GeneratorLabelMap[trackID]; }
            inline Int_t            GetPDGCode_TrackID(TrackID_t trackID)        { return sTrackID_PDGCodeMap[trackID]; }
            inline Int_t            GetAbsPDGCode_TrackID(TrackID_t trackID)     { return std::abs(sTrackID_PDGCodeMap[trackID]); }
            inline ProcessType      GetProcess_TrackID(TrackID_t trackID)        { return sTrackID_ProcessMap[trackID]; }
            inline ProcessType      GetEndProcess_TrackID(TrackID_t trackID)     { return sTrackID_EndProcessMap[trackID]; }
            inline Double_t         GetEnergy_TrackID(TrackID_t trackID)         { return sTrackID_EnergyMap[trackID];}
            inline std::vector<TrackID_t> GetDaughterTrackID_TrackID(TrackID_t trackID)  { return sTrackID_DaughterTrackIDMap[trackID]; }
            inline std::vector<TrackID_t> GetProgenyTrackID_TrackID(TrackID_t trackID)   { return sTrackID_ProgenyTrackIDMap[trackID]; }
            inline std::vector<TrackID_t> GetAncestryTrackID_TrackID(TrackID_t trackID)  { return sTrackID_AncestryTrackIDMap[trackID]; }

            inline TrackID_t        GetParentTrackID_TrackID(TrackID_t trackID)          { return sTrackID_ParentTrackIDMap[trackID]; }
            inline ParticleID_t     GetParentParticleID_TrackID(TrackID_t trackID)       { return sTrackID_ParticleIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline GeneratorLabel   GetParentGeneratorLabel_TrackID(TrackID_t trackID)   { return sTrackID_GeneratorLabelMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline Int_t            GetParentPDGCode_TrackID(TrackID_t trackID)          { return sTrackID_ParentPDGCodeMap[trackID]; }
            inline Int_t            GetParentAbsPDGCode_TrackID(TrackID_t trackID)       { return std::abs(sTrackID_ParentPDGCodeMap[trackID]); }
            inline ProcessType      GetParentProcess_TrackID(TrackID_t trackID)          { return sTrackID_ProcessMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline ProcessType      GetParentEndProcess_TrackID(TrackID_t trackID)       { return sTrackID_EndProcessMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline Double_t         GetParentEnergy_TrackID(TrackID_t trackID)           { return sTrackID_EnergyMap[sTrackID_ParentTrackIDMap[trackID]]; }

            inline TrackID_t        GetAncestorTrackID_TrackID(TrackID_t trackID)        { return sTrackID_AncestorTrackIDMap[trackID]; }
            inline Int_t            GetAncestorLevel_TrackID(TrackID_t trackID)          { return sTrackID_AncestorLevelMap[trackID]; }
            inline ParticleID_t     GetAncestorParticleID_TrackID(TrackID_t trackID)     { return sTrackID_ParticleIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline GeneratorLabel   GetAncestorGeneratorLabel_TrackID(TrackID_t trackID) { return sTrackID_GeneratorLabelMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline Int_t            GetAncestorPDGCode_TrackID(TrackID_t trackID)        { return sTrackID_AncestorPDGCodeMap[trackID]; }
            inline Int_t            GetAncestorAbsPDGCode_TrackID(TrackID_t trackID)     { return std::abs(sTrackID_AncestorPDGCodeMap[trackID]); }
            inline ProcessType      GetAncestorProcess_TrackID(TrackID_t trackID)        { return sTrackID_ProcessMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline ProcessType      GetAncestorEndProcess_TrackID(TrackID_t trackID)     { return sTrackID_EndProcessMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline Double_t         GetAncestorEnergy_TrackID(TrackID_t trackID)         { return sTrackID_EnergyMap[sTrackID_AncestorTrackIDMap[trackID]]; }
                        

            
            inline std::vector<EdepID_t> GetParticleEdep(TrackID_t trackID){ return sParticleEdepMap[trackID]; }
            inline std::vector<ProcessType> GetParticleEdepProcess(TrackID_t trackID)  { return sParticleEdepProcessMap[trackID]; }
            inline std::vector<DetSimID_t> GetParticleDetectorSimulation(TrackID_t trackID) { return sParticleDetectorSimulationMap[trackID]; }
            inline std::vector<DetSimID_t> GetRandomDetectorSimulation(TrackID_t trackID)   { return sRandomDetectorSimulationMap[trackID]; }

            void PrintParticleData(TrackID_t trackID);
            void PrintEdepData(EdepID_t edepID);
            void PrintDetSimData(DetSimID_t detsimID);

            // maps from edep to process
            inline ProcessType GetEdepProcess(EdepID_t edepID)   { return sEdepProcessMap[edepID]; }

            // helper functions for organizing data
            ProcessType DetermineEdepProcess(const sim::SimEnergyDeposit& edep);
            std::vector<EdepID_t> DetermineDetectorSimulationEdeps(const std::vector<sim::IDE>& det_ide, DetSimID_t detsim_id);

            // functions for collecting track ids
            //std::vector<Int_t> GetPrimariesByProcess(ProcessType process_type);
            std::vector<TrackID_t> GetPrimariesByGeneratorLabel(GeneratorLabel label);
            std::vector<TrackID_t> GetPrimariesByPDG(Int_t pdg);
            std::vector<TrackID_t> GetParticlesByPDG(Int_t pdg);
            std::vector<TrackID_t> GetDaughtersByPDG(TrackID_t track_id, Int_t pdg);
            std::vector<TrackID_t> GetProgenyByPDG(TrackID_t track_id, Int_t pdg);

            std::vector<TrackID_t> FilterParticlesByProcess(std::vector<TrackID_t> track_ids, ProcessType process_type);

            std::vector<EdepID_t> GetParticleAndProgenyEdeps(TrackID_t track_id);
            std::vector<EdepID_t> GetEdepsByParticles(std::vector<TrackID_t> track_ids);
            std::vector<EdepID_t> FilterEdepsByVolume(std::vector<EdepID_t> edep_ids, geometry::VolumeType volume_type);
            std::vector<EdepID_t> FilterEdepsByPDG(std::vector<EdepID_t> edep_ids, Int_t pdg);
            
            std::vector<DetSimID_t> GetDetectorSimulationByParticles(std::vector<TrackID_t> track_ids);
            std::vector<DetSimID_t> GetDetectorSimulationByParticleAndProgeny(TrackID_t track_id);
            std::vector<DetSimID_t> GetDetectorSimulationByEdeps(std::vector<EdepID_t> edep_ids);
            std::vector<DetSimID_t> GetDetectorSimulationByParticleVolume(TrackID_t track_id, geometry::VolumeType volume_type);
            std::vector<DetSimID_t> GetDetectorSimulationByParticleProgenyVolume(TrackID_t track_id, geometry::VolumeType volume_type);
            std::vector<DetSimID_t> GetDetectorSimulationByParticleAndProgenyVolume(TrackID_t track_id, geometry::VolumeType volume_type);


            // fill TTree
            void FillTTree();

        protected:
            MCData();
            ~MCData() {}

        private:
            static MCData * sInstance;
            static std::mutex sMutex;

            // ADC threshold
            Double_t sADCThreshold = {0.0};

            // Output TTree
            art::ServiceHandle<art::TFileService> sTFileService;
            TTree *sMCDataTree;

            // handles
            std::vector<art::Handle<std::vector<simb::MCTruth>>> sMCTruthHandles;
            std::vector<std::string>                        sMCTruthHandleLabels;
            art::Handle<std::vector<simb::MCTruth>>         sMCTruthHandle;
            art::Handle<std::vector<simb::MCParticle>>      sMCParticleHandle;
            art::Handle<std::vector<sim::SimEnergyDeposit>> sMCSimEnergyDepositHandle;
            art::Handle<std::vector<sim::SimChannel>>       sMCSimChannelHandle;
            art::Handle<std::vector<raw::RawDigit>>         sMCRawDigitHandle;

            std::map<std::string, GeneratorLabel> sGeneratorMap;

            // List of primary track IDs
            std::vector<TrackID_t> sPrimaries;

            // List of detector simulation structs
            std::vector<DetectorSimulation> sDetectorSimulation;
            DetectorSimulationNoise sDetectorSimulationNoise;

            // MCParticle TrackID maps
            std::map<TrackID_t, ParticleID_t>   sTrackID_ParticleIDMap;
            std::map<TrackID_t, GeneratorLabel> sTrackID_GeneratorLabelMap;
            std::map<TrackID_t, Int_t>          sTrackID_PDGCodeMap;
            std::map<TrackID_t, ProcessType>    sTrackID_ProcessMap;
            std::map<TrackID_t, ProcessType>    sTrackID_EndProcessMap;
            std::map<TrackID_t, Double_t>       sTrackID_EnergyMap;
            std::map<TrackID_t, std::vector<TrackID_t>> sTrackID_DaughterTrackIDMap;
            std::map<TrackID_t, std::vector<TrackID_t>> sTrackID_ProgenyTrackIDMap;
            std::map<TrackID_t, std::vector<TrackID_t>> sTrackID_AncestryTrackIDMap;

            std::map<TrackID_t, TrackID_t>      sTrackID_ParentTrackIDMap;
            std::map<TrackID_t, Int_t>          sTrackID_ParentPDGCodeMap;

            std::map<TrackID_t, TrackID_t>      sTrackID_AncestorTrackIDMap;
            std::map<TrackID_t, Int_t>          sTrackID_AncestorLevelMap;
            std::map<TrackID_t, Int_t>          sTrackID_AncestorPDGCodeMap;


            
            
            std::map<TrackID_t, std::vector<EdepID_t>>  sParticleEdepMap;
            
            std::map<TrackID_t, std::vector<DetSimID_t>> sParticleDetectorSimulationMap;
            std::map<TrackID_t, std::vector<DetSimID_t>> sRandomDetectorSimulationMap;

            // maps from edepID
            std::map<EdepID_t, ProcessType> sEdepProcessMap;
            std::map<EdepID_t,  std::vector<ProcessType>> sParticleEdepProcessMap;
            std::map<EdepID_t, std::vector<DetSimID_t>> sEdepDetectorSimulationMap;

            // maps from detsimID
            std::map<DetSimID_t, std::vector<EdepID_t>> sDetectorSimulationEdepMap;
            
        };
    }
}