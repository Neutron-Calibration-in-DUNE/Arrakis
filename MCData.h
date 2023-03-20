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

            art::Handle<std::vector<simb::MCTruth>>         GetMCTruth()        { return sMCTruthHandle; }
            art::Handle<std::vector<simb::MCParticle>>      GetMCParticles()    { return sMCParticleHandle; }
            art::Handle<std::vector<sim::SimEnergyDeposit>> GetSimEnergyDeposits() { return sMCSimEnergyDepositHandle; }
            art::Handle<std::vector<sim::SimChannel>>       GetSimChannels()    { return sMCSimChannelHandle; }
            art::Handle<std::vector<raw::RawDigit>>         GetRawDigits()      { return sMCRawDigitHandle; }

            const simb::MCParticle& GetMCParticle(ParticleID_t index)   { return sMCParticleHandle->at(index); }
            const sim::SimChannel& GetMCSimChannel(Int_t index)         { return sMCSimChannelHandle->at(index); }
            const raw::RawDigit& GetMCRawDigit(Int_t index)             { return sMCRawDigitHandle->at(index); }

            std::vector<DetectorSimulation> GetDetectorSimulation() { return sDetectorSimulation; }
            DetectorSimulationNoise GetDetectorSimulationNoise()    { return sDetectorSimulationNoise; }

            /**
             * Various accessors from TrackID.  The convention for the function names are
             * "Get<Value>_<Key>", so for example "GetPDGCode_TrackID" will return the PDGCode
             * for a given TrackID by accessing the map sTrackID_PDGCode.
             */
            inline ParticleID_t     GetParticleID_TrackID(TrackID_t trackID)     { return sTrackID_ParticleIDMap[trackID]; }
            inline GeneratorLabel   GetGeneratorLabel_TrackID(TrackID_t trackID) { return sTrackID_GeneratorLabelMap[trackID]; }
            inline Int_t            GetPDGCode_TrackID(TrackID_t trackID)        { return sTrackID_PDGCodeMap[trackID]; }
            inline Int_t            GetAbsPDGCode_TrackID(TrackID_t trackID)     { return std::abs(sTrackID_PDGCodeMap[trackID]); }
            inline ProcessType      GetProcess_TrackID(TrackID_t trackID)        { return sTrackID_ProcessMap[trackID]; }
            inline ProcessType      GetEndProcess_TrackID(TrackID_t trackID)     { return sTrackID_EndProcessMap[trackID]; }
            inline Double_t         GetEnergy_TrackID(TrackID_t trackID)         { return sTrackID_EnergyMap[trackID];}
            inline Double_t         GetEnergy_TrackID(TrackID_t trackID, Double_t precision){ return round(sTrackID_EnergyMap[trackID] * pow(10.0, precision)) / pow(10.0, precision); }
            inline TrackID_List     GetDaughterTrackID_TrackID(TrackID_t trackID)   { return sTrackID_DaughterTrackIDMap[trackID]; }
            inline TrackID_List     GetProgenyTrackID_TrackID(TrackID_t trackID)    { return sTrackID_ProgenyTrackIDMap[trackID]; }
            inline TrackID_List     GetDescendantTrackID_TrackID(TrackID_t trackID) { return sTrackID_DescendantTrackIDMap[trackID]; }
            inline TrackID_List     GetAncestryTrackID_TrackID(TrackID_t trackID)   { return sTrackID_AncestryTrackIDMap[trackID]; }
            inline EdepID_List      GetEdepID_TrackID(TrackID_t trackID)            { return sTrackID_EdepIDMap[trackID]; }
            inline ProcessType_List GetEdepProcess_TrackID(TrackID_t trackID)       { return sTrackID_EdepProcessMap[trackID]; }
            inline DetSimID_List    GetDetSimID_TrackID(TrackID_t trackID)          { return sTrackID_DetSimIDMap[trackID]; }
            inline DetSimID_List    GetRandomDetSimID_TrackID(TrackID_t trackID)    { return sTrackID_RandomDetSimIDMap[trackID]; }
            const simb::MCParticle& GetMCParticle_TrackID(TrackID_t trackID)        { return sMCParticleHandle->at(sTrackID_ParticleIDMap[trackID]); }

            inline TrackID_t        GetParentTrackID_TrackID(TrackID_t trackID)         { return sTrackID_ParentTrackIDMap[trackID]; }
            inline ParticleID_t     GetParentParticleID_TrackID(TrackID_t trackID)      { return sTrackID_ParticleIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline GeneratorLabel   GetParentGeneratorLabel_TrackID(TrackID_t trackID)  { return sTrackID_GeneratorLabelMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline Int_t            GetParentPDGCode_TrackID(TrackID_t trackID)         { return sTrackID_ParentPDGCodeMap[trackID]; }
            inline Int_t            GetParentAbsPDGCode_TrackID(TrackID_t trackID)      { return std::abs(sTrackID_ParentPDGCodeMap[trackID]); }
            inline ProcessType      GetParentProcess_TrackID(TrackID_t trackID)         { return sTrackID_ProcessMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline ProcessType      GetParentEndProcess_TrackID(TrackID_t trackID)      { return sTrackID_EndProcessMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline Double_t         GetParentEnergy_TrackID(TrackID_t trackID)          { return sTrackID_EnergyMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline Double_t         GetParentEnergy_TrackID(TrackID_t trackID, Double_t precision) { return GetEnergy_TrackID(sTrackID_ParentTrackIDMap[trackID] , precision); }
            inline TrackID_List     GetParentDaughterTrackID_TrackID(TrackID_t trackID) { return sTrackID_DaughterTrackIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline TrackID_List     GetParentProgenyTrackID_TrackID(TrackID_t trackID)  { return sTrackID_ProgenyTrackIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline TrackID_List     GetParentDescendantTrackID_TrackID(TrackID_t trackID) { return sTrackID_DescendantTrackIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline TrackID_List     GetParentAncestryTrackID_TrackID(TrackID_t trackID) { return sTrackID_AncestryTrackIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline EdepID_List      GetParentEdepID_TrackID(TrackID_t trackID)          { return sTrackID_EdepIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline ProcessType_List GetParentEdepProcess_TrackID(TrackID_t trackID)     { return sTrackID_EdepProcessMap[sTrackID_ParentTrackIDMap[trackID]]; }
            inline DetSimID_List    GetParentDetSimID_TrackID(TrackID_t trackID)        { return sTrackID_DetSimIDMap[sTrackID_ParentTrackIDMap[trackID]]; }
            const simb::MCParticle& GetParentMCParticle_TrackID(TrackID_t trackID)      { return sMCParticleHandle->at(sTrackID_ParticleIDMap[sTrackID_ParentTrackIDMap[trackID]]); }

            inline TrackID_t        GetAncestorTrackID_TrackID(TrackID_t trackID)       { return sTrackID_AncestorTrackIDMap[trackID]; }
            inline Int_t            GetAncestorLevel_TrackID(TrackID_t trackID)         { return sTrackID_AncestorLevelMap[trackID]; }
            inline ParticleID_t     GetAncestorParticleID_TrackID(TrackID_t trackID)    { return sTrackID_ParticleIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline GeneratorLabel   GetAncestorGeneratorLabel_TrackID(TrackID_t trackID){ return sTrackID_GeneratorLabelMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline Int_t            GetAncestorPDGCode_TrackID(TrackID_t trackID)       { return sTrackID_AncestorPDGCodeMap[trackID]; }
            inline Int_t            GetAncestorAbsPDGCode_TrackID(TrackID_t trackID)    { return std::abs(sTrackID_AncestorPDGCodeMap[trackID]); }
            inline ProcessType      GetAncestorProcess_TrackID(TrackID_t trackID)       { return sTrackID_ProcessMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline ProcessType      GetAncestorEndProcess_TrackID(TrackID_t trackID)    { return sTrackID_EndProcessMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline Double_t         GetAncestorEnergy_TrackID(TrackID_t trackID)        { return sTrackID_EnergyMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline Double_t         GetAncestorEnergy_TrackID(TrackID_t trackID, Double_t precision) { return GetEnergy_TrackID(sTrackID_AncestorTrackIDMap[trackID] , precision); }
            inline TrackID_List     GetAncestorDaughterTrackID_TrackID(TrackID_t trackID) { return sTrackID_DaughterTrackIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline TrackID_List     GetAncestorProgenyTrackID_TrackID(TrackID_t trackID)  { return sTrackID_ProgenyTrackIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline TrackID_List     GetAncestorDescendantTrackID_TrackID(TrackID_t trackID) { return sTrackID_DescendantTrackIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline TrackID_List     GetAncestorAncestryTrackID_TrackID(TrackID_t trackID) { return sTrackID_AncestryTrackIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline EdepID_List      GetAncestorEdepID_TrackID(TrackID_t trackID)          { return sTrackID_EdepIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline ProcessType_List GetAncestorEdepProcess_TrackID(TrackID_t trackID)   { return sTrackID_EdepProcessMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            inline DetSimID_List    GetAncestorDetSimID_TrackID(TrackID_t trackID)      { return sTrackID_DetSimIDMap[sTrackID_AncestorTrackIDMap[trackID]]; }
            const simb::MCParticle& GetAncestorMCParticle_TrackID(TrackID_t trackID)    { return sMCParticleHandle->at(sTrackID_ParticleIDMap[sTrackID_AncestorTrackIDMap[trackID]]); }            

            /**
             * Various functions for collecting TrackIDs from primaries, etc.
             * The output is always a vector of TrackIDs, but the input can vary.
             */
            TrackID_List GetPrimaries_GeneratorLabel(GeneratorLabel label);
            TrackID_List GetPrimaries_PDGCode(Int_t pdg);
            TrackID_List GetPrimaries_AbsPDGCode(Int_t pdg);
            TrackID_List GetPrimaries_Process(ProcessType process);
            TrackID_List GetPrimaries_EndProcess(ProcessType process);

            TrackID_List GetTrackID_GeneratorLabel(GeneratorLabel label);
            TrackID_List GetTrackID_PDGCode(Int_t pdg);
            TrackID_List GetTrackID_AbsPDGCode(Int_t pdg);
            TrackID_List GetTrackID_Process(ProcessType process);
            TrackID_List GetTrackID_EndProcess(ProcessType process);

            DetSimID_Collection GetDetSimID_TrackID(TrackID_List trackIDs);
            std::vector<DetSimID_Collection> GetDetSimID_TrackID(TrackID_Collection trackIDs);

            TrackID_Collection GetDaughterTrackID_GeneratorLabel(GeneratorLabel label);
            TrackID_Collection GetDaughterTrackID_PDGCode(Int_t pdg);
            TrackID_Collection GetDaughterTrackID_AbsPDGCode(Int_t pdg);
            TrackID_Collection GetDaughterTrackID_Process(ProcessType process);
            TrackID_Collection GetDaughterTrackID_EndProcess(ProcessType process);
            TrackID_Collection GetDaughterTrackID_TrackID(TrackID_List trackIDs);

            TrackID_List FilterTrackID_NotAbsPDGCode(TrackID_List& trackIDs, Int_t pdg);
            TrackID_List FilterTrackID_AbsPDGCode(TrackID_List& trackIDs, Int_t pdg);
            TrackID_List FilterTrackID_Process(TrackID_List& trackIDs, ProcessType process);
            TrackID_List FilterTrackID_NotProcess(TrackID_List& trackIDs, ProcessType process);

            TrackID_Collection FilterTrackID_NotAbsPDGCode(TrackID_Collection& trackIDs, Int_t pdg);
            TrackID_Collection FilterTrackID_AbsPDGCode(TrackID_Collection& trackIDs, Int_t pdg);
            TrackID_Collection FilterTrackID_Process(TrackID_Collection& trackIDs, ProcessType process);
            TrackID_Collection FilterTrackID_NotProcess(TrackID_Collection& trackIDs, ProcessType process);
            
            // DetSimID_List GetDaughterDetSimID_GeneratorLabel(GeneratorLabel label);
            // DetSimID_List GetDaughterDetSimID_PDGCode(Int_t pdg);
            // DetSimID_List GetDaughterDetSimID_AbsPDGCode(Int_t pdg);
            // DetSimID_List GetDaughterDetSimID_Process(ProcessType process);
            // DetSimID_List GetDaughterDetSimID_EndProcess(ProcessType process);
            // DetSimID_Collection GetDaughterDetSimID_TrackID(TrackID_t trackID);

            TrackID_Collection GetProgenyTrackID_GeneratorLabel(GeneratorLabel label);
            TrackID_Collection GetProgenyTrackID_PDGCode(Int_t pdg);
            TrackID_Collection GetProgenyTrackID_AbsPDGCode(Int_t pdg);
            TrackID_Collection GetProgenyTrackID_Process(ProcessType process);
            TrackID_Collection GetProgenyTrackID_EndProcess(ProcessType process);
            TrackID_Collection GetProgenyTrackID_TrackID(TrackID_List trackIDs);

            // DetSimID_List GetProgenyDetSimID_GeneratorLabel(GeneratorLabel label);
            // DetSimID_List GetProgenyDetSimID_PDGCode(Int_t pdg);
            // DetSimID_List GetProgenyDetSimID_AbsPDGCode(Int_t pdg);
            // DetSimID_List GetProgenyDetSimID_Process(ProcessType process);
            // DetSimID_List GetProgenyDetSimID_EndProcess(ProcessType process);
            // DetSimID_Collection GetProgenyDetSimID_TrackID(TrackID_t trackID);

            TrackID_Collection GetDescendantTrackID_GeneratorLabel(GeneratorLabel label);
            TrackID_Collection GetDescendantTrackID_PDGCode(Int_t pdg);
            TrackID_Collection GetDescendantTrackID_AbsPDGCode(Int_t pdg);
            TrackID_Collection GetDescendantTrackID_Process(ProcessType process);
            TrackID_Collection GetDescendantTrackID_EndProcess(ProcessType process);
            TrackID_Collection GetDescendantTrackID_TrackID(TrackID_List trackIDs);

            /**
             * Various accessors for EdepID.  Convention is the same as TrackID accessors,
             * "Get<Value>_EdepID".
             */
            inline ProcessType GetProcess_EdepID(EdepID_t edepID)               { return sEdepID_ProcessMap[edepID]; }
            inline DetSimID_List GetDetSimID_EdepID(EdepID_t edepID)  { return sEdepID_DetSimIDMap[edepID]; }
            const sim::SimEnergyDeposit& GetMCSimEnergyDeposit_EdepID(EdepID_t edepID)  { return sMCSimEnergyDepositHandle->at(edepID); }

            DetSimID_List GetAllDetSimID_TrackID(TrackID_t track_id);

            void PrintParticleData(TrackID_t trackID);
            void PrintEdepData(EdepID_t edepID);
            void PrintDetSimData(DetSimID_t detsimID);

            // helper functions for organizing data
            ProcessType DetermineEdepProcess(const sim::SimEnergyDeposit& edep);
            EdepID_List DetermineDetectorSimulationEdeps(const std::vector<sim::IDE>& det_ide, DetSimID_t detsim_id);


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
            TrackID_List sPrimaries;

            // List of detector simulation structs
            std::vector<DetectorSimulation> sDetectorSimulation;
            DetectorSimulationNoise sDetectorSimulationNoise;

            // TrackID maps
            std::map<TrackID_t, ParticleID_t>   sTrackID_ParticleIDMap;
            std::map<TrackID_t, GeneratorLabel> sTrackID_GeneratorLabelMap;
            std::map<TrackID_t, Int_t>          sTrackID_PDGCodeMap;
            std::map<TrackID_t, ProcessType>    sTrackID_ProcessMap;
            std::map<TrackID_t, ProcessType>    sTrackID_EndProcessMap;
            std::map<TrackID_t, Double_t>       sTrackID_EnergyMap;
            std::map<TrackID_t, TrackID_List>   sTrackID_DaughterTrackIDMap;
            std::map<TrackID_t, TrackID_List>   sTrackID_ProgenyTrackIDMap;
            std::map<TrackID_t, TrackID_List>   sTrackID_DescendantTrackIDMap;
            std::map<TrackID_t, TrackID_List>   sTrackID_AncestryTrackIDMap;
            std::map<TrackID_t, EdepID_List>    sTrackID_EdepIDMap;
            std::map<TrackID_t, ProcessType_List> sTrackID_EdepProcessMap;
            std::map<TrackID_t, DetSimID_List>  sTrackID_DetSimIDMap;
            std::map<TrackID_t, DetSimID_List>  sTrackID_RandomDetSimIDMap; 

            std::map<TrackID_t, TrackID_t>      sTrackID_ParentTrackIDMap;
            std::map<TrackID_t, Int_t>          sTrackID_ParentPDGCodeMap;

            std::map<TrackID_t, TrackID_t>      sTrackID_AncestorTrackIDMap;
            std::map<TrackID_t, Int_t>          sTrackID_AncestorLevelMap;
            std::map<TrackID_t, Int_t>          sTrackID_AncestorPDGCodeMap;

            // EdepID maps
            std::map<EdepID_t, ProcessType> sEdepID_ProcessMap;
            std::map<EdepID_t, DetSimID_List> sEdepID_DetSimIDMap;

            // DetSmID maps
            std::map<DetSimID_t, EdepID_List> sDetSimID_EdepIDMap;
            
        };
    }
}