/**
 * @file SimulationLabelingLogic.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#pragma once
#include <mutex>
#include <algorithm>

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
#include "SimulationWrangler.h"
#include "WirePlanePointCloud.h"

namespace arrakis
{
    class SimulationLabelingLogic
    {
    public:
        SimulationLabelingLogic(SimulationLabelingLogic &other) = delete;
        void operator=(const SimulationLabelingLogic &) = delete;

        static SimulationLabelingLogic* GetInstance();

    protected:
        SimulationLabelingLogic();
        ~SimulationLabelingLogic() {}

    public:
        void SetConfigurationParameters(const Parameters& config);

        // methods for processing event data
        void ResetEvent();
        Int_t UniqueTopology();
        Int_t UniquePhysicsMicro();
        Int_t UniquePhysicsMeso();
        Int_t UniquePhysicsMacro();
        void ProcessEvent(const Parameters& config, art::Event const& event);

        PhysicsMacroLabel DeterminePhysicsMacroLabel(TrackID_t trackID);

        void SetLabels(
            DetSimID_List detsimID, 
            EdepID_List edepID, 
            TrackID_t track_id,
            TopologyLabel topology, 
            PhysicsMicroLabel physicsMicroLabel,
            PhysicsMesoLabel physicsMesoLabel,
            PhysicsMacroLabel physicsMacroLabel,
            Int_t uniqueTopologyLabel,
            Int_t uniquePhysicsMicroLabel,
            Int_t uniquePhysicsMesoLabel,
            Int_t uniquePhysicsMacroLabel
        );
        void SetLabels(
            DetSimID_Collection 
            detsimIDs, 
            EdepID_Collection edepIDs, 
            TrackID_List trackIDList,
            TopologyLabel topology, 
            PhysicsMicroLabel physicsMicroLabel,
            PhysicsMesoLabel physicsMesoLabel,
            PhysicsMacroLabel physicsMacroLabel,
            Int_t uniqueTopologyLabel,
            Int_t uniquePhysicsMicroLabel,
            Int_t uniquePhysicsMesoLabel,
            Int_t uniquePhysicsMacroLabel
        );
        void SetLabels(
            std::vector<DetSimID_Collection> detsimIDs, 
            std::vector<EdepID_Collection> edepIDs,
            TrackID_Collection trackIDCollection,
            TopologyLabel topology, 
            PhysicsMicroLabel physicsMicroLabel,
            PhysicsMesoLabel physicsMesoLabel,
            PhysicsMacroLabel physicsMacroLabel,
            Int_t uniqueTopologyLabel,
            Int_t uniquePhysicsMicroLabel,
            Int_t uniquePhysicsMesoLabel,
            Int_t uniquePhysicsMacroLabel
        );

        void PrepareInitialPointClouds(const Parameters& config, art::Event const& event);

        void ProcessShowers(TrackID_t trackID, Int_t TopologyLabel);
        void ProcessShowers(TrackID_List trackID, Int_t TopologyLabel);
        void ProcessShowers(TrackID_Collection trackID);

        void ProcessNoise(const Parameters& config, art::Event const& event);

        void ProcessMCTruth(const Parameters& config, art::Event const& event);

        void ProcessElectrons(const Parameters& config, art::Event const& event);
        void ProcessPositrons(const Parameters& config, art::Event const& event);
        void ProcessElectronNeutrinos(const Parameters& config, art::Event const& event);
        void ProcessAntiElectronNeutrinos(const Parameters& config, art::Event const& event);

        void ProcessMuons(const Parameters& config, art::Event const& event);
        void ProcessMuonNeutrinos(const Parameters& config, art::Event const& event);
        void ProcessAntiMuonNeutrinos(const Parameters& config, art::Event const& event);

        void ProcessTauons(const Parameters& config, art::Event const& event);
        void ProcessTauonNeutrinos(const Parameters& config, art::Event const& event);
        void ProcessAntiTauonNeutrinos(const Parameters& config, art::Event const& event);

        void ProcessGammas(const Parameters& config, art::Event const& event);

        void ProcessPion0s(const Parameters& config, art::Event const& event);
        void ProcessPions(const Parameters& config, art::Event const& event);

        void ProcessKaon0s(const Parameters& config, art::Event const& event);
        void ProcessKaons(const Parameters& config, art::Event const& event);

        void ProcessProtons(const Parameters& config, art::Event const& event);
        
        void ProcessNeutrons(const Parameters& config, art::Event const& event);
        void ProcessNuclearRecoils(const Parameters& config, art::Event const& event);
        void ProcessElectronRecoils(const Parameters& config, art::Event const& event);

        void ProcessAr39(const Parameters& config, art::Event const& event);
        void ProcessAr42(const Parameters& config, art::Event const& event);
        void ProcessKr85(const Parameters& config, art::Event const& event);
        void ProcessRn222(const Parameters& config, art::Event const& event);

        void ProcessCosmics(const Parameters& config, art::Event const& event);
        void CleanUpPointClouds(const Parameters& config, art::Event const& event);

        void FillTTree();

    private:
        static SimulationLabelingLogic* sInstance;
        static std::mutex sMutex;

        // Configuration Parameters
        Int_t sInducedChannelInfluence;
        Int_t sInducedTDCInfluence;
        Double_t sShowerEnergyThreshold;

        // Output TTree
        art::ServiceHandle<art::TFileService> mTFileService;

        Int_t mUniqueTopologyLabel;
        Int_t mUniquePhysicsMicroLabel;
        Int_t mUniquePhysicsMesoLabel;
        Int_t mUniquePhysicsMacroLabel;

        std::vector<PhysicsMesoLabel> mRn222PhysicsMeso = {
            PhysicsMesoLabel::AlphaDecay, 
            PhysicsMesoLabel::AlphaDecay, PhysicsMesoLabel::BetaDecay,
            PhysicsMesoLabel::AlphaDecay, PhysicsMesoLabel::BetaDecay,
            PhysicsMesoLabel::AlphaDecay,
            PhysicsMesoLabel::BetaDecay,
            PhysicsMesoLabel::AlphaDecay, PhysicsMesoLabel::BetaDecay,
            PhysicsMesoLabel::AlphaDecay,
            PhysicsMesoLabel::BetaDecay,
            PhysicsMesoLabel::AlphaDecay, PhysicsMesoLabel::BetaDecay,            
            PhysicsMesoLabel::AlphaDecay, PhysicsMesoLabel::BetaDecay,
            PhysicsMesoLabel::AlphaDecay
        };
        std::vector<PhysicsMacroLabel> mRn222PhysicsMacro = {
            PhysicsMacroLabel::Rn222, 
            PhysicsMacroLabel::Po218a, PhysicsMacroLabel::Po218b,
            PhysicsMacroLabel::At218a, PhysicsMacroLabel::At218b,
            PhysicsMacroLabel::Rn218,
            PhysicsMacroLabel::Pb214,
            PhysicsMacroLabel::Bi214a, PhysicsMacroLabel::Bi214b,
            PhysicsMacroLabel::Po214,
            PhysicsMacroLabel::Tl210,
            PhysicsMacroLabel::Pb210a, PhysicsMacroLabel::Pb210b,
            PhysicsMacroLabel::Bi210a, PhysicsMacroLabel::Bi210b,
            PhysicsMacroLabel::Po210
        };
        std::vector<Int_t> mRn222PDGs = {
            1000020040,
            1000020040, 11,
            1000020040, 11,
            1000020040,
            11,
            1000020040, 11,
            1000020040,
            11,
            1000020040, 11,
            1000020040, 11,
            1000020040
        };
        std::vector<Double_t> mRn222Energies = {
            5.590,
            6.115, 0.294,
            6.874, 2.883,
            7.263,
            1.024,
            5.627, 3.272,
            7.833,
            5.484,
            3.792, 0.064,
            5.037, 1.163,
            5.407
        };
    };
}