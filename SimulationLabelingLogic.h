/**
 * @file SimulationLabelingLogic.h
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
#include "SimulationWrangler.h"
#include "WirePlanePointCloud.h"

namespace arrakis
{
    enum class NeutronCaptureGammaDetail
    {
        Simple = 0,
        Medium = 1,
        Full = 2,
    };
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
        Int_t IterateTopologyLabel();
        void ProcessEvent(const Parameters& config, art::Event const& event);

        SourceLabel DetermineSourceLabel(TrackID_t trackID);

        void SetLabels(
            DetSimID_List detsimID, EdepID_List edepID, TrackID_t track_id,
            TopologyLabel topology, PhysicsLabel physics,
            Int_t unique_topology
        );
        void SetLabels(
            DetSimID_Collection detsimIDs, EdepID_Collection edepIDs, TrackID_List trackIDList,
            TopologyLabel topology, PhysicsLabel physics,
            Int_t unique_topology
        );
        void SetLabels(
            std::vector<DetSimID_Collection> detsimIDs, 
            std::vector<EdepID_Collection> edepIDs,
            TrackID_Collection trackIDCollection,
            TopologyLabel topology, PhysicsLabel physics,
            Int_t unique_topology
        );

        void PrepareInitialPointClouds(const Parameters& config, art::Event const& event);

        void ProcessShowers(TrackID_t trackID, Int_t TopologyLabel);
        void ProcessShowers(TrackID_List trackID, Int_t TopologyLabel);
        void ProcessShowers(TrackID_Collection trackID);

        void ProcessNoise(const Parameters& config, art::Event const& event);

        void ProcessElectrons(const Parameters& config, art::Event const& event);
        void ProcessPositrons(const Parameters& config, art::Event const& event);
        void ProcessElectronNeutrinos(const Parameters& config, art::Event const& event);
        void ProcessAntiElectronNeutrinos(const Parameters& config, art::Event const& event);

        void ProcessMuons(const Parameters& config, art::Event const& event);
        void ProcessAntiMuons(const Parameters& config, art::Event const& event);
        void ProcessMuonNeutrinos(const Parameters& config, art::Event const& event);
        void ProcessAntiMuonNeutrinos(const Parameters& config, art::Event const& event);

        void ProcessTauons(const Parameters& config, art::Event const& event);
        void ProcessAntiTauons(const Parameters& config, art::Event const& event);
        void ProcessTauonNeutrinos(const Parameters& config, art::Event const& event);
        void ProcessAntiTauonNeutrinos(const Parameters& config, art::Event const& event);

        void ProcessGammas(const Parameters& config, art::Event const& event);

        void ProcessPion0s(const Parameters& config, art::Event const& event);
        void ProcessPionPlus(const Parameters& config, art::Event const& event);
        void ProcessPionMinus(const Parameters& config, art::Event const& event);

        void ProcessKaon0s(const Parameters& config, art::Event const& event);
        void ProcessKaonPlus(const Parameters& config, art::Event const& event);
        void ProcessKaonMinus(const Parameters& config, art::Event const& event);

        void ProcessProtons(const Parameters& config, art::Event const& event);
        
        void ProcessNeutronCaptures(const Parameters& config, art::Event const& event);
        void ProcessNuclearRecoils(const Parameters& config, art::Event const& event);
        void ProcessElectronRecoils(const Parameters& config, art::Event const& event);

        void ProcessAr39(const Parameters& config, art::Event const& event);
        void ProcessAr42(const Parameters& config, art::Event const& event);
        void ProcessK42(const Parameters& config, art::Event const& event);
        void ProcessKr85(const Parameters& config, art::Event const& event);
        void ProcessRn222(const Parameters& config, art::Event const& event);
        void ProcessPo218(const Parameters& config, art::Event const& event);
        void ProcessAt218(const Parameters& config, art::Event const& event);
        void ProcessRn218(const Parameters& config, art::Event const& event);
        void ProcessPb214(const Parameters& config, art::Event const& event);
        void ProcessBi214(const Parameters& config, art::Event const& event);
        void ProcessPo214(const Parameters& config, art::Event const& event);
        void ProcessTl210(const Parameters& config, art::Event const& event);
        void ProcessPb210(const Parameters& config, art::Event const& event);
        void ProcessBi210(const Parameters& config, art::Event const& event);
        void ProcessPo210(const Parameters& config, art::Event const& event);

        void ProcessCosmics(const Parameters& config, art::Event const& event);
        void CleanUpPointClouds(const Parameters& config, art::Event const& event);

        void FillTTree();

    private:
        static SimulationLabelingLogic* sInstance;
        static std::mutex sMutex;

        // Configuration Parameters
        NeutronCaptureGammaDetail sNeutronCaptureGammaDetail;
        Int_t sInducedChannelInfluence;
        Int_t sInducedTDCInfluence;
        Double_t sShowerEnergyThreshold;

        // Output TTree
        art::ServiceHandle<art::TFileService> mTFileService;

        Int_t mTopologyLabel;
        Int_t mParticleLabel;
    };
}