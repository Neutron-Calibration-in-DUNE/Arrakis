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
#include "DetectorPointCloud.h"
#include "DetectorSimulation.h"
#include "Logger.h"
#include "MCData.h"

namespace arrakis
{
    namespace melange
    {
        enum class FilterDetectorSimulation
        {
            TrackID = 0,
            EdepID = 1,
        };
        enum class NeutronCaptureGammaDetail
        {
            Simple = 0,
            Medium = 1,
            Full = 2,
        };
        class Melange
        {
        public:
            Melange(Melange &other) = delete;
            void operator=(const Melange &) = delete;

            static Melange* GetInstance();

        protected:
            Melange();
            ~Melange() {}

        public:
            void SetConfigurationParameters(const Parameters& config);

            // methods for processing event data
            void ResetEvent();
            Int_t IterateClusterLabel();
            void ProcessEvent(const Parameters& config, art::Event const& event);

            void SetLabels(std::vector<Int_t> detsimID, ShapeLabel shape, ParticleLabel particle);
            void SetLabels(std::vector<std::vector<Int_t>> detsimIDs, ShapeLabel shape, ParticleLabel particle);

            void PrepareInitialPointClouds(const Parameters& config, art::Event const& event);
            void ProcessShowers(Int_t trackID);
            void ProcessShowers(std::vector<Int_t> trackID);
            void ProcessShowers(std::vector<std::vector<Int_t>> trackID);
            void ProcessMuons(const Parameters& config, art::Event const& event);
            void ProcessAntiMuons(const Parameters& config, art::Event const& event);
            void ProcessPion0s(const Parameters& config, art::Event const& event);
            void ProcessPionPlus(const Parameters& config, art::Event const& event);
            void ProcessPionMinus(const Parameters& config, art::Event const& event);
            void ProcessNeutronCaptures(const Parameters& config, art::Event const& event);
            void ProcessAr39(const Parameters& config, art::Event const& event);
            void ProcessAr42(const Parameters& config, art::Event const& event);
            void ProcessKr85(const Parameters& config, art::Event const& event);
            void ProcessRn222(const Parameters& config, art::Event const& event);
            void ProcessCosmics(const Parameters& config, art::Event const& event);
            void CleanUpPointClouds(const Parameters& config, art::Event const& event);
            void SeparatePointClouds(const Parameters& config, art::Event const& event);

            void FillTTree();

        private:
            static Melange* sInstance;
            static std::mutex sMutex;

            // Configuration Parameters
            FilterDetectorSimulation sFilterDetectorSimulation;
            NeutronCaptureGammaDetail sNeutronCaptureGammaDetail;

            // Output TTree
            art::ServiceHandle<art::TFileService> mTFileService;
            TTree *mDetectorPointCloudTree;
            TTree *mDetectorView0PointCloudTree;
            TTree *mDetectorView1PointCloudTree;
            TTree *mDetectorView2PointCloudTree;

            TTree *mDetectorView0VoxelTree;
            TTree *mDetectorView1VoxelTree;
            TTree *mDetectorView2VoxelTree;

            // data products
            DetectorPointCloud mDetectorPointCloud;
            DetectorPointCloud mDetectorView0PointCloud;
            DetectorPointCloud mDetectorView1PointCloud;
            DetectorPointCloud mDetectorView2PointCloud;

            Int_t mClusterLabel;
        };
    }
}