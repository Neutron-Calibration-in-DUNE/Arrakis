/**
 * @file    Arrakis_module.cc
 * @brief   A module.or extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 *          Generated at Mon Oct 11 11:21:12 2021 using cetskelgen
 * @ingroup Arrakis
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
 *          Yashwanth Bezawada
**/
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"

// #include "larcore/Geometry/Geometry.h"
// #include "larcorealg/Geometry/GeometryCore.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"
#include <cmath>

#include "Configuration.h"
#include "DetectorGeometry.h"
#include "Generators.h"
#include "ParticleMaps.h"
#include "PrimaryData.h"
#include "SoloPointCloudGenerator.h"

namespace arrakis
{
    class Arrakis : public art::EDAnalyzer
    {
    public:
        explicit Arrakis(const Parameters& config);
        Arrakis(const Arrakis&) = delete;
        Arrakis(Arrakis&&) = delete;
        Arrakis& operator=(const Arrakis&) = delete;
        Arrakis& operator=(Arrakis&&) = delete;

        // required EDAnalyzer.unctions
        void analyze(const art::Event& event) override;
        void beginJob() override;
        void endJob() override;

    private:
        Parameters mParameters;
        // Set of configuration parameters
        bool    mSaveMeta;
        bool    mSaveGeometry;

        bool    mSaveParticleMaps;
        bool    mSavePrimaryData;
        bool    mSavePrimaryDataEdeps;
        bool    mSavePrimaryDataRawTPC;

        bool    mGenerateSoloPointCloudData;

        bool    mAr39Simulated;
        
        bool    mGenerate2DArrays;
        bool    mGenerateEvtLvlNInfo;
        bool    mGenerateNCapInfo;
        
        double  mADCThresholdUPlane;
        double  mADCThresholdVPlane;
        double  mADCThresholdZPlane;
        int     mClockTicks;

        // producer labels
        art::InputTag mLArGeantProducerLabel;
        art::InputTag mIonAndScintProducerLabel;
        art::InputTag mSimChannelProducerLabel;
        art::InputTag mSimChannelInstanceProducerLabel;
        art::InputTag mTPCInputLabel;
        art::InputTag mTPCInstanceLabel;

        std::vector<art::InputTag> mGeneratorInputTags;
        art::InputTag mAr39InputTag;

        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        /// TTrees
        TTree *mMetaTree;

        // Detector Geometry Instance
        DetectorGeometry* mGeometry;
        // Generator Labels
        Generators* mGenerators;
        // Particle Tree
        ParticleMaps* mParticleMaps;
        // Primary Data
        PrimaryData* mPrimaryData;
        // SoloPointCloudGenerator
        SoloPointCloudGenerator* mSoloPointCloudGenerator;

    };

    // constructor
    Arrakis::Arrakis(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {
        // Set various configuration parameters
        mSaveMeta = mParameters().SaveMeta();
        mSaveGeometry = mParameters().SaveGeometry();

        mSaveParticleMaps = mParameters().SaveParticleMaps();
        mSavePrimaryData =  mParameters().SavePrimaryData();
        mSavePrimaryDataEdeps =  mParameters().SavePrimaryDataEdeps();
        mSavePrimaryDataRawTPC =  mParameters().SavePrimaryDataRawTPC();

        // solo point clouds
        mGenerateSoloPointCloudData = mParameters().GenerateSoloPointCloudData();

        // module labels
        mLArGeantProducerLabel =    mParameters().LArGeantProducerLabel();
        mIonAndScintProducerLabel = mParameters().IonAndScintProducerLabel();
        mSimChannelProducerLabel =  mParameters().SimChannelProducerLabel();
        mSimChannelInstanceProducerLabel = mParameters().SimChannelInstanceProducerLabel();
        mTPCInputLabel =    mParameters().TPCInputLabel();
        mTPCInstanceLabel = mParameters().TPCInstanceLabel();

        // generator labels
        mAr39Simulated = mParameters().Ar39Simulated();
        if(mAr39Simulated) {
            mGeneratorLabels.emplace_back(mParameters().Ar39Label());
        }

        mGeometry = DetectorGeometry::GetInstance("Arrakis");

        mParticleMaps = new ParticleMaps(
            mSaveParticleMaps
        );

        mPrimaryData = new PrimaryData(
            mSavePrimaryData,
            mSavePrimaryDataEdeps,
            mSavePrimaryDataRawTPC
        );

        mSoloPointCloudGenerator = new SoloPointCloudGenerator();

        if(mSaveMeta) {
            mMetaTree = mTFileService->make<TTree>("meta", "meta");
        }
    }

    // begin job
    void Arrakis::beginJob()
    {
        if(mSaveGeometry) {
            mGeometry->FillTTree();
        }
    }

    // analyze.unction
    void Arrakis::analyze(art::Event const& event)
    {
        /**
         * @details  For each event, we will look through the various
         * available data products and send event info to the 
         * corresponding submodules that process them, starting with MCParticles
         * then SimEnergyDeposit, SimChannel and RawDigit.
         */
        detinfo::DetectorClocksData const clockData(
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        );
        art::Handle<std::vector<simb::MCParticle>>      mcParticleHandle;
        art::Handle<std::vector<sim::SimEnergyDeposit>> mcSimEnergyDepositHandle;
        art::Handle<std::vector<sim::SimChannel>>       mcSimChannelHandle;
        art::Handle<std::vector<raw::RawDigit>>         mcRawDigitHandle;

        // prepare generator labels, mc particles and sim energy deposits
        std::vector<art::ValidHandle<simb::MCTruth>> mcTruth;
        for(auto label : mGeneratorLabels) {
            mcTruth.emplace_back(
                event.getValidHandle<simb::MCTruth>(label)
            );
        }
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(
            mLArGeantProducerLabel
        );
        auto mcEnergyDeposits = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(
            mIonAndScintProducerLabel
        );
        // Add the particle maps and primary data
        // for MCParticle and SimEnergyDeposit.
        mGenerators->ProcessEvent(
            mGeneratorLabels,
            mcTruth
        );
        mParticleMaps->ProcessEvent(
            mGenerators, 
            mcParticles
        );
        mPrimaryData->ProcessEventMC(
            mParticleMaps, 
            mcParticles, 
            mcEnergyDeposits
        );
        // Check if SimChannel and RawDigit are available,
        // and then process those into primary data.
        if(
            event.getByLabel(
                art::InputTag(mSimChannelProducerLabel.label(), mSimChannelInstanceProducerLabel.label()), 
                mcSimChannelHandle
            ) &&
            event.getByLabel(
                art::InputTag(mTPCInputLabel.label(), mTPCInstanceLabel.label()), 
                mcRawDigitHandle
            )
        )
        {
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(art::InputTag(
                    mSimChannelProducerLabel.label(), 
                    mSimChannelInstanceProducerLabel.label()
                )
            );
            auto rawDigit = 
                event.getValidHandle<std::vector<raw::RawDigit>>(art::InputTag(
                    mTPCInputLabel.label(), 
                    mTPCInstanceLabel.label()
                )
            );
            mPrimaryData->ProcessEventDetectorSimulation(
                mParticleMaps, 
                clockData,
                mcSimChannels,
                rawDigit
            );
        }
        mPrimaryData->FillTTree();
        /**
         * @brief Now that everything is collected, we pass the data to 
         * the various classes which construct training data.
         */
        if(mGenerateSoloPointCloudData) {
            mSoloPointCloudGenerator->ProcessEvent(
                mParticleMaps,
                mPrimaryData
            );
        }
        
    }
    
    // end job
    void Arrakis::endJob()
    {
        if(mSaveMeta) {
            mMetaTree->Fill();
        }
    }
}
DEFINE_ART_MODULE(arrakis::Arrakis)
