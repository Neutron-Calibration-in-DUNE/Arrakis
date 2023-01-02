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
#include "Logger.h"
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
        Logger* mLogger;
        Parameters mParameters;
        // Set of configuration parameters
        bool    mSaveMeta;
        bool    mSaveGeometry;

        bool    mSaveParticleMaps;
        bool    mSavePrimaryData;
        bool    mSavePrimaryDataEdeps;
        bool    mSavePrimaryDataRawTPC;

        bool    mGenerateSoloPointCloudData;

        // producer labels
        art::InputTag mLArGeantProducerLabel;
        art::InputTag mIonAndScintProducerLabel;
        art::InputTag mSimChannelProducerLabel;
        art::InputTag mSimChannelInstanceProducerLabel;
        art::InputTag mTPCInputLabel;
        art::InputTag mTPCInstanceLabel;

        std::map<art::InputTag, GeneratorLabel> mGeneratorMap;
        art::InputTag mAr39Label;
        art::InputTag mSingleNeutronLabel;
        art::InputTag mPNSLabel;

        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        /// TTrees
        TTree *mMetaTree;

        // Detector Geometry Instance
        DetectorGeometry* mGeometry;
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
        if(mGenerateSoloPointCloudData)
        {
            mAr39Label = mParameters().Ar39Label();
            mSingleNeutronLabel = mParameters().SingleNeutronLabel();
            mPNSLabel = mParameters().PNSLabel();
            mGeneratorMap[mAr39Label] = GeneratorLabel::kAr39;
            mGeneratorMap[mSingleNeutronLabel] = GeneratorLabel::kSingleNeutron;
            mGeneratorMap[mPNSLabel] = GeneratorLabel::kPNS;
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
         * corresponding submodules that process them, starting with mc_particles
         * then SimEnergyDeposit, SimChannel and RawDigit.
         */
        detinfo::DetectorClocksData const clock_data(
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        );

        art::Handle<std::vector<simb::MCTruth>>         mc_truth_handle;
        art::Handle<std::vector<simb::MCParticle>>      mc_particle_handle;
        art::Handle<std::vector<sim::SimEnergyDeposit>> mc_sim_energy_deposit_handle;
        art::Handle<std::vector<sim::SimChannel>>       mc_sim_channel_handle;
        art::Handle<std::vector<raw::RawDigit>>         mc_raw_digit_handle;

        // prepare generator labels, mc particles and sim energy deposits
        auto mc_particles = event.getValidHandle<std::vector<simb::MCParticle>>(
            mLArGeantProducerLabel
        );
        auto mc_energy_deposits = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(
            mIonAndScintProducerLabel
        );
        // Add the particle maps and primary data
        // for MCParticle and SimEnergyDeposit.
        mParticleMaps->ProcessEvent(
            mc_particles
        );
        // Add MCTruth labels to particle maps.
        for(auto const& [key, val] : mGeneratorMap)
        {
            if(event.getByLabel(key, mc_truth_handle))
            {
                auto mc_truth = event.getValidHandle<std::vector<simb::MCTruth>>(key);
                mParticleMaps->ProcessMCTruth(val, mc_truth);
            }
            else
            {


            }
        }
        /**
         * This section processes MCTruth, mc_particles, SimEnergyDeposits,
         * SimChannel, and RawDigit into PrimaryData objects to be used later.
        */
        mPrimaryData->ProcessEventMC(
            mParticleMaps, 
            mc_particles, 
            mc_energy_deposits
        );
        ARRAKIS_TRACE("Processed MC...");

        // Check if SimChannel and RawDigit are available,
        // and then process those into primary data.
        if(
            event.getByLabel(
                art::InputTag(
                    mSimChannelProducerLabel.label(), 
                    mSimChannelInstanceProducerLabel.label()
                ), 
                mc_sim_channel_handle
            ) &&
            event.getByLabel(
                art::InputTag(mTPCInputLabel.label(), mTPCInstanceLabel.label()), 
                mc_raw_digit_handle
            )
        )
        {
            auto mc_sim_channels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(art::InputTag(
                    mSimChannelProducerLabel.label(), 
                    mSimChannelInstanceProducerLabel.label()
                )
            );
            auto mc_raw_digits = 
                event.getValidHandle<std::vector<raw::RawDigit>>(art::InputTag(
                    mTPCInputLabel.label(), 
                    mTPCInstanceLabel.label()
                )
            );
            mPrimaryData->ProcessEventDetectorSimulation(
                mParticleMaps, 
                clock_data,
                mc_sim_channels,
                mc_raw_digits
            );
        }
        std::cout << "Processed Detector Simulation..." << std::endl;

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
        mParticleMaps->FillTTree();
        mPrimaryData->FillTTree();
        
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
