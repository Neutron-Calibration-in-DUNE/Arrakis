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
#include "MCTree.h"

namespace arrakis
{
    /**
     * @brief Arrakis module for processing
     * MC data.
     */
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

        void ProcessMCTruth(
            const art::Event& event, art::InputTag input_tag, GeneratorLabel label
        );
        art::ValidHandle<std::vector<simb::MCParticle>> GetMCParticles(
            const art::Event& event
        );
        art::ValidHandle<std::vector<sim::SimEnergyDeposit>> GetSimEnergyDeposits(
            const art::Event& event
        );

    private:
        Parameters mParameters;
        // Set of configuration parameters
        bool    mSaveMeta;
        bool    mSaveGeometry;

        bool    mSaveParticleMaps;
        bool    mSavePrimaryData;
        bool    mSavePrimaryDataEdeps;
        bool    mSavePrimaryDataRawTPC;

        bool    mGeneratePointCloudData;

        Double_t mADCThreshold;

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

        art::Handle<std::vector<simb::MCTruth>>         mMCTruthHandle;
        art::Handle<std::vector<simb::MCParticle>>      mMCParticleHandle;
        art::Handle<std::vector<sim::SimEnergyDeposit>> mMCSimEnergyDepositHandle;
        art::Handle<std::vector<sim::SimChannel>>       mMCSimChannelHandle;
        art::Handle<std::vector<raw::RawDigit>>         mMCRawDigitHandle;

        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        /// TTrees
        TTree *mMetaTree;

        // Detector Geometry Instance
        geometry::DetectorGeometry* mGeometry;
        Logger* mLogger;
        MCData* mMCData;

        // Particle Tree
        ParticleMaps* mParticleMaps;

        // Primary Data
        //PrimaryData* mPrimaryData;

        // MCTree
        mctree::MCTree* mMCTree;

        // PointCloudGenerator
        //PointCloudGenerator* mPointCloudGenerator;

    };

    // constructor
    Arrakis::Arrakis(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {

        mGeometry = geometry::DetectorGeometry::GetInstance("arrakis_module");
        mLogger = Logger::GetInstance("arrakis_module");
        mMCData = mcdata::MCData::GetInstance("arrakis_module");

        mLogger->trace("initializing arrakis module");

        // // Set various configuration parameters
        // mSaveMeta = mParameters().SaveMeta();
        // mLogger->trace("setting SaveMeta = " + std::to_string(mSaveMeta));

        // mSaveGeometry = mParameters().SaveGeometry();
        // mLogger->trace("setting SaveGeometry = " + std::to_string(mSaveGeometry));

        // mSaveParticleMaps = mParameters().SaveParticleMaps();
        // mLogger->trace("setting SaveParticleMaps = " + std::to_string(mSaveParticleMaps));
        // mSavePrimaryData =  mParameters().SavePrimaryData();
        // mLogger->trace("setting SavePrimaryData = " + std::to_string(mSavePrimaryData));
        // mSavePrimaryDataEdeps =  mParameters().SavePrimaryDataEdeps();
        // mLogger->trace("setting SavePrimaryDataEdeps = " + std::to_string(mSavePrimaryDataEdeps));
        // mSavePrimaryDataRawTPC =  mParameters().SavePrimaryDataRawTPC();
        // mLogger->trace("setting SavePrimaryDataRawTPC = " + std::to_string(mSavePrimaryDataRawTPC));

        // // solo point clouds
        // mGeneratePointCloudData = mParameters().GeneratePointCloudData();
        // mLogger->trace("setting GeneratePointCloudData = " + std::to_string(mGeneratePointCloudData));

        // // module labels
        // mLArGeantProducerLabel =    mParameters().LArGeantProducerLabel();
        // mLogger->trace("setting LArGeantProducerLabel = " + mLArGeantProducerLabel.label());
        // mIonAndScintProducerLabel = mParameters().IonAndScintProducerLabel();
        // mLogger->trace("setting IonAndScintProducerLabel = " + mIonAndScintProducerLabel.label());
        // mSimChannelProducerLabel =  mParameters().SimChannelProducerLabel();
        // mLogger->trace("setting SimChannelProducerLabel = " + mSimChannelProducerLabel.label());
        // mSimChannelInstanceProducerLabel = mParameters().SimChannelInstanceProducerLabel();
        // mLogger->trace("setting SimChannelInstanceProducerLabel = " + mSimChannelInstanceProducerLabel.label());
        // mTPCInputLabel =    mParameters().TPCInputLabel();
        // mLogger->trace("setting TPCInputLabel = " + mTPCInputLabel.label());
        // mTPCInstanceLabel = mParameters().TPCInstanceLabel();
        // mLogger->trace("setting TPCInstanceLabel = " + mTPCInstanceLabel.label());

        // // generator labels
        // if(mGeneratePointCloudData)
        // {
        //     mAr39Label = mParameters().Ar39Label();
        //     mLogger->trace("setting Ar39Label = " + mAr39Label.label());
        //     mSingleNeutronLabel = mParameters().SingleNeutronLabel();
        //     mLogger->trace("setting SingleNeutronLabel = " + mSingleNeutronLabel.label());
        //     mPNSLabel = mParameters().PNSLabel();
        //     mLogger->trace("setting PNSLabel = " + mPNSLabel.label());
        //     mGeneratorMap[mAr39Label] = GeneratorLabel::kAr39;
        //     mGeneratorMap[mSingleNeutronLabel] = GeneratorLabel::kSingleNeutron;
        //     mGeneratorMap[mPNSLabel] = GeneratorLabel::kPNS;
        // }

        // mADCThreshold = mParameters().ADCThreshold();
        // mLogger->trace("setting ADCThreshold = " + std::to_string(mADCThreshold));
        

        // mParticleMaps = new ParticleMaps(
        //     mSaveParticleMaps
        // );

        // // mPrimaryData = new PrimaryData(
        // //     mSavePrimaryData,
        // //     mSavePrimaryDataEdeps,
        // //     mSavePrimaryDataRawTPC,
        // //     mADCThreshold
        // // );

        // //mPointCloudGenerator = new PointCloudGenerator();

        // if(mSaveMeta) {
        //     mMetaTree = mTFileService->make<TTree>("meta", "meta");
        // }
    }

    // begin job
    void Arrakis::beginJob()
    {
        if(mSaveGeometry) {
            mGeometry->FillTTree();
        }
    }
    void Arrakis::ProcessMCTruth(
        const art::Event& event, art::InputTag input_tag, GeneratorLabel label
    )
    {
        mLogger->trace("collecting simb::MCTruth from input_tag <" + input_tag.label() + ">");
        if(!event.getByLabel(input_tag, mMCTruthHandle))
        {
            mLogger->warning("no input_tag matching " + input_tag.label() + " for simb::MCTruth");
        }
        else 
        {
            auto mc_truth = event.getValidHandle<std::vector<simb::MCTruth>>(
                input_tag
            );
            auto particles = event.getValidHandle<std::vector<simb::MCParticle>>(
                mLArGeantProducerLabel
            );
            mParticleMaps->ProcessMCTruth(label, mc_truth, particles);
        }
    }
    art::ValidHandle<std::vector<simb::MCParticle>> Arrakis::GetMCParticles(
        const art::Event& event
    )
    {
        mLogger->trace("collecting simb::MCParticle from label <" + mLArGeantProducerLabel.label() + ">");
        if(!event.getByLabel(mLArGeantProducerLabel, mMCParticleHandle))
        {
            mLogger->error("no label matching " + mLArGeantProducerLabel.label() + " for simb::MCParticle!");
            exit(0);
        }
        else 
        {
            return event.getValidHandle<std::vector<simb::MCParticle>>(
                mLArGeantProducerLabel
            );
        }
    }
    art::ValidHandle<std::vector<sim::SimEnergyDeposit>> Arrakis::GetSimEnergyDeposits(
        const art::Event& event
    )
    {
        mLogger->trace("collecting sim::SimEnergyDeposit from label <" + mIonAndScintProducerLabel.label() + ">");
        if(!event.getByLabel(mIonAndScintProducerLabel, mMCSimEnergyDepositHandle))
        {
            mLogger->error("no label matching " + mIonAndScintProducerLabel.label() + " for sim::SimEnergyDeposit!");
            exit(0);
        }
        else 
        {
            return event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(
                mIonAndScintProducerLabel
            );
        }
    }

    // analyze.unction
    void Arrakis::analyze(art::Event const& event)
    {
        auto event_id = event.id().event();
        auto run_id = event.run();
        auto sub_run_id = event.subRun();
        mLogger->trace(
            "processing event " + 
            std::to_string(run_id) + ":" + 
            std::to_string(sub_run_id) + ":" + 
            std::to_string(event_id)
        );


        /**
         * @details  For each event, we will look through the various
         * available data products and send event info to the 
         * corresponding submodules that process them, starting with mc_particles
         * then SimEnergyDeposit, SimChannel and RawDigit.
         */
        // detinfo::DetectorClocksData const clock_data(
        //     art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        // );

        // // prepare generator labels, mc particles and sim energy deposits
        // auto mc_particles = GetMCParticles(event);
        // auto mc_energy_deposits = GetSimEnergyDeposits(event);

        // // Add the particle maps and primary data
        // // for MCParticle and SimEnergyDeposit.
        // mParticleMaps->ProcessEvent(
        //     mc_particles
        // );
        // // Add MCTruth labels to particle maps.
        // for(auto const& [key, val] : mGeneratorMap)
        // {
        //     ProcessMCTruth(event, key, val);
        // }
        // /**
        //  * This section processes MCTruth, mc_particles, SimEnergyDeposits,
        //  * SimChannel, and RawDigit into PrimaryData objects to be used later.
        // */
        // // mPrimaryData->ProcessEventMC(
        // //     mParticleMaps, 
        // //     mc_particles, 
        // //     mc_energy_deposits
        // // );
        // mMCTree->ProcessEventMC(
        //     mParticleMaps,
        //     mc_particles,
        //     mc_energy_deposits
        // );
        // mLogger->trace("processed MCParticle, SimEnergyDeposit and MCTruth products");

        // // Check if SimChannel and RawDigit are available,
        // // and then process those into primary data.
        // if(
        //     event.getByLabel(
        //         art::InputTag(
        //             mSimChannelProducerLabel.label(), 
        //             mSimChannelInstanceProducerLabel.label()
        //         ), 
        //         mMCSimChannelHandle
        //     ) &&
        //     event.getByLabel(
        //         art::InputTag(mTPCInputLabel.label(), mTPCInstanceLabel.label()), 
        //         mMCRawDigitHandle
        //     )
        // )
        // {
        //     auto mc_sim_channels = 
        //         event.getValidHandle<std::vector<sim::SimChannel>>(art::InputTag(
        //             mSimChannelProducerLabel.label(), 
        //             mSimChannelInstanceProducerLabel.label()
        //         )
        //     );
        //     auto mc_raw_digits = 
        //         event.getValidHandle<std::vector<raw::RawDigit>>(art::InputTag(
        //             mTPCInputLabel.label(), 
        //             mTPCInstanceLabel.label()
        //         )
        //     );
        //     // mPrimaryData->ProcessEventDetectorSimulation(
        //     //     mParticleMaps, 
        //     //     clock_data,
        //     //     mc_sim_channels,
        //     //     mc_raw_digits
        //     // );
        //     mLogger->trace("processed SimChannel and RawDigit products");
        // }
        

        // /**
        //  * @brief Now that everything is collected, we pass the data to 
        //  * the various classes which construct training data.
        //  */
        // // if(mGeneratePointCloudData) 
        // // {
        // //     mPointCloudGenerator->ProcessEvent(
        // //         mParticleMaps,
        // //         mPrimaryData
        // //     );
        // // }
        // mParticleMaps->FillTTree();
        //mPrimaryData->FillTTree();
        
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
