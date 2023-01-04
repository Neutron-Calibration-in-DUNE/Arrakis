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
        auto logger = Logger::GetInstance("arrakis_module");
        logger->trace("initializing arrakis module");

        // Set various configuration parameters
        mSaveMeta = mParameters().SaveMeta();
        logger->trace("setting SaveMeta = " + std::to_string(mSaveMeta));

        mSaveGeometry = mParameters().SaveGeometry();
        logger->trace("setting SaveGeometry = " + std::to_string(mSaveGeometry));

        mSaveParticleMaps = mParameters().SaveParticleMaps();
        logger->trace("setting SaveParticleMaps = " + std::to_string(mSaveParticleMaps));
        mSavePrimaryData =  mParameters().SavePrimaryData();
        logger->trace("setting SavePrimaryData = " + std::to_string(mSavePrimaryData));
        mSavePrimaryDataEdeps =  mParameters().SavePrimaryDataEdeps();
        logger->trace("setting SavePrimaryDataEdeps = " + std::to_string(mSavePrimaryDataEdeps));
        mSavePrimaryDataRawTPC =  mParameters().SavePrimaryDataRawTPC();
        logger->trace("setting SavePrimaryDataRawTPC = " + std::to_string(mSavePrimaryDataRawTPC));

        // solo point clouds
        mGenerateSoloPointCloudData = mParameters().GenerateSoloPointCloudData();
        logger->trace("setting GenerateSoloPointCloudData = " + std::to_string(mGenerateSoloPointCloudData));

        // module labels
        mLArGeantProducerLabel =    mParameters().LArGeantProducerLabel();
        logger->trace("setting LArGeantProducerLabel = " + mLArGeantProducerLabel.label());
        mIonAndScintProducerLabel = mParameters().IonAndScintProducerLabel();
        logger->trace("setting IonAndScintProducerLabel = " + mIonAndScintProducerLabel.label());
        mSimChannelProducerLabel =  mParameters().SimChannelProducerLabel();
        logger->trace("setting SimChannelProducerLabel = " + mSimChannelProducerLabel.label());
        mSimChannelInstanceProducerLabel = mParameters().SimChannelInstanceProducerLabel();
        logger->trace("setting SimChannelInstanceProducerLabel = " + mSimChannelInstanceProducerLabel.label());
        mTPCInputLabel =    mParameters().TPCInputLabel();
        logger->trace("setting TPCInputLabel = " + mTPCInputLabel.label());
        mTPCInstanceLabel = mParameters().TPCInstanceLabel();
        logger->trace("setting TPCInstanceLabel = " + mTPCInstanceLabel.label());

        // generator labels
        if(mGenerateSoloPointCloudData)
        {
            mAr39Label = mParameters().Ar39Label();
            logger->trace("setting Ar39Label = " + mAr39Label.label());
            mSingleNeutronLabel = mParameters().SingleNeutronLabel();
            logger->trace("setting SingleNeutronLabel = " + mSingleNeutronLabel.label());
            mPNSLabel = mParameters().PNSLabel();
            logger->trace("setting PNSLabel = " + mPNSLabel.label());
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
    void Arrakis::ProcessMCTruth(
        const art::Event& event, art::InputTag input_tag, GeneratorLabel label
    )
    {
        Logger::GetInstance("arrakis_module")->trace("collecting simb::MCTruth from input_tag <" + input_tag.label() + ">");
        if(!event.getByLabel(input_tag, mMCTruthHandle))
        {
            Logger::GetInstance("arrakis_module")->warning("no input_tag matching " + input_tag.label() + " for simb::MCTruth");
        }
        else 
        {
            const art::FindManyP<simb::MCTruth> mc_truth_ptrs(
                mMCParticleHandle, event, mLArGeantProducerLabel
            );
            auto mc_truth = event.getValidHandle<std::vector<simb::MCTruth>>(
                input_tag
            );
            auto particles = event.getValidHandle<std::vector<simb::MCParticle>>(
                mLArGeantProducerLabel
            );
            for(size_t ii = 0; ii < (*particles).size(); ii++)
            {
                auto const& truth = mc_truth_ptrs.at(ii);
                std::cout << (*particles)[ii].TrackId() << ", " << mc_truth << ", ";
                std::cout << mc_truth_ptrs.at(ii).size() << ":\n";
                for(size_t jj = 0; jj < truth.size(); jj++)
                {
                    std::cout << "\t " << truth[jj] << ", " << truth[jj].GetParticle(0).TrackId() << "\n";
                }
            }

            mParticleMaps->ProcessMCTruth(label, mc_truth);
        }
    }
    art::ValidHandle<std::vector<simb::MCParticle>> Arrakis::GetMCParticles(
        const art::Event& event
    )
    {
        Logger::GetInstance("arrakis_module")->trace("collecting simb::MCParticle from label <" + mLArGeantProducerLabel.label() + ">");
        if(!event.getByLabel(mLArGeantProducerLabel, mMCParticleHandle))
        {
            Logger::GetInstance("arrakis_module")->error("no label matching " + mLArGeantProducerLabel.label() + " for simb::MCParticle!");
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
        Logger::GetInstance("arrakis_module")->trace("collecting sim::SimEnergyDeposit from label <" + mIonAndScintProducerLabel.label() + ">");
        if(!event.getByLabel(mIonAndScintProducerLabel, mMCSimEnergyDepositHandle))
        {
            Logger::GetInstance("arrakis_module")->error("no label matching " + mIonAndScintProducerLabel.label() + " for sim::SimEnergyDeposit!");
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
        Logger::GetInstance("arrakis_module")->trace(
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
        detinfo::DetectorClocksData const clock_data(
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        );

        // prepare generator labels, mc particles and sim energy deposits
        auto mc_particles = GetMCParticles(event);
        auto mc_energy_deposits = GetSimEnergyDeposits(event);

        // Add the particle maps and primary data
        // for MCParticle and SimEnergyDeposit.
        mParticleMaps->ProcessEvent(
            mc_particles
        );
        // Add MCTruth labels to particle maps.
        for(auto const& [key, val] : mGeneratorMap)
        {
            ProcessMCTruth(event, key, val);
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
        Logger::GetInstance("arrakis_module")->trace("processed MCParticle, SimEnergyDeposit and MCTruth products");

        // Check if SimChannel and RawDigit are available,
        // and then process those into primary data.
        if(
            event.getByLabel(
                art::InputTag(
                    mSimChannelProducerLabel.label(), 
                    mSimChannelInstanceProducerLabel.label()
                ), 
                mMCSimChannelHandle
            ) &&
            event.getByLabel(
                art::InputTag(mTPCInputLabel.label(), mTPCInstanceLabel.label()), 
                mMCRawDigitHandle
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
            Logger::GetInstance("arrakis_module")->trace("processed SimChannel and RawDigit products");
        }
        

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
