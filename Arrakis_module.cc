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
#include "Core.h"
#include "Logger.h"
#include "MCData.h"

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

    private:
        Parameters mParameters;

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
        geometry::DetectorGeometry* mGeometry;
        mcdata::MCData* mMCData;

    };

    // constructor
    Arrakis::Arrakis(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {

        mGeometry = geometry::DetectorGeometry::GetInstance();
        mMCData = mcdata::MCData::GetInstance();

        Logger::GetInstance("arrakis_module")->trace("initializing arrakis module");

        // // generator labels
        // if(mGeneratePointCloudData)
        // {
        //     mAr39Label = mParameters().Ar39Label();
        //     Logger::GetInstance("arrakis_module")->trace("setting Ar39Label = " + mAr39Label.label());
        //     mSingleNeutronLabel = mParameters().SingleNeutronLabel();
        //     Logger::GetInstance("arrakis_module")->trace("setting SingleNeutronLabel = " + mSingleNeutronLabel.label());
        //     mPNSLabel = mParameters().PNSLabel();
        //     Logger::GetInstance("arrakis_module")->trace("setting PNSLabel = " + mPNSLabel.label());
        //     mGeneratorMap[mAr39Label] = GeneratorLabel::kAr39;
        //     mGeneratorMap[mSingleNeutronLabel] = GeneratorLabel::kSingleNeutron;
        //     mGeneratorMap[mPNSLabel] = GeneratorLabel::kPNS;
        // }
    }

    // begin job
    void Arrakis::beginJob()
    {
    }

    // analyze function
    void Arrakis::analyze(art::Event const& event)
    {
        auto event_id = event.id().event();
        auto run_id = event.run();
        auto sub_run_id = event.subRun();
        Logger::GetInstance("arrakis_module")->trace(
            "processing event [" + 
            std::to_string(run_id) + ":" + 
            std::to_string(sub_run_id) + ":" + 
            std::to_string(event_id) + "]"
        );

        mMCData->ProcessEvent(mParameters, event);
        
        // detinfo::DetectorClocksData const clock_data(
        //     art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        // );


        // Logger::GetInstance("arrakis_module")->trace("processed MCParticle, SimEnergyDeposit and MCTruth products");
        //     // mPrimaryData->ProcessEventDetectorSimulation(
        //     //     mParticleMaps, 
        //     //     clock_data,
        //     //     mc_sim_channels,
        //     //     mc_raw_digits
        //     // );
        //     Logger::GetInstance("arrakis_module")->trace("processed SimChannel and RawDigit products");
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
        // if(mSaveMeta) {
        //     mMetaTree->Fill();
        // }
    }
}
DEFINE_ART_MODULE(arrakis::Arrakis)
