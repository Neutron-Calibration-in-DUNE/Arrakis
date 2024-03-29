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
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dunecore/DuneObj/OpDetDivRec.h"

#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"
#include <cmath>

#include "Configuration.h"
#include "DataWrangler.h"
#include "DetectorGeometry.h"
#include "Core.h"
#include "Logger.h"
#include "SimulationWrangler.h"
#include "SimulationLabelingLogic.h"


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
        std::string mProcessType;
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        /// TTrees
        TTree *mMetaTree;

        // Detector Geometry Instance
        DetectorGeometry* mGeometry;
        DataWrangler* mDataWrangler;
        SimulationWrangler* mSimulationWrangler;
        SimulationLabelingLogic* mSimulationLabelingLogic;

    };

    // constructor
    Arrakis::Arrakis(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {
        
        mProcessType = config().ProcessType();
        mGeometry = DetectorGeometry::GetInstance();

        mDataWrangler = DataWrangler::GetInstance();
        mDataWrangler->SetConfigurationParameters(config);

        mSimulationWrangler = SimulationWrangler::GetInstance();
        mSimulationWrangler->SetConfigurationParameters(config);

        mSimulationLabelingLogic = SimulationLabelingLogic::GetInstance();
        mSimulationLabelingLogic->SetConfigurationParameters(config);

        Logger::GetInstance("arrakis_module")->trace("initializing arrakis module");      
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
        if (mProcessType == "data") {
            mDataWrangler->ProcessEvent(mParameters, event);
            mDataWrangler->FillTTree();
        }
        else {
            mSimulationWrangler->ProcessEvent(mParameters, event);
            mSimulationLabelingLogic->ProcessEvent(mParameters, event);
            mSimulationWrangler->FillTTree();
        }
    }
    
    // end job
    void Arrakis::endJob()
    {
        DetectorGeometry::GetInstance()->FillTTree();
        // if(mSaveMeta) {
        //     mMetaTree->Fill();
        // }
    }
}
DEFINE_ART_MODULE(arrakis::Arrakis)
