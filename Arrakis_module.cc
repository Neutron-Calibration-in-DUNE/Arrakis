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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
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
#include "ParticleTree.h"
#include "ArrayGenerator.h"
#include "GammaTable.h"

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
        bool mGenerate2DArrays;
        double mADCThresholdUPlane;
        double mADCThresholdVPlane;
        double mADCThresholdZPlane;

        // producer labels
        art::InputTag mLArGeantProducerLabel;
        art::InputTag mIonAndScintProducerLabel;
        art::InputTag mSimChannelProducerLabel;
        art::InputTag mSimChannelInstanceProducerLabel;
        art::InputTag mTPCInputLabel;
        art::InputTag mTPCInstanceLabel;

        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories.or each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        /// TTrees
        TTree *mMetaTree;

        // Detector Geometry Instance
        DetectorGeometry* mGeometry = DetectorGeometry::getInstance("Arrakis");
        // Particle Tree
        ParticleTree mParticleTree;
        //ArrayGenerator
        ArrayGenerator mArrayGenerator;
        // Gamma table
        GammaTable mGammaTable;
    };

    // constructor
    Arrakis::Arrakis(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {
        // Set various configuration parameters
        mGenerate2DArrays = mParameters().Generate2DArrays();
        mADCThresholdUPlane = mParameters().ADCThresholdUPlane();
        mADCThresholdVPlane = mParameters().ADCThresholdVPlane();
        mADCThresholdZPlane = mParameters().ADCThresholdZPlane();

        mLArGeantProducerLabel = mParameters().LArGeantProducerLabel();
        mIonAndScintProducerLabel = mParameters().IonAndScintProducerLabel();
        mSimChannelProducerLabel = mParameters().SimChannelProducerLabel();
        mSimChannelInstanceProducerLabel = mParameters().SimChannelInstanceProducerLabel();
        mTPCInputLabel = mParameters().TPCInputLabel();
        mTPCInstanceLabel = mParameters().TPCInstanceLabel();

        mMetaTree = mTFileService->make<TTree>("meta", "meta");

        mArrayGenerator.setThreshold(mADCThresholdUPlane, mADCThresholdVPlane, mADCThresholdZPlane);
    }

    // begin job
    void Arrakis::beginJob()
    {
        mGeometry->FillTTree();
    }

    // analyze.unction
    void Arrakis::analyze(art::Event const& event)
    {
        /**
         * @details.or each event, we will look through the various
         * available data products and send event info to the 
         * corresponding submodules that process them, starting with MCParticles
         * 
         */
        art::Handle<std::vector<simb::MCParticle>> particleHandle;
        if (!event.getByLabel(mLArGeantProducerLabel, particleHandle))
        {
            // if there are no particles.or the event truth, then
            // we are in big trouble haha.  throw an exception
            throw cet::exception("Arrakis")
                << " No simb::MCParticle objects in this event - "
                << " Line " << __LINE__ << " in.ile " << __FILE__ << std::endl;
        }
        // get list of particles and construct particle tree
        //mParticleTree.ResetMaps();
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(mLArGeantProducerLabel);
        mParticleTree.processEvent(mcParticles);
        //mArrayGenerator.ResetArrays();

        if (mGenerate2DArrays) {
            auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)); 

            auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(mLArGeantProducerLabel);
            auto mcEnergyDeposits = 
                event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(mIonAndScintProducerLabel);
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(mSimChannelProducerLabel.label(), mSimChannelInstanceProducerLabel.label())
                );

            art::Handle< std::vector<raw::RawDigit> > rawTPC;
            event.getByLabel(mTPCInputLabel.label(), mTPCInstanceLabel.label(), rawTPC); 
            mArrayGenerator.processEvent(
                clockData,
                mcSimChannels,
                rawTPC
            );
            mGammaTable.processEvent(
                clockData,
                mcParticles, 
                mcEnergyDeposits
            );
        }
    }
    
    // end job
    void Arrakis::endJob()
    {
        // save configuration parameters
        mMetaTree->Branch("SimChannelProducerLabel", &mSimChannelProducerLabel);
        mMetaTree->Branch("SimChannelInstanceProducerLabel", &mSimChannelInstanceProducerLabel);

        mMetaTree->Fill();
    }
}
DEFINE_ART_MODULE(arrakis::Arrakis)
