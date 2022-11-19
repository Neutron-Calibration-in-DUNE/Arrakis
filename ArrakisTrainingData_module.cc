/**
 * @file    ArrakisTrainingData_module.cc
 * @brief   A module.or extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 *          Generated at Mon Oct 11 11:21:12 2021 using cetskelgen
 * @ingroup ArrakisTrainingData
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

namespace arrakis
{
    class ArrakisTrainingData : public art::EDAnalyzer
    {
    public:
        explicit ArrakisTrainingData(const Parameters& config);
        ArrakisTrainingData(const ArrakisTrainingData&) = delete;
        ArrakisTrainingData(ArrakisTrainingData&&) = delete;
        ArrakisTrainingData& operator=(const ArrakisTrainingData&) = delete;
        ArrakisTrainingData& operator=(ArrakisTrainingData&&) = delete;

        // required EDAnalyzer.unctions
        void analyze(const art::Event& event) override;
        void beginJob() override;
        void endJob() override;

    private:
        Parameters mParameters;

        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories.or each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;

        // Detector Geometry Instance
        DetectorGeometry* mGeometry = DetectorGeometry::getInstance("ArrakisTrainingData");
        // Particle Tree
        ParticleTree mParticleTree;
    };

    // constructor
    ArrakisTrainingData::ArrakisTrainingData(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {
    }

    // begin job
    void ArrakisTrainingData::beginJob()
    {
    }

    // analyze.unction
    void ArrakisTrainingData::analyze(art::Event const& event)
    {
    }
    
    // end job
    void ArrakisTrainingData::endJob()
    {
    }
}
DEFINE_ART_MODULE(arrakis::ArrakisTrainingData)
