/**
 * @file LabelGenerator.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada
 * @author Junying Huang
 * @brief 
 * @version 0.1
 * @date 2022-07-21
 */
# pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <map>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
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
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

#include "Configuration.h"
#include "DetectorGeometry.h"
#include "ParticleTree.h"
#include "ArrayGenerator.h"
#include "GammaTable.h"

namespace arrakis
{
    struct Labels
    {
        std::vector<Int_t> u1_gamma_id;
        std::vector<Double_t> u1_gamma_type;
        std::vector<Int_t> u1_ancestor_pdg;
        std::vector<Int_t> u1_ancestor_id;

        std::vector<Int_t> v1_gamma_id;
        std::vector<Double_t> v1_gamma_type;
        std::vector<Int_t> v1_ancestor_pdg;
        std::vector<Int_t> v1_ancestor_id;

        std::vector<Int_t> z1_gamma_id;
        std::vector<Double_t> z1_gamma_type;
        std::vector<Int_t> z1_ancestor_pdg;
        std::vector<Int_t> z1_ancestor_id;

        std::vector<Int_t> u2_gamma_id;
        std::vector<Double_t> u2_gamma_type;
        std::vector<Int_t> u2_ancestor_pdg;
        std::vector<Int_t> u2_ancestor_id;

        std::vector<Int_t> v2_gamma_id;
        std::vector<Double_t> v2_gamma_type;
        std::vector<Int_t> v2_ancestor_pdg;
        std::vector<Int_t> v2_ancestor_id;

        std::vector<Int_t> z2_gamma_id;
        std::vector<Double_t> z2_gamma_type;
        std::vector<Int_t> z2_ancestor_pdg;
        std::vector<Int_t> z2_ancestor_id;
    };
    
    class LabelGenerator
    {
    public:
        LabelGenerator();
        ~LabelGenerator();

        void ResetLabels();

        int getLargestEnergyTrackID(std::vector<int> trackID_vec, std::vector<double> energy_vec);

        void processEvent(
            ParticleTree particleTree,
            GammaTable gammaTable,
            ArrayGenerator arrayGenerator
        );

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mLabelTree;

        Labels mLabels;
    };
}