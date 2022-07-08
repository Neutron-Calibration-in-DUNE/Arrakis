/**
 * @file ArrayGenerator.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Junying Huang
 * @author Yashwanth Bezawada
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 */
#pragma once
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

#include "larsim/Utils/TruthMatchUtils.h"

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
#include "ParticleTree.h"
#include "Configuration.h"
#include "DetectorGeometry.h"

namespace arrakis
{

    struct EventArray
    {
        std::vector<Int_t> u1_tdc;
        std::vector<Int_t> u1_channel;
        std::vector<Int_t> u1_adc;
        std::vector<std::vector<Int_t>> u1_track_ids;
        std::vector<std::vector<Double_t>> u1_energy;

        std::vector<Int_t> v1_tdc;
        std::vector<Int_t> v1_channel;
        std::vector<Int_t> v1_adc;
        std::vector<std::vector<Int_t>> v1_track_ids;
        std::vector<std::vector<Double_t>> v1_energy;

        std::vector<Int_t> z1_tdc;
        std::vector<Int_t> z1_channel;
        std::vector<Int_t> z1_adc;
        std::vector<std::vector<Int_t>> z1_track_ids;
        std::vector<std::vector<Double_t>> z1_energy;

        std::vector<Int_t> u2_tdc;
        std::vector<Int_t> u2_channel;
        std::vector<Int_t> u2_adc;
        std::vector<std::vector<Int_t>> u2_track_ids;
        std::vector<std::vector<Double_t>> u2_energy;

        std::vector<Int_t> v2_tdc;
        std::vector<Int_t> v2_channel;
        std::vector<Int_t> v2_adc;
        std::vector<std::vector<Int_t>> v2_track_ids;
        std::vector<std::vector<Double_t>> v2_energy;

        std::vector<Int_t> z2_tdc;
        std::vector<Int_t> z2_channel;
        std::vector<Int_t> z2_adc;
        std::vector<std::vector<Int_t>> z2_track_ids;
        std::vector<std::vector<Double_t>> z2_energy;
    };

    class ArrayGenerator
    {
    public:
        ArrayGenerator();
        ~ArrayGenerator();

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            //arrakis::ParticleTree const& ParticleMaps,
            // const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<raw::RawDigit>>& rawTPC
        );

        void setThreshold(Double_t Threshold) { fThreshold = Threshold; }

    private:
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fArrayTTree;

        // struct for holding event information
        EventArray fEventArray;

        geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());

        // pdg codes to construct
        VolumeType fBoundingBoxType;

        //Threshold
        Double_t fThreshold;

        // struct for holding event information
        // RawTrainingSetStruct fRawTrainingSetStruct;

        // TPC // Number of channels in each planes
        unsigned int fNUCh;
        unsigned int fNVCh;
        unsigned int fNZCh;

        // find channel boundaries for each view
        unsigned int fUChanMin;
        unsigned int fUChanMax;
        unsigned int fVChanMin;
        unsigned int fVChanMax;
        unsigned int fZChanMin;
        unsigned int fZChanMax;

        unsigned int fNofAPA; //Number of APAs
        unsigned int fChansPerAPA; //Number of channels in each APA

        // define nADC counts for uncompressed vs compressed
        unsigned int nADC_uncompPed;
    };
}