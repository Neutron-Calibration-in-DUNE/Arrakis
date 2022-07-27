/**
 * @file EvtLvlNeutronInfo.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * Stores total event level information (Like total summed ADC of neutrons, total neutrons captured, etc.)
 * @version 0.1
 * @date 2022-07-26
 */
#pragma once
#include <string>
#include <vector>
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

#include "DetectorGeometry.h"

namespace arrakis
{

    struct EventStatistics
    {
        Int_t total_neutrons_captured = 0;  // Total neutrons captured in the TPC in an event

        std::vector<Int_t> u_summed_adc;    // ADC summed across all channels for a TDC value
        Int_t u_total_summed_adc = 0;       // Total summed ADC in an event across all channels and TDCs
        Double_t u_total_summed_energy = 0; // Total evergy deposited in the TPC in an event

        std::vector<Int_t> v_summed_adc;
        Int_t v_total_summed_adc = 0;
        Double_t v_total_summed_energy = 0;

        std::vector<Int_t> z_summed_adc;
        Int_t z_total_summed_adc = 0;
        Double_t z_total_summed_energy = 0;
    };

    class EvtLvlNeutronInfo
    {
    public:
        EvtLvlNeutronInfo();
        ~EvtLvlNeutronInfo();

        void ResetArrays();

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::Handle<std::vector<raw::RawDigit>>& rawTPC
        );

        void setThreshold(Double_t U_Threshold, Double_t V_Threshold, Double_t Z_Threshold) { 
            fUPlaneThreshold = U_Threshold;
            fVPlaneThreshold = V_Threshold;
            fZPlaneThreshold = Z_Threshold;
        }

        void setClockTicks(int clockTicks) {
            fClockTicks = clockTicks;
        }

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mEventStatisticsTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("EvtLvlNeutronInfo");

        geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());

        //Threshold
        Double_t fUPlaneThreshold;
        Double_t fVPlaneThreshold;
        Double_t fZPlaneThreshold;

        EventStatistics mEventStatistics;

        // max total clock ticks per event
        int fClockTicks;

        // number of clock ticks as unsigned int
        unsigned int clock_ticks;
    };
}