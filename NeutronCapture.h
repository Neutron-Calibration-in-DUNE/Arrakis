/**
 * @file NeutronCapture.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-09-21
 */
#pragma once
#include <string>
#include <vector>
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

#include "DetectorGeometry.h"

namespace arrakis
{
    struct NCapture{

        std::vector<Int_t> mPDGCode;
        std::vector<std::string> mProcess;
        std::vector<std::string> mEndProcess;

    };

    class NeutronCapture
    {
    public:
        NeutronCapture();
        ~NeutronCapture();
        void ResetArrays();

        void setBoundingBoxType(std::string volumeType);

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
            //const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
        );

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mNeutronCaptureTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("NeutronCapture");

        NCapture mNCapture;
        // std::map<Int_t, Int_t> mNeutronCaptureIndex;
        // std::vector<Int_t> mNeutronCaptureTrackIDs;
    };
}