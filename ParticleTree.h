/**
 * @file ParticleTree.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
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

namespace arrakis
{
    class ParticleTree
    {
    public:
        ParticleTree();
        ~ParticleTree();

        void processEvent(const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles);

        Int_t GetParentTrackID(Int_t trackID) const { return mParentTrackIDMap[trackID]; }
        Int_t GetAncestorTrackID(Int_t trackID) const { return mAncestorTrackIDMap[trackID]; }

        Int_t GetPDGCode(Int_t trackID) const { return mPDGMap[trackID]; }
        Int_t GetParentPDG(Int_t trackID) const { return mParentPDGMap[trackID]; }
        Int_t GetAncestorPDG(Int_t trackID) const { return mAncestorPDGMap[trackID]; }

        Int_t GetAncestorLevel(Int_t trackID) const { return mAncestorLevelMap[trackID]; }

    private:
        std::map<Int_t, Int_t> mParentTrackIDMap;
        std::map<Int_t, Int_t> mAncestorTrackIDMap;

        std::map<Int_t, Int_t> mPDGMap;
        std::map<Int_t, Int_t> mParentPDGMap;
        std::map<Int_t, Int_t> mAncestorPDGMap;
        
        std::map<Int_t, Int_t> mAncestorLevelMap;
    };
}