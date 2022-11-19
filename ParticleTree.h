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

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>

namespace arrakis
{
    class ParticleTree
    {
    public:
        ParticleTree();
        ~ParticleTree();

        void ResetMaps();
        void processEvent(const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles);
        
        inline Int_t GetPDGCode(Int_t trackID)          { return mPDGMap[trackID]; }
        inline Int_t GetParentPDG(Int_t trackID)        { return mParentPDGMap[trackID]; }
        inline Int_t GetParentTrackID(Int_t trackID)    { return mParentTrackIDMap[trackID]; }
        inline Int_t GetParticleEnergy(Int_t trackID)   { return mParticleEnergyMap[trackID];}
        inline Int_t GetAncestorPDG(Int_t trackID)      { return mAncestorPDGMap[trackID]; }
        inline Int_t GetAncestorTrackID(Int_t trackID)  { return mAncestorTrackIDMap[trackID]; }
        inline Int_t GetAncestorLevel(Int_t trackID)    { return mAncestorLevelMap[trackID]; }
        inline Int_t GetAncestorEnergy(Int_t trackID)   { return mAncestorEnergyMap[trackID]; }

    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mMapTTree;
        
        std::map<Int_t, Int_t> mPDGMap;
        std::map<Int_t, Int_t> mParentPDGMap;
        std::map<Int_t, Int_t> mParentTrackIDMap;
        std::map<Int_t, Double_t> mParticleEnergyMap;
        std::map<Int_t, Int_t> mAncestorPDGMap;
        std::map<Int_t, Int_t> mAncestorTrackIDMap;
        std::map<Int_t, Int_t> mAncestorLevelMap;
        std::map<Int_t, Double_t> mAncestorEnergyMap;

    };
}