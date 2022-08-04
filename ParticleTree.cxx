/**
 * @file ParticleTree.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 */
#include "ParticleTree.h"

namespace arrakis
{
    ParticleTree::ParticleTree()
    {
        fMapTTree = fTFileService->make<TTree>("Map_Tree", "Map_Tree");
        fMapTTree->Branch("ParentTrackIDMap", &mParentTrackIDMap);
        fMapTTree->Branch("AncestorTrackIDMap", &mAncestorTrackIDMap);
        fMapTTree->Branch("PDGMap", &mPDGMap);
        fMapTTree->Branch("ParentPDGMap", &mParentPDGMap);
        fMapTTree->Branch("AncestorPDGMap", &mAncestorPDGMap);
        fMapTTree->Branch("AncestorLevelMap", &mAncestorLevelMap);
    }

    ParticleTree::~ParticleTree()
    {
    }

    void ParticleTree::ResetMaps()
    {
        mParentTrackIDMap.clear();
        mAncestorTrackIDMap.clear();
        mPDGMap.clear();
        mParentPDGMap.clear();
        mAncestorPDGMap.clear();
        mAncestorLevelMap.clear();
    }

    void ParticleTree::processEvent(const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles)
    {
        if (!mcParticles.isValid()) {
            return;
        }
        
        ResetMaps();

        for (auto particle : *mcParticles)
        {
            mParentTrackIDMap[particle.TrackId()] = particle.Mother();
            mPDGMap[particle.TrackId()] = particle.PdgCode();
        }
        for (auto particle : *mcParticles)
        {
            Int_t mother = particle.Mother();
            Int_t track_id = particle.TrackId();
            //Int_t prev_track_id = 0;
            Int_t level = 0;
            while (mother != 0)
            {
                level += 1;
                track_id = mother;
                mother = mParentTrackIDMap[track_id];
            }
            mAncestorLevelMap[particle.TrackId()] = level;
            if (level == 0) {
                mParentPDGMap[particle.TrackId()] = 0;
                mAncestorPDGMap[particle.TrackId()] = 0;
                mAncestorTrackIDMap[particle.TrackId()] = 0;
            }
            else {
                mParentPDGMap[particle.TrackId()] = mPDGMap[particle.Mother()];
                mAncestorPDGMap[particle.TrackId()] = mPDGMap[track_id];
                mAncestorTrackIDMap[particle.TrackId()] = track_id;
            }
        }
        fMapTTree->Fill();
    }
}