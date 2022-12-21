/**
 * @file ParticleMaps.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 */
#include "ParticleMaps.h"

namespace arrakis
{
    ParticleMaps::ParticleMaps(bool SaveParticleMaps)
    : mSaveParticleMaps(SaveParticleMaps)
    {
        if(mSaveParticleMaps) 
        {
            mTTree = mTFileService->make<TTree>("particle_map", "particle_map");
            mTTree->Branch("parent_track_id_map", &mParentTrackIDMap);
            mTTree->Branch("ancestor_track_id_map", &mAncestorTrackIDMap);
            mTTree->Branch("pdg_map", &mPDGMap);
            mTTree->Branch("parent_pdg_map", &mParentPDGMap);
            mTTree->Branch("ancestor_pdg_map", &mAncestorPDGMap);
            mTTree->Branch("ancestor_level_map", &mAncestorLevelMap);
            mTTree->Branch("ancestor_energy_map", &mAncestorEnergyMap);
            mTTree->Branch("parent_energy_map", &mParticleEnergyMap);
        }
    }

    ParticleMaps::~ParticleMaps()
    {
    }

    void ParticleMaps::ResetEvent()
    {
        mParentTrackIDMap.clear();
        mAncestorTrackIDMap.clear();
        mPDGMap.clear();
        mParentPDGMap.clear();
        mAncestorPDGMap.clear();
        mAncestorLevelMap.clear();
        mAncestorEnergyMap.clear();
        mParticleEnergyMap.clear();
    }

    void ParticleMaps::ProcessEvent(
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
    )
    {
        if (!mcParticles.isValid()) {
            return;
        }
        
        ResetEvent();

        for (auto particle : *mcParticles)
        {
            mParentTrackIDMap[particle.TrackId()] = particle.Mother();
            mPDGMap[particle.TrackId()] = particle.PdgCode();

            mParticleEnergyMap[particle.TrackId()] = round(particle.E()*10e6)/10e6;
        }
        for (auto particle : *mcParticles)
        {
            Int_t mother = particle.Mother();
            Int_t track_id = particle.TrackId();
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
                mAncestorEnergyMap[particle.TrackId()] = mParticleEnergyMap[track_id];
            }
            else {
                mParentPDGMap[particle.TrackId()] = mPDGMap[particle.Mother()];
                mAncestorPDGMap[particle.TrackId()] = mPDGMap[track_id];
                mAncestorTrackIDMap[particle.TrackId()] = track_id;
                mAncestorEnergyMap[particle.TrackId()] = mParticleEnergyMap[track_id];
            }
        }
        if(mSaveParticleMaps) {
            mTTree->Fill();
        }
    }
}