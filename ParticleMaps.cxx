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
            mTTree = mTFileService->make<TTree>("particle_maps", "particle_maps");
            mTTree->Branch("origin", &mOriginMap);
            mTTree->Branch("generator", &mGeneratorMap);
            mTTree->Branch("pdg_map", &mPDGMap);
            mTTree->Branch("parent_pdg_map", &mParentPDGMap);
            mTTree->Branch("parent_track_id_map", &mParentTrackIDMap);
            mTTree->Branch("ancestor_pdg_map", &mAncestorPDGMap);
            mTTree->Branch("ancestor_track_id_map", &mAncestorTrackIDMap);
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
        mOriginMap.clear();
        mPDGMap.clear();
        mParentPDGMap.clear();
        mParentTrackIDMap.clear();
        mAncestorPDGMap.clear();
        mAncestorTrackIDMap.clear();
        mAncestorLevelMap.clear();
        mAncestorEnergyMap.clear();
        mParticleEnergyMap.clear();
    }

    void ParticleMaps::ProcessEvent(
        const art::FindManyP<simb::MCTruth>& mcTruth,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
    )
    {
        ResetEvent();
        for(size_t ii = 0; ii < (*mcParticles).size(); ii++)
        {
            auto const& truth = mcTruth.at(ii);
            const simb::MCGeneratorInfo& generator_info = truth[0]->GeneratorInfo();
            mOriginMap[truth[0]->GetParticle(0).TrackId()] = truth[0]->Origin();
            mGeneratorMap[truth[0]->GetParticle(0).TrackId()] = generator_info.generator;
        }
        for (auto particle : *mcParticles)
        {
            mPDGMap[particle.TrackId()] = particle.PdgCode();
            mParentTrackIDMap[particle.TrackId()] = particle.Mother();
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