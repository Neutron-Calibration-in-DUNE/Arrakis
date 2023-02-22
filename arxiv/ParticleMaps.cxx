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
        auto logger = Logger::GetInstance("particle_maps");
        if(mSaveParticleMaps) 
        {
            logger->trace("setting up TTree.");
            mTTree = mTFileService->make<TTree>("particle_maps", "particle_maps");
            mTTree->Branch("generator_label_map", &mGeneratorLabelMap);
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
        mGeneratorLabelMap.clear();
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
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
    )
    {
        ResetEvent();
        Logger::GetInstance("particle_maps")->trace("constructing particle maps for " + std::to_string((*mcParticles).size()) + " particles");
        for (auto particle : *mcParticles)
        {
            mGeneratorLabelMap[particle.TrackId()] = GeneratorLabel::kNone;
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
    }


    void ParticleMaps::ProcessMCTruth(
        GeneratorLabel label,
        const art::ValidHandle<std::vector<simb::MCTruth>>& mcTruth,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
    )
    {
        Logger::GetInstance("particle_maps")->trace(
            "adding MCTruth data to particle maps for label " + 
            std::to_string(label)
        );
        for(auto truth : *mcTruth)
        {
            /**
             * MCTruth stores MCParticles starting with trackID = 0,
             * rather than Geant4 which starts with trackID = 1.
            */
            Logger::GetInstance("particle_maps")->trace(
                "adding labels of type " + 
                std::to_string(label) + 
                " for " + std::to_string(truth.NParticles()) + 
                " particles starting with track ID = " + 
                std::to_string(truth.GetParticle(0).TrackId()+1)
            );
            for(Int_t ii = 0; ii < truth.NParticles(); ii++)
            {
                //Double_t init_t = truth.GetParticle(ii).T();
                Double_t init_x = truth.GetParticle(ii).Vx();
                Double_t init_y = truth.GetParticle(ii).Vy();
                Double_t init_z = truth.GetParticle(ii).Vz();
                //Double_t energy = truth.GetParticle(ii).E();
                Int_t pdg_code = truth.GetParticle(ii).PdgCode();
                if(truth.GetParticle(ii).Process() == "primary")
                {
                    for(auto particle : *mcParticles)
                    {
                        if(
                            //particle.T() == init_t &&
                            particle.Vx() == init_x &&
                            particle.Vy() == init_y &&
                            particle.Vz() == init_z &&
                            //particle.E() == energy &&
                            particle.PdgCode() == pdg_code
                        )
                        {
                            mGeneratorLabelMap[particle.TrackId()] = label;
                        }
                    }
                }
                
            }
        }
    }

    void ParticleMaps::FillTTree()
    {
        if(mSaveParticleMaps) {
            mTTree->Fill();
        }
    }
}