/**
 * @file PrimaryData.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-21
 */
#include "PrimaryData.h"

namespace arrakis
{
    PrimaryData::PrimaryData()
    {
        mTTree = mTFileService->make<TTree>("primary_data", "primary_data");
        mTTree->Branch("track_id", &mPrimary.track_id);
        mTTree->Branch("pdg", &mPrimary.pdg);
        mTTree->Branch("init_process", &mPrimary.init_process);
        mTTree->Branch("init_energy", &mPrimary.init_energy);
        mTTree->Branch("init_x", &mPrimary.init_x);
        mTTree->Branch("init_y", &mPrimary.init_y);
        mTTree->Branch("init_y", &mPrimary.init_z);
        mTTree->Branch("end_process", &mPrimary.end_process);
        mTTree->Branch("end_energy", &mPrimary.end_energy);
        mTTree->Branch("end_x", &mPrimary.end_x);
        mTTree->Branch("end_y", &mPrimary.end_y);
        mTTree->Branch("end_y", &mPrimary.end_z);
        mTTree->Branch("daughter_ids", &mPrimary.daughter_ids);
        mTTree->Branch("daughter_level", &mPrimary.daughter_level);
        mTTree->Branch("daughter_init_process", &mPrimary.daughter_init_process);
        mTTree->Branch("daughter_init_energy", &mPrimary.daughter_init_energy);
        mTTree->Branch("daughter_init_x", &mPrimary.daughter_init_x);
        mTTree->Branch("daughter_init_y", &mPrimary.daughter_init_y);
        mTTree->Branch("daughter_init_y", &mPrimary.daughter_init_z);
        mTTree->Branch("daughter_end_process", &mPrimary.daughter_end_process);
        mTTree->Branch("daughter_end_energy", &mPrimary.daughter_end_energy);
        mTTree->Branch("daughter_end_x", &mPrimary.daughter_end_x);
        mTTree->Branch("daughter_end_y", &mPrimary.daughter_end_y);
        mTTree->Branch("daughter_end_y", &mPrimary.daughter_end_z);
    }

    PrimaryData::~PrimaryData()
    {
    }

    Int_t PrimaryData::FindPrimary(Int_t track_id)
    {
        for(size_t ii = 0; ii < mPrimaries.size(); ii++)
        {
            if(mPrimaries[ii].track_id == track_id) {
                return ii;
            }
        }
        return -1;
    }

    void PrimaryData::ResetEvent()
    {
        mPrimaries.clear();
    }

    void PrimaryData::ProcessEvent(
        ParticleMaps particle_maps,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::EnergyDeposit>>& mcEnergyDeposits
    )
    {
        if (!mcParticles.isValid()) {
            return;
        }
        ResetEvent();

        for (auto particle : *mcParticles)
        {
            if(particle.Mother() == 0) 
            {

            }
        }
        for(size_t ii = 0; ii < mPrimaries.size(); ii++)
        {
            mPrimary = mPrimaries[ii];
            mTTree->Fill();
        }
    }
}