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
        mTTree->Branch("", );
    }

    PrimaryData::~PrimaryData()
    {
    }

    void PrimaryData::ResetEvent()
    {

    }

    void PrimaryData::ProcessEvent(
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
        }
        
        mTTree->Fill();
    }
}