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
    PrimaryData::PrimaryData(
        bool SavePrimaryData, bool SavePrimaryDataEdeps,
        bool SavePrimaryDataSimChannel, bool SavePrimaryDataRawTPC
    )
    : mSavePrimaryData(SavePrimaryData)
    , mSavePrimaryDataEdeps(SavePrimaryDataEdeps)
    , mSavePrimaryDataSimChannel(SavePrimaryDataSimChannel)
    , mSavePrimaryDataRawTPC(SavePrimaryDataRawTPC)
    {
        if(SavePrimaryData) {
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

            if(mSavePrimaryDataEdeps)
            {
                mTTree->Branch("edep_energy", &mPrimary.edep_energy);
                mTTree->Branch("edep_x", &mPrimary.edep_x);
                mTTree->Branch("edep_y", &mPrimary.edep_y);
                mTTree->Branch("edep_y", &mPrimary.edep_z);
                
                mTTree->Branch("daughter_edep_ids", &mPrimary.daughter_edep_ids);
                mTTree->Branch("daughter_edep_energy", &mPrimary.daughter_edep_energy);
                mTTree->Branch("daughter_edep_x", &mPrimary.daughter_edep_x);
                mTTree->Branch("daughter_edep_y", &mPrimary.daughter_edep_y);
                mTTree->Branch("daughter_edep_y", &mPrimary.daughter_edep_z);
            }
        }
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
        detinfo::DetectorClocksData const& clockData,
        const std::vector<simb::MCParticle>& mcParticles,
        const std::vector<sim::SimEnergyDeposit>& mcEnergyDeposits,
        const std::vector<sim::SimChannel>& mcChannels,
        const std::vector<raw::RawDigit>& rawTPC
    )
    {
        ResetEvent();

        for (auto particle : mcParticles)
        {
            if(particle.Mother() == 0) 
            {
                mPrimaries.emplace_back(Primary(
                    particle.TrackId(),
                    particle.PdgCode(),
                    particle.Process(),
                    particle.E(),
                    particle.Vx(),
                    particle.Vy(),
                    particle.Vz(),
                    particle.EndProcess(),
                    particle.EndE(),
                    particle.EndX(),
                    particle.EndY(),
                    particle.EndZ()
                ));
            }
            else
            {
                Int_t primary_index = FindPrimary(
                    particle_maps.GetAncestorTrackID(particle.TrackId())
                );
                if(primary_index == -1) {
                    continue;
                }
                mPrimaries[primary_index].AddDaughter(
                    particle.TrackId(),
                    particle_maps.GetAncestorLevel(particle.TrackId()),
                    particle.Process(),
                    particle.E(),
                    particle.Vx(),
                    particle.Vy(),
                    particle.Vz(),
                    particle.EndProcess(),
                    particle.EndE(),
                    particle.EndX(),
                    particle.EndY(),
                    particle.EndZ()
                );
            }
        }
        for(auto edep : mcEnergyDeposits)
        {
            Int_t primary_index = FindPrimary(
                edep.TrackID()
            );
            if(primary_index != -1) {
                mPrimaries[primary_index].AddEdep(
                    edep.Energy(),
                    edep.MidPointX(),
                    edep.MidPointY(),
                    edep.MidPointZ()
                );
            }
            else 
            {
                primary_index = FindPrimary(
                    particle_maps.GetAncestorTrackID(edep.TrackID())
                );
                mPrimaries[primary_index].AddDaughterEdep(
                    edep.TrackID(),
                    edep.Energy(),
                    edep.MidPointX(),
                    edep.MidPointY(),
                    edep.MidPointZ()
                );
            }
        }
        if(mSavePrimaryData)
        {
            for(size_t ii = 0; ii < mPrimaries.size(); ii++)
            {
                mPrimary = mPrimaries[ii];
                mTTree->Fill();
            }
        }
    }
}