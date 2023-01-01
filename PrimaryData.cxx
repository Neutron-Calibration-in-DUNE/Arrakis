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
        bool SavePrimaryDataRawTPC
    )
    : mSavePrimaryData(SavePrimaryData)
    , mSavePrimaryDataEdeps(SavePrimaryDataEdeps)
    , mSavePrimaryDataRawTPC(SavePrimaryDataRawTPC)
    {
        if(SavePrimaryData) {
            mTTree = mTFileService->make<TTree>("primary_data", "primary_data");
            mTTree->Branch("track_id", &mPrimary.track_id);
            mTTree->Branch("pdg", &mPrimary.pdg);
            mTTree->Branch("init_process", &mPrimary.init_process);
            mTTree->Branch("init_energy", &mPrimary.init_energy);
            mTTree->Branch("init_t", &mPrimary.init_t);
            mTTree->Branch("init_x", &mPrimary.init_x);
            mTTree->Branch("init_y", &mPrimary.init_y);
            mTTree->Branch("init_z", &mPrimary.init_z);
            mTTree->Branch("end_process", &mPrimary.end_process);
            mTTree->Branch("end_energy", &mPrimary.end_energy);
            mTTree->Branch("end_t", &mPrimary.end_t);
            mTTree->Branch("end_x", &mPrimary.end_x);
            mTTree->Branch("end_y", &mPrimary.end_y);
            mTTree->Branch("end_z", &mPrimary.end_z);

            mTTree->Branch("daughter_ids", &mPrimary.daughter_ids);
            mTTree->Branch("daughter_pdgs", &mPrimary.daughter_pdgs);
            mTTree->Branch("daughter_level", &mPrimary.daughter_level);
            mTTree->Branch("daughter_init_process", &mPrimary.daughter_init_process);
            mTTree->Branch("daughter_init_energy", &mPrimary.daughter_init_energy);
            mTTree->Branch("daughter_init_t", &mPrimary.daughter_init_t);
            mTTree->Branch("daughter_init_x", &mPrimary.daughter_init_x);
            mTTree->Branch("daughter_init_y", &mPrimary.daughter_init_y);
            mTTree->Branch("daughter_init_z", &mPrimary.daughter_init_z);
            mTTree->Branch("daughter_end_process", &mPrimary.daughter_end_process);
            mTTree->Branch("daughter_end_energy", &mPrimary.daughter_end_energy);
            mTTree->Branch("daughter_end_t", &mPrimary.daughter_end_t);
            mTTree->Branch("daughter_end_x", &mPrimary.daughter_end_x);
            mTTree->Branch("daughter_end_y", &mPrimary.daughter_end_y);
            mTTree->Branch("daughter_end_z", &mPrimary.daughter_end_z);

            if(mSavePrimaryDataEdeps)
            {
                mTTree->Branch("total_edep_energy", &mPrimary.total_edep_energy);
                mTTree->Branch("edep_energy", &mPrimary.edep_energy);
                mTTree->Branch("edep_volume", &mPrimary.edep_volume);
                mTTree->Branch("edep_material", &mPrimary.edep_material);
                mTTree->Branch("edep_process", &mPrimary.edep_process);
                mTTree->Branch("edep_t", &mPrimary.edep_t);
                mTTree->Branch("edep_x", &mPrimary.edep_x);
                mTTree->Branch("edep_y", &mPrimary.edep_y);
                mTTree->Branch("edep_z", &mPrimary.edep_z);
                
                mTTree->Branch("total_daughter_edep_energy", &mPrimary.total_daughter_edep_energy);
                mTTree->Branch("daughter_edep_ids", &mPrimary.daughter_edep_ids);
                mTTree->Branch("daughter_edep_energy", &mPrimary.daughter_edep_energy);
                mTTree->Branch("daughter_edep_volume", &mPrimary.daughter_edep_volume);
                mTTree->Branch("daughter_edep_material", &mPrimary.daughter_edep_material);
                mTTree->Branch("daughter_edep_process", &mPrimary.daughter_edep_process);
                mTTree->Branch("daughter_edep_t", &mPrimary.daughter_edep_t);
                mTTree->Branch("daughter_edep_x", &mPrimary.daughter_edep_x);
                mTTree->Branch("daughter_edep_y", &mPrimary.daughter_edep_y);
                mTTree->Branch("daughter_edep_z", &mPrimary.daughter_edep_z);
            }

            if(mSavePrimaryDataRawTPC)
            {
                mTTree->Branch("det_energy_fraction", &mPrimary.det_energy_fraction);
                mTTree->Branch("det_energy", &mPrimary.det_energy);
                mTTree->Branch("det_channel", &mPrimary.det_channel);
                mTTree->Branch("det_tick", &mPrimary.det_tick);
                mTTree->Branch("det_tdc", &mPrimary.det_channel);
                mTTree->Branch("det_adc", &mPrimary.det_adc);
                mTTree->Branch("det_edep", &mPrimary.det_edep);
                mTTree->Branch("det_process", &mPrimary.det_process);

                mTTree->Branch("daughter_det_track_id", &mPrimary.daughter_det_track_id);
                mTTree->Branch("daughter_det_energy_fraction", &mPrimary.daughter_det_energy_fraction);
                mTTree->Branch("daughter_det_energy", &mPrimary.daughter_det_energy);
                mTTree->Branch("daughter_det_channel", &mPrimary.daughter_det_channel);
                mTTree->Branch("daughter_det_tick", &mPrimary.daughter_det_tick);
                mTTree->Branch("daughter_det_tdc", &mPrimary.daughter_det_channel);
                mTTree->Branch("daughter_det_adc", &mPrimary.daughter_det_adc);
                mTTree->Branch("daughter_det_edep", &mPrimary.daughter_det_edep);
                mTTree->Branch("daughter_det_process", &mPrimary.daughter_det_process);
            }
        }
    }

    PrimaryData::~PrimaryData()
    {
    }

    void PrintPrimaryEnergyDepositions()
    {
        for(size_t ii = 0; ii < mPrimaries.size(); ii++)
        {
            mPrimaries[ii].PrintPrimaryEnergyDepositions();
        }
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

    void PrimaryData::ProcessEventMC(
        ParticleMaps* particle_maps,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        if (!mcParticles.isValid()) {
            return;
        }
        ResetEvent();

        // Loop through every particle and save
        // information to the primary that was
        // generated.
        for (auto particle : *mcParticles)
        {
            // If the particle is a primary, make
            // a new entry in mPrimaries.
            if(particle.Mother() == 0) 
            {
                mPrimaries.emplace_back(Primary(
                    particle_maps->GetGeneratorLabel(particle.TrackId()),
                    particle
                ));
            }
            // Otherwise, find the associated primary
            // using the ancestor map.
            else
            {
                Int_t primary_index = FindPrimary(
                    particle_maps->GetAncestorTrackID(particle.TrackId())
                );
                if(primary_index == -1) {
                    continue;
                }
                mPrimaries[primary_index].AddDaughter(
                    particle, particle_maps->GetAncestorLevel(particle.TrackId())
                );
            }
        }
        // Now loop through the SimEnergyDeposits.
        for(auto edep : *mcEnergyDeposits)
        {
            Int_t primary_index = FindPrimary(
                edep.TrackID()
            );
            // If the primary is found, find the process that 
            // caused the energy deposit, and then save the information
            // associated with it.
            if(primary_index != -1) 
            {
                mPrimaries[primary_index].AddEdep(
                    edep
                );
            }
            else 
            {
                primary_index = FindPrimary(
                    particle_maps->GetAncestorTrackID(edep.TrackID())
                );
                mPrimaries[primary_index].AddDaughterEdep(
                    edep
                );
            }
        }
    }

    void PrimaryData::ProcessEventDetectorSimulation(
        ParticleMaps* particle_maps,
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<raw::RawDigit>>& rawTPC
    )
    {
        for(auto digit : *rawTPC)
        {
            // Get the channel number for this digit, number of samples,
            // and the pedestal value so that we can uncompress and
            // remove the pedestal.
            raw::ChannelID_t channel = digit.Channel();
            int num_samples = digit.Samples();
            int pedestal = (int)digit.GetPedestal();
            
            // uncompress the digits and remove the pedestal
            std::vector<short> uncompressed(num_samples);
            raw::Uncompress(
                digit.ADCs(), uncompressed, 
                pedestal, digit.Compression()
            );
            for (auto uncomp : uncompressed) {
                uncomp -= pedestal;
            }
            sim::SimChannel truth_channel = (*mcChannels)[channel]; 

            for(int l=0; l < num_samples; l++) 
            {
                auto const& trackIDsAndEnergy = truth_channel.TrackIDsAndEnergies(l, l);
                if (trackIDsAndEnergy.size() == 0) { 
                    continue; 
                }
                Double_t total_energy = 0;
                for(auto track : trackIDsAndEnergy) {
                    total_energy += track.energy;
                }
                for(auto track : trackIDsAndEnergy)
                {
                    Int_t primary_index = FindPrimary(
                        track.trackID
                    );
                    if(primary_index != -1) 
                    {
                        mPrimaries[primary_index].AddPrimaryDetectorSimulation(
                            clockData,
                            track,
                            total_energy,
                            l,
                            channel,
                            (Int_t) (std::abs(uncompressed[l]))
                        );
                    }  
                    else
                    {
                        Int_t primary_index = FindPrimary(
                            particle_maps->GetAncestorTrackID(track.trackID)
                        );
                        if(primary_index != -1) 
                        {
                            mPrimaries[primary_index].AddDaughterDetectorSimulation(
                                clockData,
                                track,
                                total_energy,
                                l,
                                channel,
                                (Int_t) (std::abs(uncompressed[l]))
                            );
                        }  
                    }
                }
            }
        }
    }

    void PrimaryData::FillTTree()
    {
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