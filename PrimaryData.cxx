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
        bool SavePrimaryDataRawTPC,
        Double_t ADCThreshold
    )
    : mSavePrimaryData(SavePrimaryData)
    , mSavePrimaryDataEdeps(SavePrimaryDataEdeps)
    , mSavePrimaryDataRawTPC(SavePrimaryDataRawTPC)
    , mADCThreshold(ADCThreshold)
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
                mTTree->Branch("det_view", &mPrimary.det_view);
                mTTree->Branch("det_channel", &mPrimary.det_channel);
                mTTree->Branch("det_tick", &mPrimary.det_tick);
                mTTree->Branch("det_tdc", &mPrimary.det_channel);
                mTTree->Branch("det_adc", &mPrimary.det_adc);
                mTTree->Branch("det_edep", &mPrimary.det_edep);
                mTTree->Branch("det_process", &mPrimary.det_process);

                mTTree->Branch("daughter_det_track_id", &mPrimary.daughter_det_track_id);
                mTTree->Branch("daughter_det_energy_fraction", &mPrimary.daughter_det_energy_fraction);
                mTTree->Branch("daughter_det_energy", &mPrimary.daughter_det_energy);
                mTTree->Branch("daughter_det_view", &mPrimary.daughter_det_view);
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

    void PrimaryData::PrintPrimaryEnergyDepositions()
    {
        for(size_t ii = 0; ii < mPrimaries.size(); ii++)
        {
            mPrimaries[ii].PrintPrimaryEnergyDepositions();
        }
    }
    void PrimaryData::PrintDaughterEnergyDepositions()
    {
        for(size_t ii = 0; ii < mPrimaries.size(); ii++)
        {
            mPrimaries[ii].PrintDaughterEnergyDepositions();
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
            Logger::GetInstance("primary_data")->error("MCParticles handle is not valid!");
            return;
        }
        if (!mcEnergyDeposits.isValid()) {
            Logger::GetInstance("primary_data")->error("SimEnergyDeposits handle is not valid!");
            return;
        }
        ResetEvent();
        Logger::GetInstance("primary_data")->trace("processing " + std::to_string((*mcParticles).size()) + " MCParticles");
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
                if(primary_index == -1) 
                {
                    Logger::GetInstance("primary_data")->warning(
                        "could not find primary with track id " + 
                        std::to_string(particle_maps->GetAncestorTrackID(particle.TrackId())) + 
                        " from particle ancestor with track id " + 
                        std::to_string(particle.TrackId())
                    );
                    continue;
                }
                mPrimaries[primary_index].AddDaughter(
                    particle, particle_maps->GetAncestorLevel(particle.TrackId())
                );
            }
        }
        Logger::GetInstance("primary_data")->trace("processing " + std::to_string((*mcEnergyDeposits).size()) + " SimEnergyDeposits");
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
                if(primary_index == -1) 
                {
                    Logger::GetInstance("primary_data")->warning(
                        "could not find primary with track id " + 
                        std::to_string(particle_maps->GetAncestorTrackID(edep.TrackID())) + 
                        " from energy deposit ancestor with track id " + 
                        std::to_string(edep.TrackID())
                    );
                    continue;
                }
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
        Logger::GetInstance("primary_data")->trace("processing " + std::to_string((*rawTPC).size()) + " raw digits");
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
                /**
                 * Check wether the raw digit at this step passes 
                 * the threshold cut.
                 */
                if(std::abs(uncompressed[l]) < mADCThreshold) {
                    continue;
                }
                auto const& trackIDsAndEnergy = truth_channel.TrackIDsAndEnergies(l, l);
                /**
                 * This step distinguishes noise from true MC particles.
                 * If the input is noise, pass the detector output to a
                 * noise variable, otherwise, attach the output to the
                 * associated primary.
                 */
                if (trackIDsAndEnergy.size() == 0) { 
                    mJunk.AddJunkDetectorSimulation(
                        clockData,
                        l,
                        channel,
                        (Int_t) (std::abs(uncompressed[l]))
                    ); 
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
                        else
                        {
                            Logger::GetInstance("primary_data")->warning(
                                "could not find track id = " + 
                                std::to_string(particle_maps->GetAncestorTrackID(track.trackID)) + 
                                " associated to track in sim channel " + 
                                std::to_string(channel) + 
                                " with track id = " + 
                                std::to_string(track.trackID)
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