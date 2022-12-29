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
            mTTree->Branch("init_y", &mPrimary.init_z);
            mTTree->Branch("end_process", &mPrimary.end_process);
            mTTree->Branch("end_energy", &mPrimary.end_energy);
            mTTree->Branch("end_t", &mPrimary.end_t);
            mTTree->Branch("end_x", &mPrimary.end_x);
            mTTree->Branch("end_y", &mPrimary.end_y);
            mTTree->Branch("end_y", &mPrimary.end_z);

            mTTree->Branch("daughter_ids", &mPrimary.daughter_ids);
            mTTree->Branch("daughter_level", &mPrimary.daughter_level);
            //mTTree->Branch("daughter_init_process", &mPrimary.daughter_init_process);
            mTTree->Branch("daughter_init_energy", &mPrimary.daughter_init_energy);
            mTTree->Branch("daughter_init_t", &mPrimary.daughter_init_t);
            mTTree->Branch("daughter_init_x", &mPrimary.daughter_init_x);
            mTTree->Branch("daughter_init_y", &mPrimary.daughter_init_y);
            mTTree->Branch("daughter_init_y", &mPrimary.daughter_init_z);
            //mTTree->Branch("daughter_end_process", &mPrimary.daughter_end_process);
            mTTree->Branch("daughter_end_energy", &mPrimary.daughter_end_energy);
            mTTree->Branch("daughter_end_t", &mPrimary.daughter_end_t);
            mTTree->Branch("daughter_end_x", &mPrimary.daughter_end_x);
            mTTree->Branch("daughter_end_y", &mPrimary.daughter_end_y);
            mTTree->Branch("daughter_end_y", &mPrimary.daughter_end_z);

            if(mSavePrimaryDataEdeps)
            {
                mTTree->Branch("total_edep_energy", &mPrimary.total_edep_energy);
                mTTree->Branch("edep_energy", &mPrimary.edep_energy);
                //mTTree->Branch("edep_volume", &mPrimary.edep_volume);
                //mTTree->Branch("edep_material", &mPrimary.edep_material);
                //mTTree->Branch("edep_process", &mPrimary.edep_process);
                mTTree->Branch("edep_t", &mPrimary.edep_t);
                mTTree->Branch("edep_x", &mPrimary.edep_x);
                mTTree->Branch("edep_y", &mPrimary.edep_y);
                mTTree->Branch("edep_y", &mPrimary.edep_z);
                
                mTTree->Branch("total_daughter_edep_energy", &mPrimary.total_daughter_edep_energy);
                mTTree->Branch("daughter_edep_ids", &mPrimary.daughter_edep_ids);
                mTTree->Branch("daughter_edep_energy", &mPrimary.daughter_edep_energy);
                //mTTree->Branch("daughter_edep_volume", &mPrimary.daughter_edep_volume);
                //mTTree->Branch("daughter_edep_material", &mPrimary.daughter_edep_material);
                //mTTree->Branch("daughter_edep_process", &mPrimary.daughter_edep_process);
                mTTree->Branch("daughter_edep_t", &mPrimary.daughter_edep_t);
                mTTree->Branch("daughter_edep_x", &mPrimary.daughter_edep_x);
                mTTree->Branch("daughter_edep_y", &mPrimary.daughter_edep_y);
                mTTree->Branch("daughter_edep_y", &mPrimary.daughter_edep_z);
            }

            if(mSavePrimaryDataRawTPC)
            {
                mTTree->Branch("det_track_id", &mPrimary.det_track_id);
                mTTree->Branch("det_energy_fraction", &mPrimary.det_energy_fraction);
                mTTree->Branch("det_energy", &mPrimary.det_energy);
                mTTree->Branch("det_channel", &mPrimary.det_channel);
                mTTree->Branch("det_tdc", &mPrimary.det_channel);
                mTTree->Branch("det_adc", &mPrimary.det_adc);
                mTTree->Branch("det_edep", &mPrimary.det_edep);
                //mTTree->Branch("det_process", &mPrimary.det_process);
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
    /**
     * Here we are trying to match up the point at which
     * the energy deposition was created with the MC Particle
     * process that caused the energy deposition.  We do this
     * by checking if the energy values and local time of the 
     * event match.
    */
    std::string PrimaryData::FindEnergyDepositionProcess(
        Int_t primary_index, Int_t track_id,
        Double_t energy, Double_t t
    )
    {
        std::string process = "not_found";

        return process;    
    }

    void PrimaryData::FindDetectorProcess(
        detinfo::DetectorClocksData const& clockData,
        Int_t primary_index, Int_t track_id, 
        Double_t energy, unsigned int tdc
    )
    {
        Int_t edep_index = -1;
        std::string process = "not_found";
        if(mPrimaries[primary_index].track_id == track_id)
        {
            for(size_t ii = 0; ii < mPrimaries[primary_index].edep_energy.size(); ii++)
            {
                if(
                    mPrimaries[primary_index].edep_energy[ii] == energy &&
                    clockData.TPCG4Time2TDC(mPrimaries[primary_index].edep_t[ii]) == tdc
                ) 
                {
                    edep_index = ii;
                    process = mPrimaries[primary_index].edep_process[ii];
                    break;
                }
            }
        }
        else
        {
            for(size_t ii = 0; ii < mPrimaries[primary_index].daughter_edep_energy.size(); ii++)
            {
                if(
                    mPrimaries[primary_index].daughter_edep_ids[ii] == track_id &&
                    mPrimaries[primary_index].daughter_edep_energy[ii] == energy &&
                    clockData.TPCG4Time2TDC(mPrimaries[primary_index].daughter_edep_t[ii]) == tdc
                ) 
                {
                    edep_index = ii;
                    process = mPrimaries[primary_index].daughter_edep_process[ii];
                    break;
                }
            }
        }
        mPrimaries[primary_index].det_edep.emplace_back(edep_index);
        mPrimaries[primary_index].det_process.emplace_back(
            process
        );
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
                    particle.TrackId(),
                    particle.PdgCode(),
                    particle.Process(),
                    particle.E(),
                    particle.T(),
                    particle.Vx(),
                    particle.Vy(),
                    particle.Vz(),
                    particle.EndProcess(),
                    particle.EndE(),
                    particle.EndT(),
                    particle.EndX(),
                    particle.EndY(),
                    particle.EndZ()
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
                    particle.TrackId(),
                    particle_maps->GetAncestorLevel(particle.TrackId()),
                    particle.Process(),
                    particle.E(),
                    particle.T(),
                    particle.Vx(),
                    particle.Vy(),
                    particle.Vz(),
                    particle.EndProcess(),
                    particle.EndE(),
                    particle.EndT(),
                    particle.EndX(),
                    particle.EndY(),
                    particle.EndZ()
                );
            }
        }
        // Now loop through the SimEnergyDeposits.
        for(auto edep : *mcEnergyDeposits)
        {
            Int_t primary_index = FindPrimary(
                edep.TrackID()
            );
            // Get the volume information for the energy deposit.
            auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                edep.MidPointX(),
                edep.MidPointY(),
                edep.MidPointZ()
            );
            // If the primary is found, find the process that 
            // caused the energy deposit, and then save the information
            // associated with it.
            if(primary_index != -1) 
            {
                std::string process = FindEnergyDepositionProcess(
                    primary_index, edep.TrackID(),
                    edep.Energy(), edep.Time()
                );
                mPrimaries[primary_index].AddEdep(
                    edep.Energy(),
                    process,
                    volume.volume_name,
                    volume.material_name,
                    edep.Time(),
                    edep.MidPointX(),
                    edep.MidPointY(),
                    edep.MidPointZ()
                );
            }
            else 
            {
                primary_index = FindPrimary(
                    particle_maps->GetAncestorTrackID(edep.TrackID())
                );
                std::string process = FindEnergyDepositionProcess(
                    primary_index, edep.TrackID(),
                    edep.Energy(), edep.Time()
                );
                mPrimaries[primary_index].AddDaughterEdep(
                    edep.TrackID(),
                    edep.Energy(),
                    process,
                    volume.volume_name,
                    volume.material_name,
                    edep.Time(),
                    edep.MidPointX(),
                    edep.MidPointY(),
                    edep.MidPointZ()
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
                auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                if (trackIDs.size() == 0) { 
                    continue; 
                }

                for (auto track : trackIDs)
                {
                    Int_t primary_index = FindPrimary(
                        particle_maps->GetAncestorTrackID(track.trackID)
                    );
                    if(primary_index != -1) 
                    {
                        mPrimaries[primary_index].AddDetectorSimulation(
                            track.trackID,
                            track.energyFrac,
                            track.energy,
                            l,
                            channel,
                            (Int_t) (std::abs(uncompressed[l]))
                        );
                        FindDetectorProcess(
                            clockData, 
                            primary_index, track.trackID, 
                            track.energy, l
                        ); 
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