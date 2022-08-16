/**
 * @file SingleNeutronCalibration.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-16
 */
#include "SingleNeutronCalibration.h"

namespace arrakis
{
    SingleNeutronCalibration::SingleNeutronCalibration()
    {
        mSingleNeutronTree = mTFileService->make<TTree>("single_neutron", "single_neutron");
        mSingleNeutronTree->Branch("edep_x", &mSingleNeutron.edep_x);
        mSingleNeutronTree->Branch("edep_y", &mSingleNeutron.edep_y);
        mSingleNeutronTree->Branch("edep_z", &mSingleNeutron.edep_z);
        mSingleNeutronTree->Branch("edep_energy", &mSingleNeutron.edep_energy);
        mSingleNeutronTree->Branch("edep_num_electrons", &mSingleNeutron.edep_num_electrons);
        mSingleNeutronTree->Branch("edep_gamma_id", &mSingleNeutron.edep_gamma_id);
        mSingleNeutronTree->Branch("edep_total_energy", &mSingleNeutron.edep_total_energy);
        mSingleNeutronTree->Branch("complete_apa", &mSingleNeutron.complete_apa);
        mSingleNeutronTree->Branch("complete_capture", &mSingleNeutron.complete_capture);
    }

    SingleNeutronCalibration::~SingleNeutronCalibration()
    {
    }

    void SingleNeutronCalibration::ResetSingleNeutron()
    {
        mSingleNeutron.edep_x.clear();
        mSingleNeutron.edep_y.clear();
        mSingleNeutron.edep_z.clear();
        mSingleNeutron.edep_energy.clear();
        mSingleNeutron.edep_num_electrons.clear();
        mSingleNeutron.edep_gamma_id.clear();
        mSingleNeutron.edep_total_energy = 0.0;

        mSingleNeutron.complete_apa = true;
        mSingleNeutron.complete_capture = true;
    }

    void SingleNeutronCalibration::processEvent(
        ParticleTree particleTree,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        ResetSingleNeutron();

        bool complete_apa = true;
        bool positive = false;
        bool negative = false;
        for (auto energyDeposit : *mcEnergyDeposits)
        {
            if(
                particleTree.GetAncestorPDG(energyDeposit.TrackID()) == 2112 &&
                particleTree.GetParentPDG(energyDeposit.TrackID()) != 2112
            )
            {
                mSingleNeutron.edep_x.emplace_back(energyDeposit.StartX());
                if (energyDeposit.StartX() <= 0.0) {
                    negative = true;
                }
                else {
                    positive = true;
                }
                mSingleNeutron.edep_y.emplace_back(energyDeposit.StartY());
                mSingleNeutron.edep_z.emplace_back(energyDeposit.StartZ());
                mSingleNeutron.edep_energy.emplace_back(energyDeposit.Energy());
                mSingleNeutron.edep_num_electrons.emplace_back(energyDeposit.NumElectrons());
                mSingleNeutron.edep_total_energy += energyDeposit.Energy();

                Int_t pdg = particleTree.GetPDGCode(energyDeposit.TrackID());
                Int_t track_id = energyDeposit.TrackID();
                while (pdg != 22)
                {
                    Int_t temp_track_id = particleTree.GetParentTrackID(track_id);
                    pdg = particleTree.GetPDGCode(temp_track_id);
                    track_id = temp_track_id;
                }
                mSingleNeutron.edep_gamma_id.emplace_back(track_id);
            }
        }
        if (positive && negative) {
            complete_apa = false;
        }
        mSingleNeutron.complete_apa = complete_apa;
        if (round(mSingleNeutron.edep_total_energy) != 6.1) { 
            mSingleNeutron.complete_capture = false;
        }
        mSingleNeutronTree->Fill();
    }
}