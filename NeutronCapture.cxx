/**
 * @file NeutronCapture.cxx
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-09-21
 */
#include "NeutronCapture.h"

namespace arrakis
{
    NeutronCapture::NeutronCapture()
    {
        // mNeutronCaptureTree = mTFileService->make<TTree>("NeutronCapture", "NeutronCapture");
        // mNeutronCaptureTree->Branch("PDGCode", &mNCapture.mPDGCode);
        // mNeutronCaptureTree->Branch("Process", &mNCapture.mProcess);
        // mNeutronCaptureTree->Branch("EndProcess", &mNCapture.mEndProcess);
    }

    NeutronCapture::~NeutronCapture()
    {}

    void NeutronCapture::setBoundingBoxType(std::string volumeType)
    {
        if (volumeType == "TPC" or volumeType == "tpc") { 
            fBoundingBoxType = VolumeType::TPC;
        }
        else if (volumeType == "Cryo" or volumeType == "cryo") {
            fBoundingBoxType = VolumeType::Cryostat;
        }
        else {
            fBoundingBoxType = VolumeType::World;
        }
    }

    // void NeutronCapture::ResetArrays(){
    //     mNCapture.mPDGCode.clear();
    //     mNCapture.mProcess.clear();
    //     mNCapture.mEndProcess.clear();
    // }

    bool NeutronCapture::processEvent(
        ParticleTree particleTree,
        // const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        // ResetArrays();

        if (mcEnergyDeposits.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */

            //Loop over particles
            // for (auto particle : *mcParticles)
            // {
            //     mNCapture.mPDGCode.emplace_back(particle.PdgCode());
            //     mNCapture.mProcess.emplace_back(particle.Process());
            //     mNCapture.mEndProcess.emplace_back(particle.EndProcess());
            // }

            Double_t total_energy = 0;
            bool complete_apa = true;
            bool complete_capture = true;
            bool positive = false;
            bool negative = false;
            //Loop over all the energy deposits
            for (auto energyDeposit : *mcEnergyDeposits)
            {
                if(
                    particleTree.GetAncestorPDG(energyDeposit.TrackID()) == 2112 &&
                    particleTree.GetParentPDG(energyDeposit.TrackID()) != 2112 &&
                    particleTree.GetPDGCode(energyDeposit.TrackID()) == 11
                )
                {
                    if (energyDeposit.StartX() <= 0.0) {
                        negative = true;
                    }
                    else {
                        positive = true;
                    }
                    total_energy += energyDeposit.Energy();
                }
            }
            if (positive && negative) {
                complete_apa = false;
            }
            
            if (round(total_energy) != 6.1) { 
                complete_capture = false;
            }

            return (complete_apa && complete_capture);
        }

        // mNeutronCaptureTree->Fill();
    }
}