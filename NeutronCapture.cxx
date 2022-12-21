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
        mNeutronCaptureTree = mTFileService->make<TTree>("NeutronCapture", "NeutronCapture");
        // mNeutronCaptureTree->Branch("PDGCode", &mNCapture.mPDGCode);
        // mNeutronCaptureTree->Branch("Process", &mNCapture.mProcess);
        // mNeutronCaptureTree->Branch("EndProcess", &mNCapture.mEndProcess);
        mNeutronCaptureTree->Branch("complete_apa", &mNCapture.complete_apa);
        mNeutronCaptureTree->Branch("complete_capture", &mNCapture.complete_capture);
        mNeutronCaptureTree->Branch("total_energy", &mNCapture.total_energy);
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

    void NeutronCapture::ResetArrays(){
        // mNCapture.mPDGCode.clear();
        // mNCapture.mProcess.clear();
        // mNCapture.mEndProcess.clear();
        mNCapture.complete_apa.clear();
        mNCapture.complete_capture.clear();
        mNCapture.total_energy.clear();
    }

    bool NeutronCapture::processEvent(
        ParticleMaps particle_maps,
        // const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        ResetArrays();
        Double_t total_energy = 0;
        bool complete_apa = true;
        bool complete_capture = true;
        bool positive = false;
        bool negative = false;
            
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

            //Loop over all the energy deposits
            for (auto energyDeposit : *mcEnergyDeposits)
            {
                if(
                    // particle_maps.GetAncestorPDG(energyDeposit.TrackID()) == 2112 &&
                    // particle_maps.GetParentPDG(energyDeposit.TrackID()) != 2112 &&
                    std::abs( particle_maps.GetPDGCode(energyDeposit.TrackID()) ) == 11 ||
                    particle_maps.GetPDGCode(energyDeposit.TrackID()) == 22
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
            
            if (round(total_energy*10)/10 != 6.1) { 
                complete_capture = false;
            }
        
            mNCapture.complete_apa.emplace_back(complete_apa);
            mNCapture.complete_capture.emplace_back(complete_capture);
            mNCapture.total_energy.emplace_back(total_energy);
            mNeutronCaptureTree->Fill();
        }

        bool storeEvent = (complete_apa && complete_capture);
        return storeEvent;
    }
}