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
        mNeutronCaptureTree->Branch("PDGCode", &mNCapture.mPDGCode);
        mNeutronCaptureTree->Branch("Process", &mNCapture.mProcess);
        mNeutronCaptureTree->Branch("EndProcess", &mNCapture.mEndProcess);
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
        mNCapture.mPDGCode.clear();
        mNCapture.mProcess.clear();
        mNCapture.mEndProcess.clear();
    }

    void NeutronCapture::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
        //const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        ResetArrays();

        if (mcParticles.isValid())
        {
            // std::map<int, int> gamma_map;
            // std::vector<int> gamma_map_track_ids;
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            // std::map<Int_t, Int_t> parentDaughterMap;
            // std::map<Int_t, Int_t> particlePDGMap;

            // std::vector<int> neutron_captures;

            //Loop over particles
            for (auto particle : *mcParticles)
            {
                // parentDaughterMap[particle.TrackId()] = particle.Mother();
                // particlePDGMap[particle.TrackId()] = particle.PdgCode();

                mNCapture.mPDGCode.emplace_back(particle.PdgCode());
                mNCapture.mProcess.emplace_back(particle.Process());
                mNCapture.mEndProcess.emplace_back(particle.EndProcess());

                // check if the particle is a neutron
                // if (particle.PdgCode() == 2112)
                // {
                    
                //     if (particle.EndProcess() == "nCapture")
                //     {
                //         DetectorVolume ending_volume = fGeometry->getVolume(
                //             particle.EndX(), particle.EndY(), particle.EndZ()
                //         );
                //         // if the capture is in LAr and inside the cryostat
                //         if (ending_volume.material_name == "LAr" and ending_volume.volume_type == 2) {
                //             neutron_captures.emplace_back(particle.TrackId());
                //         }
                //     }
                // }
            }
        }

        mNeutronCaptureTree->Fill();
    }
}