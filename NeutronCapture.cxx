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
        mNeutronCaptureTree->Branch("PDGCode", &mPDGCode);
        mNeutronCaptureTree->Branch("Process", &mProcess);
        mNeutronCaptureTree->Branch("EndProcess", &mEndProcess);
    }

    NeutronCapture::~NeutronCapture()
    {}


    void NeutronCapture::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
        //const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
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

                mPDGCode.emplace_back(particle.PdgCode());
                mProcess.emplace_back(particle.Process());
                mEndProcess.emplace_back(particle.EndProcess());

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
    }
}