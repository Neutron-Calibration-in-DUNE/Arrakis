/**
 * @file SoloPointCloudGenerator.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12/21
 */
#include "SoloPointCloudGenerator.h"

namespace arrakis
{
    SoloPointCloudGenerator::SoloPointCloudGenerator()
    {
        mTTree = mTFileService->make<TTree>("solo_point_cloud", "solo_point_cloud");   
        mTTree->Branch("event_id", &mSoloPointCloud.event_id);
        mTTree->Branch("point_cloud_id", &mSoloPointCloud.point_cloud_id);
        mTTree->Branch("channel", &mSoloPointCloud.channel);
        mTTree->Branch("tdc", &mSoloPointCloud.tdc);
        mTTree->Branch("adc", &mSoloPointCloud.adc);
        mTTree->Branch("energy", &mSoloPointCloud.energy);
        mTTree->Branch("label", &mSoloPointCloud.label);
        mTTree->Branch("all_deposited", &mSoloPointCloud.all_deposited);
        mTTree->Branch("all_lar", &mSoloPointCloud.all_lar);
        mTTree->Branch("same_apa", &mSoloPointCloud.same_apa);
        mTTree->Branch("lar_edep_fraction", &mSoloPointCloud.lar_edep_fraction);
    }

    SoloPointCloudGenerator::~SoloPointCloudGenerator()
    {
    }

    void SoloPointCloudGenerator::ProcessEvent(
        ParticleMaps particle_maps,
        PrimaryData primary_data
    )
    {
        for(auto primary : primary_data)
        {
            /**
             * @brief Depending on the type of particle,
             * we will attach labels to the data that is stored from 
             * here.  
             */
            if(primary.pdg == 2112) {
                ProcessNeutron(primary);
            }
            else if(primary.pdg == 11) {
                ProcessGamma(primary);
            }
            /**
             * @brief Determine the fraction of initial energy which was
             * deposited in the detector, as well as whether it was 
             * deposited all in the active volume LAr, and whether
             * it was deposited in the same TPC.  
             */
            SoloPointCloud solo_point_cloud;
            Double_t edep_energy_fraction = 
                1.0 - (primary.total_edep_energy + primary.total_daughter_edep_energy) / primary.init_energy;
            if(edep_energy_fraction < mEdepEnergyThreshold) {
                solo_point_cloud.all_deposited = true;
            }
            Int_t num_lar_edeps = 0;
            for(auto material : primary.edep_material)
            {
                if(material == "LAr") {
                    num_lar_edeps += 1;
                }
            }
            for(auto material : primary.daughter_edep_material)
            {
                if(material == "LAr") {
                    num_lar_edeps += 1;
                }
            }
            if(num_lar_edeps == (primary.edep_material.size() + primary.daughter_edep_material.size())) {
                solo_point_cloud.all_lar = true;
            }
            solo_point_cloud.lar_edep_fraction = num_lar_edeps / (primary.edep_material.size() + primary.daughter_edep_material.size());
        }
    }

    void SoloPointCloudGenerator::ProcessNeutron(
        PrimaryData neutron
    )
    {
        /**
         * @brief Neutrons produce different types of energy depositions.
         * (1) - elastic scatter
         * (2) - inelastic scatter
         * (3) - capture <-- the one we are most interested in
         * Depending on the type, we will generate different point
         * clouds for the elastic scatters, the total neutron capture,
         * and the individual gammas which come out.
         * For elastic scatters, this requires assuming that the
         * energy deposits made by neutrons are all elastic scatters?
         */

        // First check the elastic scatters

        // Now the gammas
        // Collect the unique gamma ids from the collection
        /**
         * @brief Getting the gammas amounts to finding all
         * of the edeps from particles that are daughters of the gamma,
         * plus the gamma itself and associating those track IDs
         * with the channel/adc values.  Each gamma will have
         * a list of track ids that we will use to build the point cloud.
         */
        std::vector<Int_t> gamma_ids;
        for(auto daughter : neutron.daughter_ids)
        {
            if(particle_maps.GetPDGCode[daughter] == 22) {
                gamma_ids.emplace_back(daughter);
            }
        }


        // Now the entire capture
        /**
         * @brief Once we have all the gamma track id lists, its
         * the same for the neutron except we just add them all up.
         */

    }
    void SoloPointCloudGenerator::ProcessGamma(
        PrimaryData primary
    )
    {
        SoloPointCloud solo_point_cloud;
    }
}