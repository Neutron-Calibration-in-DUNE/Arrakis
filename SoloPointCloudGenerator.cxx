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
        mTTree->Branch("label_id", &mSoloPointCloud.label_id);
        mTTree->Branch("all_deposited", &mSoloPointCloud.all_deposited);
        mTTree->Branch("all_lar", &mSoloPointCloud.all_lar);
        mTTree->Branch("same_apa", &mSoloPointCloud.same_apa);
        mTTree->Branch("lar_edep_fraction", &mSoloPointCloud.lar_edep_fraction);

        mGeneratorLabelNameMap = 
        {
            {kNone, "none"},
            {kAr39, "ar39"},
            {kSingleNeutron, "neutron"}
            {kPNS, "pns"},
            {kNeutronCaptureGamma4_75, "neutron_capture_gamma_4_75"},
            {kNeutronCaptureGamma1_18, "neutron_capture_gamma_1_18"},
        };
        mGeneratorLabelIDMap = 
        {
            {kNone, -1},
            {kAr39, 0},
            {kSingleNeutron, 1}
            {kPNS, 2},
            {kNeutronCaptureGamma4_75, 3},
            {kNeutronCaptureGamma1_18, 4},
        };
    }

    SoloPointCloudGenerator::~SoloPointCloudGenerator()
    {
    }

    void SoloPointCloudGenerator::ProcessEvent(
        ParticleMaps* particle_maps,
        PrimaryData* primary_data
    )
    {
        mPointCloudID = 0;
        for(auto primary : primary_data->GetPrimaries())
        {
            /**
             * @brief Depending on the type of particle,
             * we will attach labels to the data that is stored from 
             * here.  
             */
            if(primary.generator_label == kAr39) {
                ProcessAr39(primary, particle_maps);
            }
            else if(primary.generator_label == kPNS) {
                ProcessPNS(primary, particle_maps);
            }
            else {
                ProcessLES(primary, particle_maps);
            }
            mPointCloudID += 1;
        }
    }

    void SoloPointCloudGenerator::CollectStatistics(
        Primary primary, SoloPointCloud& solo_point_cloud
    )
    {
        /**
         * @brief Determine the fraction of initial energy which was
         * deposited in the detector, as well as whether it was 
         * deposited all in the active volume LAr, and whether
         * it was deposited in the same TPC.  
         */
        Double_t edep_energy_fraction = 
            (primary.total_edep_energy + primary.total_daughter_edep_energy) 
            / primary.init_energy;
        if(edep_energy_fraction < mEdepEnergyThreshold) {
            solo_point_cloud.all_deposited = true;
        }
        size_t num_lar_edeps = 0;
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
        if(num_lar_edeps == 0) {
            solo_point_cloud.lar_edep_fraction = 0;
        }
        else {
            solo_point_cloud.lar_edep_fraction = num_lar_edeps / (primary.edep_material.size() + primary.daughter_edep_material.size());
            if(num_lar_edeps == (primary.edep_material.size() + primary.daughter_edep_material.size())) {
                solo_point_cloud.all_lar = true;
            }
        }
    }

    void SoloPointCloudGenerator::ProcessAr39(
        Primary ar39, ParticleMaps* particle_maps
    )
    {
        SoloPointCloud solo_point_cloud;
        solo_point_cloud.point_cloud_id = mPointCloudID;
        std::cout << "HERE" << std::endl;
        solo_point_cloud.label = mGeneratorLabelNameMap[ar39.generator_label];
        solo_point_cloud.label_id = mGeneratorLabelIDMap[ar39.generator_label];
        std::cout << solo_point_cloud.label << std::endl;
        std::cout << solo_point_cloud.label_id << std::endl;

        // copy channels
        solo_point_cloud.channel.insert(
            solo_point_cloud.channel.end(),
            ar39.det_channel.begin(), 
            ar39.det_channel.end()
        );
        solo_point_cloud.channel.insert(
            solo_point_cloud.channel.end(),
            ar39.daughter_det_channel.begin(), 
            ar39.daughter_det_channel.end()
        );
        // copy tdc
        solo_point_cloud.tdc.insert(
            solo_point_cloud.tdc.end(),
            ar39.det_tdc.begin(), 
            ar39.det_tdc.end()
        );
        solo_point_cloud.tdc.insert(
            solo_point_cloud.tdc.end(),
            ar39.daughter_det_tdc.begin(), 
            ar39.daughter_det_tdc.end()
        );
        // copy adc
        solo_point_cloud.adc.insert(
            solo_point_cloud.adc.end(),
            ar39.det_adc.begin(), 
            ar39.det_adc.end()
        );
        solo_point_cloud.adc.insert(
            solo_point_cloud.adc.end(),
            ar39.daughter_det_adc.begin(), 
            ar39.daughter_det_adc.end()
        );
        CollectStatistics(ar39, solo_point_cloud);
        mSoloPointCloud = solo_point_cloud;
        mTTree->Fill();
    }

    void SoloPointCloudGenerator::ProcessPNS(
        Primary neutron, ParticleMaps* particle_maps
    )
    {
        SoloPointCloud solo_point_cloud;
        CollectStatistics(neutron, solo_point_cloud);
        mSoloPointCloud = solo_point_cloud;
        mTTree->Fill();
    }

    void SoloPointCloudGenerator::ProcessLES(
        Primary primary, ParticleMaps* particle_maps
    )
    {
        SoloPointCloud solo_point_cloud;
        CollectStatistics(primary, solo_point_cloud);
        mSoloPointCloud = solo_point_cloud;
        mTTree->Fill();
    }

    void SoloPointCloudGenerator::ProcessNeutron(
        Primary neutron, ParticleMaps* particle_maps
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
            if(particle_maps->GetPDGCode(daughter) == 22) {
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
        Primary primary, ParticleMaps* particle_maps
    )
    {
        SoloPointCloud solo_point_cloud;
    }
}