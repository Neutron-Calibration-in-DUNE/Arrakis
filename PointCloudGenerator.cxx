/**
 * @file PointCloudGenerator.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12/21
 */
#include "PointCloudGenerator.h"

namespace arrakis
{
    PointCloudGenerator::PointCloudGenerator()
    {
        mSoloPointCloudTTree = mTFileService->make<TTree>("solo_point_cloud", "solo_point_cloud");   
        mSoloPointCloudTTree->Branch("event_id", &mSoloPointCloud.event_id);
        mSoloPointCloudTTree->Branch("point_cloud_id", &mSoloPointCloud.point_cloud_id);
        mSoloPointCloudTTree->Branch("view", &mSoloPointCloud.view);
        mSoloPointCloudTTree->Branch("channel", &mSoloPointCloud.channel);
        mSoloPointCloudTTree->Branch("wire", &mSoloPointCloud.wire);
        mSoloPointCloudTTree->Branch("tdc", &mSoloPointCloud.tdc);
        mSoloPointCloudTTree->Branch("adc", &mSoloPointCloud.adc);
        mSoloPointCloudTTree->Branch("total_energy", &mSoloPointCloud.total_energy);
        mSoloPointCloudTTree->Branch("energy", &mSoloPointCloud.energy);
        mSoloPointCloudTTree->Branch("group_label", &mSoloPointCloud.group_label);
        mSoloPointCloudTTree->Branch("group_label_id", &mSoloPointCloud.group_label_id);
        mSoloPointCloudTTree->Branch("label", &mSoloPointCloud.label);
        mSoloPointCloudTTree->Branch("label_id", &mSoloPointCloud.label_id);
        mSoloPointCloudTTree->Branch("all_deposited", &mSoloPointCloud.all_deposited);
        mSoloPointCloudTTree->Branch("all_lar", &mSoloPointCloud.all_lar);
        mSoloPointCloudTTree->Branch("same_apa", &mSoloPointCloud.same_apa);
        mSoloPointCloudTTree->Branch("lar_edep_fraction", &mSoloPointCloud.lar_edep_fraction);

        mDetectorPointCloudTTree = mTFileService->make<TTree>("detector_point_cloud", "detector_point_cloud");   
        mDetectorPointCloudTTree->Branch("event_id", &mDetectorPointCloud.event_id);
        mDetectorPointCloudTTree->Branch("view", &mDetectorPointCloud.view);
        mDetectorPointCloudTTree->Branch("channel", &mDetectorPointCloud.channel);
        mDetectorPointCloudTTree->Branch("wire", &mDetectorPointCloud.wire);
        mDetectorPointCloudTTree->Branch("tdc", &mDetectorPointCloud.tdc);
        mDetectorPointCloudTTree->Branch("adc", &mDetectorPointCloud.adc);
        mDetectorPointCloudTTree->Branch("total_energy", &mDetectorPointCloud.total_energy);
        mDetectorPointCloudTTree->Branch("energy", &mDetectorPointCloud.energy);
        mDetectorPointCloudTTree->Branch("group_label", &mDetectorPointCloud.group_label);
        mDetectorPointCloudTTree->Branch("group_label_id", &mDetectorPointCloud.group_label_id);
        mDetectorPointCloudTTree->Branch("label", &mDetectorPointCloud.label);
        mDetectorPointCloudTTree->Branch("label_id", &mDetectorPointCloud.label_id);

        mGeneratorLabelNameMap = 
        {
            {kNone, "none"},
            {kAr39, "ar39"},
            {kSingleNeutron, "neutron"},
            {kPNS, "pns"},
            {kNeutronCaptureGamma4_75, "neutron_capture_gamma_4_75"},
            {kNeutronCaptureGamma1_18, "neutron_capture_gamma_1_18"},
        };
        mGeneratorLabelIDMap = 
        {
            {kNone, -1},
            {kAr39, 0},
            {kSingleNeutron, 1},
            {kPNS, 2},
            {kNeutronCaptureGamma4_75, 3},
            {kNeutronCaptureGamma1_18, 4},
        };
    }

    PointCloudGenerator::~PointCloudGenerator()
    {
    }

    void PointCloudGenerator::ProcessEvent(
        ParticleMaps* particle_maps,
        PrimaryData* primary_data
    )
    {
        Logger::GetInstance("solo_point_cloud_generator")->trace(
            "generating point cloud data for " + 
            std::to_string(primary_data->GetPrimaries().size()) + 
            " primaries"
        );
        mPointCloudID = 0;
        for(auto primary : primary_data->GetPrimaries())
        {
            // only process primaries that have tdc values
            if(primary.det_tdc.empty() && primary.daughter_det_tdc.empty()) { 
                continue;
            }
            /**
             * @brief Depending on the type of particle,
             * we will attach labels to the data that is stored from 
             * here.  
             */
            if(primary.generator_label == kAr39) {
                ProcessAr39(primary, particle_maps);
            }
            else if(primary.generator_label == kSingleNeutron) {
                ProcessSingleNeutron(primary, particle_maps);
            }
            else if(primary.generator_label == kPNS) {
                ProcessPNS(primary, particle_maps);
            }
            else {
                ProcessLES(primary, particle_maps);
            }
        }
        auto junk = primary_data->GetJunk();
        if(!junk.det_tdc.empty()) {
            ProcessJunk(junk);
        }
        for(auto point_cloud : mSoloPointClouds)
        {
            mDetectorPointCloud.AddPointCloud(point_cloud);
            mSoloPointCloud = point_cloud;
            mSoloPointCloudTTree->Fill();
        }
        mDetectorPointCloudTTree->Fill();
    }

    void PointCloudGenerator::CollectStatistics(
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

    void PointCloudGenerator::ProcessAr39(
        Primary ar39, ParticleMaps* particle_maps
    )
    {
        SoloPointCloud ar39_point_cloud;
        ar39_point_cloud.point_cloud_id = mPointCloudID;
        ar39_point_cloud.group_label = "ar39";
        ar39_point_cloud.group_label_id = 4;
        ar39_point_cloud.label = mGeneratorLabelNameMap[ar39.generator_label];
        ar39_point_cloud.label_id = mGeneratorLabelIDMap[ar39.generator_label];
        ar39_point_cloud.total_energy = ar39.init_energy;

        // copy views
        ar39_point_cloud.view.insert(ar39_point_cloud.view.end(),ar39.det_view.begin(),ar39.det_view.end());
        ar39_point_cloud.view.insert(ar39_point_cloud.view.end(),ar39.daughter_det_view.begin(), ar39.daughter_det_view.end());
        // copy channels
        ar39_point_cloud.channel.insert(ar39_point_cloud.channel.end(),ar39.det_channel.begin(), ar39.det_channel.end());
        ar39_point_cloud.channel.insert(ar39_point_cloud.channel.end(),ar39.daughter_det_channel.begin(), ar39.daughter_det_channel.end());
        // copy tdc
        ar39_point_cloud.tick.insert(ar39_point_cloud.tick.end(),ar39.det_tick.begin(), ar39.det_tick.end());
        ar39_point_cloud.tick.insert(ar39_point_cloud.tick.end(),ar39.daughter_det_tick.begin(), ar39.daughter_det_tick.end());
        // copy tdc
        ar39_point_cloud.tdc.insert(ar39_point_cloud.tdc.end(),ar39.det_tdc.begin(), ar39.det_tdc.end());
        ar39_point_cloud.tdc.insert(ar39_point_cloud.tdc.end(),ar39.daughter_det_tdc.begin(), ar39.daughter_det_tdc.end());
        // copy adc
        ar39_point_cloud.adc.insert(ar39_point_cloud.adc.end(),ar39.det_adc.begin(), ar39.det_adc.end());
        ar39_point_cloud.adc.insert(ar39_point_cloud.adc.end(),ar39.daughter_det_adc.begin(), ar39.daughter_det_adc.end());

        CollectStatistics(ar39, ar39_point_cloud);
        if(ar39_point_cloud.tdc.size() > 0)
        {
            mSoloPointClouds.emplace_back(ar39_point_cloud);
            mPointCloudID += 1;
        }
    }

    void PointCloudGenerator::ProcessSingleNeutron(
        Primary neutron, ParticleMaps* particle_maps
    )
    {
        /**
         * @brief We are interested in a few different types of
         * neutron interactions.  The main one is nCapture,
         * which produces a certain spectrum of gammas totaling
         * 6.1MeV.  We want to generate a point cloud if there is
         * a capture for the entire spectrum as well as the individual
         * gammas.
         */
        // First, find all the gammas whose process is nCapture
        std::vector<Int_t> capture_gamma_ids;
        std::vector<Double_t> capture_gamma_energy;
        std::vector<std::vector<Int_t>> capture_gamma_daughters;
        for(size_t ii = 0; ii < neutron.daughter_ids.size(); ii++) 
        {
            if(neutron.daughter_init_process[ii] == "nCapture" && abs(neutron.daughter_pdgs[ii]) == 22) 
            {
                capture_gamma_ids.emplace_back(neutron.daughter_ids[ii]);
                capture_gamma_energy.emplace_back(neutron.daughter_init_energy[ii]);
            }
        }

        // Now, for each capture_gamma_id, find all the daughters which
        // originated from them.
        for(size_t ii = 0; ii < capture_gamma_ids.size(); ii++)
        {
            std::vector<Int_t> gamma_daughters;
            for(size_t jj = 0; jj < neutron.daughter_ids.size(); jj++)
            {
                Int_t daughter_parent = particle_maps->GetParentTrackID(neutron.daughter_ids[jj]);
                while(daughter_parent != 0)
                {
                    if(daughter_parent == capture_gamma_ids[ii])
                    {
                        gamma_daughters.emplace_back(neutron.daughter_ids[jj]);
                        break;
                    }
                    Int_t new_parent = particle_maps->GetParentTrackID(daughter_parent);
                    daughter_parent = new_parent;
                }
            }
            capture_gamma_daughters.emplace_back(gamma_daughters);
        }

        // Now, gather all (tdc,channel,adc) values for each of 
        // the gammas.
        std::vector<Int_t> capture_view;
        std::vector<Double_t> capture_channel;
        std::vector<Double_t> capture_wire;
        std::vector<Double_t> capture_tick;
        std::vector<Double_t> capture_adc;
        std::vector<Double_t> capture_tdc;
        Double_t capture_total_energy = 0;
        for(size_t ii = 0; ii < capture_gamma_ids.size(); ii++)
        {
            std::vector<Int_t> gamma_view;
            std::vector<Double_t> gamma_channel;
            std::vector<Double_t> gamma_wire;
            std::vector<Double_t> gamma_tick;
            std::vector<Double_t> gamma_adc;
            std::vector<Double_t> gamma_tdc;
            for(size_t jj = 0; jj < neutron.daughter_det_channel.size(); jj++)
            {
                if(neutron.daughter_det_track_id[jj] == capture_gamma_ids[ii])
                {
                    gamma_view.emplace_back(neutron.daughter_det_view[jj]);
                    gamma_channel.emplace_back(neutron.daughter_det_channel[jj]);
                    gamma_wire.emplace_back(neutron.daughter_det_wire[jj]);
                    gamma_tick.emplace_back(neutron.daughter_det_tick[jj]);
                    gamma_adc.emplace_back(neutron.daughter_det_adc[jj]);
                    gamma_tdc.emplace_back(neutron.daughter_det_tdc[jj]);
                }
                for(size_t kk = 0; kk < capture_gamma_daughters[ii].size(); kk++)
                {
                    if(neutron.daughter_det_track_id[jj] == capture_gamma_daughters[ii][kk])
                    {
                        gamma_view.emplace_back(neutron.daughter_det_view[jj]);
                        gamma_channel.emplace_back(neutron.daughter_det_channel[jj]);
                        gamma_wire.emplace_back(neutron.daughter_det_wire[jj]);
                        gamma_tick.emplace_back(neutron.daughter_det_tick[jj]);
                        gamma_adc.emplace_back(neutron.daughter_det_adc[jj]);
                        gamma_tdc.emplace_back(neutron.daughter_det_tdc[jj]);
                    }
                }
            }
            // If we collected actual points for this gamma,
            // then create a point cloud for it and add
            // the total to the neutron capture.
            if(gamma_channel.size() != 0)
            {
                SoloPointCloud gamma_point_cloud;
                Int_t gamma_energy = std::floor(capture_gamma_energy[ii] * 10e5 + 0.5);
                /**
                 * Group label will be used in BLIP to distinguish between gammas
                 * of different types.  If the gamma is 4.75 or 1.18, then those get 
                 * a special label, otherwise the label is just "gamma_neutron_other".
                 * 
                 */
                if(gamma_energy == 4745) {
                    gamma_point_cloud.group_label = "gamma_neutron_4745";
                    gamma_point_cloud.group_label_id = 1;
                }
                else if(gamma_energy == 1187) {
                    gamma_point_cloud.group_label = "gamma_neutron_1187";
                    gamma_point_cloud.group_label_id = 2;
                }
                else {
                    gamma_point_cloud.group_label = "gamma_neutron_other";
                    gamma_point_cloud.group_label_id = 3;
                }
                gamma_point_cloud.label = "gamma_" + std::to_string(gamma_energy);
                gamma_point_cloud.view = gamma_view;
                gamma_point_cloud.channel = gamma_channel;
                gamma_point_cloud.wire = gamma_wire;
                gamma_point_cloud.tick = gamma_tick;
                gamma_point_cloud.adc = gamma_adc;
                gamma_point_cloud.tdc = gamma_tdc;
                gamma_point_cloud.point_cloud_id = mPointCloudID;
                gamma_point_cloud.total_energy = capture_gamma_energy[ii];

                capture_view.insert(
                    capture_view.end(),
                    gamma_view.begin(),
                    gamma_view.end()
                );
                capture_channel.insert(
                    capture_channel.end(),
                    gamma_channel.begin(),
                    gamma_channel.end()
                );
                capture_wire.insert(
                    capture_wire.end(),
                    gamma_wire.begin(),
                    gamma_wire.end()
                );
                capture_tick.insert(
                    capture_tick.end(),
                    gamma_tick.begin(),
                    gamma_tick.end()
                );
                capture_adc.insert(
                    capture_adc.end(),
                    gamma_adc.begin(),
                    gamma_adc.end()
                );
                capture_tdc.insert(
                    capture_tdc.end(),
                    gamma_tdc.begin(),
                    gamma_tdc.end()
                );
                capture_total_energy += gamma_point_cloud.total_energy;
                mSoloPointClouds.emplace_back(gamma_point_cloud);
                mPointCloudID += 1;
            }
        }
        if(capture_channel.size() != 0 && capture_total_energy >= 0.006)
        {
            SoloPointCloud capture_point_cloud;
            capture_point_cloud.group_label = "neutron";
            capture_point_cloud.group_label_id = 0;
            capture_point_cloud.label = "capture";
            capture_point_cloud.view = capture_view;
            capture_point_cloud.channel = capture_channel;
            capture_point_cloud.wire = capture_wire;
            capture_point_cloud.tick = capture_tick;
            capture_point_cloud.adc = capture_adc;
            capture_point_cloud.tdc = capture_tdc;
            capture_point_cloud.point_cloud_id = mPointCloudID;
            capture_point_cloud.total_energy = capture_total_energy;
            mPointCloudID += 1;
            mSoloPointClouds.emplace_back(capture_point_cloud);
        }
    }

    void PointCloudGenerator::ProcessPNS(
        Primary neutron, ParticleMaps* particle_maps
    )
    {
        SoloPointCloud solo_point_cloud;
        CollectStatistics(neutron, solo_point_cloud);
    }

    void PointCloudGenerator::ProcessLES(
        Primary primary, ParticleMaps* particle_maps
    )
    {
        SoloPointCloud solo_point_cloud;
        CollectStatistics(primary, solo_point_cloud);
    }

    void PointCloudGenerator::ProcessNeutron(
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
    void PointCloudGenerator::ProcessGamma(
        Primary primary, ParticleMaps* particle_maps
    )
    {
    }

    void PointCloudGenerator::ProcessJunk(
        Junk junk
    )
    {
    }
}