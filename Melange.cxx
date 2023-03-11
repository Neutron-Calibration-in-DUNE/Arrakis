/**
 * @file Melange.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#include "Melange.h"

namespace arrakis
{
    namespace melange
    {
        Melange::Melange()
        {
            Logger::GetInstance("melange")->trace(
                "setting up melange trees."
            );
            mDetectorPointCloudTree = mTFileService->make<TTree>("det_point_cloud", "det_point_cloud");

            mDetectorView0PointCloudTree = mTFileService->make<TTree>("det_view0_point_cloud", "det_view0_point_cloud");
            mDetectorView0PointCloudTree->Branch("channel", &mDetectorView0PointCloud.channel);
            mDetectorView0PointCloudTree->Branch("tdc",     &mDetectorView0PointCloud.tdc);
            mDetectorView0PointCloudTree->Branch("adc",     &mDetectorView0PointCloud.adc);
            mDetectorView0PointCloudTree->Branch("shape_label",     &mDetectorView0PointCloud.shape_label);
            mDetectorView0PointCloudTree->Branch("particle_label",  &mDetectorView0PointCloud.particle_label);

            mDetectorView1PointCloudTree = mTFileService->make<TTree>("det_view1_point_cloud", "det_view1_point_cloud");
            mDetectorView1PointCloudTree->Branch("channel", &mDetectorView1PointCloud.channel);
            mDetectorView1PointCloudTree->Branch("tdc",     &mDetectorView1PointCloud.tdc);
            mDetectorView1PointCloudTree->Branch("adc",     &mDetectorView1PointCloud.adc);
            mDetectorView1PointCloudTree->Branch("shape_label",     &mDetectorView1PointCloud.shape_label);
            mDetectorView1PointCloudTree->Branch("particle_label",  &mDetectorView1PointCloud.particle_label);

            mDetectorView2PointCloudTree = mTFileService->make<TTree>("det_view2_point_cloud", "det_view2_point_cloud");
            mDetectorView2PointCloudTree->Branch("channel", &mDetectorView2PointCloud.channel);
            mDetectorView2PointCloudTree->Branch("tdc",     &mDetectorView2PointCloud.tdc);
            mDetectorView2PointCloudTree->Branch("adc",     &mDetectorView2PointCloud.adc);
            mDetectorView2PointCloudTree->Branch("shape_label",     &mDetectorView2PointCloud.shape_label);
            mDetectorView2PointCloudTree->Branch("particle_label",  &mDetectorView2PointCloud.particle_label);

            mDetectorView0VoxelTree = mTFileService->make<TTree>("det_view0_voxel", "det_view0_voxel");
            mDetectorView1VoxelTree = mTFileService->make<TTree>("det_view1_voxel", "det_view1_voxel");
            mDetectorView2VoxelTree = mTFileService->make<TTree>("det_view2_voxel", "det_view2_voxel");
        }
        Melange::~Melange()
        {
        }

        void Melange::ResetEvent()
        {
            mDetectorPointCloud.clear();
            mDetectorView0PointCloud.clear();
            mDetectorView1PointCloud.clear();
            mDetectorView2PointCloud.clear();
        }

        void Melange::FillTTree()
        {
            mDetectorView0PointCloudTree->Fill();
            mDetectorView1PointCloudTree->Fill();
            mDetectorView2PointCloudTree->Fill();
        }

        void Melange::ProcessEvent(
            const Parameters& config, art::Event const& event
        )
        {
            ResetEvent();
            PrepareInitialPointClouds(config, event);
            ProcessMuons(config, event);
            ProcessAntiMuons(config, event);
            ProcessPion0s(config, event);
            ProcessPionPlus(config, event);
            ProcessPionMinus(config, event);
            ProcessNeutronCaptures(config, event);
            ProcessAr39(config, event);
            CleanUpPointClouds(config, event);
            SeparatePointClouds(config, event);
            FillTTree();
        }
        
        void Melange::PrepareInitialPointClouds(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto det_sim = mc_data->GetDetectorSimulation();
            auto det_sim_noise = mc_data->GetDetectorSimulationNoise();
            for(size_t ii = 0; ii < det_sim.size(); ii++)
            {
                mDetectorPointCloud.channel.emplace_back(det_sim[ii].channel);
                mDetectorPointCloud.tdc.emplace_back(det_sim[ii].tdc);
                mDetectorPointCloud.adc.emplace_back(det_sim[ii].adc);
                mDetectorPointCloud.view.emplace_back(det_sim[ii].view);
                mDetectorPointCloud.shape_label.emplace_back(LabelCast(ShapeLabel::Undefined));
                mDetectorPointCloud.particle_label.emplace_back(LabelCast(ParticleLabel::Undefined));
            }
            for(size_t ii = 0; ii < det_sim_noise.channel.size(); ii++)
            {
                mDetectorPointCloud.channel.emplace_back(det_sim_noise.channel[ii]);
                mDetectorPointCloud.tdc.emplace_back(det_sim_noise.tdc[ii]);
                mDetectorPointCloud.adc.emplace_back(det_sim_noise.adc[ii]);
                mDetectorPointCloud.view.emplace_back(det_sim_noise.view[ii]);
                mDetectorPointCloud.shape_label.emplace_back(LabelCast(ShapeLabel::Noise));
                mDetectorPointCloud.particle_label.emplace_back(LabelCast(ParticleLabel::Noise));
            }
        }

        void Melange::CleanUpPointClouds(
            const Parameters& config, art::Event const& event
        )
        {
        }

        void Melange::SeparatePointClouds(
            const Parameters& config, art::Event const& event
        )
        {
            for(size_t ii = 0; ii < mDetectorPointCloud.channel.size(); ii++)
            {
                if(mDetectorPointCloud.view[ii] == 0) 
                {
                    mDetectorView0PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
                    mDetectorView0PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
                    mDetectorView0PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
                    mDetectorView0PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
                    mDetectorView0PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
                }
                else if(mDetectorPointCloud.view[ii] == 1) 
                {
                    mDetectorView1PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
                    mDetectorView1PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
                    mDetectorView1PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
                    mDetectorView1PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
                    mDetectorView1PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
                }
                else
                {
                    mDetectorView2PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
                    mDetectorView2PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
                    mDetectorView2PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
                    mDetectorView2PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
                    mDetectorView2PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
                }
            }
        }
        void Melange::ProcessMuons(
            const Parameters& config, art::Event const& event
        )
        {
        }
        void Melange::ProcessAntiMuons(
            const Parameters& config, art::Event const& event
        )
        {
        }
        void Melange::ProcessPion0s(
            const Parameters& config, art::Event const& event
        )
        {
        }
        void Melange::ProcessPionPlus(
            const Parameters& config, art::Event const& event
        )
        {
        }
        void Melange::ProcessPionMinus(
            const Parameters& config, art::Event const& event
        )
        {
        }
        void Melange::ProcessNeutronCaptures(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto neutrons = mc_data->GetParticlesByPDG(2112);
            for(auto neutron : neutrons)
            {
                auto gammas = mc_data->GetProgenyByPDG(neutron, 22);
                auto capture_gammas = mc_data->FilterParticlesByProcess(gammas, ProcessType::NeutronCapture);
                for(auto gamma : capture_gammas)
                {
                    auto gamma_edeps = mc_data->GetParticleAndProgenyEdeps(gamma);
                    PrintIndices("gamma_edeps", gamma_edeps);
                    auto tpc_gamma_edeps = mc_data->FilterEdepsByVolume(gamma_edeps, geometry::VolumeType::TPC);
                    auto tpc_gamma_det_sim = mc_data->GetDetectorSimulationByEdeps(tpc_gamma_edeps);
                    PrintIndices("tpc_gamma_det_sim", tpc_gamma_det_sim);
                    ParticleLabel particle_label = ParticleLabel::NeutronCaptureGammaOther;
                    Double_t gamma_energy = mc_data->RoundParticleEnergy(gamma, 1);
                    if(gamma_energy == 4.7) {
                        particle_label = ParticleLabel::NeutronCaptureGamma475;
                    }
                    else if(gamma_energy == 1.8) {
                        particle_label = ParticleLabel::NeutronCaptureGamma181;
                    }
                    for(auto detsim : tpc_gamma_det_sim) 
                    {
                        mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::NeutronCapture);
                        mDetectorPointCloud.particle_label[detsim] = LabelCast(particle_label);
                    }
                }
            }
        }
        void Melange::ProcessAr39(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto ar39 = mc_data->GetPrimariesByGeneratorLabel(GeneratorLabel::kAr39);
            for(auto elec : ar39)
            {
                auto ar39_edeps = mc_data->GetParticleAndProgenyEdeps(elec);
                auto tpc_ar39_edeps = mc_data->FilterEdepsByVolume(ar39_edeps, geometry::VolumeType::TPC);
                auto tpc_ar39_detsim = mc_data->GetDetectorSimulationByEdeps(tpc_ar39_edeps);
                for(auto detsim : tpc_ar39_detsim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Blip);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::Ar39);
                }
            }
        }
    }
}