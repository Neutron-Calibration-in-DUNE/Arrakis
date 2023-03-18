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
        std::map<std::string, FilterDetectorSimulation> FilterDetectorSimulationMap = 
        {
            {"TrackID", FilterDetectorSimulation::TrackID},
            {"EdepID",  FilterDetectorSimulation::EdepID}
        };
        std::map<std::string, NeutronCaptureGammaDetail> NeutronCaptureGammaDetailMap = 
        {
            {"Simple", NeutronCaptureGammaDetail::Simple},
            {"Medium", NeutronCaptureGammaDetail::Medium},
            {"Full", NeutronCaptureGammaDetail::Full}
        };

        Melange* Melange::sInstance{nullptr};
        std::mutex Melange::sMutex;

        Melange *Melange::GetInstance()
        {
            std::lock_guard<std::mutex> lock(sMutex);
            if (sInstance == nullptr)
            {
                sInstance = new Melange();
            }
            return sInstance;
        }
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
            mDetectorView0PointCloudTree->Branch("unique_label",  &mDetectorView0PointCloud.unique_label);

            mDetectorView1PointCloudTree = mTFileService->make<TTree>("det_view1_point_cloud", "det_view1_point_cloud");
            mDetectorView1PointCloudTree->Branch("channel", &mDetectorView1PointCloud.channel);
            mDetectorView1PointCloudTree->Branch("tdc",     &mDetectorView1PointCloud.tdc);
            mDetectorView1PointCloudTree->Branch("adc",     &mDetectorView1PointCloud.adc);
            mDetectorView1PointCloudTree->Branch("shape_label",     &mDetectorView1PointCloud.shape_label);
            mDetectorView1PointCloudTree->Branch("particle_label",  &mDetectorView1PointCloud.particle_label);
            mDetectorView1PointCloudTree->Branch("unique_label",  &mDetectorView1PointCloud.unique_label);

            mDetectorView2PointCloudTree = mTFileService->make<TTree>("det_view2_point_cloud", "det_view2_point_cloud");
            mDetectorView2PointCloudTree->Branch("channel", &mDetectorView2PointCloud.channel);
            mDetectorView2PointCloudTree->Branch("tdc",     &mDetectorView2PointCloud.tdc);
            mDetectorView2PointCloudTree->Branch("adc",     &mDetectorView2PointCloud.adc);
            mDetectorView2PointCloudTree->Branch("shape_label",     &mDetectorView2PointCloud.shape_label);
            mDetectorView2PointCloudTree->Branch("particle_label",  &mDetectorView2PointCloud.particle_label);
            mDetectorView2PointCloudTree->Branch("unique_label",  &mDetectorView2PointCloud.unique_label);

            mDetectorView0VoxelTree = mTFileService->make<TTree>("det_view0_voxel", "det_view0_voxel");
            mDetectorView1VoxelTree = mTFileService->make<TTree>("det_view1_voxel", "det_view1_voxel");
            mDetectorView2VoxelTree = mTFileService->make<TTree>("det_view2_voxel", "det_view2_voxel");
        }

        void Melange::SetConfigurationParameters(const Parameters& config)
        {
            Logger::GetInstance("melange")->trace(
                "setting up configuration parameters."
            );
            fhicl::ParameterSet const& melange_params = config().melange_parameters.get_PSet();
            std::map<std::string, std::string> melange_tags;
            for(std::string const& name : melange_params.get_names()) {
                melange_tags[name] = melange_params.get<std::string>(name);
            }
            sFilterDetectorSimulation = FilterDetectorSimulationMap[melange_tags["FilterDetectorSimulation"]];
            Logger::GetInstance("melange")->trace(
                "setting filter detector simulation to: " + melange_tags["FilterDetectorSimulation"]
            );
            sNeutronCaptureGammaDetail = NeutronCaptureGammaDetailMap[melange_tags["NeutronCaptureGammaDetail"]];
            Logger::GetInstance("melange")->trace(
                "setting neutron capture gamma detail to: " + melange_tags["NeutronCaptureGammaDetail"]
            );
        }

        void Melange::ResetEvent()
        {
            mDetectorPointCloud.clear();
            mDetectorView0PointCloud.clear();
            mDetectorView1PointCloud.clear();
            mDetectorView2PointCloud.clear();
            mClusterLabel = 0;
        }
        Int_t Melange::IterateClusterLabel()
        {
            mClusterLabel += 1;
            return mClusterLabel;
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
            ProcessAr42(config, event);
            ProcessKr85(config, event);
            ProcessRn222(config, event);
            ProcessCosmics(config, event);
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
                mDetectorPointCloud.unique_label.emplace_back(-1);
            }
            for(size_t ii = 0; ii < det_sim_noise.channel.size(); ii++)
            {
                mDetectorPointCloud.channel.emplace_back(det_sim_noise.channel[ii]);
                mDetectorPointCloud.tdc.emplace_back(det_sim_noise.tdc[ii]);
                mDetectorPointCloud.adc.emplace_back(det_sim_noise.adc[ii]);
                mDetectorPointCloud.view.emplace_back(det_sim_noise.view[ii]);
                mDetectorPointCloud.shape_label.emplace_back(LabelCast(ShapeLabel::Noise));
                mDetectorPointCloud.particle_label.emplace_back(LabelCast(ParticleLabel::Noise));
                mDetectorPointCloud.unique_label.emplace_back(-1);
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
                    mDetectorView0PointCloud.unique_label.emplace_back(mDetectorPointCloud.unique_label[ii]);
                }
                else if(mDetectorPointCloud.view[ii] == 1) 
                {
                    mDetectorView1PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
                    mDetectorView1PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
                    mDetectorView1PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
                    mDetectorView1PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
                    mDetectorView1PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
                    mDetectorView1PointCloud.unique_label.emplace_back(mDetectorPointCloud.unique_label[ii]);
                }
                else
                {
                    mDetectorView2PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
                    mDetectorView2PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
                    mDetectorView2PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
                    mDetectorView2PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
                    mDetectorView2PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
                    mDetectorView2PointCloud.unique_label.emplace_back(mDetectorPointCloud.unique_label[ii]);
                }
            }
        }
        void Melange::ProcessShowers(Int_t trackID)
        {
            auto mc_data = mcdata::MCData::GetInstance();
            std::vector<Int_t> particle_det_sim;
            if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(trackID, geometry::VolumeType::TPC);
            }
            else {
                particle_det_sim = mc_data->GetDetSimID_TrackID(trackID);
            }
            for(auto detsim : particle_det_sim)
            {
                mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Shower);
                if(std::abs(mc_data->GetPDGCode_TrackID(trackID)) == 11) {
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::ElectronShower);
                }
                else if(std::abs(mc_data->GetPDGCode_TrackID(trackID)) == 22) {
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::PhotonShower);
                }
            }
            auto progeny = mc_data->GetProgenyTrackID_TrackID(trackID);
            for(auto particle : progeny)
            {
                if(mc_data->GetPDGCode_TrackID(particle) == 11 || mc_data->GetPDGCode_TrackID(particle) == 22) {
                    ProcessShowers(particle);
                }
            }
        }
        void Melange::ProcessMuons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto muons = mc_data->GetParticlesByPDG(13);
            std::vector<Int_t> particle_det_sim;
            for(auto muon : muons)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(muon, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(muon);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::Muon);
                }
                auto muon_progeny = mc_data->GetProgenyTrackID_TrackID(muon);
                /**
                 * Now go through and grab all electrons which are direct descendants of
                 * muons.  These can be classified into several types, depending on how
                 * they were generated.  Electrons which are direct descendants (daughters)
                 * of muons are technically delta-rays, and should be labelled as tracks.
                 * Unless of course the electron comes from a particular decay, then
                 * it is a special electron called a Michel electron.
                 */
                for(auto particle : muon_progeny)
                {
                    if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                        particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(particle, geometry::VolumeType::TPC);
                    }
                    else {
                        particle_det_sim = mc_data->GetDetSimID_TrackID(particle);
                    }
                    if(std::abs(mc_data->GetPDGCode_TrackID(particle)) == 11 && mc_data->GetParentTrackID_TrackID(particle) == muon)
                    {
                        for(auto detsim : particle_det_sim)
                        {
                            mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                            if(mc_data->GetProcess_TrackID(particle) == ProcessType::Decay) {
                                mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::MichelElectron);
                            }
                            else {
                                mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::DeltaElectron);    
                            }
                        }
                    }
                    else if(
                        mc_data->GetParentTrackID_TrackID(particle) != muon && 
                        (std::abs(mc_data->GetPDGCode_TrackID(particle)) == 11 || std::abs(mc_data->GetPDGCode_TrackID(particle)) == 22)
                    )
                    {
                        ProcessShowers(particle);
                    }
                    else
                    {
                    }
                }
            }
        }
        void Melange::ProcessAntiMuons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto muons = mc_data->GetParticlesByPDG(-13);
            std::vector<Int_t> particle_det_sim;
            for(auto muon : muons)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(muon, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(muon);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::AntiMuon);
                }
                auto muon_progeny = mc_data->GetProgenyTrackID_TrackID(muon);
                /**
                 * Now go through and grab all electrons which are direct descendants of
                 * muons.  These can be classified into several types, depending on how
                 * they were generated.  Electrons which are direct descendants (daughters)
                 * of muons are technically delta-rays, and should be labelled as tracks.
                 * Unless of course the electron comes from a particular decay, then
                 * it is a special electron called a Michel electron.
                 */
                for(auto particle : muon_progeny)
                {
                    if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                        particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(particle, geometry::VolumeType::TPC);
                    }
                    else {
                        particle_det_sim = mc_data->GetDetSimID_TrackID(particle);
                    }
                    if(std::abs(mc_data->GetPDGCode_TrackID(particle)) == 11 && mc_data->GetParentTrackID_TrackID(particle) == muon)
                    {
                        for(auto detsim : particle_det_sim)
                        {
                            mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                            if(mc_data->GetProcess_TrackID(particle) == ProcessType::Decay) {
                                mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::MichelElectron);
                            }
                            else {
                                mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::DeltaElectron);    
                            }
                        }
                    }
                    else if(
                        mc_data->GetParentTrackID_TrackID(particle) != muon && 
                        (std::abs(mc_data->GetPDGCode_TrackID(particle)) == 11 || std::abs(mc_data->GetPDGCode_TrackID(particle)) == 22)
                    )
                    {
                        ProcessShowers(particle);
                    }
                    else
                    {
                    }
                }
            }
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
            auto mc_data = mcdata::MCData::GetInstance();
            auto pipluses = mc_data->GetParticlesByPDG(211);
            std::vector<Int_t> particle_det_sim;
            for(auto piplus : pipluses)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(piplus, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(piplus);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::PionPlus);
                }
                auto piplus_progeny = mc_data->GetProgenyTrackID_TrackID(piplus);
                /**
                 * Now go through and grab all electrons which are direct descendants of
                 * pipluses.
                 */
                for(auto particle : piplus_progeny)
                {
                    if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                        particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(particle, geometry::VolumeType::TPC);
                    }
                    else {
                        particle_det_sim = mc_data->GetDetSimID_TrackID(particle);
                    }
                    if(std::abs(mc_data->GetPDGCode_TrackID(particle)) == 11 && mc_data->GetParentTrackID_TrackID(particle) == piplus)
                    {
                        for(auto detsim : particle_det_sim)
                        {
                            mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                            mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::PionPlus);    
                        }
                    }
                    else {
                        ProcessShowers(particle);
                    }
                }
            }
        }
        void Melange::ProcessPionMinus(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto piminuses = mc_data->GetParticlesByPDG(-211);
            std::vector<Int_t> particle_det_sim;
            for(auto piminus : piminuses)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(piminus, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(piminus);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::PionMinus);
                }
                auto piminus_progeny = mc_data->GetProgenyTrackID_TrackID(piminus);
                /**
                 * Now go through and grab all electrons which are direct descendants of
                 * piminuses.
                 */
                for(auto particle : piminus_progeny)
                {
                    if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                        particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(particle, geometry::VolumeType::TPC);
                    }
                    else {
                        particle_det_sim = mc_data->GetDetSimID_TrackID(particle);
                    }
                    if(std::abs(mc_data->GetPDGCode_TrackID(particle)) == 11 && mc_data->GetParentTrackID_TrackID(particle) == piminus)
                    {
                        for(auto detsim : particle_det_sim)
                        {
                            mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Track);
                            mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::PionMinus);    
                        }
                    }
                    else {
                        ProcessShowers(particle);
                    }
                }
            }
        }
        void Melange::ProcessNeutronCaptures(
            const Parameters& config, art::Event const& event
        )
        {
            /**
             * Neutron captures produce a standard candle of 6.1 MeV
             * gammas, which are generated according to a cascade: 
             * 
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto neutrons = mc_data->GetParticlesByPDG(2112);
            std::vector<Int_t> particle_det_sim;
            for(auto neutron : neutrons)
            {
                auto gammas = mc_data->GetProgenyByPDG(neutron, 22);
                auto capture_gammas = mc_data->FilterParticlesByProcess(gammas, ProcessType::NeutronCapture);
                for(auto gamma : capture_gammas)
                {
                    if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                        particle_det_sim = mc_data->GetDetectorSimulationByParticleAndProgenyVolume(gamma, geometry::VolumeType::TPC);
                    }
                    else {
                        particle_det_sim = mc_data->GetDetectorSimulationByParticleAndProgeny(gamma);
                    }
                    ParticleLabel particle_label = ParticleLabel::NeutronCaptureGammaOther;
                    Double_t gamma_energy = mc_data->GetEnergyTrackID(gamma, 1);
                    if(gamma_energy == 4.7) {
                        particle_label = ParticleLabel::NeutronCaptureGamma475;
                    }
                    else if(gamma_energy == 1.8) {
                        particle_label = ParticleLabel::NeutronCaptureGamma181;
                    }
                    for(auto detsim : particle_det_sim) 
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
            /**
             * Argon-39 decays via beta decay into Potassium-39,
             * with a Q-value of 565 keV: http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=180039.
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto ar39 = mc_data->GetPrimariesByGeneratorLabel(GeneratorLabel::Ar39);
            std::vector<Int_t> particle_det_sim;
            for(auto elec : ar39)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(elec, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(elec);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Blip);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::Ar39);
                }
            }
        }
        void Melange::ProcessAr42(
            const Parameters& config, art::Event const& event
        )
        {
            /**
             * Argon-42 decays via beta decay into Potassium-42,
             * with a Q-value of 599 keV: http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=180042.
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto ar42 = mc_data->GetPrimariesByGeneratorLabel(GeneratorLabel::Ar42);
            std::vector<Int_t> particle_det_sim;
            for(auto elec : ar42)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(elec, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(elec);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Blip);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::Ar42);
                }
            }
        }
        void Melange::ProcessKr85(
            const Parameters& config, art::Event const& event
        )
        {
            /**
             * Krypton-85 decays via beta decay into Rubidium 85
             * with two prominent betas with energies of 687 keV (99.56 %) and
             * 173 keV (.43 %): http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=360085. 
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto kr85 = mc_data->GetPrimariesByGeneratorLabel(GeneratorLabel::Kr85);
            std::vector<Int_t> particle_det_sim;
            for(auto elec : kr85)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(elec, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(elec);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Blip);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::Kr85);
                }
            }
        }
        void Melange::ProcessRn222(
            const Parameters& config, art::Event const& event
        )
        {
            /**
             * Radon-222 decays via alpha decay through a chain that ends in lead
             * (https://en.wikipedia.org/wiki/Radon-222).  The alpha has an energy of
             * 5.5904 MeV, which bounces around locally in Argon, but quickly thermalizes
             * due to the short scattering length.  The CSDA range of a 5.5 MeV alpha in 
             * Argon is about 7.5e-3 g/cm^2.  Using a density of 1.3954 g/cm^3, the 
             * scattering length is (~0.005 cm) or (~50 um).
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto rn222 = mc_data->GetPrimariesByGeneratorLabel(GeneratorLabel::Rn222);
            std::vector<Int_t> particle_det_sim;
            for(auto alpha : rn222)
            {
                if(sFilterDetectorSimulation == FilterDetectorSimulation::EdepID) {
                    particle_det_sim = mc_data->GetDetectorSimulationByParticleVolume(alpha, geometry::VolumeType::TPC);
                }
                else {
                    particle_det_sim = mc_data->GetDetSimID_TrackID(alpha);
                }
                for(auto detsim : particle_det_sim)
                {
                    mDetectorPointCloud.shape_label[detsim] = LabelCast(ShapeLabel::Blip);
                    mDetectorPointCloud.particle_label[detsim] = LabelCast(ParticleLabel::Rn222);
                }
            }
        }
        void Melange::ProcessCosmics(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto cosmics = mc_data->GetPrimariesByGeneratorLabel(GeneratorLabel::Cosmics);
            for(auto cosmic : cosmics)
            {
                mc_data->PrintParticleData(cosmic);
            }
        }
    }
}