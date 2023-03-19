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
            mDetectorView0PointCloudTree->Branch("unique_shape",    &mDetectorView0PointCloud.unique_shape);
            mDetectorView0PointCloudTree->Branch("unique_particle", &mDetectorView0PointCloud.unique_particle);

            mDetectorView1PointCloudTree = mTFileService->make<TTree>("det_view1_point_cloud", "det_view1_point_cloud");
            mDetectorView1PointCloudTree->Branch("channel", &mDetectorView1PointCloud.channel);
            mDetectorView1PointCloudTree->Branch("tdc",     &mDetectorView1PointCloud.tdc);
            mDetectorView1PointCloudTree->Branch("adc",     &mDetectorView1PointCloud.adc);
            mDetectorView1PointCloudTree->Branch("shape_label",     &mDetectorView1PointCloud.shape_label);
            mDetectorView1PointCloudTree->Branch("particle_label",  &mDetectorView1PointCloud.particle_label);
            mDetectorView1PointCloudTree->Branch("unique_shape",    &mDetectorView1PointCloud.unique_shape);
            mDetectorView1PointCloudTree->Branch("unique_particle", &mDetectorView1PointCloud.unique_particle);

            mDetectorView2PointCloudTree = mTFileService->make<TTree>("det_view2_point_cloud", "det_view2_point_cloud");
            mDetectorView2PointCloudTree->Branch("channel", &mDetectorView2PointCloud.channel);
            mDetectorView2PointCloudTree->Branch("tdc",     &mDetectorView2PointCloud.tdc);
            mDetectorView2PointCloudTree->Branch("adc",     &mDetectorView2PointCloud.adc);
            mDetectorView2PointCloudTree->Branch("shape_label",     &mDetectorView2PointCloud.shape_label);
            mDetectorView2PointCloudTree->Branch("particle_label",  &mDetectorView2PointCloud.particle_label);
            mDetectorView2PointCloudTree->Branch("unique_shape",    &mDetectorView2PointCloud.unique_shape);
            mDetectorView2PointCloudTree->Branch("unique_particle", &mDetectorView2PointCloud.unique_particle);

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
                mDetectorPointCloud.unique_shape.emplace_back(-1);
                mDetectorPointCloud.unique_particle.emplace_back(-1);
            }
            for(size_t ii = 0; ii < det_sim_noise.channel.size(); ii++)
            {
                mDetectorPointCloud.channel.emplace_back(det_sim_noise.channel[ii]);
                mDetectorPointCloud.tdc.emplace_back(det_sim_noise.tdc[ii]);
                mDetectorPointCloud.adc.emplace_back(det_sim_noise.adc[ii]);
                mDetectorPointCloud.view.emplace_back(det_sim_noise.view[ii]);
                mDetectorPointCloud.shape_label.emplace_back(LabelCast(ShapeLabel::Noise));
                mDetectorPointCloud.particle_label.emplace_back(LabelCast(ParticleLabel::Noise));
                mDetectorPointCloud.unique_shape.emplace_back(-1);
                mDetectorPointCloud.unique_particle.emplace_back(-1);
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
                    mDetectorView0PointCloud.unique_shape.emplace_back(mDetectorPointCloud.unique_shape[ii]);
                    mDetectorView0PointCloud.unique_particle.emplace_back(mDetectorPointCloud.unique_particle[ii]);
                }
                else if(mDetectorPointCloud.view[ii] == 1) 
                {
                    mDetectorView1PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
                    mDetectorView1PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
                    mDetectorView1PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
                    mDetectorView1PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
                    mDetectorView1PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
                    mDetectorView1PointCloud.unique_shape.emplace_back(mDetectorPointCloud.unique_shape[ii]);
                    mDetectorView1PointCloud.unique_particle.emplace_back(mDetectorPointCloud.unique_particle[ii]);
                }
                else
                {
                    mDetectorView2PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
                    mDetectorView2PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
                    mDetectorView2PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
                    mDetectorView2PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
                    mDetectorView2PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
                    mDetectorView2PointCloud.unique_shape.emplace_back(mDetectorPointCloud.unique_shape[ii]);
                    mDetectorView2PointCloud.unique_particle.emplace_back(mDetectorPointCloud.unique_particle[ii]);
                }
            }
        }
        void Melange::SetLabels(
            DetSimID_List detSimIDList, 
            ShapeLabel shape, ParticleLabel particle
        )
        {
            for(auto detsim : detSimIDList)
            {
                mDetectorPointCloud.shape_label[detsim] = LabelCast(shape);
                mDetectorPointCloud.particle_label[detsim] = LabelCast(particle);
            }
        }
        void Melange::SetLabels(
            DetSimID_Collection detSimIDCollection, 
            ShapeLabel shape, ParticleLabel particle
        )
        {
            for(auto detsimIDList : detSimIDCollection) {
                SetLabels(detsimIDList, shape, particle);
            }
        }
        void Melange::SetLabels(
            std::vector<DetSimID_Collection> detSimIDCollection, 
            ShapeLabel shape, ParticleLabel particle
        )
        {
            for(auto detsimIDList : detSimIDCollection) {
                SetLabels(detsimIDList, shape, particle);
            }
        }
        void Melange::ProcessShowers(TrackID_t trackID)
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto particle_det_sim = mc_data->GetDetSimID_TrackID(trackID);
            if(mc_data->GetAbsPDGCode_TrackID(trackID) == 11) {
                SetLabels(particle_det_sim, ShapeLabel::Shower, ParticleLabel::ElectronShower);
            }
            else if(mc_data->GetAbsPDGCode_TrackID(trackID) == 22) {
                SetLabels(particle_det_sim, ShapeLabel::Shower, ParticleLabel::PhotonShower);
            }
            auto daughters = mc_data->GetDaughterTrackID_TrackID(trackID);
            auto progeny = mc_data->GetProgenyTrackID_TrackID(trackID);
            auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(daughters, 11);
            auto photon_daughters = mc_data->FilterTrackID_AbsPDGCode(daughters, 22);
            ProcessShowers(elec_daughters);
            ProcessShowers(photon_daughters);
        }
        void Melange::ProcessShowers(TrackID_List trackIDList)
        {
            for(auto track_id : trackIDList) {
                ProcessShowers(track_id);
            }
        }
        void Melange::ProcessShowers(TrackID_Collection trackIDCollection)
        {
            for(auto track_id : trackIDCollection) {
                ProcessShowers(track_id);
            }
        }
        void Melange::ProcessMuons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto muons = mc_data->GetTrackID_PDGCode(13);
            auto muon_daughters = mc_data->GetDaughterTrackID_TrackID(muons);
            auto muon_progeny = mc_data->GetProgenyTrackID_TrackID(muons);
            auto muon_det_sim = mc_data->GetDetSimID_TrackID(muons);
            // Set muon detsim labels to Track:Muon
            SetLabels(muon_det_sim, ShapeLabel::Track, ParticleLabel::Muon);
            // Filter for Michel and Delta Electrons
            auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(muon_daughters, 11);
            auto michel_daughters = mc_data->FilterTrackID_Process(elec_daughters, ProcessType::Decay);
            auto delta_daughters = mc_data->FilterTrackID_NotProcess(elec_daughters, ProcessType::Decay);
            auto michel_det_sim = mc_data->GetDetSimID_TrackID(michel_daughters);
            auto delta_det_sim = mc_data->GetDetSimID_TrackID(delta_daughters);
            SetLabels(michel_det_sim, ShapeLabel::Track, ParticleLabel::MichelElectron);
            SetLabels(delta_det_sim, ShapeLabel::Track, ParticleLabel::DeltaElectron);
            // Process progeny as electron and photon showers
            ProcessShowers(muon_progeny);
        }
        void Melange::ProcessAntiMuons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto muons = mc_data->GetTrackID_PDGCode(-13);
            auto muon_daughters = mc_data->GetDaughterTrackID_TrackID(muons);
            auto muon_progeny = mc_data->GetProgenyTrackID_TrackID(muons);
            auto muon_det_sim = mc_data->GetDetSimID_TrackID(muons);
            // Set muon detsim labels to Track:Muon
            SetLabels(muon_det_sim, ShapeLabel::Track, ParticleLabel::Muon);
            // Filter for Michel and Delta Electrons
            auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(muon_daughters, 11);
            auto michel_daughters = mc_data->FilterTrackID_Process(elec_daughters, ProcessType::Decay);
            auto delta_daughters = mc_data->FilterTrackID_NotProcess(elec_daughters, ProcessType::Decay);
            auto michel_det_sim = mc_data->GetDetSimID_TrackID(michel_daughters);
            auto delta_det_sim = mc_data->GetDetSimID_TrackID(delta_daughters);
            SetLabels(michel_det_sim, ShapeLabel::Track, ParticleLabel::MichelElectron);
            SetLabels(delta_det_sim, ShapeLabel::Track, ParticleLabel::DeltaElectron);
            // Process progeny as electron and photon showers
            ProcessShowers(muon_progeny);
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
            auto pipluses = mc_data->GetTrackID_PDGCode(211);
            auto piplus_daughters = mc_data->GetDaughterTrackID_TrackID(pipluses);
            auto piplus_progeny = mc_data->GetProgenyTrackID_TrackID(pipluses);
            auto piplus_det_sim = mc_data->GetDetSimID_TrackID(pipluses);
            // Set piplus detsim labels to Track:PionPlus
            SetLabels(piplus_det_sim, ShapeLabel::Track, ParticleLabel::PionPlus);
            // Filter for electron daughters
            auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(piplus_daughters, 11);
            auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
            SetLabels(elec_det_sim, ShapeLabel::Track, ParticleLabel::PionPlus);
            ProcessShowers(piplus_progeny);
        }
        void Melange::ProcessPionMinus(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto pipluses = mc_data->GetTrackID_PDGCode(-211);
            auto piplus_daughters = mc_data->GetDaughterTrackID_TrackID(pipluses);
            auto piplus_progeny = mc_data->GetProgenyTrackID_TrackID(pipluses);
            auto piplus_det_sim = mc_data->GetDetSimID_TrackID(pipluses);
            // Set piplus detsim labels to Track:PionPlus
            SetLabels(piplus_det_sim, ShapeLabel::Track, ParticleLabel::PionPlus);
            // Filter for electron daughters
            auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(piplus_daughters, 11);
            auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
            SetLabels(elec_det_sim, ShapeLabel::Track, ParticleLabel::PionPlus);
            ProcessShowers(piplus_progeny);
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
            auto neutrons = mc_data->GetTrackID_PDGCode(2112);
            auto neutron_daughters = mc_data->GetDaughterTrackID_TrackID(neutrons);
            auto gamma_daughters = mc_data->FilterTrackID_AbsPDGCode(neutron_daughters, 22);
            auto capture_daughters = mc_data->FilterTrackID_Process(gamma_daughters, ProcessType::NeutronCapture);
            for(auto capture : capture_daughters)
            {
                std::cout << "neutron: " << std::endl;
                for(auto gamma : capture)
                {
                    std::cout << "   gamma: " << gamma << std::endl;
                    Double_t gamma_energy = mc_data->GetEnergy_TrackID(gamma, 1000);
                    std::cout << "   energy: " << mc_data->GetEnergy_TrackID(gamma) << std::endl;
                    auto gamma_det_sim = mc_data->GetAllDetSimID_TrackID(gamma);
                    std::cout << "   num: " << gamma_det_sim.size() << std::endl;
                    if(gamma_energy == 0.0047) {
                        SetLabels(gamma_det_sim, ShapeLabel::NeutronCapture, ParticleLabel::NeutronCaptureGamma475);
                    }
                    else if(gamma_energy == 0.0018) {
                        SetLabels(gamma_det_sim, ShapeLabel::NeutronCapture, ParticleLabel::NeutronCaptureGamma181);
                    }
                    else {
                        SetLabels(gamma_det_sim, ShapeLabel::NeutronCapture, ParticleLabel::NeutronCaptureGammaOther);
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
            auto ar39 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Ar39);
            auto ar39_det_sim = mc_data->GetDetSimID_TrackID(ar39);
            SetLabels(ar39_det_sim, ShapeLabel::Blip, ParticleLabel::Ar39);
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
            auto ar42 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Ar42);
            auto ar42_det_sim = mc_data->GetDetSimID_TrackID(ar42);
            SetLabels(ar42_det_sim, ShapeLabel::Blip, ParticleLabel::Ar42);
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
            auto kr85 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Kr85);
            auto kr85_det_sim = mc_data->GetDetSimID_TrackID(kr85);
            SetLabels(kr85_det_sim, ShapeLabel::Blip, ParticleLabel::Kr85);
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
            auto rn222 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Rn222);
            auto rn222_det_sim = mc_data->GetDetSimID_TrackID(rn222);
            SetLabels(rn222_det_sim, ShapeLabel::Blip, ParticleLabel::Rn222);
        }
        void Melange::ProcessCosmics(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto cosmics = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Cosmics);
            for(auto cosmic : cosmics)
            {
                mc_data->PrintParticleData(cosmic);
            }
        }
    }
}