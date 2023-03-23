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
            // Logger::GetInstance("melange")->trace(
            //     "setting up melange trees."
            // );
            // mDetectorPointCloudTree = mTFileService->make<TTree>("det_point_cloud", "det_point_cloud");

            // mDetectorView0PointCloudTree = mTFileService->make<TTree>("det_view0_point_cloud", "det_view0_point_cloud");
            // mDetectorView0PointCloudTree->Branch("channel", &mDetectorView0PointCloud.channel);
            // mDetectorView0PointCloudTree->Branch("tdc",     &mDetectorView0PointCloud.tdc);
            // mDetectorView0PointCloudTree->Branch("adc",     &mDetectorView0PointCloud.adc);
            // mDetectorView0PointCloudTree->Branch("shape_label",     &mDetectorView0PointCloud.shape_label);
            // mDetectorView0PointCloudTree->Branch("particle_label",  &mDetectorView0PointCloud.particle_label);
            // mDetectorView0PointCloudTree->Branch("unique_shape",    &mDetectorView0PointCloud.unique_shape);
            // mDetectorView0PointCloudTree->Branch("unique_particle", &mDetectorView0PointCloud.unique_particle);

            // mDetectorView1PointCloudTree = mTFileService->make<TTree>("det_view1_point_cloud", "det_view1_point_cloud");
            // mDetectorView1PointCloudTree->Branch("channel", &mDetectorView1PointCloud.channel);
            // mDetectorView1PointCloudTree->Branch("tdc",     &mDetectorView1PointCloud.tdc);
            // mDetectorView1PointCloudTree->Branch("adc",     &mDetectorView1PointCloud.adc);
            // mDetectorView1PointCloudTree->Branch("shape_label",     &mDetectorView1PointCloud.shape_label);
            // mDetectorView1PointCloudTree->Branch("particle_label",  &mDetectorView1PointCloud.particle_label);
            // mDetectorView1PointCloudTree->Branch("unique_shape",    &mDetectorView1PointCloud.unique_shape);
            // mDetectorView1PointCloudTree->Branch("unique_particle", &mDetectorView1PointCloud.unique_particle);

            // mDetectorView2PointCloudTree = mTFileService->make<TTree>("det_view2_point_cloud", "det_view2_point_cloud");
            // mDetectorView2PointCloudTree->Branch("channel", &mDetectorView2PointCloud.channel);
            // mDetectorView2PointCloudTree->Branch("tdc",     &mDetectorView2PointCloud.tdc);
            // mDetectorView2PointCloudTree->Branch("adc",     &mDetectorView2PointCloud.adc);
            // mDetectorView2PointCloudTree->Branch("shape_label",     &mDetectorView2PointCloud.shape_label);
            // mDetectorView2PointCloudTree->Branch("particle_label",  &mDetectorView2PointCloud.particle_label);
            // mDetectorView2PointCloudTree->Branch("unique_shape",    &mDetectorView2PointCloud.unique_shape);
            // mDetectorView2PointCloudTree->Branch("unique_particle", &mDetectorView2PointCloud.unique_particle);

            // mDetectorView0VoxelTree = mTFileService->make<TTree>("det_view0_voxel", "det_view0_voxel");
            // mDetectorView1VoxelTree = mTFileService->make<TTree>("det_view1_voxel", "det_view1_voxel");
            // mDetectorView2VoxelTree = mTFileService->make<TTree>("det_view2_voxel", "det_view2_voxel");
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
            // mDetectorPointCloud.clear();
            // mDetectorView0PointCloud.clear();
            // mDetectorView1PointCloud.clear();
            // mDetectorView2PointCloud.clear();
            mShapeLabel = 0;
            mParticleLabel = 0;
        }
        Int_t Melange::IterateShapeLabel()
        {
            mShapeLabel += 1;
            return mShapeLabel;
        }
        Int_t Melange::IterateParticleLabel()
        {
            mParticleLabel += 1;
            return mParticleLabel;
        }
        void Melange::FillTTree()
        {
            // mDetectorView0PointCloudTree->Fill();
            // mDetectorView1PointCloudTree->Fill();
            // mDetectorView2PointCloudTree->Fill();
        }

        void Melange::ProcessEvent(
            const Parameters& config, art::Event const& event
        )
        {
            ResetEvent();
            PrepareInitialPointClouds(config, event);
            ProcessElectrons(config, event);
            ProcessPositrons(config, event);
            ProcessGammas(config, event);
            ProcessMuons(config, event);
            ProcessAntiMuons(config, event);
            ProcessPion0s(config, event);
            ProcessPionPlus(config, event);
            ProcessPionMinus(config, event);
            ProcessProtons(config, event);
            ProcessNeutronCaptures(config, event);
            ProcessNuclearRecoils(config, event);
            ProcessElectronRecoils(config, event);
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
            // auto mc_data = mcdata::MCData::GetInstance();
            // auto det_sim = mc_data->GetDetectorSimulation();
            // auto det_sim_noise = mc_data->GetDetectorSimulationNoise();
            // for(size_t ii = 0; ii < det_sim.size(); ii++)
            // {
            //     mDetectorPointCloud.channel.emplace_back(det_sim[ii].channel);
            //     mDetectorPointCloud.tdc.emplace_back(det_sim[ii].tdc);
            //     mDetectorPointCloud.adc.emplace_back(det_sim[ii].adc);
            //     mDetectorPointCloud.view.emplace_back(det_sim[ii].view);
            //     mDetectorPointCloud.track_id.emplace_back(det_sim[ii].largest_energy_track_id());
            //     mDetectorPointCloud.shape_label.emplace_back(LabelCast(ShapeLabel::Undefined));
            //     mDetectorPointCloud.particle_label.emplace_back(LabelCast(ParticleLabel::Undefined));
            //     mDetectorPointCloud.unique_shape.emplace_back(-1);
            //     mDetectorPointCloud.unique_particle.emplace_back(-1);
            // }
            // for(size_t ii = 0; ii < det_sim_noise.channel.size(); ii++)
            // {
            //     mDetectorPointCloud.channel.emplace_back(det_sim_noise.channel[ii]);
            //     mDetectorPointCloud.tdc.emplace_back(det_sim_noise.tdc[ii]);
            //     mDetectorPointCloud.adc.emplace_back(det_sim_noise.adc[ii]);
            //     mDetectorPointCloud.view.emplace_back(det_sim_noise.view[ii]);
            //     mDetectorPointCloud.track_id.emplace_back(-1);
            //     mDetectorPointCloud.shape_label.emplace_back(LabelCast(ShapeLabel::Noise));
            //     mDetectorPointCloud.particle_label.emplace_back(LabelCast(ParticleLabel::Noise));
            //     mDetectorPointCloud.unique_shape.emplace_back(-1);
            //     mDetectorPointCloud.unique_particle.emplace_back(-1);
            // }
        }

        void Melange::CleanUpPointClouds(
            const Parameters& config, art::Event const& event
        )
        {
            // auto mc_data = mcdata::MCData::GetInstance();
            // for(size_t ii = 0; ii < mDetectorPointCloud.channel.size(); ii++)
            // {
            //     if(mDetectorPointCloud.particle_label[ii] == LabelCast(ParticleLabel::Undefined))
            //     {
            //         Logger::GetInstance("melange")->warning(
            //             "undefined point: " + std::to_string(ii) + " - trackid: " + 
            //             std::to_string(mDetectorPointCloud.track_id[ii]) + " - pdg: " +
            //             std::to_string(mc_data->GetPDGCode_TrackID(mDetectorPointCloud.track_id[ii])) + " - process: " +
            //             std::to_string(ProcessTypeInt(mc_data->GetProcess_TrackID(mDetectorPointCloud.track_id[ii])))
            //         );
            //         // auto ancestry = mc_data->GetAncestryTrackID_TrackID(mDetectorPointCloud.track_id[ii]);
            //         // for(auto ancestor : ancestry)
            //         // {
            //         //     std::cout << "\tancestor: " << ancestor << " - pdg: ";
            //         //     std::cout << mc_data->GetPDGCode_TrackID(ancestor) << " - process: ";
            //         //     std::cout << ProcessTypeInt(mc_data->GetProcess_TrackID(ancestor)) << " - parent: ";
            //         //     std::cout << mc_data->GetParentTrackID_TrackID(ancestor) << std::endl;
            //         // }
            //     }
            // }
        }

        void Melange::SeparatePointClouds(
            const Parameters& config, art::Event const& event
        )
        {
            // for(size_t ii = 0; ii < mDetectorPointCloud.channel.size(); ii++)
            // {
            //     if(mDetectorPointCloud.view[ii] == 0) 
            //     {
            //         mDetectorView0PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
            //         mDetectorView0PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
            //         mDetectorView0PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
            //         mDetectorView0PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
            //         mDetectorView0PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
            //         mDetectorView0PointCloud.unique_shape.emplace_back(mDetectorPointCloud.unique_shape[ii]);
            //         mDetectorView0PointCloud.unique_particle.emplace_back(mDetectorPointCloud.unique_particle[ii]);
            //     }
            //     else if(mDetectorPointCloud.view[ii] == 1) 
            //     {
            //         mDetectorView1PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
            //         mDetectorView1PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
            //         mDetectorView1PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
            //         mDetectorView1PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
            //         mDetectorView1PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
            //         mDetectorView1PointCloud.unique_shape.emplace_back(mDetectorPointCloud.unique_shape[ii]);
            //         mDetectorView1PointCloud.unique_particle.emplace_back(mDetectorPointCloud.unique_particle[ii]);
            //     }
            //     else
            //     {
            //         mDetectorView2PointCloud.channel.emplace_back(mDetectorPointCloud.channel[ii]);
            //         mDetectorView2PointCloud.tdc.emplace_back(mDetectorPointCloud.tdc[ii]);
            //         mDetectorView2PointCloud.adc.emplace_back(mDetectorPointCloud.adc[ii]);
            //         mDetectorView2PointCloud.shape_label.emplace_back(mDetectorPointCloud.shape_label[ii]);
            //         mDetectorView2PointCloud.particle_label.emplace_back(mDetectorPointCloud.particle_label[ii]);
            //         mDetectorView2PointCloud.unique_shape.emplace_back(mDetectorPointCloud.unique_shape[ii]);
            //         mDetectorView2PointCloud.unique_particle.emplace_back(mDetectorPointCloud.unique_particle[ii]);
            //     }
            // }
        }
        void Melange::SetLabels(
            DetSimID_List detSimIDList, 
            ShapeLabel shape, ParticleLabel particle,
            Int_t shape_label, Int_t particle_label
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            for(auto detsim : detSimIDList)
            {
                mc_data->SetWirePlanePointCloudLabels(
                    detsim, LabelCast(shape), LabelCast(particle), shape_label, particle_label
                );

                // mDetectorPointCloud.shape_label[detsim] = LabelCast(shape);
                // mDetectorPointCloud.particle_label[detsim] = LabelCast(particle);
                // mDetectorPointCloud.unique_shape[detsim] = shape_label;
                // mDetectorPointCloud.unique_particle[detsim] = particle_label;
            }
        }
        void Melange::SetLabels(
            DetSimID_Collection detSimIDCollection, 
            ShapeLabel shape, ParticleLabel particle,
            Int_t shape_label, Int_t particle_label
        )
        {
            for(auto detsimIDList : detSimIDCollection) {
                SetLabels(detsimIDList, shape, particle, shape_label, particle_label);
            }
        }
        void Melange::SetLabels(
            std::vector<DetSimID_Collection> detSimIDCollection, 
            ShapeLabel shape, ParticleLabel particle,
            Int_t shape_label, Int_t particle_label
        )
        {
            for(auto detsimIDList : detSimIDCollection) {
                SetLabels(detsimIDList, shape, particle, shape_label, particle_label);
            }
        }
        void Melange::ProcessShowers(TrackID_t trackID, Int_t shapeLabel)
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto particle_det_sim = mc_data->GetDetSimID_TrackID(trackID);
            if(mc_data->GetPDGCode_TrackID(trackID) == 11) {
                SetLabels(particle_det_sim, ShapeLabel::Shower, ParticleLabel::ElectronShower, shapeLabel, IterateParticleLabel());
            }
            else if(mc_data->GetPDGCode_TrackID(trackID) == -11) {
                SetLabels(particle_det_sim, ShapeLabel::Shower, ParticleLabel::PositronShower, shapeLabel, IterateParticleLabel());
            }
            else if(mc_data->GetAbsPDGCode_TrackID(trackID) == 22) {
                SetLabels(particle_det_sim, ShapeLabel::Shower, ParticleLabel::PhotonShower, shapeLabel, IterateParticleLabel());
            }
            auto daughters = mc_data->GetDaughterTrackID_TrackID(trackID);
            auto progeny = mc_data->GetProgenyTrackID_TrackID(trackID);
            auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(daughters, 11);
            auto photon_daughters = mc_data->FilterTrackID_AbsPDGCode(daughters, 22);
            ProcessShowers(elec_daughters, shapeLabel);
            ProcessShowers(photon_daughters, shapeLabel);
        }
        void Melange::ProcessShowers(TrackID_List trackIDList, Int_t shapeLabel)
        {
            for(auto track_id : trackIDList) {
                ProcessShowers(track_id, shapeLabel);
            }
        }
        void Melange::ProcessShowers(TrackID_Collection trackIDCollection)
        {
            for(auto track_id : trackIDCollection) {
                ProcessShowers(track_id, IterateShapeLabel());
            }
        }
        void Melange::ProcessElectrons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto electrons = mc_data->GetPrimaries_PDGCode(11);
            for(auto electron : electrons)
            {
                auto electron_daughters = mc_data->GetDaughterTrackID_TrackID(electron);
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(electron_daughters, 11);
                auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(electron_daughters, 11);
                auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
                auto electron_progeny = mc_data->GetProgenyTrackID_TrackID(electron);
                auto electron_det_sim = mc_data->GetDetSimID_TrackID(electron);
                // Set electron detsim labels to Shower::ElectronShower
                Int_t shower_label = IterateShapeLabel();
                Int_t electron_label = IterateParticleLabel();
                SetLabels(electron_det_sim, ShapeLabel::Shower, ParticleLabel::ElectronShower, shower_label, electron_label);
                SetLabels(elec_det_sim, ShapeLabel::Shower, ParticleLabel::ElectronShower, shower_label, electron_label);            
                ProcessShowers(electron_progeny, shower_label);
                ProcessShowers(other_daughters, IterateShapeLabel());
            }
        }
        void Melange::ProcessPositrons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto positrons = mc_data->GetPrimaries_PDGCode(-11);
            for(auto positron : positrons)
            {
                auto positron_daughters = mc_data->GetDaughterTrackID_TrackID(positron);
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(positron_daughters, 11);
                auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(positron_daughters, 11);
                auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
                auto positron_progeny = mc_data->GetProgenyTrackID_TrackID(positron);
                auto positron_det_sim = mc_data->GetDetSimID_TrackID(positron);
                // Set positron detsim labels to Shower::positronShower
                Int_t shower_label = IterateShapeLabel();
                Int_t positron_label = IterateParticleLabel();
                SetLabels(positron_det_sim, ShapeLabel::Shower, ParticleLabel::PositronShower, shower_label, positron_label);
                SetLabels(elec_det_sim, ShapeLabel::Shower, ParticleLabel::PositronShower, shower_label, positron_label);            
                ProcessShowers(positron_progeny, shower_label);
                ProcessShowers(other_daughters, IterateShapeLabel());
            }
        }
        void Melange::ProcessGammas(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto gammas = mc_data->GetPrimaries_AbsPDGCode(22);
            for(auto gamma : gammas)
            {
                auto gamma_daughters = mc_data->GetDaughterTrackID_TrackID(gamma);
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(gamma_daughters, 11);
                auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
                auto gamma_progeny = mc_data->GetProgenyTrackID_TrackID(gamma);
                auto gamma_det_sim = mc_data->GetDetSimID_TrackID(gamma);
                // Set gamma detsim labels to Shower::PhotonShower
                Int_t shower_label = IterateShapeLabel();
                Int_t gamma_label = IterateParticleLabel();
                SetLabels(gamma_det_sim, ShapeLabel::Shower, ParticleLabel::PhotonShower, shower_label, gamma_label);
                SetLabels(elec_det_sim, ShapeLabel::Shower, ParticleLabel::PhotonShower, shower_label, gamma_label);            
                ProcessShowers(gamma_progeny, shower_label);
            }
        }
        void Melange::ProcessMuons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto muons = mc_data->GetTrackID_PDGCode(13);
            for(auto muon : muons)
            {
                auto muon_daughters = mc_data->GetDaughterTrackID_TrackID(muon);
                auto muon_progeny = mc_data->GetProgenyTrackID_TrackID(muon);
                auto muon_det_sim = mc_data->GetDetSimID_TrackID(muon);
                Int_t muon_label = IterateShapeLabel();
                Int_t particle_label = IterateParticleLabel();
                // Set muon detsim labels to Track:Muon
                SetLabels(muon_det_sim, ShapeLabel::Track, ParticleLabel::Muon, muon_label, particle_label);
                // Filter for Michel and Delta Electrons
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(muon_daughters, 11);
                auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(muon_daughters, 11);
                auto michel_daughters = mc_data->FilterTrackID_Process(elec_daughters, ProcessType::Decay);
                auto delta_daughters = mc_data->FilterTrackID_NotProcess(elec_daughters, ProcessType::Decay);
                auto michel_det_sim = mc_data->GetDetSimID_TrackID(michel_daughters);
                auto delta_det_sim = mc_data->GetDetSimID_TrackID(delta_daughters);
                SetLabels(michel_det_sim, ShapeLabel::Track, ParticleLabel::MichelElectron, IterateShapeLabel(), IterateParticleLabel());
                SetLabels(delta_det_sim, ShapeLabel::Track, ParticleLabel::DeltaElectron, IterateShapeLabel(), IterateParticleLabel());

                auto delta_elec_daughters = mc_data->GetDaughterTrackID_TrackID(delta_daughters);
                auto michel_elec_daughters = mc_data->GetDaughterTrackID_TrackID(michel_daughters);
                // Process progeny as electron and photon showers
                ProcessShowers(delta_elec_daughters);
                ProcessShowers(michel_elec_daughters);
                ProcessShowers(other_daughters, IterateShapeLabel());
                ProcessShowers(muon_progeny, IterateShapeLabel());
            }
        }
        void Melange::ProcessAntiMuons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto muons = mc_data->GetTrackID_PDGCode(-13);
            for(auto muon : muons)
            {
                auto muon_daughters = mc_data->GetDaughterTrackID_TrackID(muon);
                auto muon_progeny = mc_data->GetProgenyTrackID_TrackID(muon);
                auto muon_det_sim = mc_data->GetDetSimID_TrackID(muon);
                Int_t muon_label = IterateShapeLabel();
                Int_t particle_label = IterateParticleLabel();
                // Set muon detsim labels to Track:Muon
                SetLabels(muon_det_sim, ShapeLabel::Track, ParticleLabel::AntiMuon, muon_label, particle_label);
                // Filter for Michel and Delta Electrons
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(muon_daughters, 11);
                auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(muon_daughters, 11);
                auto michel_daughters = mc_data->FilterTrackID_Process(elec_daughters, ProcessType::Decay);
                auto delta_daughters = mc_data->FilterTrackID_NotProcess(elec_daughters, ProcessType::Decay);
                auto michel_det_sim = mc_data->GetDetSimID_TrackID(michel_daughters);
                auto delta_det_sim = mc_data->GetDetSimID_TrackID(delta_daughters);
                SetLabels(michel_det_sim, ShapeLabel::Track, ParticleLabel::MichelElectron, IterateShapeLabel(), IterateParticleLabel());
                SetLabels(delta_det_sim, ShapeLabel::Track, ParticleLabel::DeltaElectron, IterateShapeLabel(), IterateParticleLabel());
                
                auto delta_elec_daughters = mc_data->GetDaughterTrackID_TrackID(delta_daughters);
                auto michel_elec_daughters = mc_data->GetDaughterTrackID_TrackID(michel_daughters);
                // Process progeny as electron and photon showers
                ProcessShowers(delta_elec_daughters);
                ProcessShowers(michel_elec_daughters);
                ProcessShowers(other_daughters, IterateShapeLabel());
                ProcessShowers(muon_progeny, IterateShapeLabel());
            }
        }
        void Melange::ProcessPion0s(
            const Parameters& config, art::Event const& event
        )
        {
            /** (From Wikipedia: https://en.wikipedia.org/wiki/Pion)
             * The dominant π0 decay mode, with a branching ratio of BR2γ = 0.98823, 
             * is into two photons:
             *      π0 → 2γ.
             * The decay π0 → 3γ (as well as decays into any odd number of photons) is 
             * forbidden by the C-symmetry of the electromagnetic interaction: The intrinsic 
             * C-parity of the π0 is +1, while the C-parity of a system of n photons is (−1)n.
             * The second largest π0 decay mode ( BRγee = 0.01174 ) is the Dalitz decay 
             * (named after Richard Dalitz), which is a two-photon decay with an internal 
             * photon conversion resulting a photon and an electron-positron pair in the final state:
             *      π0 → γ + e− + e+.
             * The third largest established decay mode ( BR2e2e = 3.34×10−5 ) is the double-Dalitz 
             * decay, with both photons undergoing internal conversion which leads to further 
             * suppression of the rate:
             *      π0 → e− + e+ + e− +	e+.
             * The fourth largest established decay mode is the loop-induced and therefore suppressed 
             * (and additionally helicity-suppressed) leptonic decay mode ( BRee = 6.46×10−8 ):
             *      π0 → e− + e+.
             * The neutral pion has also been observed to decay into positronium with a branching 
             * fraction on the order of 10−9. No other decay modes have been established experimentally. 
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto pi0s = mc_data->GetTrackID_PDGCode(111);
            for(auto pi0 : pi0s)
            {
                auto pi0_daughters = mc_data->GetDaughterTrackID_TrackID(pi0);
                auto pi0_det_sim = mc_data->GetAllDetSimID_TrackID(pi0);
                SetLabels(pi0_det_sim, ShapeLabel::Shower, ParticleLabel::Pion0, IterateShapeLabel(), IterateParticleLabel());
            }
        }
        void Melange::ProcessPionPlus(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto pipluses = mc_data->GetTrackID_PDGCode(211);
            for(auto piplus : pipluses)
            {
                auto piplus_daughters = mc_data->GetDaughterTrackID_TrackID(piplus);
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(piplus_daughters, 11);
                auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(piplus_daughters, 11);
                auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
                auto piplus_progeny = mc_data->GetProgenyTrackID_TrackID(piplus);
                auto piplus_det_sim = mc_data->GetDetSimID_TrackID(piplus);
                // Set piplus detsim labels to Track:PionMinus
                Int_t track_label = IterateShapeLabel();
                Int_t piplus_label = IterateParticleLabel();
                SetLabels(piplus_det_sim, ShapeLabel::Track, ParticleLabel::PionPlus, track_label, piplus_label);
                SetLabels(elec_det_sim, ShapeLabel::Track, ParticleLabel::PionPlus, track_label, piplus_label);   
                ProcessShowers(other_daughters, IterateShapeLabel());        
                ProcessShowers(piplus_progeny, IterateShapeLabel());
            }
        }
        void Melange::ProcessPionMinus(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto piminuses = mc_data->GetTrackID_PDGCode(-211);
            for(auto piminus : piminuses)
            {
                auto piminus_daughters = mc_data->GetDaughterTrackID_TrackID(piminus);
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(piminus_daughters, 11);
                auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(piminus_daughters, 11);
                auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
                auto piminus_progeny = mc_data->GetProgenyTrackID_TrackID(piminus);
                auto piminus_det_sim = mc_data->GetDetSimID_TrackID(piminus);
                // Set piminus detsim labels to Track:PionMinus
                Int_t track_label = IterateShapeLabel();
                Int_t piminus_label = IterateParticleLabel();
                SetLabels(piminus_det_sim, ShapeLabel::Track, ParticleLabel::PionMinus, track_label, piminus_label);
                SetLabels(elec_det_sim, ShapeLabel::Track, ParticleLabel::PionMinus, track_label, piminus_label);  
                ProcessShowers(other_daughters, IterateShapeLabel());          
                ProcessShowers(piminus_progeny, IterateShapeLabel());
            }
        }
        void Melange::ProcessProtons(
            const Parameters& config, art::Event const& event
        )
        {
            auto mc_data = mcdata::MCData::GetInstance();
            auto protons = mc_data->GetTrackID_PDGCode(2212);
            for(auto proton : protons)
            {
                auto proton_daughters = mc_data->GetDaughterTrackID_TrackID(proton);
                auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(proton_daughters, 11);
                auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(proton_daughters, 11);
                auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
                auto proton_progeny = mc_data->GetProgenyTrackID_TrackID(proton);
                auto proton_det_sim = mc_data->GetDetSimID_TrackID(proton);
                // Set proton detsim labels to Track:PionMinus
                Int_t track_label = IterateShapeLabel();
                Int_t proton_label = IterateParticleLabel();
                SetLabels(proton_det_sim, ShapeLabel::Track, ParticleLabel::Proton, track_label, proton_label);
                SetLabels(elec_det_sim, ShapeLabel::Track, ParticleLabel::Proton, track_label, proton_label); 
                ProcessShowers(other_daughters, IterateShapeLabel());           
                ProcessShowers(proton_progeny, IterateShapeLabel());
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
            auto neutrons = mc_data->GetTrackID_PDGCode(2112);
            auto neutron_daughters = mc_data->GetDaughterTrackID_TrackID(neutrons);
            auto gamma_daughters = mc_data->FilterTrackID_AbsPDGCode(neutron_daughters, 22);
            auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(neutron_daughters, 22);
            auto capture_daughters = mc_data->FilterTrackID_Process(gamma_daughters, ProcessType::NeutronCapture);
            auto other_gammas = mc_data->FilterTrackID_NotProcess(gamma_daughters, ProcessType::NeutronCapture);
            for(auto capture : capture_daughters)
            {
                Int_t neutron_label = IterateShapeLabel();
                for(auto gamma : capture)
                {
                    Double_t gamma_energy = mc_data->GetEnergy_TrackID(gamma, 5);
                    auto gamma_det_sim = mc_data->GetAllDetSimID_TrackID(gamma);
                    Int_t gamma_label = IterateParticleLabel();
                    auto particle_label = ParticleLabel::NeutronCaptureGammaOther;
                    if(gamma_energy == 0.00474 || gamma_energy == 0.00475) {
                        particle_label = ParticleLabel::NeutronCaptureGamma474;
                    }
                    else if(gamma_energy == 0.00336 || gamma_energy == 0.00337) {
                        particle_label = ParticleLabel::NeutronCaptureGamma336;
                    }
                    else if(gamma_energy == 0.00256 || gamma_energy == 0.00257) {
                        particle_label = ParticleLabel::NeutronCaptureGamma256;
                    }
                    else if(gamma_energy == 0.00118 || gamma_energy == 0.00119) {
                        particle_label = ParticleLabel::NeutronCaptureGamma118;
                    }
                    else if(gamma_energy == 0.00083 || gamma_energy == 0.00084) {
                        particle_label = ParticleLabel::NeutronCaptureGamma083;
                    }
                    else if(gamma_energy == 0.00051 || gamma_energy == 0.00052) {
                        particle_label = ParticleLabel::NeutronCaptureGamma051;
                    }
                    else if(gamma_energy == 0.00016 || gamma_energy == 0.00017) {
                        particle_label = ParticleLabel::NeutronCaptureGamma016;
                    }
                    SetLabels(
                        gamma_det_sim, 
                        ShapeLabel::NeutronCapture, particle_label,
                        neutron_label, gamma_label
                    );
                }
            }
            ProcessShowers(other_daughters);
            ProcessShowers(other_gammas);
        }
        void Melange::ProcessNuclearRecoils(
            const Parameters& config, art::Event const& event
        )
        {
            /**
             * Nuclear recoils can come from many things, but are essentially
             * edeps created by Ar40, Ar39 and Ar36 particles (with PDGCodes 
             * 1000180400, 1000180390 and 1000180360).
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto ar40 = mc_data->GetTrackID_PDGCode(1000180400);
            auto ar39 = mc_data->GetTrackID_PDGCode(1000180390);
            auto ar38 = mc_data->GetTrackID_PDGCode(1000180380);
            auto ar36 = mc_data->GetTrackID_PDGCode(1000180360);
            auto ar40_daughters = mc_data->GetDaughterTrackID_TrackID(ar40);
            auto ar39_daughters = mc_data->GetDaughterTrackID_TrackID(ar39);
            auto ar38_daughters = mc_data->GetDaughterTrackID_TrackID(ar38);
            auto ar36_daughters = mc_data->GetDaughterTrackID_TrackID(ar36);
            for(auto ar : ar40)
            {
                auto ar40_det_sim = mc_data->GetDetSimID_TrackID(ar);
                SetLabels(
                    ar40_det_sim,
                    ShapeLabel::Blip, ParticleLabel::NuclearRecoil,
                    IterateShapeLabel(), IterateParticleLabel()
                );
            }
            for(auto ar : ar39)
            {
                auto ar39_det_sim = mc_data->GetDetSimID_TrackID(ar);
                SetLabels(
                    ar39_det_sim,
                    ShapeLabel::Blip, ParticleLabel::NuclearRecoil,
                    IterateShapeLabel(), IterateParticleLabel()
                );
            }
            for(auto ar : ar38)
            {
                auto ar38_det_sim = mc_data->GetDetSimID_TrackID(ar);
                SetLabels(
                    ar38_det_sim,
                    ShapeLabel::Blip, ParticleLabel::NuclearRecoil,
                    IterateShapeLabel(), IterateParticleLabel()
                );
            }
            for(auto ar : ar36)
            {
                auto ar36_det_sim = mc_data->GetDetSimID_TrackID(ar);
                SetLabels(
                    ar36_det_sim,
                    ShapeLabel::Blip, ParticleLabel::NuclearRecoil,
                    IterateShapeLabel(), IterateParticleLabel()
                );
            }
            ProcessShowers(ar40_daughters);
            ProcessShowers(ar39_daughters);
            ProcessShowers(ar38_daughters);
            ProcessShowers(ar36_daughters);
        }
        void Melange::ProcessElectronRecoils(
            const Parameters& config, art::Event const& event
        )
        {
            /**
             * Electron recoils can come from many things, such as
             * edeps created by deuterons/tritons coming out of 
             * neutron inelastic scatters.
             */
            auto mc_data = mcdata::MCData::GetInstance();
            auto deuterons = mc_data->GetTrackID_PDGCode(1000010020);
            auto tritons = mc_data->GetTrackID_PDGCode(1000010030);
            auto deuteron_daughters = mc_data->GetDaughterTrackID_TrackID(deuterons);
            auto triton_daughters = mc_data->GetDaughterTrackID_TrackID(tritons);
            for(auto deuteron : deuterons)
            {
                auto deuteron_det_sim = mc_data->GetDetSimID_TrackID(deuteron);
                SetLabels(
                    deuteron_det_sim,
                    ShapeLabel::Blip, ParticleLabel::ElectronRecoil,
                    IterateShapeLabel(), IterateParticleLabel()
                );
            }
            for(auto triton : tritons)
            {
                auto triton_det_sim = mc_data->GetDetSimID_TrackID(triton);
                SetLabels(
                    triton_det_sim,
                    ShapeLabel::Blip, ParticleLabel::ElectronRecoil,
                    IterateShapeLabel(), IterateParticleLabel()
                );
            }
            ProcessShowers(deuteron_daughters);
            ProcessShowers(triton_daughters);
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
            auto ar39_daughters = mc_data->GetDaughterTrackID_TrackID(ar39);
            SetLabels(ar39_det_sim, ShapeLabel::Blip, ParticleLabel::Ar39, IterateShapeLabel(), IterateParticleLabel());
            ProcessShowers(ar39_daughters);
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
            auto ar42_daughters = mc_data->GetDaughterTrackID_TrackID(ar42);
            SetLabels(ar42_det_sim, ShapeLabel::Blip, ParticleLabel::Ar42, IterateShapeLabel(), IterateParticleLabel());
            ProcessShowers(ar42_daughters);
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
            auto kr85_daughters = mc_data->GetDaughterTrackID_TrackID(kr85);
            SetLabels(kr85_det_sim, ShapeLabel::Blip, ParticleLabel::Kr85, IterateShapeLabel(), IterateParticleLabel());
            ProcessShowers(kr85_daughters);
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
            auto rn222_daughters = mc_data->GetDaughterTrackID_TrackID(rn222);
            SetLabels(rn222_det_sim, ShapeLabel::Blip, ParticleLabel::Rn222, IterateShapeLabel(), IterateParticleLabel());
            ProcessShowers(rn222_daughters);
        }
        void Melange::ProcessCosmics(
            const Parameters& config, art::Event const& event
        )
        {
            // auto mc_data = mcdata::MCData::GetInstance();
            // auto cosmics = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Cosmics);
            // for(auto cosmic : cosmics)
            // {
            // }
        }
    }
}