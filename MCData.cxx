/**
 * @file MCData.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#include "MCData.h"

namespace arrakis
{
    std::map<ProcessType, std::string> ProcessTypeToString
    {
        {ProcessType::NotDefined,           "NotDefined"},
        {ProcessType::Unknown,              "Unknown"},
        {ProcessType::Primary,              "Primary"},
        {ProcessType::HadronElastic,        "HadronElastic"},
        {ProcessType::PiMinusInelastic,     "PiMinusInelastic"},
        {ProcessType::PiPlusInelastic,      "PiPlusInelastic"},
        {ProcessType::KaonMinusInelastic,   "KaonMinusInelastic"},
        {ProcessType::KaonPlusInelastic,    "KaonPlusInelastic"},
        {ProcessType::ProtonInelastic,      "ProtonInelastic"},
        {ProcessType::NeutronInelastic,     "NeutronInelastic"},
        {ProcessType::CoulombScatter,       "CoulombScatter"},
        {ProcessType::NeutronCapture,       "NeutronCapture"},
        {ProcessType::Transportation,       "Transportation"},
    };
    std::map<std::string, ProcessType> StringToProcessType
    {
        {"NotDefined",          ProcessType::NotDefined},
        {"Unknown",             ProcessType::Unknown},
        {"Primary",             ProcessType::Primary},
        {"HadronElastic",       ProcessType::HadronElastic},
        {"PiMinusInelastic",    ProcessType::PiMinusInelastic},
        {"PiPlusInelastic",     ProcessType::PiPlusInelastic},
        {"KaonMinusInelastic",  ProcessType::KaonMinusInelastic},
        {"KaonPlusInelastic",   ProcessType::KaonPlusInelastic},
        {"ProtonInelastic",     ProcessType::ProtonInelastic},
        {"NeutronInelastic",    ProcessType::NeutronInelastic},
        {"CoulombScatter",      ProcessType::CoulombScatter},
        {"NeutronCapture",      ProcessType::NeutronCapture},
        {"Transportation",      ProcessType::Transportation},
    };
    std::map<std::string, ProcessType> TrajectoryStringToProcessType
    {
        {"NotDefined",      ProcessType::NotDefined},
        {"Unknown",         ProcessType::Unknown},
        {"primary",         ProcessType::Primary},
        {"hadElastic",      ProcessType::HadronElastic},
        {"pi-Inelastic",    ProcessType::PiMinusInelastic},
        {"pi+Inelastic",    ProcessType::PiPlusInelastic},
        {"kaon-Inelastic",  ProcessType::KaonMinusInelastic},
        {"kaon+Inelastic",  ProcessType::KaonPlusInelastic},
        {"protonInelastic", ProcessType::ProtonInelastic},
        {"neutronInelastic",ProcessType::NeutronInelastic},
        {"CoulombScatter",  ProcessType::CoulombScatter},
        {"nCapture",        ProcessType::NeutronCapture},
        {"Transportation",  ProcessType::Transportation},
    };
    std::map<ProcessType, std::string> TrajectoryProcessTypeToString
    {
        {ProcessType::NotDefined,           "NotDefined"},
        {ProcessType::Unknown,              "Unknown"},
        {ProcessType::Primary,              "primary"},
        {ProcessType::HadronElastic,        "hadElastic"},
        {ProcessType::PiMinusInelastic,     "pi-Inelastic"},
        {ProcessType::PiPlusInelastic,      "pi+Inelastic"},
        {ProcessType::KaonMinusInelastic,   "kaon-Inelastic"},
        {ProcessType::KaonPlusInelastic,    "kaon+Inelastic"},
        {ProcessType::ProtonInelastic,      "protonInelastic"},
        {ProcessType::NeutronInelastic,     "neutronInelastic"},
        {ProcessType::CoulombScatter,       "CoulombScatter"},
        {ProcessType::NeutronCapture,       "nCapture"},
        {ProcessType::Transportation,       "Transportation"},
    };

    namespace mcdata
    {
        MCData* MCData::sInstance{nullptr};
        std::mutex MCData::sMutex;

        MCData *MCData::GetInstance()
        {
            std::lock_guard<std::mutex> lock(sMutex);
            if (sInstance == nullptr)
            {
                sInstance = new MCData();
            }
            return sInstance;
        }
        MCData::MCData()
        {
            Logger::GetInstance("mcdata")->trace(
                "setting up mcdata tree."
            );
            sMCDataTree = sTFileService->make<TTree>("mcdata", "mcdata");
            sMCDataTree->Branch("generator_label_map", &sGeneratorLabelMap);
            sMCDataTree->Branch("generator_map", &sGeneratorMap);
            sMCDataTree->Branch("pdg_map", &sPDGMap);
            sMCDataTree->Branch("parent_pdg_map", &sParentPDGMap);
            sMCDataTree->Branch("parent_track_id_map", &sParentTrackIDMap);
            sMCDataTree->Branch("particle_energy_map", &sParticleEnergyMap);
            sMCDataTree->Branch("ancestor_pdg_map", &sAncestorPDGMap);
            sMCDataTree->Branch("ancestor_track_id_map", &sAncestorTrackIDMap);
            sMCDataTree->Branch("ancestor_level_map", &sAncestorLevelMap);
            sMCDataTree->Branch("ancestor_energy_map", &sAncestorEnergyMap);
            sMCDataTree->Branch("daughter_map", &sDaughterMap);
            sMCDataTree->Branch("progeny_map", &sProgenyMap);
            sMCDataTree->Branch("ancestry_map", &sAncestryMap);
            sMCDataTree->Branch("particle_edep_map", &sParticleEdepMap);
            sMCDataTree->Branch("particle_edep_process_map", &sParticleEdepProcessMap);
            sMCDataTree->Branch("particle_detsim_map", &sParticleDetectorSimulationMap);
            sMCDataTree->Branch("random_detsim_map", &sRandomDetectorSimulationMap);
            sMCDataTree->Branch("edep_process_map", &sEdepProcessMap);
            sMCDataTree->Branch("edep_detsim_map", &sEdepDetectorSimulationMap);
            sMCDataTree->Branch("detsim_edep_map", &sDetectorSimulationEdepMap);
  
            sGeneratorMap["Ar39Label"] = GeneratorLabel::Ar39;
            sGeneratorMap["Ar42Label"] = GeneratorLabel::Ar42;
            sGeneratorMap["Kr85Label"] = GeneratorLabel::Kr85;
            sGeneratorMap["Rn222Label"] = GeneratorLabel::Rn222;
            sGeneratorMap["CosmicsLabel"] = GeneratorLabel::Cosmics;
            sGeneratorMap["SingleNeutronLabel"] = GeneratorLabel::SingleNeutron;
            sGeneratorMap["PNSLabel"] = GeneratorLabel::PNS;
        }
        void MCData::SetADCThreshold(Double_t ADCThreshold)
        {
            Logger::GetInstance("mcdata")->trace(
                "setting ADC threshold to " + std::to_string(ADCThreshold)
            );
            sADCThreshold = ADCThreshold;
        }
        void MCData::ResetEvent()
        {
            sMCTruthHandles.clear();
            sPrimaries.clear();
            sDetectorSimulation.clear();
            sDetectorSimulationNoise.clear();

            sGeneratorLabelMap.clear();
            sPDGMap.clear();
            sParentPDGMap.clear();
            sParentTrackIDMap.clear();
            sParticleEnergyMap.clear();
            sAncestorPDGMap.clear();
            sAncestorTrackIDMap.clear();
            sAncestorLevelMap.clear();
            sAncestorEnergyMap.clear();
            sDaughterMap.clear();
            sProgenyMap.clear();
            sAncestryMap.clear();
            sParticleEdepMap.clear();
            sParticleEdepProcessMap.clear();
            sParticleDetectorSimulationMap.clear();
            sRandomDetectorSimulationMap.clear();

            sEdepProcessMap.clear();
            sEdepDetectorSimulationMap.clear();

            sDetectorSimulationEdepMap.clear();
        }
        void MCData::PrintParticleData(TrackID_t trackID)
        {
            auto particle = (*sMCParticleHandle)[sParticleMap[trackID]];
            std::cout << "## MCParticle #################################\n";
            std::cout << "## TrackID:                 [" << std::setw(10) << std::setfill('.') << trackID << "]\n";
            std::cout << "## PDG:                     [" << std::setw(10) << std::setfill('.') << particle.PdgCode() << "]\n";
            std::cout << "## Energy [MeV]:            [" << std::setw(10) << std::setfill('.') << particle.E() << "]\n";
            std::cout << "## Process:                 [" << std::setw(10) << std::setfill('.') << particle.Process() << "]\n";
            std::cout << "## Parent TrackID:          [" << std::setw(10) << std::setfill('.') << particle.Mother() << "]\n";
            std::cout << "## Parent PDG:              [" << std::setw(10) << std::setfill('.') << sParentPDGMap[trackID] << "]\n";
            std::cout << "## Ancestor TrackID:        [" << std::setw(10) << std::setfill('.') << sAncestorTrackIDMap[trackID] << "]\n";
            std::cout << "## Ancestor PDG:            [" << std::setw(10) << std::setfill('.') << sAncestorPDGMap[trackID] << "]\n";
            std::cout << "## Ancestor level:          [" << std::setw(10) << std::setfill('.') << sAncestorLevelMap[trackID] << "]\n";
            std::cout << "## Progeny  [level, trackID, PDG]\n";
            auto progeny = sProgenyMap[trackID];
            auto particle_level = sAncestorLevelMap[trackID];
            for(auto progeny_track_id : progeny)
            {
                std::cout << "##         [" << sAncestorLevelMap[progeny_track_id] - particle_level << ", ";
                std::cout << progeny_track_id << ", " << sPDGMap[progeny_track_id] << "]\n";
            }
            std::cout << "###############################################" << std::endl;
        }
        void MCData::PrintEdepData(EdepID_t edepID)
        {
            auto edep = (*sMCSimEnergyDepositHandle)[edepID];
            std::cout << "## MCSimEnergyDeposit #########################\n";
            std::cout << "## TrackID:                " << edep.TrackID() << "\n";
            std::cout << "## PDG:                    " << sPDGMap[edep.TrackID()] << "\n";
            std::cout << "## Energy:                 " << edep.Energy() << "\n";
            std::cout << "## MidPoint [x, y, z]:     [" << edep.MidPointX() << ", " << edep.MidPointY() << ", " << edep.MidPointZ() << "]\n"; 
            std::cout << "###############################################" << std::endl;
        }
        void MCData::PrintDetSimData(DetSimID_t detsimID)
        {

        }
        void MCData::ProcessEvent(
            const Parameters& config, art::Event const& event
        )
        {
            ResetEvent();
            ProcessMCTruth(event, config().labels.get_PSet());
            ProcessMCParticles(event, config().LArGeantProducerLabel());
            ProcessSimEnergyDeposits(event, 
                config().IonAndScintProducerLabel(), config().IonAndScintInstanceLabel()
            );
            ProcessSimChannels(event, 
                config().SimChannelProducerLabel(), config().SimChannelInstanceLabel()
            );
            ProcessRawDigits(event,
                config().RawDigitProducerLabel(), config().RawDigitInstanceLabel()
            );
        }
        void MCData::ProcessMCTruth(
            art::Event const& event, fhicl::ParameterSet const& generator_labels
        )
        {
            std::map<std::string, art::InputTag> generator_tags;
            for(std::string const& name : generator_labels.get_names()) {
                generator_tags[name] = generator_labels.get<art::InputTag>(name);
            }
            for(auto const& [key, tag] : generator_tags)
            {
                Logger::GetInstance("mcdata")->trace(
                    "collecting simb::MCTruth from input_tag <" + 
                    tag.label() + ">"
                );
                if(!event.getByLabel(tag, sMCTruthHandle))
                {
                    Logger::GetInstance("mcdata")->warning(
                        "no input_tag matching " + tag.label() + 
                        " for simb::MCTruth"
                    );
                }
                else 
                {
                    sMCTruthHandles.emplace_back(event.getHandle<std::vector<simb::MCTruth>>(
                        tag
                    ));
                    sMCTruthHandleLabels.emplace_back(key);
                    if(!sMCTruthHandles.back().isValid()) 
                    {
                        Logger::GetInstance("mcdata")->error(
                            "data product " + tag.label() + 
                            " for simb::MCTruth is invalid!"
                        );
                        exit(0);
                    }
                }
            }
        }
        void MCData::ProcessMCParticles(
            art::Event const& event, art::InputTag input_tag
        )
        {
            Logger::GetInstance("mcdata")->trace(
                "collecting simb::MCParticle from label <" + 
                input_tag.label() + ">"
            );
            if(!event.getByLabel(input_tag, sMCParticleHandle))
            {
                Logger::GetInstance("mcdata")->error(
                    "no label matching " + input_tag.label() + 
                    " for simb::MCParticle!"
                );
                exit(0);
            }
            else 
            {
                sMCParticleHandle = event.getHandle<std::vector<simb::MCParticle>>(
                    input_tag
                );
                if(!sMCParticleHandle.isValid()) 
                {
                    Logger::GetInstance("mcdata")->error(
                        "data product " + input_tag.label() + 
                        " for simb::MCParticle is invalid!"
                    );
                    exit(0);
                }
            }
            Int_t particle_index = 0;
            Logger::GetInstance("mcdata")->trace(
                "creating particle track ID maps for " +
                std::to_string((*sMCParticleHandle).size()) + 
                " <simb::MCParticle>s."
            );
            for (auto particle : *sMCParticleHandle)
            {
                sParticleMap[particle.TrackId()] = particle_index;
                particle_index += 1;

                sGeneratorLabelMap[particle.TrackId()] = GeneratorLabel::None;
                sPDGMap[particle.TrackId()] = particle.PdgCode();
                sParentTrackIDMap[particle.TrackId()] = particle.Mother();
                sParticleEnergyMap[particle.TrackId()] = round(particle.E()*10e6)/10e6;
                // Construct daughter map
                std::vector<Int_t> daughters = {};
                for(auto ii = 0; ii < particle.NumberDaughters(); ii++) {
                    daughters.emplace_back(particle.Daughter(ii));
                }
                sDaughterMap[particle.TrackId()] = daughters;
                // construct progeny map
                sProgenyMap[particle.TrackId()] = daughters;

                // construct ancestry map
                std::vector<Int_t> ancestry = {};
                Int_t mother = particle.Mother();
                Int_t track_id = particle.TrackId();
                Int_t level = 0;
                while (mother != 0)
                {
                    level += 1;
                    track_id = mother;
                    ancestry.emplace_back(mother);
                    mother = sParentTrackIDMap[track_id];
                    if(level > 1) {
                        sProgenyMap[mother].emplace_back(particle.TrackId());
                    }
                }
                sAncestorLevelMap[particle.TrackId()] = level;
                sAncestryMap[particle.TrackId()] = ancestry;

                if (level == 0) {
                    sPrimaries.emplace_back(particle.TrackId());
                    sParentPDGMap[particle.TrackId()] = 0;
                    sAncestorPDGMap[particle.TrackId()] = 0;
                    sAncestorTrackIDMap[particle.TrackId()] = 0;
                    sAncestorEnergyMap[particle.TrackId()] = sParticleEnergyMap[track_id];
                }
                else {
                    sParentPDGMap[particle.TrackId()] = sPDGMap[particle.Mother()];
                    sAncestorPDGMap[particle.TrackId()] = sPDGMap[track_id];
                    sAncestorTrackIDMap[particle.TrackId()] = track_id;
                    sAncestorEnergyMap[particle.TrackId()] = sParticleEnergyMap[track_id];
                }

                // initialize edep maps
                sParticleEdepMap[particle.TrackId()] = {};
                sParticleEdepProcessMap[particle.TrackId()] = {};
                sParticleDetectorSimulationMap[particle.TrackId()] = {};
            }
            for(size_t jj = 0; jj < sMCTruthHandles.size(); jj++)
            {
                for(auto truth : *sMCTruthHandles[jj])
                {
                    /**
                     * MCTruth stores MCParticles starting with trackID = 0,
                     * rather than Geant4 which starts with trackID = 1.
                    */
                    if(truth.NParticles() == 0)
                    {
                        Logger::GetInstance("mcdata")->trace(
                            "MCTruth for " + sMCTruthHandleLabels[jj] + 
                            " contains no simulated particles."
                        );
                        continue;
                    }
                    Logger::GetInstance("mcdata")->trace(
                        "adding labels of type " + 
                        sMCTruthHandleLabels[jj] + 
                        " for " + std::to_string(truth.NParticles()) + 
                        " particles starting with track ID = " + 
                        std::to_string(truth.GetParticle(0).TrackId()+1)
                    );
                    for(Int_t ii = 0; ii < truth.NParticles(); ii++)
                    {
                        if(truth.GetParticle(ii).Process() == "primary")
                        {
                            Double_t init_x = truth.GetParticle(ii).Vx();
                            Double_t init_y = truth.GetParticle(ii).Vy();
                            Double_t init_z = truth.GetParticle(ii).Vz();
                            Int_t pdg_code = truth.GetParticle(ii).PdgCode();
                            for(auto particle : *sMCParticleHandle)
                            {
                                if(
                                    particle.Vx() == init_x &&
                                    particle.Vy() == init_y &&
                                    particle.Vz() == init_z &&
                                    particle.PdgCode() == pdg_code
                                )
                                {
                                    sGeneratorLabelMap[particle.TrackId()] = sGeneratorMap[sMCTruthHandleLabels[jj]];
                                }
                            }
                        }
                    }
                }
            }
        }
        void MCData::ProcessSimEnergyDeposits(art::Event const& event, 
            art::InputTag producer_label, art::InputTag instance_label
        )
        {
            Logger::GetInstance("mcdata")->trace(
                "collecting sim::SimEnergyDeposit from label <" + 
                producer_label.label() + ":" + instance_label.label() + ">"
            );
            if(!event.getByLabel(
                art::InputTag(producer_label.label(), instance_label.label()), 
                sMCSimEnergyDepositHandle
            ))
            {
                Logger::GetInstance("mcdata")->error(
                    "no label matching " + producer_label.label() + ":" + 
                    instance_label.label() + " for sim::SimEnergyDeposit!"
                );
                exit(0);
            }
            else 
            {
                sMCSimEnergyDepositHandle = event.getHandle<std::vector<sim::SimEnergyDeposit>>(
                    art::InputTag(
                        producer_label.label(), instance_label.label()
                    )
                );
                if(!sMCSimEnergyDepositHandle.isValid()) 
                {
                    Logger::GetInstance("mcdata")->error(
                        "data product " + producer_label.label() + ":" + 
                        instance_label.label() + " for simb::SimEnergyDeposit is invalid!"
                    );
                    exit(0);
                }
            }
            Int_t edep_index = 0;
            Logger::GetInstance("mcdata")->trace(
                "creating particle edep ID maps for " +
                std::to_string((*sMCSimEnergyDepositHandle).size()) + 
                " <sim::SimEnergyDeposit>s."
            );
            for(auto edep : *sMCSimEnergyDepositHandle)
            {
                sParticleEdepMap[edep.TrackID()].emplace_back(edep_index);
                ProcessType process = DetermineEdepProcess(edep);
                sParticleEdepProcessMap[edep.TrackID()].emplace_back(process);
                sEdepDetectorSimulationMap[edep_index] = {};
                edep_index += 1;
            }
        }
        void MCData::ProcessSimChannels(art::Event const& event,
            art::InputTag producer_label, art::InputTag instance_label
        )
        {
            Logger::GetInstance("mcdata")->trace(
                "collecting sim::SimChannel from label <" + 
                producer_label.label() + ":" + instance_label.label() + ">"
            );
            if(!event.getByLabel(
                art::InputTag(producer_label.label(), instance_label.label()),
                sMCSimChannelHandle
            ))
            {
                Logger::GetInstance("mcdata")->error(
                    "no label matching " + producer_label.label() + ":" + 
                    instance_label.label() + " for sim::SimChannel!"
                );
                exit(0);
            }
            else 
            {
                sMCSimChannelHandle = event.getHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(
                        producer_label.label(), instance_label.label()
                    )
                );
                if(!sMCSimChannelHandle.isValid()) 
                {
                    Logger::GetInstance("mcdata")->error(
                        "data product " + producer_label.label() + ":" + 
                        instance_label.label() + " for sim::SimChannel is invalid!"
                    );
                    exit(0);
                }
            }
        }
        void MCData::ProcessRawDigits(art::Event const& event,
            art::InputTag producer_label, art::InputTag instance_label
        )
        {
            Logger::GetInstance("mcdata")->trace(
                "collecting raw::RawDigit from label <" + 
                producer_label.label() + ":" + instance_label.label() + ">"
            );
            if(!event.getByLabel(
                art::InputTag(producer_label.label(),instance_label.label()),
                sMCRawDigitHandle
            ))
            {
                Logger::GetInstance("mcdata")->error(
                    "no label matching " + producer_label.label() + ":" + 
                    instance_label.label() + " for raw::RawDigit!"
                );
                exit(0);
            }
            else 
            {
                sMCRawDigitHandle = event.getHandle<std::vector<raw::RawDigit>>(
                    art::InputTag(
                        producer_label.label(), instance_label.label()
                    )
                );
                if(!sMCRawDigitHandle.isValid()) 
                {
                    Logger::GetInstance("mcdata")->error(
                        "data product " + producer_label.label() + ":" + 
                        instance_label.label() + " for raw::RawDigit is invalid!"
                    );
                    exit(0);
                }
            }
            detinfo::DetectorClocksData const clock_data(
                art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
            ); 
            Int_t digit_index = 0;
            Logger::GetInstance("mcdata")->trace(
                "creating detector simulation and particle ID maps for " +
                std::to_string((*sMCRawDigitHandle).size()) + 
                " <raw::RawDigit>s."
            );
            for(auto digit : *sMCRawDigitHandle)
            {
                // Get the channel number for this digit, number of samples,
                // and the pedestal value so that we can uncompress and
                // remove the pedestal.
                raw::ChannelID_t channel = digit.Channel();
                int num_samples = digit.Samples();
                int pedestal = (int)digit.GetPedestal();
                
                // uncompress the digits and remove the pedestal
                std::vector<short> uncompressed(num_samples);
                raw::Uncompress(
                    digit.ADCs(), uncompressed, 
                    pedestal, digit.Compression()
                );
                for (int ii = 0; ii < num_samples; ii++) {
                    uncompressed[ii] -= pedestal;
                }
                sim::SimChannel truth_channel = (*sMCSimChannelHandle)[channel]; 

                // iterate over each tdc value
                for(int l=0; l < num_samples; l++) 
                {
                    // check if threshold has been met
                    if(std::abs(uncompressed[l]) < sADCThreshold) {
                        continue;
                    }
                    auto const& trackIDsAndEnergy = truth_channel.TrackIDsAndEnergies(l, l);
                    /**
                     * This step distinguishes noise from true MC particles.
                     * If the input is noise, pass the detector output to a
                     * noise variable, otherwise, attach the output to the
                     * associated primary.
                     */
                    if(trackIDsAndEnergy.size() == 0) 
                    {
                        sDetectorSimulationNoise.AddNoise(
                            clock_data,
                            l,
                            channel,
                            (Int_t) (std::abs(uncompressed[l]))
                        );
                        continue;
                    }
                    sDetectorSimulation.emplace_back(
                        DetectorSimulation(
                            clock_data,
                            trackIDsAndEnergy,
                            l,
                            channel,
                            (Int_t) (std::abs(uncompressed[l]))
                        )
                    );
                    // associate this detector simulation with a particle particle
                    for(auto track : trackIDsAndEnergy)
                    {
                        if(track.trackID > 0) {
                            sParticleDetectorSimulationMap[track.trackID].emplace_back(digit_index);
                        }
                        else {
                            sRandomDetectorSimulationMap[track.trackID] = {digit_index};
                        }
                    }
                    // determine the edeps associated with this detector simulation
                    sDetectorSimulationEdepMap[digit_index] = DetermineDetectorSimulationEdeps(
                        trackIDsAndEnergy,
                        digit_index
                    );
                    digit_index += 1;
                }
            }
        }
        ProcessType MCData::DetermineEdepProcess(const sim::SimEnergyDeposit& edep)
        {
            std::string process = "NotDefined";
            auto mc_particle = GetMCParticleTrackID(edep.TrackID());
            simb::MCTrajectory trajectory = mc_particle.Trajectory();
            auto trajectory_processes = trajectory.TrajectoryProcesses();
            for(size_t ii = 0; ii < mc_particle.NumberTrajectoryPoints(); ii++)
            {
                if(
                    edep.StartT() <= trajectory.T(ii) &&
                    edep.EndT() >= trajectory.T(ii)
                )
                for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
                {
                    if(trajectory_processes[jj].first == ii) {
                        process = trajectory_processes[jj].second;
                    }
                }
            }
            return TrajectoryStringToProcessType[process];
        }
        std::vector<Int_t> MCData::DetermineDetectorSimulationEdeps(
            const std::vector<sim::IDE>& det_ide,
            Int_t detsim_id
        )
        {
            std::vector<Int_t> edep_ids;
            //std::cout << "DETSIM: " << detsim_id << std::endl;
            for(auto track : det_ide)
            {
                if(track.trackID > 0)
                {
                    //std::cout << "track_id: " << track.trackID << " - (x,y,z) (";
                    //std::cout << track.x << "," << track.y << "," << track.z << ")" << std::endl;
                    std::vector<Int_t> candidate_edeps = sParticleEdepMap[track.trackID];
                    for(auto edep_id : candidate_edeps)
                    {
                        //std::cout << "\tedep_id: " << edep_id << " - (x,y,z) (";
                        auto edep = GetMCSimEnergyDeposit(edep_id);
                        //std::cout << edep.MidPointX() << "," << edep.MidPointY() << "," << edep.MidPointZ() << ")" << std::endl;
                        if(
                            edep.MidPointX() == track.x && 
                            edep.MidPointY() == track.y &&
                            edep.MidPointZ() == track.z
                        ) {
                            //std::cout << "\t\tMATCH" << std::endl;
                            edep_ids.emplace_back(edep_id);
                            sEdepDetectorSimulationMap[edep_id].emplace_back(detsim_id);
                        }
                    }
                }
            }
            return edep_ids;
        }
        void MCData::FillTTree()
        {
            Logger::GetInstance("mcdata")->trace(
                "saving mcdata info to root file."
            );
            // add mcdata info
            sMCDataTree->Fill();
        }    
        std::vector<TrackID_t> MCData::GetPrimariesByGeneratorLabel(GeneratorLabel label)
        {
            std::vector<TrackID_t> primaries;
            for(auto primary : sPrimaries)
            {
                if(sGeneratorLabelMap[primary] == label) {
                    primaries.emplace_back(primary);
                }
            }
            return primaries;
        }
        std::vector<TrackID_t> MCData::GetPrimariesByPDG(Int_t pdg)
        {
            std::vector<TrackID_t> primaries;
            for(auto primary : sPrimaries)
            {
                if (sPDGMap[primary] == pdg) {
                    primaries.emplace_back(primary);
                }
            }
            return primaries;
        }
        std::vector<TrackID_t> MCData::GetParticlesByPDG(Int_t pdg)
        {
            std::vector<TrackID_t> particles;
            for(auto& [key, value]: sPDGMap)
            {
                if (value == pdg) {
                    particles.emplace_back(key);
                }
            }
            return particles;
        }
        std::vector<TrackID_t> MCData::GetDaughtersByPDG(TrackID_t track_id, Int_t pdg)
        {
            std::vector<TrackID_t> daughters;
            for(auto daughter : sDaughterMap[track_id])
            {
                if(sPDGMap[daughter] == pdg) {
                    daughters.emplace_back(daughter);
                }
            }
            return daughters;
        }
        std::vector<TrackID_t> MCData::GetProgenyByPDG(TrackID_t track_id, Int_t pdg)
        {
            std::vector<TrackID_t> progeny;
            for(auto daughter : sProgenyMap[track_id])
            {
                if(sPDGMap[daughter] == pdg) {
                    progeny.emplace_back(daughter);
                }
            }
            return progeny;
        }
        std::vector<TrackID_t> MCData::FilterParticlesByProcess(std::vector<TrackID_t> track_ids, ProcessType process_type)
        {
            std::vector<TrackID_t> particles;
            for(auto track_id : track_ids)
            {
                if(sMCParticleHandle->at(sParticleMap[track_id]).Process() == TrajectoryProcessTypeToString[process_type]) {
                    particles.emplace_back(track_id);
                }
            }
            return particles;
        }    
        std::vector<EdepID_t> MCData::GetParticleAndProgenyEdeps(TrackID_t track_id)
        {
            std::vector<EdepID_t> edeps = GetParticleEdep(track_id);
            auto track_ids = GetProgeny(track_id);
            for(auto track_id : track_ids)
            {
                auto edep_ids = sParticleEdepMap[track_id];
                edeps.insert(edeps.end(), edep_ids.begin(), edep_ids.end());
            }
            return edeps;
        }       
        std::vector<EdepID_t> MCData::GetEdepsByParticles(std::vector<TrackID_t> track_ids)
        {
            std::vector<EdepID_t> edeps;
            for(auto track_id : track_ids)
            {
                auto edep_ids = sParticleEdepMap[track_id];
                edeps.insert(edeps.end(), edep_ids.begin(), edep_ids.end());
            }
            return edeps;
        }
        std::vector<EdepID_t> MCData::FilterEdepsByVolume(std::vector<EdepID_t> edep_ids, geometry::VolumeType volume_type)
        {
            std::vector<EdepID_t> edeps;
            for(auto edep_id : edep_ids)
            {
                geometry::DetectorVolume volume = geometry::DetectorGeometry::GetInstance()->GetVolume(
                    sMCSimEnergyDepositHandle->at(edep_id).MidPointX(),
                    sMCSimEnergyDepositHandle->at(edep_id).MidPointX(),
                    sMCSimEnergyDepositHandle->at(edep_id).MidPointX()
                );
                if(volume.volume_type == volume_type) {
                    edeps.emplace_back(edep_id);
                }
            }
            return edeps;
        }
        std::vector<EdepID_t> MCData::FilterEdepsByPDG(std::vector<EdepID_t> edep_ids, Int_t pdg)
        {
            std::vector<EdepID_t> edeps;
            for(auto edep_id : edep_ids)
            {
                if(sPDGMap[sMCSimEnergyDepositHandle->at(edep_id).TrackID()] == pdg) {
                    edeps.emplace_back(edep_id);
                }
            }
            return edeps;
        }
        std::vector<EdepID_t> MCData::GetDetectorSimulationByEdeps(std::vector<EdepID_t> edep_ids)
        {
            std::vector<DetSimID_t> detsim;
            for(auto edep_id : edep_ids)
            {
                auto detsim_ids = sEdepDetectorSimulationMap[edep_id];
                detsim.insert(detsim.end(), detsim_ids.begin(), detsim_ids.end());
            }
            return detsim;
        }
    }
}