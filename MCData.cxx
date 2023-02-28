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
        void MCData::ResetEvent()
        {
            sMCTruthHandles.clear();
            sPrimaries.clear();

            sGeneratorLabelMap.clear();
            sGeneratorMap.clear();
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

            sEdepProcessMap.clear();
        }

        void MCData::ProcessEvent(
            const Parameters& config, art::Event const& event
        )
        {
            ResetEvent();
            ProcessMCTruth(event, config().labels.get_PSet());
            ProcessMCParticles(event, config().LArGeantProducerLabel());
            ProcessSimEnergyDeposits(event, config().IonAndScintProducerLabel());
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
            for (auto particle : *sMCParticleHandle)
            {
                sParticleMap[particle.TrackId()] = particle_index;
                particle_index += 1;

                sGeneratorLabelMap[particle.TrackId()] = GeneratorLabel::kNone;
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
                    sProgenyMap[track_id].emplace_back(particle.TrackId());
                }

                // initialize edep maps
                sParticleEdepMap[particle.TrackId()] = {};
                sParticleEdepProcessMap[particle.TrackId()] = {};
            }
            for(auto const& [key, val] : sGeneratorMap)
            {
                for(auto mcTruth : sMCTruthHandles)
                {
                    for(auto truth : *mcTruth)
                    {
                        /**
                         * MCTruth stores MCParticles starting with trackID = 0,
                         * rather than Geant4 which starts with trackID = 1.
                        */
                        Logger::GetInstance("mcdata")->trace(
                            "adding labels of type " + 
                            std::to_string(val) + 
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
                                        sGeneratorLabelMap[particle.TrackId()] = val;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        void MCData::ProcessSimEnergyDeposits(
            art::Event const& event, art::InputTag input_tag
        )
        {
            Logger::GetInstance("mcdata")->trace(
                "collecting sim::SimEnergyDeposit from label <" + 
                input_tag.label() + ">"
            );
            if(!event.getByLabel(input_tag, sMCSimEnergyDepositHandle))
            {
                Logger::GetInstance("mcdata")->error(
                    "no label matching " + input_tag.label() + 
                    " for sim::SimEnergyDeposit!"
                );
                exit(0);
            }
            else 
            {
                sMCSimEnergyDepositHandle = event.getHandle<std::vector<sim::SimEnergyDeposit>>(
                    input_tag
                );
                if(!sMCSimEnergyDepositHandle.isValid()) 
                {
                    Logger::GetInstance("mcdata")->error(
                        "data product " + input_tag.label() + 
                        " for simb::SimEnergyDeposit is invalid!"
                    );
                    exit(0);
                }
            }
            Int_t edep_index = 0;
            for(auto edep : *sMCSimEnergyDepositHandle)
            {
                sParticleEdepMap[edep.TrackID()].emplace_back(edep_index);
                ProcessType process = DetermineEdepProcess(edep);
                sParticleEdepProcessMap[edep.TrackID()].emplace_back(process);
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
                        std::cout << trajectory_processes[jj].second << std::endl;
                        process = trajectory_processes[jj].second;
                    }
                }
            }
            return TrajectoryStringToProcessType[process];
        }
    }
}