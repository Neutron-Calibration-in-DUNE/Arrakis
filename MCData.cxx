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

        void MCData::ProcessEvent(
            const Parameters& config, art::Event const& event
        )
        {
            std::map<std::string, art::InputTag> generator_tags;
            fhicl::ParameterSet const& generator_labels = config().labels().get_PSet();
            for(std::string const& name : generator_labels.get_names()) {
                generator_tags[name] = generator_labels.get<art::InputTag>(name);
            }
            ProcessMCTruth(event, generator_tags);
            ProcessMCParticle(event, config().LArGeantProducerLabel());
            ProcessSimEnergyDeposit(event, config().IonAndScintProducerLabel());
        }
        void MCData::ProcessMCTruth(
            art::Event const& event, std::map<std::string, art::InputTag> input_tags
        )
        {
            for(auto const& [key, tag] : input_tags)
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
                    // if(sGeneratePointCloudData)
                    // {
                    //     mAr39Label = mParameters().Ar39Label();
                    //     Logger::GetInstance("arrakis_module")->trace("setting Ar39Label = " + mAr39Label.label());
                    //     mSingleNeutronLabel = mParameters().SingleNeutronLabel();
                    //     Logger::GetInstance("arrakis_module")->trace("setting SingleNeutronLabel = " + mSingleNeutronLabel.label());
                    //     mPNSLabel = mParameters().PNSLabel();
                    //     Logger::GetInstance("arrakis_module")->trace("setting PNSLabel = " + mPNSLabel.label());
                    //     sGeneratorMap[mAr39Label] = GeneratorLabel::kAr39;
                    //     sGeneratorMap[mSingleNeutronLabel] = GeneratorLabel::kSingleNeutron;
                    //     sGeneratorMap[mPNSLabel] = GeneratorLabel::kPNS;
                    // }

                }
            }
        }
        void MCData::ProcessMCParticle(
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
            for (auto particle : *sMCParticleHandle)
            {
                sGeneratorLabelMap[particle.TrackId()] = GeneratorLabel::kNone;
                sPDGMap[particle.TrackId()] = particle.PdgCode();
                sParentTrackIDMap[particle.TrackId()] = particle.Mother();
                sParticleEnergyMap[particle.TrackId()] = round(particle.E()*10e6)/10e6;

                Int_t mother = particle.Mother();
                Int_t track_id = particle.TrackId();
                Int_t level = 0;
                while (mother != 0)
                {
                    level += 1;
                    track_id = mother;
                    mother = sParentTrackIDMap[track_id];
                }
                sAncestorLevelMap[particle.TrackId()] = level;
                if (level == 0) {
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
        void MCData::ProcessSimEnergyDeposit(
            art::Event const& event, art::InputTag input_tag
        )
        {
            Logger::GetInstance("mc_data")->trace(
                "collecting sim::SimEnergyDeposit from label <" + 
                input_tag.label() + ">"
            );
            if(!event.getByLabel(input_tag, sMCSimEnergyDepositHandle))
            {
                Logger::GetInstance("mc_data")->error(
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
        }
    }
}