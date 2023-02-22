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
            //ProcessMCTruth(event);
            ProcessMCParticle(event, config().LArGeantProducerLabel());
            //ProcessSimEnergyDeposit(event, );
        }
        void MCData::ProcessMCTruth(
            art::Event const& event, art::InputTag input_tag
        )
        {
        }
        void MCData::ProcessMCParticle(
            art::Event const& event, art::InputTag input_tag
        )
        {
            mLogger->trace(
                "collecting simb::MCParticle from label <" + 
                input_tag.label() + ">"
            );
            if(!event.getByLabel(input_tag, sMCParticleHandle))
            {
                mLogger->error(
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
        }
        void MCData::ProcessSimEnergyDeposit(
            art::Event const& event, art::InputTag input_tag
        )
        {

        }
    }
}