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
            if(!event.getByLabel(input_tag, mMCParticleHandle))
            {
                mLogger->error(
                    "no label matching " + input_tag.label() + 
                    " for simb::MCParticle!"
                );
                exit(0);
            }
            else 
            {
                mMCParticleHandle = event.getHandle<std::vector<simb::MCParticle>>(
                    input_tag
                );
            }
        }
        void MCData::ProcessSimEnergyDeposit(
            art::Event const& event, art::InputTag input_tag
        )
        {

        }
    }
}