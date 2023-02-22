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
                mLogger = Logger::GetInstance("mcdata");
            }
            return sInstance;
        }

        void MCData::ProcessEvent(
            const Parameters& params, art::Event const& event
        )
        {

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

        }
        void MCData::ProcessSimEnergyDeposit(
            art::Event const& event, art::InputTag input_tag
        )
        {

        }
    }
}