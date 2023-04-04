/**
 * @file DataWrangler.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#include "DataWrangler.h"

namespace arrakis
{
    DataWrangler* DataWrangler::sInstance{nullptr};
    std::mutex DataWrangler::sMutex;

    DataWrangler *DataWrangler::GetInstance()
    {
        std::lock_guard<std::mutex> lock(sMutex);
        if (sInstance == nullptr)
        {
            sInstance = new DataWrangler();
        }
        return sInstance;
    }
    DataWrangler::DataWrangler()
    {
        // Logger::GetInstance("DataWrangler")->trace(
        //     "setting up DataWrangler tree."
        // );
        // sDataWranglerTree = sTFileService->make<TTree>(
        //     "edep_point_cloud", "edep_point_cloud"
        // );
    }
    void DataWrangler::SetConfigurationParameters(const Parameters& config)
    {
        Logger::GetInstance("DataWrangler")->trace(
            "setting up configuration parameters."
        );
    }
    
    void DataWrangler::ResetEvent()
    {
    }
    void DataWrangler::ProcessEvent(
        const Parameters& config, art::Event const& event
    )
    {
        ResetEvent();
    }
    void DataWrangler::FillTTree()
    {
    }
}