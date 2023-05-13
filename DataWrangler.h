/**
 * @file DataWrangler.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#pragma once
#include <mutex>
#include <iostream>
#include <iomanip>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "Configuration.h"
#include "Core.h"
#include "Logger.h"
#include "WirePlanePointCloud.h"

namespace arrakis
{
    class DataWrangler
    {
    public:
        // this singleton cannot be cloned
        DataWrangler(DataWrangler &other) = delete;
        // singleton should also not be assignable
        void operator=(const DataWrangler &) = delete;

        // static method that controls access to 
        // the singleton instance
        static DataWrangler* GetInstance();

        void SetConfigurationParameters(const Parameters& config);

        // methods for processing event data
        void ResetEvent();
        void ProcessEvent(const Parameters& config, art::Event const& event);
        void ProcessSimChannels(art::Event const& event,
            art::InputTag producer_label, art::InputTag instance_label
        );
        void ProcessRawDigits(art::Event const& event,
            art::InputTag producer_label, art::InputTag instance_label
        );

        // fill TTree
        void FillTTree();

    protected:
        DataWrangler();
        ~DataWrangler() {}

    private:
        static DataWrangler * sInstance;
        static std::mutex sMutex;

        // Output TTree
        art::ServiceHandle<art::TFileService> sTFileService;   

        // service handles
        art::Handle<std::vector<sim::SimChannel>>       sMCSimChannelHandle;
        art::Handle<std::vector<raw::RawDigit>>         sMCRawDigitHandle;

        WirePlanePointCloud sWirePlanePointCloud;    
    };
}