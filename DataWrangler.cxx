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
        Logger::GetInstance("DataWrangler")->trace(
            "setting up DataWrangler tree."
        );
        sWirePlanePointCloudTree = sTFileService->make<TTree>(
            "data_wire_plane_point_cloud", "data_wire_plane_point_cloud"
        );
        sWirePlanePointCloudTree->Branch("channel", &sWirePlanePointCloud.channel);
        sWirePlanePointCloudTree->Branch("wire",    &sWirePlanePointCloud.wire);
        sWirePlanePointCloudTree->Branch("tick",    &sWirePlanePointCloud.tick);
        sWirePlanePointCloudTree->Branch("tdc",     &sWirePlanePointCloud.tdc);
        sWirePlanePointCloudTree->Branch("adc",     &sWirePlanePointCloud.adc);
        sWirePlanePointCloudTree->Branch("view",    &sWirePlanePointCloud.view);
    }
    void DataWrangler::SetConfigurationParameters(const Parameters& config)
    {
        Logger::GetInstance("DataWrangler")->trace(
            "setting up configuration parameters."
        );
    }
    
    void DataWrangler::ResetEvent()
    {
        sWirePlanePointCloud.clear();
    }
    void DataWrangler::ProcessEvent(
        const Parameters& config, art::Event const& event
    )
    {
        ResetEvent();
        if(config().ProcessRawDigits())
        {
            Logger::GetInstance("DataWrangler")->trace(
                "processing RawDigits"
            );
            ProcessRawDigits(event,
                config().RawDigitProducerLabel(), config().RawDigitInstanceLabel()
            );
        }
    }
    void DataWrangler::ProcessSimChannels(art::Event const& event,
        art::InputTag producer_label, art::InputTag instance_label
    )
    {
        Logger::GetInstance("DataWrangler")->trace(
            "collecting sim::SimChannel from label <" + 
            producer_label.label() + ":" + instance_label.label() + ">"
        );
        if(!event.getByLabel(
            art::InputTag(producer_label.label(), instance_label.label()),
            sMCSimChannelHandle
        ))
        {
            Logger::GetInstance("DataWrangler")->error(
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
                Logger::GetInstance("DataWrangler")->error(
                    "data product " + producer_label.label() + ":" + 
                    instance_label.label() + " for sim::SimChannel is invalid!"
                );
                exit(0);
            }
        }
    }
    void DataWrangler::ProcessRawDigits(art::Event const& event,
        art::InputTag producer_label, art::InputTag instance_label
    )
    {
        Logger::GetInstance("DataWrangler")->trace(
            "collecting raw::RawDigit from label <" + 
            producer_label.label() + ":" + instance_label.label() + ">"
        );
        if(!event.getByLabel(
            art::InputTag(producer_label.label(),instance_label.label()),
            sMCRawDigitHandle
        ))
        {
            Logger::GetInstance("DataWrangler")->error(
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
                Logger::GetInstance("DataWrangler")->error(
                    "data product " + producer_label.label() + ":" + 
                    instance_label.label() + " for raw::RawDigit is invalid!"
                );
                exit(0);
            }
        }
        detinfo::DetectorClocksData const clock_data(
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        ); 
        Logger::GetInstance("DataWrangler")->trace(
            "creating detector Data and particle ID maps for " +
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
                /**
                 */
                sWirePlanePointCloud.AddPoint(
                    clock_data,
                    l,
                    channel,
                    (Int_t) (std::abs(uncompressed[l]))
                );
            }
        }
    }
    void DataWrangler::FillTTree()
    {
        Logger::GetInstance("DataWrangler")->trace(
            "saving wire plane point cloud data to root file."
        );
        sWirePlanePointCloudTree->Fill();
    }
}