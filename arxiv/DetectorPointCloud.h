/**
 * @file DetectorPointCloud.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12/21
 */
#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <map>

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "Logger.h"
#include "SoloPointCloud.h"

namespace arrakis
{
    struct DetectorPointCloud
    {
        Int_t event_id = {-1};
        std::vector<Int_t> view = {};
        std::vector<Double_t> wire = {};
        std::vector<Double_t> channel = {};
        std::vector<Double_t> tick = {};
        std::vector<Double_t> tdc = {};
        std::vector<Double_t> adc = {};
        std::vector<Double_t> energy = {};
        std::vector<std::string> group_label = {};
        std::vector<Int_t> group_label_id = {};
        std::vector<std::string> label = {};
        std::vector<Int_t> label_id = {};

        DetectorPointCloud() {}

        void AddPointCloud(
            SoloPointCloud point_cloud
        )
        {
            view.insert(view.end(),point_cloud.view.begin(),point_cloud.view.end());
            wire.insert(wire.end(),point_cloud.wire.begin(),point_cloud.wire.end());
            channel.insert(channel.end(),point_cloud.channel.begin(),point_cloud.channel.end());
            tick.insert(tick.end(),point_cloud.tick.begin(),point_cloud.tick.end());
            tdc.insert(tdc.end(),point_cloud.tdc.begin(),point_cloud.tdc.end());
            adc.insert(adc.end(),point_cloud.adc.begin(),point_cloud.adc.end());
            energy.insert(energy.end(),point_cloud.energy.begin(),point_cloud.energy.end());
            for(size_t ii = 0; ii < point_cloud.view.size(); ii++)
            {
                group_label.emplace_back(point_cloud.group_label);
                group_label_id.emplace_back(point_cloud.group_label_id);
                label.emplace_back(point_cloud.label);
                label_id.emplace_back(point_cloud.label_id);
            }
        }
    };
}