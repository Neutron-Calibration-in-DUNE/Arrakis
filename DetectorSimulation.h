/**
 * @file DetectorSimulation.h
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
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/sim.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "DetectorGeometry.h"

namespace arrakis
{
    namespace mcdata
    {
        struct DetectorSimulation
        {
            Int_t view = {0};
            Int_t channel = {0};
            Int_t wire = {0};
            Int_t tick = {0};
            Int_t adc = {0};
            Int_t tdc = {0};

            std::vector<Int_t> track_ids = {};
            std::vector<Double_t> energies = {};

            DetectorSimulation(){}
            DetectorSimulation(
                detinfo::DetectorClocksData const clock_data,
                std::vector<sim::IDE> det_ide,
                Int_t det_tick,
                Int_t det_channel,
                Int_t det_adc
            )
            {
                auto wires = geometry::DetectorGeometry::GetInstance()->ChannelToWire(channel);
                auto det_view = geometry::DetectorGeometry::GetInstance()->View(channel);
                Double_t wire_multiple = geometry::DetectorGeometry::GetInstance()->GetWirePitch(det_view);

                view = det_view;
                channel = det_channel;
                wire = wires[view].Wire * wire_multiple;
                tick = det_tick;
                adc = det_adc;
                tdc = clock_data.TPCTick2TDC(tick);

                for(auto ide : det_ide) {
                    track_ids.emplace_back(ide.trackID);
                    energies.emplace_back(ide.energy);
                }
            }
        };
    }
}