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
            std::vector<Double_t> x = {};
            std::vector<Double_t> y = {};
            std::vector<Double_t> z = {};

            DetectorSimulation(){}
            DetectorSimulation(
                detinfo::DetectorClocksData const clock_data,
                std::vector<sim::IDE> det_ide,
                Int_t det_tick,
                Int_t det_channel,
                Int_t det_adc
            )
            {
                auto wires = geometry::DetectorGeometry::GetInstance()->ChannelToWire(det_channel);
                auto det_view = geometry::DetectorGeometry::GetInstance()->View(det_channel);
                Double_t wire_multiple = geometry::DetectorGeometry::GetInstance()->GetWirePitch(det_view);

                view = det_view;
                channel = det_channel;
                wire = wires[det_view].Wire * wire_multiple;
                tick = det_tick;
                adc = det_adc;
                tdc = clock_data.TPCTick2TDC(det_tick);

                for(auto ide : det_ide) {
                    track_ids.emplace_back(ide.trackID);
                    energies.emplace_back(ide.energy);
                    x.emplace_back(ide.x);
                    y.emplace_back(ide.y);
                    z.emplace_back(ide.z);
                }
            }

            Double_t largest_energy()
            {
                Double_t energy = 0;
                for(size_t ii = 0; ii < energies.size(); ii++)
                {
                    if(energies[ii] > energy) {
                        energy = energies[ii];
                    }
                }
                return energy;
            }
            Int_t largest_energy_track_id()
            {
                Int_t index = 0;
                Double_t energy = 0;
                for(size_t ii = 0; ii < energies.size(); ii++)
                {
                    if(energies[ii] > energy) {
                        index = ii; energy = energies[ii];
                    }
                }
                return track_ids[index];
            }
            Double_t largest_energy_x()
            {
                Int_t index = 0;
                Double_t energy = 0;
                for(size_t ii = 0; ii < energies.size(); ii++)
                {
                    if(energies[ii] > energy) {
                        index = ii; energy = energies[ii];
                    }
                }
                return x[index];
            }
            Double_t largest_energy_y()
            {
                Int_t index = 0;
                Double_t energy = 0;
                for(size_t ii = 0; ii < energies.size(); ii++)
                {
                    if(energies[ii] > energy) {
                        index = ii; energy = energies[ii];
                    }
                }
                return y[index];
            }
            Double_t largest_energy_z()
            {
                Int_t index = 0;
                Double_t energy = 0;
                for(size_t ii = 0; ii < energies.size(); ii++)
                {
                    if(energies[ii] > energy) {
                        index = ii; energy = energies[ii];
                    }
                }
                return z[index];
            }
        };
        
        struct DetectorSimulationNoise
        {
            std::vector<Int_t> view = {};
            std::vector<Int_t> channel = {};
            std::vector<Int_t> wire = {};
            std::vector<Int_t> tick = {};
            std::vector<Int_t> adc = {};
            std::vector<Int_t> tdc = {};

            DetectorSimulationNoise(){}
            void AddNoise(
                detinfo::DetectorClocksData const clock_data,
                Int_t det_tick,
                Int_t det_channel,
                Int_t det_adc
            )
            {
                auto wires = geometry::DetectorGeometry::GetInstance()->ChannelToWire(det_channel);
                auto det_view = geometry::DetectorGeometry::GetInstance()->View(det_channel);
                Double_t wire_multiple = geometry::DetectorGeometry::GetInstance()->GetWirePitch(det_view);

                view.emplace_back(det_view);
                channel.emplace_back(det_channel);
                wire.emplace_back(wires[det_view].Wire * wire_multiple);
                tick.emplace_back(det_tick);
                adc.emplace_back(det_adc);
                tdc.emplace_back(clock_data.TPCTick2TDC(det_tick));
            }

            void clear()
            {
                view.clear();
                channel.clear();
                wire.clear();
                tick.clear();
                adc.clear();
                tdc.clear();
            }
        };
    }
}