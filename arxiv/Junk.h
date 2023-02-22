/**
 * @file Junk.h
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
#include "ParticleMaps.h"
#include "Trajectory.h"

namespace arrakis
{
    struct Junk
    {
        GeneratorLabel generator_label = {kJunk};
        
        // Raw digit information.
        std::vector<Int_t> det_view = {};
        std::vector<Int_t> det_channel = {};
        std::vector<Int_t> det_wire = {};
        std::vector<Int_t> det_tick = {};
        std::vector<Int_t> det_adc = {};
        std::vector<Double_t> det_tdc = {};

        Junk(){}

        /**
         * @brief Construct a new Junk object with an input label
         * and a simb::MCParticle.
         * 
         * @param label 
         * @param particle 
         */
        Junk(GeneratorLabel label)
        {
            generator_label = label;
        }

        void AddJunkDetectorSimulation(
            detinfo::DetectorClocksData const& clockData,
            Int_t tick, 
            Int_t channel,
            Int_t adc
        )
        {
            auto view = DetectorGeometry::GetInstance("junk")->View(channel);
            auto wires = DetectorGeometry::GetInstance("junk")->ChannelToWire(channel);
            Double_t wire_multiple = DetectorGeometry::GetInstance("junk")->GetWirePitch(view);

            det_view.emplace_back(view);
            det_channel.emplace_back(channel);
            det_wire.emplace_back(wires[view].Wire * wire_multiple);
            det_tick.emplace_back(tick);
            det_adc.emplace_back(adc);
            det_tdc.emplace_back(clockData.TPCTick2TDC(tick));     
        }
    };
}