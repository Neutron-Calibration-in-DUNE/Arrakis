/**
 * @file WirePlaneTrackTopology.h
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

#include "Core.h"
#include "DetectorGeometry.h"

namespace arrakis
{
    struct WirePlaneTrackTopology
    {
        /**
         * WirePlaneTrackTopology concerns topological information
         * for tracks in the WirePlane configuration space.
        */
        std::vector<Int_t> track_begin_channel = {};
        std::vector<Int_t> track_begin_wire = {};
        std::vector<Int_t> track_begin_tick = {};
        std::vector<Int_t> track_begin_tdc = {};
        std::vector<Int_t> track_begin_adc = {};
        std::vector<Int_t> track_begin_view = {};
        
        std::vector<Int_t> track_end_channel = {};
        std::vector<Int_t> track_end_wire = {};
        std::vector<Int_t> track_end_tick = {};
        std::vector<Int_t> track_end_tdc = {};
        std::vector<Int_t> track_end_adc = {};
        std::vector<Int_t> track_end_view = {};

        std::vector<Int_t> vertex_channel = {};
        std::vector<Int_t> vertex_wire = {};
        std::vector<Int_t> vertex_tick = {};
        std::vector<Int_t> vertex_tdc = {};
        std::vector<Int_t> vertex_adc = {};
        std::vector<Int_t> vertex_view = {};
    }
}