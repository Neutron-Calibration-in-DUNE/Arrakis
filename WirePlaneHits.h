/**
 * @file WirePlaneHits.h
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
    struct WirePlaneHits
    {
        /**
         * WirePlaneHits attempts to reconstruct hits from 
         * raw digits.  Hits are constructed by assumed that 
         * each energy deposition is a small clump of electrons
         * which drifts to the wire planes.  Since the influence
         * of the electric field goes like 1/r^2, the shape of the
         * recorded pulses are expected to be Gaussian, with the center of
         * charge corresponding to the mean, which defines where the 
         * crossing occurs in the collection plane.  In the induction 
         * planes, the signal is bipolar so that the crossing point
         * occurs where adc = 0.
        */
        std::vector<Int_t> hit_channel = {};
        std::vector<Int_t> hit_wire = {};
        std::vector<Int_t> hit_tick = {};
        std::vector<Int_t> hit_tdc = {};
        std::vector<Int_t> hit_adc = {};
        std::vector<Int_t> hit_view = {};
        
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