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
         * 
         * The variables for hits are functions of the four vectors of
         * the particles at the energy deposition location, and the
         * drift parameters theta,
         *      
         *      channel = f_channel(t,x,y,z,edep,theta)
         *      tdc = f_tdc(t,x,y,z,edep,theta)
         *      adc = f_adc(t,x,y,z,edep,theta)
        */
        std::vector<Int_t> hit_channel = {};
        std::vector<Int_t> hit_wire = {};
        std::vector<Int_t> hit_tick = {};
        std::vector<Int_t> hit_tdc = {};
        std::vector<Int_t> hit_adc = {};
        std::vector<Int_t> hit_view = {};

        std::vector<Int_t> hit_mean = {};
        std::vector<Double_t> hit_rms = {};
        std::vector<Double_t> hit_amplitude = {};
        std::vector<Double_t> hit_charge = {};

        void clear()
        {
            hit_channel.clear();
            hit_wire.clear();
            hit_tick.clear();
            hit_tdc.clear();
            hit_adc.clear();
            hit_view.clear();
            hit_mean.clear();
            hit_rms.clear();
            hit_amplitude.clear();
            hit_charge.clear();
        }
        void CreateHit(
            detinfo::DetectorClocksData const& clock_data,
            Int_t det_tick,
            Int_t det_channel,
            Int_t det_adc
        )
        {
            auto wires = DetectorGeometry::GetInstance()->ChannelToWire(det_channel);
            auto det_view = DetectorGeometry::GetInstance()->View(det_channel);
            Double_t wire_multiple = DetectorGeometry::GetInstance()->GetWirePitch(det_view);

            hit_channel.emplace_back(det_channel);
            hit_wire.emplace_back(wires[det_view].Wire * wire_multiple);
            hit_tick.emplace_back(det_tick);
            hit_tdc.emplace_back(clock_data.TPCTick2TDC(det_tick));
            hit_adc.emplace_back(det_adc);
            hit_view.emplace_back(det_view);

            hit_mean.emplace_back(-1);
            hit_rms.emplace_back(-1);
            hit_amplitude.emplace_back(-1);
            hit_charge.emplace_back(-1);
        }

        void AddHit(
            DetSimID_t detsim,
            Double_t mean,
            Double_t rms,
            Double_t amplitude,
            Double_t charge
        )
        {
            hit_mean[detsim] = mean;
            hit_rms[detsim] = rms;
            hit_amplitude[detsim] = amplitude;
            hit_charge[detsim] = charge;
        }
    };
}