/**
 * @file WirePlanePointCloud.h
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
    enum class ShapeLabel
    {
        Undefined = -1,
        Noise = 0,
        Blip = 1,
        Track = 2,
        Shower = 3,
        NeutronCapture = 4,
    };
    using ShapeLabelInt = std::underlying_type<ShapeLabel>::type;
    inline Int_t LabelCast(ShapeLabel label) 
    { 
        return static_cast<ShapeLabelInt>(label);
    }
    enum class ParticleLabel
    {
        Undefined = -1,
        Noise = 0,
        Muon = 1,
        AntiMuon = 2,
        Pion0 = 3,
        PionPlus = 4,
        PionMinus = 5,
        KaonPlus = 6,
        KaonMinus = 7,
        Proton = 8,
        DeltaElectron = 9,
        MichelElectron = 10,
        ElectronShower = 11,
        PositronShower = 12,
        PhotonShower = 13,
        NeutronCaptureGamma = 14,
        NeutronCaptureGamma474 = 15,
        NeutronCaptureGamma336 = 16,
        NeutronCaptureGamma256 = 17,
        NeutronCaptureGamma118 = 18,
        NeutronCaptureGamma083 = 19,
        NeutronCaptureGamma051 = 20,
        NeutronCaptureGamma016 = 21,
        NeutronCaptureGammaOther = 22,
        Ar39 = 23,
        Ar42 = 24,
        Kr85 = 25,
        Rn222 = 26,
        NuclearRecoil = 27,
        ElectronRecoil = 28,
    };
    using ParticleLabelInt = std::underlying_type<ParticleLabel>::type;
    inline Int_t LabelCast(ParticleLabel label) 
    { 
        return static_cast<ParticleLabelInt>(label);
    }

    struct WirePlanePointCloud
    {
        std::vector<Int_t> channel = {};
        std::vector<Int_t> wire = {};
        std::vector<Int_t> tick = {};
        std::vector<Int_t> tdc = {};
        std::vector<Int_t> adc = {};
        std::vector<Int_t> view = {};
        std::vector<Double_t> energy = {};

        std::vector<std::vector<Int_t>> track_ids = {};
        std::vector<std::vector<Double_t>> energies = {};
        std::vector<std::vector<Double_t>> x = {};
        std::vector<std::vector<Double_t>> y = {};
        std::vector<std::vector<Double_t>> z = {};

        std::vector<ShapeLabelInt> shape_label = {};
        std::vector<ParticleLabelInt> particle_label = {};

        std::vector<Int_t> unique_shape = {};
        std::vector<Int_t> unique_particle = {};

        detinfo::DetectorClocksData clock_data;

        WirePlanePointCloud()
        {
        }

        WirePlanePointCloud(
            detinfo::DetectorClocksData const _clock_data
        )
        : clock_data(_clock_data)
        {
        }

        void clear()
        {
            channel.clear();
            wire.clear();
            tick.clear();
            tdc.clear();
            adc.clear();
            view.clear();
            energy.clear();

            track_ids.clear();
            energies.clear();
            
            shape_label.clear();
            particle_label.clear();
            unique_shape.clear();
            unique_particle.clear();
        }

        void AddPoint(
            std::vector<sim::IDE> det_ide,
            Int_t det_tick,
            Int_t det_channel,
            Int_t det_adc,
            bool det_noise
        )
        {
            auto wires = geometry::DetectorGeometry::GetInstance()->ChannelToWire(det_channel);
            auto det_view = geometry::DetectorGeometry::GetInstance()->View(det_channel);
            Double_t wire_multiple = geometry::DetectorGeometry::GetInstance()->GetWirePitch(det_view);

            channel.emplace_back(det_channel);
            wire.emplace_back(wires[det_view].Wire * wire_multiple);
            tick.emplace_back(det_tick);
            tdc.emplace_back(clock_data.TPCTick2TDC(det_tick));
            adc.emplace_back(det_adc);
            view.emplace_back(det_view);

            std::vector<Int_t> det_track_ids;
            std::vector<Double_t> det_energies;
            std::vector<Double_t> det_x;
            std::vector<Double_t> det_y;
            std::vector<Double_t> det_z;
            Double_t det_energy = 0.0;
            for(auto ide : det_ide)
            {
                det_track_ids.emplace_back(ide.trackID);
                det_energies.emplace_back(ide.energy);
                det_x.emplace_back(ide.x);
                det_y.emplace_back(ide.y);
                det_z.emplace_back(ide.z);
                det_energy += ide.energy;
            }
            energy.emplace_back(det_energy);
            track_ids.emplace_back(det_track_ids);
            energies.emplace_back(det_energies);
            x.emplace_back(det_x);
            y.emplace_back(det_y);
            z.emplace_back(det_z);

            if(det_noise)
            {
                shape_label.emplace_back(LabelCast(ShapeLabel::Noise));
                particle_label.emplace_back(LabelCast(ParticleLabel::Noise));
            }
            else
            {
                shape_label.emplace_back(LabelCast(ShapeLabel::Undefined));
                particle_label.emplace_back(LabelCast(ParticleLabel::Undefined));
            }
            unique_shape.emplace_back(-1);
            unique_particle.emplace_back(-1);
        }
    };
}