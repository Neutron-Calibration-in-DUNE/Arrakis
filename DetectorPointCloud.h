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

    struct DetectorPointCloud
    {
        std::vector<Double_t> channel = {};
        std::vector<Double_t> tdc = {};
        std::vector<Double_t> adc = {};
        std::vector<Int_t> view = {};
        std::vector<Int_t> track_id = {};
        std::vector<ShapeLabelInt> shape_label = {};
        std::vector<ParticleLabelInt> particle_label = {};
        std::vector<Int_t> unique_shape = {};
        std::vector<Int_t> unique_particle = {};

        void clear()
        {
            channel.clear();
            tdc.clear();
            adc.clear();
            view.clear();
            track_id.clear();
            shape_label.clear();
            particle_label.clear();
            unique_shape.clear();
            unique_particle.clear();
        }
    };
}