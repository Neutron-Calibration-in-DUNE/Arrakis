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
        DeltaElectron = 8,
        MichelElectron = 9,
        ElectronShower = 10,
        NeutronCapture = 11,
        NeutronCaptureGamma475 = 12,
        NeutronCaptureGamma181 = 13,
        Ar39 = 14,
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
        std::vector<ShapeLabelInt> shape_label = {};
        std::vector<ParticleLabelInt> particle_label = {};
        std::vector<Int_t> unique_label = {};

        void clear()
        {
            channel.clear();
            tdc.clear();
            adc.clear();
            view.clear();
            shape_label.clear();
            particle_label.clear();
            unique_label.clear();
        }
    };
}