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
    enum class DetectorLabel
    {
        Undefined = -1,
        Noise = 0,
        Muon = 1,
        NeutronCapture = 2,
        NeutronCaptureGamma475 = 3,
        NeutronCaptureGamma181 = 4,
        Ar39 = 5,
        LowEnergyScatter = 6,
        MichelElectron = 7,
        DeltaElectron = 8,
        ElectronShower = 9,
        PionPlus = 10,
        PionMinus = 11,
        KaonPlus = 12,
        KaonMinus = 13,
    };
    using DetectorLabelInt = std::underlying_type<DetectorLabel>::type;
    inline Int_t Label(DetectorLabel label) 
    { 
        return static_cast<DetectorLabelInt>(label);
    }

    struct DetectorPointCloud
    {
        std::vector<Double_t> channel = {};
        std::vector<Double_t> tdc = {};
        std::vector<Double_t> adc = {};
        std::vector<Int_t> view = {};
        std::vector<DetectorLabelInt> label = {};
        std::vector<Int_t> unique_label = {};

        void clear()
        {
            channel.clear();
            tdc.clear();
            adc.clear();
            view.clear();
            label.clear();
            unique_label.clear();
        }
    };
}