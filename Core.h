/**
 * @file Core.h
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
    enum GeneratorLabel
    {
        Junk = -1,
        None = 0,
        Ar39 = 1,
        Ar42 = 2,
        Kr85 = 3,
        Rn222 = 4,
        Cosmics = 5,
        HEPevt = 6,
        PNS = 7,
    };
    using GeneratorLabelInt = std::underlying_type<GeneratorLabel>::type;
    inline Int_t Generator(GeneratorLabel generator_label)
    {
        return static_cast<GeneratorLabelInt>(generator_label);
    }
    
    enum class ProcessType
    {
        NotDefined =           -1,
        Unknown =               0,
        Primary =               1,
        HadronElastic =         2,
        PiMinusInelastic =      3,
        PiPlusInelastic =       4,
        KaonMinusInelastic =    5,
        KaonPlusInelastic =     6,
        ProtonInelastic =       7,
        NeutronInelastic =      8,
        CoulombScatter =        9,
        NeutronCapture =        10,
        Transportation =        11,
        Decay =                 12,
        ComptonScatter =        13,
        PhotoelectricEffect =   14,
        ElectronBremsstrahlung = 15,
        ElectronIonization =    16,
        PositronAnnihilation =  17,
        MuonIonization =        18,
        GammaConversion =       19,
        IonIonization =         20,
        MuonCaptureAtRest =     21,
    };
    using ProcessTypeInt = std::underlying_type<ProcessType>::type;
    inline Int_t Process(ProcessType process) 
    { 
        return static_cast<ProcessTypeInt>(process);
    }

    inline void PrintIndices(std::string label, std::vector<Int_t> indices)
    {
        std::cout << label << ": [ ";
        for(auto index : indices)
        {
            std::cout << index << " ";
        }
        std::cout << "]" << std::endl;
    }
}