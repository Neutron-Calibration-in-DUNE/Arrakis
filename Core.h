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
    /**
     * To reduce confusion throughout MCData code, we will use different names
     * for integers that mean different things.  TrackID_t refers to the actual 
     * track id of a particle, while ParticleID_t refers to the index of a particle
     * in the simb::MCParticle vector.  EdepID_t and DetSimID_t are also the indices
     * of the sim::SimEnergyDeposit and arrakis::DetectorSimulation vectors respectively.
     */
    using TrackID_t = Int_t;
    using EdepID_t = Int_t;
    using ParticleID_t = Int_t;
    using DetSimID_t = Int_t;

    using TrackID_List = std::vector<TrackID_t>;
    using EdepID_List = std::vector<EdepID_t>;
    using ParticleID_List = std::vector<ParticleID_t>;
    using DetSimID_List = std::vector<DetSimID_t>;

    using TrackID_Collection = std::vector<TrackID_List>;
    using EdepID_Collection = std::vector<EdepID_List>;
    using ParticleID_Collection = std::vector<ParticleID_List>;
    using DetSimID_Collection = std::vector<DetSimID_List>;

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
    using ProcessType_List = std::vector<ProcessType>;

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