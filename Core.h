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
    using OpDetChannelID_t = Int_t;
    using OpDetBacktrackerID_t = Int_t;

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
        Beam = 5,
        Cosmics = 6,
        HEPevt = 7,
        PNS = 8,
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

    /**
     * 
    */
    enum class SourceLabel
    {
        Undefined = -1,
        Noise = 0,
        Cosmics = 1,
        Beam = 2,
        Radiological = 3,
        PulsedNeutronSource = 4,
        HEPevt = 5,
        Mixed = 6
    };
    using SourceLabelInt = std::underlying_type<SourceLabel>::type;
    inline Int_t LabelCast(SourceLabel label)
    {
        return static_cast<SourceLabelInt>(label);
    }
    /**
     * 
    */
    enum class ShapeLabel
    {
        Undefined = -1,
        Noise = 0,
        Blip = 1,
        Track = 2,
        Shower = 3,
        NeutronCapture = 4,
        Mixed = 5
    };
    using ShapeLabelInt = std::underlying_type<ShapeLabel>::type;
    inline Int_t LabelCast(ShapeLabel label) 
    { 
        return static_cast<ShapeLabelInt>(label);
    }
    /**
     * 
    */
    enum class ParticleLabel
    {
        Undefined = -1,
        Noise = 0,
        Muon = 1,
        AntiMuon = 2,
        Pion0 = 3,
        PionPlus = 4,
        PionMinus = 5,
        Kaon0 = 6,
        KaonPlus = 7,
        KaonMinus = 8,
        Proton = 9,
        DeltaElectron = 10,
        MichelElectron = 11,
        ElectronShower = 12,
        PositronShower = 13,
        PhotonShower = 14,
        NeutronCaptureGamma = 15,
        NeutronCaptureGamma474 = 16,
        NeutronCaptureGamma336 = 17,
        NeutronCaptureGamma256 = 18,
        NeutronCaptureGamma118 = 19,
        NeutronCaptureGamma083 = 20,
        NeutronCaptureGamma051 = 21,
        NeutronCaptureGamma016 = 22,
        NeutronCaptureGammaOther = 23,
        Ar39 = 24,
        Ar42 = 25,
        Kr85 = 26,
        Rn222 = 27,
        NuclearRecoil = 28,
        ElectronRecoil = 29,
        Mixed = 30
    };
    using ParticleLabelInt = std::underlying_type<ParticleLabel>::type;
    inline Int_t LabelCast(ParticleLabel label) 
    { 
        return static_cast<ParticleLabelInt>(label);
    }
}