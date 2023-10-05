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
    /**
     * Geant4 breaks up processes into seven major categories
     * (https://agenda.infn.it/event/22667/sessions/17787/attachments/71916/90483/Physics1_Belgrado2019.pdf)
     * :
     * 
     *  1. electromagnetic
     *      a. AdjointAlongStepWeightCorrection: "ContinuousWeightCorrection"
     *      b. ContinuousGainOfEnergy:           "EnergyGain"
     *      c. PolarizedAnnihilation:            "pol-annihil"
     *      d. PolarizedBremsstrahlung:          "pol-eBrem"
     *      e. PolarizedIonization:              "pol-eIoni"
     *      f. IonParameterizedLossModel:        "ParamICRU73"
     *      g. IonIonization:                    "ionIoni"
     *      h. ElectronBremsstrahlung:           "eBrem"
     *      i. ElectronPositronAnnihilation:     "annihil"
     *      h. ElectronIonization:               "eIoni"
     *      j. CoulombScattering:                "CoulombScat"/"CoulombScatter"
     *      k. HadronIonization:                 "hIoni"
     *      l. MuonIonization:                   "muIoni"
     *      m. ElectronElectronToHadrons:        "ee2hadr"
     *      n. MPLIonization:                    "mplIoni"
     *      o. HadronHadronIonization:           "hhIoni"
     *      p. EnergyLoss:                       "EnergyLoss"
     *      q. MultipleScattering:               "msc"
     *  2. hadronic
     *      a. TheoFSGenerator:                  "TheoFSGenerator"
     *      b. NeutrinoNucleusModel:             "neutrino-nucleus"
     *      c. NuElNucleusNCModel:               "NuElNuclNcModel"
     *      d. NuMuNucleusNCModel:               "NuMuNuclNcModel"
     *      e. NuTauNucleusNCModel:              "NuTauNuclNcModel"
     *      f. AntiNuElNucleusNCModel:           "ANuElNuclNcModel"
     *      g. AntiNuMuNucleusNCModel:           "ANuMuNuclNcModel"
     *      h. AntiNuTauNucleusNCModel:          "ANuTauNuclNcModel"
     *      i. NuElNucleusCCModel:               "NuElNuclCcModel"
     *      j. NuMuNucleusCCModel:               "NuMuNuclCcModel"     
     *      k. NuTauNucleusCCModel:              "NuTauNuclCcModel"
     *      l. AntiNuElNucleusCCModel:           "ANuElNuclCcModel"
     *      m. AntiNuMuNucleusCCModel:           "ANuMuNuclCcModel"
     *      n. AntiNuTauNucleusCCModel:          "ANuTauNuclCcModel"
     *      o. CoherentNeutronElectronElModel:   "n-e-elastic"     
     *      p. CoherentHadronElastic:            "hElasticLHEP"
     *      q. CoherentNeutrinoElectronNcModel:  "nu-e-elastic"
     *      r. CoherentElasticHadronNucleusHE:   "hElasticGlauber"
     *      s. Fission:                          "G4LFission"
     *      t. NeutronHPInelastic:               "NeutronHPInelastic"
     *      u. HadronCaptureAtRest:              "hadronCaptureAtRest"
     *      v. HadronElastic:                    "hadElastic"
     *      w. PiMinusInelastic:                 "pi-Inelastic"
     *      x. PiPlusInelastic:                  "pi+Inelastic"
     *      y. KaonMinusInelastic:               "kaon-Inelastic"
     *      z. KaonPlusInelastic:                "kaon+Inelastic"
     *      aa. ProtonInelastic:                 "protonInelastic"
     *      ab. NeutronInelastic:                "neutronInelastic"
     *      ac. NeutronCapture:                  "nCapture"
     *  3. decay
     *      a. Decay:               "Decay"
     *  4. photolepton-hadron
     *  5. optical
     *      a. GammaConversion:     "conv"
     *      b. ComptonScatter:      "compt"
     *      c. PhotoelectricEffect: "phot"
     *  6. parameterization
     *  7. transportation
     *      a. Transportation:  "Transportation"
     * 
    */
    enum class ProcessType
    {
        NotDefined =           -1,
        Unknown =               0,
        Primary =               1,
        // Electromagnetic
        AdjointAlongStepWeightCorrection =  2,
        ContinuousGainOfEnergy =            3,
        PolarizedAnnihilation =             4,
        PolarizedBremsstrahlung =           5,
        PolarizedIonization =               6,
        IonParameteizedLossModel =          7,
        IonIonization =                     8,
        ElectronBremsstrahlung =            9,
        ElectronPositronAnnihilation =      10,
        ElectronIonization =                11,
        CoulombScattering =                 12,
        HadronIonization =                  13,
        MuonIonization =                    14,
        ElectronElectronToHadrons =         15,
        MPLIonization =                     16,
        HadronHadronIonization =            17,
        EnergyLoss =                        18,
        MultipleScattering =                19,
        // Hadronic
        TheoFSGenerator =                   20,
        NeutrinoNucleusModel =              21,
        NuElNucleusNCModel =                22,
        NuMuNucleusNCModel =                23,
        NuTauNucleusNCModel =               24,
        AntiNuElNucleusNCModel =            25,
        AntiNuMuNucleusNCModel =            26,
        AntiNuTauNucleusNCModel =           27,
        NuElNucleusCCModel =                28,
        NuMuNucleusCCModel =                29,
        NuTauNucleusCCModel =               30,
        AntiNuElNucleusCCModel =            31,
        AntiNuMuNucleusCCModel =            32,
        AntiNuTauNucleusCCModel =           33,
        CoherentNeutronElectronElModel =    34,
        CoherentHadronElastic =             35,
        CoherentNeutrinoElectronNCModel =   36,
        CoherentElasticHadronNucleusHE =    37,
        Fission =                           38,
        HadronElastic =                     39,
        HadronInelastic =                   40,
        HadronCaptureAtRest =               41,
        PiMinusElastic =                    42,
        PiMinusInelastic =                  43,
        PiMinusCaptureAtRest =              44,
        PiPlusElastic =                     45,
        PiPlusInelastic =                   46,
        PiPlusCaptureAtRest =               47,
        KaonMinusElastic =                  48,
        KaonMinusInelastic =                49,
        KaonMinusCaptureAtRest =            50,
        KaonPlusElastic =                   51,
        KaonPlusInelastic =                 52,
        KaonPlusCaptureAtRest =             53,
        ProtonElastic =                     54,
        ProtonInelastic =                   55,
        ProtonCaptureAtRest =               56,
        NeutronElastic =                    57,
        NeutronInelastic =                  58,
        NeutronCapture =                    59,
        NeutronHPElastic =                  60,
        NeutronHPInelastic =                61,
        NeutronHPCapture =                  62,
        MuonCaptureAtRest =                 63,
        AntiMuonCaptureAtRest =             64,
        MuonPairProduction =                65,
        // Decay
        Decay =                             66,
        // Optical
        GammaConversion =                   67,
        ComptonScatter =                    68,
        PhotoelectricEffect =               69,
        // Transportation
        Transportation =                    70,
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
        HEPevt = 5
    };
    using SourceLabelInt = std::underlying_type<SourceLabel>::type;
    inline Int_t LabelCast(SourceLabel label)
    {
        return static_cast<SourceLabelInt>(label);
    }
    /**
     * 
    */
    enum class TopologyLabel
    {
        Undefined = -1,
        Noise = 0,
        Blip = 1,
        Track = 2,
        Shower = 3
    };
    using TopologyLabelInt = std::underlying_type<TopologyLabel>::type;
    inline Int_t LabelCast(TopologyLabel label) 
    { 
        return static_cast<TopologyLabelInt>(label);
    }
    /**
     * 
    */
    enum class ParticleLabel
    {
        Undefined = -1,
        Noise = 0,
        // Particle labels are simply the PDG codes
        Electron = 11,
        Positron = -11,
        ElectronNeutrino = 12,
        AntiElectronNeutrino = -12,
        Muon = 13,
        AntiMuon = -13,
        MuonNeutrino = 14,
        AntiMuonNeutrino = -14,
        Tauon = 15,
        AntiTauon = -15,
        TauonNeutrino = 16,
        AntiTauonNeutrino = -16,
        Gamma = 22,
        Pion0 = 111,
        PionPlus = 211,
        PionMinus = -211,
        Kaon0 = 311,
        KaonPlus = 321,
        KaonMinus = -321,
        Neutron = 2112,
        AntiNeutron = -2112,
        Proton = 2212,
        AntiProton = -2212,
        Deuteron = 1000010020,
        Triton = 1000010030,
        Alpha = 1000020040
    };
    using ParticleLabelInt = std::underlying_type<ParticleLabel>::type;
    inline Int_t LabelCast(ParticleLabel label) 
    { 
        return static_cast<ParticleLabelInt>(label);
    }
    enum class PhysicsLabel
    {
        Undefined = -1,
        Noise = 0,
        // Track-like objects
        MIPIonization = 1,
        HIPIonization = 2,
        DeltaElectron = 3,
        MichelElectron = 4,
        // Shower-like objects
        ElectronShower = 5,
        PositronShower = 6,
        PhotonShower = 7,
        // Blip-like objects
        NeutronCaptureGamma474 = 8,
        NeutronCaptureGamma336 = 9,
        NeutronCaptureGamma256 = 10,
        NeutronCaptureGamma118 = 11,
        NeutronCaptureGamma083 = 12,
        NeutronCaptureGamma051 = 13,
        NeutronCaptureGamma016 = 14,
        NeutronCaptureGammaOther = 15,
        Ar39 = 16,
        Ar42 = 17,
        K42 = 18,
        Kr85 = 19,
        Rn222 = 20,
        Po218 = 21,
        At218 = 22,
        Rn218 = 23,
        Pb214 = 24,
        Bi214 = 25,
        Po214 = 26,
        Tl210 = 27,
        Pb210 = 28,
        Bi210 = 29,
        Po210 = 30,
        NuclearRecoil = 31,
        ElectronRecoil = 32
    };
    using PhysicsLabelInt = std::underlying_type<PhysicsLabel>::type;
    inline Int_t LabelCast(PhysicsLabel label) 
    { 
        return static_cast<PhysicsLabelInt>(label);
    }
    enum class LabelingFunction
    {
        Undefined = -1,
        None = 0,
        ProcessElectrons = 1,
        ProcessPositrons = 2,
        ProcessElectronNeutrinos = 3,
        ProcessAntiElectronNeutrinos = 4,
        ProcessMuons = 5,
        ProcessAntiMuons = 6,
        ProcessMuonNeutrinos = 7,
        ProcessAntiMuonNeutrinos = 8,
        ProcessTauons = 9,
        ProcessAntiTauons = 10,
        ProcessTauonNeutrinos = 11,
        ProcessAntiTauonNeutrinos = 12,
        ProcessGammas = 13,
        ProcessPion0s = 14,
        ProcessPionPlus = 15,
        ProcessPionMinus = 16,
        ProcessKaon0s = 17,
        ProcessKaonPlus = 18,
        ProcessKaonMinus = 19,
        ProcessProtons = 20,
        ProcessNeutronCaptures = 21,
        ProcessNuclearRecoils = 22,
        ProcessElectronRecoils = 23,
        ProcessAr39 = 24,
        ProcessAr42 = 25,
        ProcessK42 = 26,
        ProcessKr85 = 27,
        ProcessRn222 = 28,
        ProcessPo218 = 29,
        ProcessAt218 = 30,
        ProcessRn218 = 31,
        ProcessPb214 = 32,
        ProcessBi214 = 33,
        ProcessPo214 = 34,
        ProcessTl210 = 35,
        ProcessPb210 = 36,
        ProcessBi210 = 37,
        ProcessPo210 = 38,
        ProcessCosmics = 39,
        ProcessShowers = 40,
        ProcessNoise = 41,
    };
    using LabelingFunctionInt = std::underlying_type<LabelingFunction>::type;
    inline Int_t LabelCast(LabelingFunction label) 
    { 
        return static_cast<LabelingFunctionInt>(label);
    }
}