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
        Alpha = 1000020040,
        Sulfur32 = 1000160320,
        Sulfur33 = 1000160330,
        Sulfur34 = 1000160340,
        Sulfur35 = 1000160350,
        Sulfur36 = 1000160360,
        Chlorine35 = 1000170350,
        Chlorine36 = 1000170360,
        Chlorine37 = 1000170370,
        Chlorine38 = 1000170380,
        Chlorine39 = 1000170390,
        Chlorine40 = 1000170400,
        Argon36 = 1000180360,
        Argon37 = 1000180370,
        Argon38 = 1000180380,
        Argon39 = 1000180390,
        Argon40 = 1000180400,
        Argon41 = 1000180410,
        Ion = 1000000000
    };
    constexpr std::array<ParticleLabel, 46> ParticleLabels = {
        ParticleLabel::Undefined,
        ParticleLabel::Noise,
        // Particle labels are simply the PDG codes
        ParticleLabel::Electron,
        ParticleLabel::Positron,
        ParticleLabel::ElectronNeutrino,
        ParticleLabel::AntiElectronNeutrino,
        ParticleLabel::Muon,
        ParticleLabel::AntiMuon,
        ParticleLabel::MuonNeutrino,
        ParticleLabel::AntiMuonNeutrino,
        ParticleLabel::Tauon,
        ParticleLabel::AntiTauon,
        ParticleLabel::TauonNeutrino,
        ParticleLabel::AntiTauonNeutrino,
        ParticleLabel::Gamma,
        ParticleLabel::Pion0,
        ParticleLabel::PionPlus,
        ParticleLabel::PionMinus,
        ParticleLabel::Kaon0,
        ParticleLabel::KaonPlus,
        ParticleLabel::KaonMinus,
        ParticleLabel::Neutron,
        ParticleLabel::AntiNeutron,
        ParticleLabel::Proton,
        ParticleLabel::AntiProton,
        ParticleLabel::Deuteron,
        ParticleLabel::Triton,
        ParticleLabel::Alpha,
        ParticleLabel::Sulfur32,
        ParticleLabel::Sulfur33,
        ParticleLabel::Sulfur34,
        ParticleLabel::Sulfur35,
        ParticleLabel::Sulfur36,
        ParticleLabel::Chlorine35,
        ParticleLabel::Chlorine36,
        ParticleLabel::Chlorine37,
        ParticleLabel::Chlorine38,
        ParticleLabel::Chlorine39,
        ParticleLabel::Chlorine40,
        ParticleLabel::Argon36,
        ParticleLabel::Argon37,
        ParticleLabel::Argon38,
        ParticleLabel::Argon39,
        ParticleLabel::Argon40,
        ParticleLabel::Argon41,
        ParticleLabel::Ion
    };
    using ParticleLabelInt = std::underlying_type<ParticleLabel>::type;
    inline Int_t LabelCast(ParticleLabel label) 
    { 
        return static_cast<ParticleLabelInt>(label);
    }
    inline bool IsParticleLabel(int particleLabelInt)
    {
        for (auto currentEnum : ParticleLabels)
        {
            if (static_cast<int>(currentEnum) == particleLabelInt) {
                return true;
            }
        }
        return false;
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
        Shower = 3,
    };
    using TopologyLabelInt = std::underlying_type<TopologyLabel>::type;
    inline Int_t LabelCast(TopologyLabel label) 
    { 
        return static_cast<TopologyLabelInt>(label);
    }

    enum class PhysicsMicroLabel
    {
        Undefined = -1,
        Noise = 0,
        MIPIonization = 1,
        HIPIonization = 2,
        ElectronIonization = 3,
        Bremsstrahlung = 4,
        Annihilation = 5,
        PhotoElectric = 6,
        GammaCompton = 7,
        GammaConversion = 8,
        HadronElastic = 9,
        HadronInelastic = 10
    };
    using PhysicsMicroLabelInt = std::underlying_type<PhysicsMicroLabel>::type;
    inline Int_t LabelCast(PhysicsMicroLabel label) 
    { 
        return static_cast<PhysicsMicroLabelInt>(label);
    }

    /**
     * 
    */
    enum class PhysicsMesoLabel
    {
        Undefined = -1,
        Noise = 0,
        MIP = 1,
        HIP = 2,
        DeltaElectron = 3,
        MichelElectron = 4,
        ElectronShower = 5,
        PositronShower = 6,
        PhotonShower = 7,
        LowEnergyIonization = 8,
        NeutronCaptureGamma474 = 9,
        NeutronCaptureGamma336 = 10,
        NeutronCaptureGamma256 = 11,
        NeutronCaptureGamma118 = 12,
        NeutronCaptureGamma083 = 13,
        NeutronCaptureGamma051 = 14,
        NeutronCaptureGamma016 = 15,
        NeutronCaptureGammaOther = 16,
        Pi0Decay = 17,
        AlphaDecay = 18,
        BetaDecay = 19,
        GammaDecay = 20,
        NuclearRecoil = 21,
        ElectronRecoil = 22
    };
    using PhysicsMesoLabelInt = std::underlying_type<PhysicsMesoLabel>::type;
    inline Int_t LabelCast(PhysicsMesoLabel label) 
    { 
        return static_cast<PhysicsMesoLabelInt>(label);
    }

    /**
     * 
    */
    enum class PhysicsMacroLabel
    {
        Undefined = -1,
        Noise = 0,

        // Neutrino interactions
        CCNue = 1,
        CCNuMu = 2,
        NC = 3,

        Cosmics = 5,

        // Radiological interactions
        Ar39 = 6,
        Ar42 = 7,
        K42 = 8,
        Kr85 = 9,
        Rn222 = 10,
        Po218a = 11,
        Po218b = 12,
        At218a = 13,
        At218b = 14,
        Rn218 = 15,
        Pb214 = 16,
        Bi214a = 17,
        Bi214b = 18,
        Po214 = 19,
        Tl210 = 20,
        Pb210a = 21,
        Pb210b = 22,
        Bi210a = 23,
        Bi210b = 24,
        Po210 = 25,

    };
    using PhysicsMacroLabelInt = std::underlying_type<PhysicsMacroLabel>::type;
    inline Int_t LabelCast(PhysicsMacroLabel label)
    {
        return static_cast<PhysicsMacroLabelInt>(label);
    }

    /**
     * 
    */
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