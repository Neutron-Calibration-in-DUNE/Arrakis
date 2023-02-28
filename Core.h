/**
 * @file Generators.h
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

namespace arrakis
{
    enum GeneratorLabel
    {
        kJunk = -1,
        kNone = 0,
        kAr39 = 1,
        kSingleNeutron = 2,
        kPNS = 3,
        kNeutronCaptureGamma4_75 = 4,
        kNeutronCaptureGamma1_18 = 5,
    };

    enum class ProcessType
    {
        NotDefined = -1,
        Unknown = 0,
        HadronElastic = 1,
        PiMinusInelastic = 2,
        PiPlusInelastic = 3,
        KaonMinusInelastic = 4,
        KaonPlusInelastic = 5,
        ProtonInelastic = 6,
        NeutronInelastic = 7,
        CoulombScatter = 8,
        NeutronCapture = 9,
        Transportation = 10,
    };
    using ProcessTypeInt = std::underlying_type<ProcessType>::type;
    inline Int_t Process(ProcessType process) 
    { 
        return static_cast<ProcessTypeInt>(process);
    }

    const std::map<ProcessType, std::string> ProcessTypeToString
    {
        {ProcessType::NotDefined, "NotDefined"},
        {ProcessType::Unknown, "Unknown"},
        {ProcessType::HadronElastic, "HadronElastic"},
        {ProcessType::PiMinusInelastic, "PiMinusInelastic"},
        {ProcessType::PiPlusInelastic, "PiPlusInelastic"},
        {ProcessType::KaonMinusInelastic, "KaonMinusInelastic"},
        {ProcessType::KaonPlusInelastic, "KaonPlusInelastic"},
        {ProcessType::ProtonInelastic, "ProtonInelastic"},
        {ProcessType::NeutronInelastic, "NeutronInelastic"},
        {ProcessType::CoulombScatter, "CoulombScatter"},
        {ProcessType::NeutronCapture, "NeutronCapture"},
        {ProcessType::Transportation, "Transportation"},
    };
    const std::map<std::string, ProcessType> StringToProcessType
    {
        {"NotDefined", ProcessType::NotDefined},
        {"Unknown", ProcessType::Unknown},
        {"HadronElastic", ProcessType::HadronElastic},
        {"PiMinusInelastic", ProcessType::PiMinusInelastic},
        {"PiPlusInelastic", ProcessType::PiPlusInelastic},
        {"KaonMinusInelastic", ProcessType::KaonMinusInelastic},
        {"KaonPlusInelastic", ProcessType::KaonPlusInelastic},
        {"ProtonInelastic", ProcessType::ProtonInelastic},
        {"NeutronInelastic", ProcessType::NeutronInelastic},
        {"CoulombScatter", ProcessType::CoulombScatter},
        {"NeutronCapture", ProcessType::NeutronCapture},
        {"Transportation", ProcessType::Transportation},
    };

}