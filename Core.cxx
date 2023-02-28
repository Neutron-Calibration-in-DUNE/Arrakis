/**
 * @file Core.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12/21
 */
#include "Core.h"

namespace arrakis
{
    std::map<ProcessType, std::string> ProcessTypeToString
    {
        {ProcessType::NotDefined,           "NotDefined"},
        {ProcessType::Unknown,              "Unknown"},
        {ProcessType::Primary,              "Primary"},
        {ProcessType::HadronElastic,        "HadronElastic"},
        {ProcessType::PiMinusInelastic,     "PiMinusInelastic"},
        {ProcessType::PiPlusInelastic,      "PiPlusInelastic"},
        {ProcessType::KaonMinusInelastic,   "KaonMinusInelastic"},
        {ProcessType::KaonPlusInelastic,    "KaonPlusInelastic"},
        {ProcessType::ProtonInelastic,      "ProtonInelastic"},
        {ProcessType::NeutronInelastic,     "NeutronInelastic"},
        {ProcessType::CoulombScatter,       "CoulombScatter"},
        {ProcessType::NeutronCapture,       "NeutronCapture"},
        {ProcessType::Transportation,       "Transportation"},
    };
    std::map<std::string, ProcessType> StringToProcessType
    {
        {"NotDefined",          ProcessType::NotDefined},
        {"Unknown",             ProcessType::Unknown},
        {"Primary",             ProcessType::Primary},
        {"HadronElastic",       ProcessType::HadronElastic},
        {"PiMinusInelastic",    ProcessType::PiMinusInelastic},
        {"PiPlusInelastic",     ProcessType::PiPlusInelastic},
        {"KaonMinusInelastic",  ProcessType::KaonMinusInelastic},
        {"KaonPlusInelastic",   ProcessType::KaonPlusInelastic},
        {"ProtonInelastic",     ProcessType::ProtonInelastic},
        {"NeutronInelastic",    ProcessType::NeutronInelastic},
        {"CoulombScatter",      ProcessType::CoulombScatter},
        {"NeutronCapture",      ProcessType::NeutronCapture},
        {"Transportation",      ProcessType::Transportation},
    };
    std::map<std::string, ProcessType> TrajectoryStringToProcessType
    {
        {"NotDefined",      ProcessType::NotDefined},
        {"Unknown",         ProcessType::Unknown},
        {"primary",         ProcessType::Primary},
        {"hadElastic",      ProcessType::HadronElastic},
        {"pi-Inelastic",    ProcessType::PiMinusInelastic},
        {"pi+Inelastic",    ProcessType::PiPlusInelastic},
        {"kaon-Inelastic",  ProcessType::KaonMinusInelastic},
        {"kaon+Inelastic",  ProcessType::KaonPlusInelastic},
        {"protonInelastic", ProcessType::ProtonInelastic},
        {"neutronInelastic",ProcessType::NeutronInelastic},
        {"CoulombScatter",  ProcessType::CoulombScatter},
        {"nCapture",        ProcessType::NeutronCapture},
        {"Transportation",  ProcessType::Transportation},
    };
}