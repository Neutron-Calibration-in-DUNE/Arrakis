/**
 * @file SimulationWrangler.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-22
 */
#include "SimulationWrangler.h"

namespace arrakis
{
    std::map<ProcessType, std::string> ProcessTypeToString
    {
        {ProcessType::NotDefined,           "NotDefined"},
        {ProcessType::Unknown,              "Unknown"},
        {ProcessType::Primary,              "Primary"},
        {ProcessType::AdjointAlongStepWeightCorrection, "AdjointAlongStepWeightCorrection"},
        {ProcessType::ContinuousGainOfEnergy,           "ContinuousGainOfEnergy"},
        {ProcessType::PolarizedAnnihilation,            "PolarizedAnnihilation"},
        {ProcessType::PolarizedBremsstrahlung,          "PolarizedBremsstrahlung"},
        {ProcessType::PolarizedIonization,              "PolarizedIonization"},
        {ProcessType::IonParameteizedLossModel,         "IonParameteizedLossModel"},
        {ProcessType::IonIonization,                    "IonIonization"},
        {ProcessType::ElectronBremsstrahlung,           "ElectronBremsstrahlung"},
        {ProcessType::ElectronPositronAnnihilation,     "ElectronPositronAnnihilation"},
        {ProcessType::ElectronIonization,               "ElectronIonization"},
        {ProcessType::CoulombScattering,                "CoulombScattering"},
        {ProcessType::HadronIonization,                 "HadronIonization"},
        {ProcessType::MuonIonization,                   "MuonIonization"},
        {ProcessType::ElectronElectronToHadrons,        "ElectronElectronToHadrons"},
        {ProcessType::MPLIonization,                    "MPLIonization"},
        {ProcessType::HadronHadronIonization,           "HadronHadronIonization"},
        {ProcessType::EnergyLoss,                       "EnergyLoss"},
        {ProcessType::MultipleScattering,               "MultipleScattering"},
        {ProcessType::TheoFSGenerator,                  "TheoFSGenerator"},
        {ProcessType::NeutrinoNucleusModel,             "NeutrinoNucleusModel"},
        {ProcessType::NuElNucleusNCModel,               "NuElNucleusNCModel"},
        {ProcessType::NuMuNucleusNCModel,               "NuMuNucleusNCModel"},
        {ProcessType::NuTauNucleusNCModel,              "NuTauNucleusNCModel"},
        {ProcessType::AntiNuElNucleusNCModel,           "AntiNuElNucleusNCModel"},
        {ProcessType::AntiNuMuNucleusNCModel,           "AntiNuMuNucleusNCModel"},
        {ProcessType::AntiNuTauNucleusNCModel,          "AntiNuTauNucleusNCModel"},
        {ProcessType::NuElNucleusCCModel,               "NuElNucleusCCModel"},
        {ProcessType::NuMuNucleusCCModel,               "NuMuNucleusCCModel"},
        {ProcessType::NuTauNucleusCCModel,              "NuTauNucleusCCModel"},
        {ProcessType::AntiNuElNucleusCCModel,           "AntiNuElNucleusCCModel"},
        {ProcessType::AntiNuMuNucleusCCModel,           "AntiNuMuNucleusCCModel"},
        {ProcessType::AntiNuTauNucleusCCModel,          "AntiNuTauNucleusCCModel"},
        {ProcessType::CoherentNeutronElectronElModel,   "CoherentNeutronElectronElModel"},
        {ProcessType::CoherentHadronElastic,            "CoherentHadronElastic"},
        {ProcessType::CoherentNeutrinoElectronNCModel,  "CoherentNeutrinoElectronNCModel"},
        {ProcessType::CoherentElasticHadronNucleusHE,   "CoherentElasticHadronNucleusHE"},
        {ProcessType::Fission,                          "Fission"},
        {ProcessType::HadronElastic,                    "HadronElastic"},
        {ProcessType::HadronInelastic,                  "HadronInelastic"},
        {ProcessType::HadronCaptureAtRest,              "HadronCaptureAtRest"},
        {ProcessType::PiMinusElastic,                   "PiMinusElastic"},
        {ProcessType::PiMinusInelastic,                 "PiMinusInelastic"},
        {ProcessType::PiMinusCaptureAtRest,             "PiMinusCaptureAtRest"},
        {ProcessType::PiPlusElastic,                    "PiPlusElastic"},
        {ProcessType::PiPlusInelastic,                  "PiPlusInelastic"},
        {ProcessType::PiPlusCaptureAtRest,              "PiPlusCaptureAtRest"},
        {ProcessType::KaonMinusElastic,                 "KaonMinusElastic"},
        {ProcessType::KaonMinusInelastic,               "KaonMinusInelastic"},
        {ProcessType::KaonMinusCaptureAtRest,           "KaonMinusCaptureAtRest"},
        {ProcessType::KaonPlusElastic,                  "KaonPlusElastic"},
        {ProcessType::KaonPlusInelastic,                "KaonPlusInelastic"},
        {ProcessType::KaonPlusCaptureAtRest,            "KaonPlusCaptureAtRest"},
        {ProcessType::ProtonElastic,                    "ProtonElastic"},
        {ProcessType::ProtonInelastic,                  "ProtonInelastic"},
        {ProcessType::ProtonCaptureAtRest,              "ProtonCaptureAtRest"},
        {ProcessType::NeutronElastic,                   "NeutronElastic"},
        {ProcessType::NeutronInelastic,                 "NeutronInelastic"},
        {ProcessType::NeutronCapture,                   "NeutronCapture"},
        {ProcessType::NeutronHPElastic,                 "NeutronHPElastic"},
        {ProcessType::NeutronHPInelastic,               "NeutronHPInelastic"},
        {ProcessType::NeutronHPCapture,                 "NeutronHPCapture"},
        {ProcessType::MuonCaptureAtRest,                "MuonCaptureAtRest"},
        {ProcessType::AntiMuonCaptureAtRest,            "AntiMuonCaptureAtRest"},
        {ProcessType::Decay,                            "Decay"},
        {ProcessType::GammaConversion,                  "GammaConversion"},
        {ProcessType::ComptonScatter,                   "ComptonScatter"},
        {ProcessType::PhotoelectricEffect,              "PhotoelectricEffect"},
        {ProcessType::Transportation,                   "Transportation"},
    };
    std::map<std::string, ProcessType> StringToProcessType
    {
        {"NotDefined",         ProcessType::NotDefined},
        {"Unknown",            ProcessType::Unknown},
        {"Primary",            ProcessType::Primary},
        {"AdjointAlongStepWeightCorrection",ProcessType::AdjointAlongStepWeightCorrection},
        {"ContinuousGainOfEnergy",         ProcessType::ContinuousGainOfEnergy},
        {"PolarizedAnnihilation",          ProcessType::PolarizedAnnihilation},
        {"PolarizedBremsstrahlung",        ProcessType::PolarizedBremsstrahlung},
        {"PolarizedIonization",            ProcessType::PolarizedIonization},
        {"IonParameteizedLossModel",       ProcessType::IonParameteizedLossModel},
        {"IonIonization",                  ProcessType::IonIonization},
        {"ElectronBremsstrahlung",         ProcessType::ElectronBremsstrahlung},
        {"ElectronPositronAnnihilation",   ProcessType::ElectronPositronAnnihilation},
        {"ElectronIonization",             ProcessType::ElectronIonization},
        {"CoulombScattering",              ProcessType::CoulombScattering},
        {"HadronIonization",               ProcessType::HadronIonization},
        {"MuonIonization",                 ProcessType::MuonIonization},
        {"ElectronElectronToHadrons",      ProcessType::ElectronElectronToHadrons},
        {"MPLIonization",                  ProcessType::MPLIonization},
        {"HadronHadronIonization",         ProcessType::HadronHadronIonization},
        {"EnergyLoss",                     ProcessType::EnergyLoss},
        {"MultipleScattering",             ProcessType::MultipleScattering},
        {"TheoFSGenerator",                ProcessType::TheoFSGenerator},
        {"NeutrinoNucleusModel",           ProcessType::NeutrinoNucleusModel},
        {"NuElNucleusNCModel",             ProcessType::NuElNucleusNCModel},
        {"NuMuNucleusNCModel",             ProcessType::NuMuNucleusNCModel},
        {"NuTauNucleusNCModel",            ProcessType::NuTauNucleusNCModel},
        {"AntiNuElNucleusNCModel",         ProcessType::AntiNuElNucleusNCModel},
        {"AntiNuMuNucleusNCModel",         ProcessType::AntiNuMuNucleusNCModel},
        {"AntiNuTauNucleusNCModel",        ProcessType::AntiNuTauNucleusNCModel},
        {"NuElNucleusCCModel",             ProcessType::NuElNucleusCCModel},
        {"NuMuNucleusCCModel",             ProcessType::NuMuNucleusCCModel},
        {"NuTauNucleusCCModel",            ProcessType::NuTauNucleusCCModel},
        {"AntiNuElNucleusCCModel",         ProcessType::AntiNuElNucleusCCModel},
        {"AntiNuMuNucleusCCModel",         ProcessType::AntiNuMuNucleusCCModel},
        {"AntiNuTauNucleusCCModel",        ProcessType::AntiNuTauNucleusCCModel},
        {"CoherentNeutronElectronElModel", ProcessType::CoherentNeutronElectronElModel},
        {"CoherentHadronElastic",          ProcessType::CoherentHadronElastic},
        {"CoherentNeutrinoElectronNCModel",ProcessType::CoherentNeutrinoElectronNCModel},
        {"CoherentElasticHadronNucleusHE", ProcessType::CoherentElasticHadronNucleusHE},
        {"Fission",                        ProcessType::Fission},
        {"HadronElastic",                  ProcessType::HadronElastic},
        {"HadronInelastic",                ProcessType::HadronInelastic},
        {"HadronCaptureAtRest",            ProcessType::HadronCaptureAtRest},
        {"PiMinusElastic",                 ProcessType::PiMinusElastic},
        {"PiMinusInelastic",               ProcessType::PiMinusInelastic},
        {"PiMinusCaptureAtRest",           ProcessType::PiMinusCaptureAtRest},
        {"PiPlusElastic",                  ProcessType::PiPlusElastic},
        {"PiPlusInelastic",                ProcessType::PiPlusInelastic},
        {"PiPlusCaptureAtRest",            ProcessType::PiPlusCaptureAtRest},
        {"KaonMinusElastic",               ProcessType::KaonMinusElastic},
        {"KaonMinusInelastic",             ProcessType::KaonMinusInelastic},
        {"KaonMinusCaptureAtRest",         ProcessType::KaonMinusCaptureAtRest},
        {"KaonPlusElastic",                ProcessType::KaonPlusElastic},
        {"KaonPlusInelastic",              ProcessType::KaonPlusInelastic},
        {"KaonPlusCaptureAtRest",          ProcessType::KaonPlusCaptureAtRest},
        {"ProtonElastic",                  ProcessType::ProtonElastic},
        {"ProtonInelastic",                ProcessType::ProtonInelastic},
        {"ProtonCaptureAtRest",            ProcessType::ProtonCaptureAtRest},
        {"NeutronElastic",                 ProcessType::NeutronElastic},
        {"NeutronInelastic",               ProcessType::NeutronInelastic},
        {"NeutronCapture",                 ProcessType::NeutronCapture},
        {"NeutronHPElastic",               ProcessType::NeutronHPElastic},
        {"NeutronHPInelastic",             ProcessType::NeutronHPInelastic},
        {"NeutronHPCapture",               ProcessType::NeutronHPCapture},
        {"MuonCaptureAtRest",              ProcessType::MuonCaptureAtRest},
        {"AntiMuonCaptureAtRest",          ProcessType::AntiMuonCaptureAtRest},
        {"Decay",                          ProcessType::Decay},
        {"GammaConversion",                ProcessType::GammaConversion},
        {"ComptonScatter",                 ProcessType::ComptonScatter},
        {"PhotoelectricEffect",            ProcessType::PhotoelectricEffect},
        {"Transportation",                 ProcessType::Transportation},
    };
    std::map<std::string, ProcessType> TrajectoryStringToProcessType
    {
        {"NotDefined",         ProcessType::NotDefined},
        {"Unknown",            ProcessType::Unknown},
        {"primary",            ProcessType::Primary},
        {"ContinuousWeightCorrection",  ProcessType::AdjointAlongStepWeightCorrection},
        {"EnergyGain",                  ProcessType::ContinuousGainOfEnergy},
        {"pol-annihil",                 ProcessType::PolarizedAnnihilation},
        {"pol-eBrem",                   ProcessType::PolarizedBremsstrahlung},
        {"pol-eIoni",                   ProcessType::PolarizedIonization},
        {"ParamICRU73",                 ProcessType::IonParameteizedLossModel},
        {"ionIoni",                     ProcessType::IonIonization},
        {"eBrem",                       ProcessType::ElectronBremsstrahlung},
        {"annihil",                     ProcessType::ElectronPositronAnnihilation},
        {"eIoni",                       ProcessType::ElectronIonization},
        {"CoulombScat",                 ProcessType::CoulombScattering},
        {"CoulombScatter",              ProcessType::CoulombScattering},
        {"hIoni",                       ProcessType::HadronIonization},
        {"muIoni",                      ProcessType::MuonIonization},
        {"ee2hadr",                     ProcessType::ElectronElectronToHadrons},
        {"mplIoni",                     ProcessType::MPLIonization},
        {"hhIoni",                      ProcessType::HadronHadronIonization},
        {"EnergyLoss",                  ProcessType::EnergyLoss},
        {"msc",                         ProcessType::MultipleScattering},
        {"TheoFSGenerator",             ProcessType::TheoFSGenerator},
        {"neutrino-nucleus",            ProcessType::NeutrinoNucleusModel},
        {"NuElNuclNcModel",             ProcessType::NuElNucleusNCModel},
        {"NuMuNuclNcModel",             ProcessType::NuMuNucleusNCModel},
        {"NuTauNuclNcModel",            ProcessType::NuTauNucleusNCModel},
        {"ANuElNuclNcModel",            ProcessType::AntiNuElNucleusNCModel},
        {"ANuMuNuclNcModel",            ProcessType::AntiNuMuNucleusNCModel},
        {"ANuTauNuclNcModel",           ProcessType::AntiNuTauNucleusNCModel},
        {"NuElNuclCcModel",             ProcessType::NuElNucleusCCModel},
        {"NuMuNuclCcModel",             ProcessType::NuMuNucleusCCModel},
        {"NuTauNuclCcModel",            ProcessType::NuTauNucleusCCModel},
        {"ANuElNuclCcModel",            ProcessType::AntiNuElNucleusCCModel},
        {"ANuMuNuclCcModel",            ProcessType::AntiNuMuNucleusCCModel},
        {"ANuTauNuclCcModel",           ProcessType::AntiNuTauNucleusCCModel},
        {"n-e-elastic",                 ProcessType::CoherentNeutronElectronElModel},
        {"hElasticLHEP",                ProcessType::CoherentHadronElastic},
        {"nu-e-elastic",                ProcessType::CoherentNeutrinoElectronNCModel},
        {"hElasticGlauber",             ProcessType::CoherentElasticHadronNucleusHE},
        {"G4LFission",                  ProcessType::Fission},
        {"hadElastic",                  ProcessType::HadronElastic},
        {"hadronInelastic",             ProcessType::HadronInelastic},
        {"hadronCaptureAtRest",         ProcessType::HadronCaptureAtRest},
        {"pi-Elastic",                  ProcessType::PiMinusElastic},
        {"pi-Inelastic",                ProcessType::PiMinusInelastic},
        {"pi-CaptureAtRest",            ProcessType::PiMinusCaptureAtRest},
        {"pi+Elastic",                  ProcessType::PiPlusElastic},
        {"pi+Inelastic",                ProcessType::PiPlusInelastic},
        {"pi+CaptureAtRest",            ProcessType::PiPlusCaptureAtRest},
        {"kaon-Elastic",                ProcessType::KaonMinusElastic},
        {"kaon-Inelastic",              ProcessType::KaonMinusInelastic},
        {"kaon-CaptureAtRest",          ProcessType::KaonMinusCaptureAtRest},
        {"kaon+Elastic",                ProcessType::KaonPlusElastic},
        {"kaon+Inelastic",              ProcessType::KaonPlusInelastic},
        {"kaon+CaptureAtRest",          ProcessType::KaonPlusCaptureAtRest},
        {"protonElastic",               ProcessType::ProtonElastic},
        {"protonInelastic",             ProcessType::ProtonInelastic},
        {"protonCaptureAtRest",         ProcessType::ProtonCaptureAtRest},
        {"neutronElastic",              ProcessType::NeutronElastic},
        {"neutronInelastic",            ProcessType::NeutronInelastic},
        {"nCapture",                    ProcessType::NeutronCapture},
        {"NeutronHPElastic",            ProcessType::NeutronHPElastic},
        {"NeutronHPInelastic",          ProcessType::NeutronHPInelastic},
        {"NeutronHPCapture",            ProcessType::NeutronHPCapture},
        {"muMinusCaptureAtRest",        ProcessType::MuonCaptureAtRest},
        {"muPlusCaptureAtRest",         ProcessType::AntiMuonCaptureAtRest},
        {"Decay",                       ProcessType::Decay},
        {"conv",                        ProcessType::GammaConversion},
        {"compt",                       ProcessType::ComptonScatter},
        {"phot",                        ProcessType::PhotoelectricEffect},
        {"Transportation",              ProcessType::Transportation},
    };
    std::map<ProcessType, std::string> TrajectoryProcessTypeToString
    {
        {ProcessType::NotDefined,           "NotDefined"},
        {ProcessType::Unknown,              "Unknown"},
        {ProcessType::Primary,              "primary"},
        {ProcessType::AdjointAlongStepWeightCorrection, "ContinuousWeightCorrection"},
        {ProcessType::ContinuousGainOfEnergy,           "EnergyGain"},
        {ProcessType::PolarizedAnnihilation,            "pol-annihil"},
        {ProcessType::PolarizedBremsstrahlung,          "pol-eBrem"},
        {ProcessType::PolarizedIonization,              "pol-eIoni"},
        {ProcessType::IonParameteizedLossModel,         "ParamICRU73"},
        {ProcessType::IonIonization,                    "ionIoni"},
        {ProcessType::ElectronBremsstrahlung,           "eBrem"},
        {ProcessType::ElectronPositronAnnihilation,     "annihil"},
        {ProcessType::ElectronIonization,               "eIoni"},
        {ProcessType::CoulombScattering,                "CoulombScat"},
        {ProcessType::CoulombScattering,                "CoulombScatter"},
        {ProcessType::HadronIonization,                 "hIoni"},
        {ProcessType::MuonIonization,                   "muIoni"},
        {ProcessType::ElectronElectronToHadrons,        "ee2hadr"},
        {ProcessType::MPLIonization,                    "mplIoni"},
        {ProcessType::HadronHadronIonization,           "hhIoni"},
        {ProcessType::EnergyLoss,                       "EnergyLoss"},
        {ProcessType::MultipleScattering,               "msc"},
        {ProcessType::TheoFSGenerator,                  "TheoFSGenerator"},
        {ProcessType::NeutrinoNucleusModel,             "neutrino-nucleus"},
        {ProcessType::NuElNucleusNCModel,               "NuElNuclNcModel"},
        {ProcessType::NuMuNucleusNCModel,               "NuMuNuclNcModel"},
        {ProcessType::NuTauNucleusNCModel,              "NuTauNuclNcModel"},
        {ProcessType::AntiNuElNucleusNCModel,           "ANuElNuclNcModel"},
        {ProcessType::AntiNuMuNucleusNCModel,           "ANuMuNuclNcModel"},
        {ProcessType::AntiNuTauNucleusNCModel,          "ANuTauNuclNcModel"},
        {ProcessType::NuElNucleusCCModel,               "NuElNuclCcModel"},
        {ProcessType::NuMuNucleusCCModel,               "NuMuNuclCcModel"},
        {ProcessType::NuTauNucleusCCModel,              "NuTauNuclCcModel"},
        {ProcessType::AntiNuElNucleusCCModel,           "ANuElNuclCcModel"},
        {ProcessType::AntiNuMuNucleusCCModel,           "ANuMuNuclCcModel"},
        {ProcessType::AntiNuTauNucleusCCModel,          "ANuTauNuclCcModel"},
        {ProcessType::CoherentNeutronElectronElModel,   "n-e-elastic"},
        {ProcessType::CoherentHadronElastic,            "hElasticLHEP"},
        {ProcessType::CoherentNeutrinoElectronNCModel,  "nu-e-elastic"},
        {ProcessType::CoherentElasticHadronNucleusHE,   "hElasticGlauber"},
        {ProcessType::Fission,                          "G4LFission"},
        {ProcessType::HadronElastic,                    "hadElastic"},
        {ProcessType::HadronInelastic,                  "hadronInelastic"},
        {ProcessType::HadronCaptureAtRest,              "hadronCaptureAtRest"},
        {ProcessType::PiMinusElastic,                   "pi-Elastic"},
        {ProcessType::PiMinusInelastic,                 "pi-Inelastic"},
        {ProcessType::PiMinusCaptureAtRest,             "pi-CaptureAtRest"},
        {ProcessType::PiPlusElastic,                    "pi+Elastic"},
        {ProcessType::PiPlusInelastic,                  "pi+Inelastic"},
        {ProcessType::PiPlusCaptureAtRest,              "pi+CaptureAtRest"},
        {ProcessType::KaonMinusElastic,                 "kaon-Elastic"},
        {ProcessType::KaonMinusInelastic,               "kaon-Inelastic"},
        {ProcessType::KaonMinusCaptureAtRest,           "kaon-CaptureAtRest"},
        {ProcessType::KaonPlusElastic,                  "kaon+Elastic"},
        {ProcessType::KaonPlusInelastic,                "kaon+Inelastic"},
        {ProcessType::KaonPlusCaptureAtRest,            "kaon+CaptureAtRest"},
        {ProcessType::ProtonElastic,                    "protonElastic"},
        {ProcessType::ProtonInelastic,                  "protonInelastic"},
        {ProcessType::ProtonCaptureAtRest,              "protonCaptureAtRest"},
        {ProcessType::NeutronElastic,                   "neutronElastic"},
        {ProcessType::NeutronInelastic,                 "neutronInelastic"},
        {ProcessType::NeutronCapture,                   "nCapture"},
        {ProcessType::NeutronHPElastic,                 "NeutronHPElastic"},
        {ProcessType::NeutronHPInelastic,               "NeutronHPInelastic"},
        {ProcessType::NeutronHPCapture,                 "NeutronHPCapture"},
        {ProcessType::MuonCaptureAtRest,                "muMinusCaptureAtRest"},
        {ProcessType::AntiMuonCaptureAtRest,            "muPlusCaptureAtRest"},
        {ProcessType::Decay,                            "Decay"},
        {ProcessType::GammaConversion,                  "conv"},
        {ProcessType::ComptonScatter,                   "compt"},
        {ProcessType::PhotoelectricEffect,              "phot"},
        {ProcessType::Transportation,                   "Transportation"},
    };
    SimulationWrangler* SimulationWrangler::sInstance{nullptr};
    std::mutex SimulationWrangler::sMutex;

    SimulationWrangler *SimulationWrangler::GetInstance()
    {
        std::lock_guard<std::mutex> lock(sMutex);
        if (sInstance == nullptr)
        {
            sInstance = new SimulationWrangler();
        }
        return sInstance;
    }
    SimulationWrangler::SimulationWrangler()
    {
        if (sSaveEnergyDepositPointCloud)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "setting up EnergyDepositPointCloud tree."
            );
            sEnergyDepositPointCloudTree = sTFileService->make<TTree>(
                "mc_edep_point_cloud", "mc_edep_point_cloud"
            );
            sEnergyDepositPointCloudTree->Branch("edep_t",   &sEnergyDepositPointCloud.edep_t);
            sEnergyDepositPointCloudTree->Branch("edep_x",   &sEnergyDepositPointCloud.edep_x);
            sEnergyDepositPointCloudTree->Branch("edep_y",   &sEnergyDepositPointCloud.edep_y);
            sEnergyDepositPointCloudTree->Branch("edep_z",   &sEnergyDepositPointCloud.edep_z);
            sEnergyDepositPointCloudTree->Branch("edep_energy",         &sEnergyDepositPointCloud.edep_energy);
            sEnergyDepositPointCloudTree->Branch("edep_num_photons",    &sEnergyDepositPointCloud.edep_num_photons);
            sEnergyDepositPointCloudTree->Branch("edep_num_electrons",  &sEnergyDepositPointCloud.edep_num_electrons);
            sEnergyDepositPointCloudTree->Branch("edep_track_id",       &sEnergyDepositPointCloud.edep_track_id);
            sEnergyDepositPointCloudTree->Branch("topology_label",      &sEnergyDepositPointCloud.topology_label);
            sEnergyDepositPointCloudTree->Branch("particle_label",      &sEnergyDepositPointCloud.particle_label);
            sEnergyDepositPointCloudTree->Branch("physics_micro_label", &sEnergyDepositPointCloud.physics_micro_label);
            sEnergyDepositPointCloudTree->Branch("physics_meso_label",  &sEnergyDepositPointCloud.physics_meso_label);
            sEnergyDepositPointCloudTree->Branch("physics_macro_label", &sEnergyDepositPointCloud.physics_macro_label);
            sEnergyDepositPointCloudTree->Branch("unique_topology",     &sEnergyDepositPointCloud.unique_topology_label);
            sEnergyDepositPointCloudTree->Branch("unique_particle",     &sEnergyDepositPointCloud.unique_particle_label);
            sEnergyDepositPointCloudTree->Branch("unique_physics_micro_label",  &sEnergyDepositPointCloud.unique_physics_micro_label);
            sEnergyDepositPointCloudTree->Branch("unique_physics_meso_label",   &sEnergyDepositPointCloud.unique_physics_meso_label);
            sEnergyDepositPointCloudTree->Branch("unique_physics_macro_label",  &sEnergyDepositPointCloud.unique_physics_macro_label);
            sEnergyDepositPointCloudTree->Branch("edep_detsim_id",      &sEnergyDepositPointCloud.edep_detsim_id);
        }

        if (sSaveSimulationWrangler)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "setting up SimulationWrangler tree."
            );
            sSimulationWranglerTree = sTFileService->make<TTree>("mc_maps", "mc_maps");
            sSimulationWranglerTree->Branch("pdg_code_map",         &sTrackID_PDGCodeMap);
            sSimulationWranglerTree->Branch("generator_map",        &sGeneratorMap);
            sSimulationWranglerTree->Branch("generator_label_map",  &sTrackID_GeneratorLabelMap);
            sSimulationWranglerTree->Branch("particle_id_map",      &sTrackID_ParticleIDMap);
            sSimulationWranglerTree->Branch("particle_energy_map",  &sTrackID_EnergyMap);       
            sSimulationWranglerTree->Branch("parent_track_id_map",  &sTrackID_ParentTrackIDMap);
            sSimulationWranglerTree->Branch("parent_pdg_code_map",  &sTrackID_ParentPDGCodeMap);
            sSimulationWranglerTree->Branch("ancestor_track_id_map",    &sTrackID_AncestorTrackIDMap);
            sSimulationWranglerTree->Branch("ancestor_level_map",       &sTrackID_AncestorLevelMap);
            sSimulationWranglerTree->Branch("ancestor_pdg_code_map",    &sTrackID_AncestorPDGCodeMap);
            sSimulationWranglerTree->Branch("daughter_track_id_map",    &sTrackID_DaughterTrackIDMap);
            sSimulationWranglerTree->Branch("progeny_track_id_map",     &sTrackID_ProgenyTrackIDMap);
            sSimulationWranglerTree->Branch("ancestry_track_id_map",    &sTrackID_AncestryTrackIDMap);
            sSimulationWranglerTree->Branch("edep_id_map",              &sTrackID_EdepIDMap);
            sSimulationWranglerTree->Branch("edep_process_map",         &sTrackID_EdepProcessMap);
            sSimulationWranglerTree->Branch("detsim_map",               &sTrackID_DetSimIDMap);
            sSimulationWranglerTree->Branch("random_detsim_map",        &sTrackID_RandomDetSimIDMap);
            sSimulationWranglerTree->Branch("edep_detsim_map",          &sEdepID_DetSimIDMap);
            sSimulationWranglerTree->Branch("detsim_edep_map",          &sDetSimID_EdepIDMap);
        }

        if (sSaveWirePlanePointCloud)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "setting up WirePlanePointCloud tree."
            );
            sWirePlanePointCloudTree = sTFileService->make<TTree>(
                "mc_wire_plane_point_cloud", "mc_wire_plane_point_cloud"
            );
            sWirePlanePointCloudTree->Branch("channel", &sWirePlanePointCloud.channel);
            sWirePlanePointCloudTree->Branch("wire",    &sWirePlanePointCloud.wire);
            sWirePlanePointCloudTree->Branch("tick",    &sWirePlanePointCloud.tick);
            sWirePlanePointCloudTree->Branch("tdc",     &sWirePlanePointCloud.tdc);
            sWirePlanePointCloudTree->Branch("adc",     &sWirePlanePointCloud.adc);
            sWirePlanePointCloudTree->Branch("view",    &sWirePlanePointCloud.view);
            sWirePlanePointCloudTree->Branch("energy",  &sWirePlanePointCloud.energy);
            sWirePlanePointCloudTree->Branch("topology_label",      &sWirePlanePointCloud.topology_label);
            sWirePlanePointCloudTree->Branch("particle_label",      &sWirePlanePointCloud.particle_label);
            sWirePlanePointCloudTree->Branch("physics_micro_label", &sWirePlanePointCloud.physics_micro_label);
            sWirePlanePointCloudTree->Branch("physics_meso_label",  &sWirePlanePointCloud.physics_meso_label);
            sWirePlanePointCloudTree->Branch("physics_macro_label", &sWirePlanePointCloud.physics_macro_label);
            sWirePlanePointCloudTree->Branch("unique_topology",     &sWirePlanePointCloud.unique_topology_label);
            sWirePlanePointCloudTree->Branch("unique_particle",     &sWirePlanePointCloud.unique_particle_label);
            sWirePlanePointCloudTree->Branch("unique_physics_micro_label",  &sWirePlanePointCloud.unique_physics_micro_label);
            sWirePlanePointCloudTree->Branch("unique_physics_meso_label",   &sWirePlanePointCloud.unique_physics_meso_label);
            sWirePlanePointCloudTree->Branch("unique_physics_macro_label",  &sWirePlanePointCloud.unique_physics_macro_label);
            if (sSaveWirePlaneInductionFlag) {
                sWirePlanePointCloudTree->Branch("induction_flag",  &sWirePlanePointCloud.induction_flag);
            }
        }

        if (sSaveWirePlaneHits)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "setting up hits in WirePlanePointCloud tree."
            );
            sWirePlanePointCloudTree->Branch("hit_mean",    &sWirePlanePointCloud.hit_mean);
            sWirePlanePointCloudTree->Branch("hit_rms",     &sWirePlanePointCloud.hit_rms);
            sWirePlanePointCloudTree->Branch("hit_amplitude",     &sWirePlanePointCloud.hit_amplitude);
            sWirePlanePointCloudTree->Branch("hit_charge",        &sWirePlanePointCloud.hit_charge);
        }

        if (sSaveWirePlaneTrackTopology)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "setting up WirePlaneTrackTopology tree."
            );
            sWirePlaneTrackTopologyTree = sTFileService->make<TTree>(
                "mc_wire_plane_track_topology", "mc_wire_plane_track_topology"
            );
            sWirePlaneTrackTopologyTree->Branch("track_begin_channel", &sWirePlaneTrackTopology.track_begin_channel);
            sWirePlaneTrackTopologyTree->Branch("track_begin_wire",    &sWirePlaneTrackTopology.track_begin_wire);
            sWirePlaneTrackTopologyTree->Branch("track_begin_tick",    &sWirePlaneTrackTopology.track_begin_tick);
            sWirePlaneTrackTopologyTree->Branch("track_begin_tdc",     &sWirePlaneTrackTopology.track_begin_tdc);
            sWirePlaneTrackTopologyTree->Branch("track_begin_adc",     &sWirePlaneTrackTopology.track_begin_adc);
            sWirePlaneTrackTopologyTree->Branch("track_begin_view",    &sWirePlaneTrackTopology.track_begin_view);
            sWirePlaneTrackTopologyTree->Branch("track_end_channel", &sWirePlaneTrackTopology.track_end_channel);
            sWirePlaneTrackTopologyTree->Branch("track_end_wire",    &sWirePlaneTrackTopology.track_end_wire);
            sWirePlaneTrackTopologyTree->Branch("track_end_tick",    &sWirePlaneTrackTopology.track_end_tick);
            sWirePlaneTrackTopologyTree->Branch("track_end_tdc",     &sWirePlaneTrackTopology.track_end_tdc);
            sWirePlaneTrackTopologyTree->Branch("track_end_adc",     &sWirePlaneTrackTopology.track_end_adc);
            sWirePlaneTrackTopologyTree->Branch("track_end_view",    &sWirePlaneTrackTopology.track_end_view);
            sWirePlaneTrackTopologyTree->Branch("vertex_channel", &sWirePlaneTrackTopology.vertex_channel);
            sWirePlaneTrackTopologyTree->Branch("vertex_wire",    &sWirePlaneTrackTopology.vertex_wire);
            sWirePlaneTrackTopologyTree->Branch("vertex_tick",    &sWirePlaneTrackTopology.vertex_tick);
            sWirePlaneTrackTopologyTree->Branch("vertex_tdc",     &sWirePlaneTrackTopology.vertex_tdc);
            sWirePlaneTrackTopologyTree->Branch("vertex_adc",     &sWirePlaneTrackTopology.vertex_adc);
            sWirePlaneTrackTopologyTree->Branch("vertex_view",    &sWirePlaneTrackTopology.vertex_view);
        }

        if (sSaveOpDetPointCloud)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "setting up OpDetPointCloud tree."
            );
            sOpDetPointCloudTree = sTFileService->make<TTree>(
                "mc_op_det_point_cloud", "mc_op_det_point_cloud"
            );
            sOpDetPointCloudTree->Branch("channel", &sOpDetPointCloud.channel);
            sOpDetPointCloudTree->Branch("tick",    &sOpDetPointCloud.tick);
            sOpDetPointCloudTree->Branch("adc",     &sOpDetPointCloud.adc);
            sOpDetPointCloudTree->Branch("energy",  &sOpDetPointCloud.energy);
            sOpDetPointCloudTree->Branch("topology_label",      &sOpDetPointCloud.topology_label);
            sOpDetPointCloudTree->Branch("particle_label",      &sOpDetPointCloud.particle_label);
            sOpDetPointCloudTree->Branch("physics_micro_label", &sOpDetPointCloud.physics_micro_label);
            sOpDetPointCloudTree->Branch("physics_meso_label",  &sOpDetPointCloud.physics_meso_label);
            sOpDetPointCloudTree->Branch("physics_macro_label", &sOpDetPointCloud.physics_macro_label);
            sOpDetPointCloudTree->Branch("unique_topology",     &sOpDetPointCloud.unique_topology_label);
            sOpDetPointCloudTree->Branch("unique_particle",     &sOpDetPointCloud.unique_particle_label);
            sOpDetPointCloudTree->Branch("unique_physics_micro_label",  &sOpDetPointCloud.unique_physics_micro_label);
            sOpDetPointCloudTree->Branch("unique_physics_meso_label",   &sOpDetPointCloud.unique_physics_meso_label);
            sOpDetPointCloudTree->Branch("unique_physics_macro_label",  &sOpDetPointCloud.unique_physics_macro_label);
        }

        sGeneratorMap["Ar39Label"] =    GeneratorLabel::Ar39;
        sGeneratorMap["Ar42Label"] =    GeneratorLabel::Ar42;
        sGeneratorMap["Kr85Label"] =    GeneratorLabel::Kr85;
        sGeneratorMap["Rn222Label"] =   GeneratorLabel::Rn222;
        sGeneratorMap["BeamLabel"] =    GeneratorLabel::Beam;
        sGeneratorMap["CosmicsLabel"] = GeneratorLabel::Cosmics;
        sGeneratorMap["HEPevtLabel"] =  GeneratorLabel::HEPevt;
        sGeneratorMap["PNSLabel"] =     GeneratorLabel::PNS;
    }
    void SimulationWrangler::SetConfigurationParameters(const Parameters& config)
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting up configuration parameters."
        );
        sProcessMCTruth = config().ProcessMCTruth();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessMCTruth: " + std::to_string(sProcessMCTruth)
        );
        sProcessMCParticles = config().ProcessMCParticles();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessMCParticles: " + std::to_string(sProcessMCParticles)
        );
        sProcessSimEnergyDeposits = config().ProcessSimEnergyDeposits();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessSimEnergyDeposits: " + std::to_string(sProcessSimEnergyDeposits)
        );
        sProcessSimChannels = config().ProcessSimChannels();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessSimChannels: " + std::to_string(sProcessSimChannels)
        );
        sProcessRawDigits = config().ProcessRawDigits();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessRawDigits: " + std::to_string(sProcessRawDigits)
        );
        sProcessHits = config().ProcessHits();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessHits: " + std::to_string(sProcessHits)
        );
        sProcessOpDetBacktrackerRecords = config().ProcessOpDetBacktrackerRecords();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessOpDetBacktrackerRecords: " + std::to_string(sProcessOpDetBacktrackerRecords)
        );
        sProcessOpDetWaveforms = config().ProcessOpDetWaveforms();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ProcessOpDetWaveforms: " + std::to_string(sProcessOpDetWaveforms)
        );
        sSaveSimulationWrangler = config().SaveSimulationWrangler();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveSimulationWrangler: " + std::to_string(sSaveSimulationWrangler)
        );
        sSaveWirePlaneHits = config().SaveWirePlaneHits();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveWirePlaneHits: " + std::to_string(sSaveWirePlaneHits)
        );
        sSaveWirePlanePointCloud = config().SaveWirePlanePointCloud();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveWirePlanePointCloud: " + std::to_string(sSaveWirePlanePointCloud)
        );
        sSaveWirePlaneTrackTopology = config().SaveWirePlaneTrackTopology();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveWirePlaneTrackTopology: " + std::to_string(sSaveWirePlaneTrackTopology)
        );
        sSaveOpDetPointCloud = config().SaveOpDetPointCloud();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveOpDetPointCloud: " + std::to_string(sSaveOpDetPointCloud)
        );
        sSaveWirePlaneInductionFlag = config().SaveWirePlaneInductionFlag();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveWirePlaneInductionFlag: " + std::to_string(sSaveWirePlaneInductionFlag)
        );
        sADCThreshold = config().ADCThreshold();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ADCThreshold: " + std::to_string(sADCThreshold)
        );
        sVoxelizeEnergyDeposits = config().VoxelizeEnergyDeposits();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting VoxelizeEnergyDeposits: " + std::to_string(sVoxelizeEnergyDeposits)
        );
    }
    void SimulationWrangler::SetWirePlanePointCloudLabels(
        DetSimID_t detSimID, 
        TrackID_t trackID,
        TopologyLabelInt topologyLabel, 
        ParticleLabelInt particleLabel, 
        PhysicsMicroLabelInt physicsMicroLabel,
        PhysicsMesoLabelInt physicsMesoLabel,
        Int_t uniqueTopologyLabel,
        Int_t uniquePhysicsMicroLabel,
        Int_t uniquePhysicsMesoLabel,
        Bool_t inductionFlag
    )
    {
        auto track_index = sWirePlanePointCloud.GetIndex_TrackID(detSimID, trackID);
        auto physics_macro = LabelCast(sTrackID_PhysicsMacroLabelMap[trackID]);
        auto unique_physics_macro = sTrackID_UniquePhysicsMacroLabelMap[trackID];
        if(track_index != -1)
        {
            sWirePlanePointCloud.topology_labels[detSimID][track_index] = topologyLabel;
            sWirePlanePointCloud.particle_labels[detSimID][track_index] = particleLabel;
            if (!IsParticleLabel(particleLabel))
            {
                if (particleLabel > 2212) {
                    sWirePlanePointCloud.particle_labels[detSimID][track_index] = 1000000000;
                }
                else {
                    sWirePlanePointCloud.particle_labels[detSimID][track_index] = -1;
                }
            }
            sWirePlanePointCloud.physics_micro_labels[detSimID][track_index] = physicsMicroLabel;
            sWirePlanePointCloud.physics_meso_labels[detSimID][track_index] = physicsMesoLabel;
            sWirePlanePointCloud.physics_macro_labels[detSimID][track_index] = physics_macro;
            sWirePlanePointCloud.unique_particle_labels[detSimID][track_index] = trackID;
            sWirePlanePointCloud.unique_topology_labels[detSimID][track_index] = uniqueTopologyLabel;
            sWirePlanePointCloud.unique_physics_micro_labels[detSimID][track_index] = uniquePhysicsMicroLabel;
            sWirePlanePointCloud.unique_physics_meso_labels[detSimID][track_index] = uniquePhysicsMesoLabel;
            sWirePlanePointCloud.unique_physics_macro_labels[detSimID][track_index] = unique_physics_macro;
        }
        sWirePlanePointCloud.topology_label[detSimID] = topologyLabel;
        sWirePlanePointCloud.particle_label[detSimID] = particleLabel;
        if (!IsParticleLabel(particleLabel))
        {
            if (particleLabel > 2212) {
                sWirePlanePointCloud.particle_label[detSimID] = 1000000000;
            }
            else {
                sWirePlanePointCloud.particle_label[detSimID] = -1;
            }
        }
        sWirePlanePointCloud.physics_micro_label[detSimID] = physicsMicroLabel;
        sWirePlanePointCloud.physics_meso_label[detSimID] = physicsMesoLabel;
        sWirePlanePointCloud.physics_macro_label[detSimID]= physics_macro;
        sWirePlanePointCloud.unique_particle_label[detSimID] = trackID;
        sWirePlanePointCloud.unique_topology_label[detSimID] = uniqueTopologyLabel;
        sWirePlanePointCloud.unique_physics_micro_label[detSimID] = uniquePhysicsMicroLabel;
        sWirePlanePointCloud.unique_physics_meso_label[detSimID] = uniquePhysicsMesoLabel;
        sWirePlanePointCloud.unique_physics_macro_label[detSimID] = unique_physics_macro;
        sWirePlanePointCloud.induction_flag[detSimID] = inductionFlag;
    }
    void SimulationWrangler::SetEnergyDepositPointCloudLabels(
        EdepID_t edepID, 
        TrackID_t trackID,
        TopologyLabelInt topologyLabel, 
        ParticleLabelInt particleLabel, 
        PhysicsMicroLabelInt physicsMicroLabel,
        PhysicsMesoLabelInt physicsMesoLabel,
        Int_t uniqueTopologyLabel,
        Int_t uniquePhysicsMicroLabel,
        Int_t uniquePhysicsMesoLabel
    )
    {
        sEnergyDepositPointCloud.topology_label[edepID] = topologyLabel;
        sEnergyDepositPointCloud.particle_label[edepID] = particleLabel;
        sEnergyDepositPointCloud.physics_micro_label[edepID] = physicsMicroLabel;
        sEnergyDepositPointCloud.physics_meso_label[edepID] = physicsMesoLabel;
        sEnergyDepositPointCloud.unique_particle_label[edepID] = trackID;
        sEnergyDepositPointCloud.unique_topology_label[edepID] = uniqueTopologyLabel;
        sEnergyDepositPointCloud.unique_physics_micro_label[edepID] = uniquePhysicsMicroLabel;
        sEnergyDepositPointCloud.unique_physics_meso_label[edepID] = uniquePhysicsMesoLabel;
    }
    void SimulationWrangler::PrintParticleData(TrackID_t trackID)
    {
        auto particle = (*sMCParticleHandle)[sTrackID_ParticleIDMap[trackID]];
        std::cout << "## MCParticle #######################################\n";
        std::cout << "## TrackID:            [" << std::setw(25) << std::setfill('.') << trackID << "] ##\n";
        std::cout << "## PDG:                [" << std::setw(25) << std::setfill('.') << particle.PdgCode() << "] ##\n";
        std::cout << "## Energy [MeV]:       [" << std::setw(25) << std::setfill('.') << particle.E() << "] ##\n";
        std::cout << "## Process:            [" << std::setw(25) << std::setfill('.') << particle.Process() << "] ##\n";
        std::cout << "## Parent TrackID:     [" << std::setw(25) << std::setfill('.') << particle.Mother() << "] ##\n";
        std::cout << "## Parent PDG:         [" << std::setw(25) << std::setfill('.') << sTrackID_ParentPDGCodeMap[trackID] << "] ##\n";
        std::cout << "## Ancestor TrackID:   [" << std::setw(25) << std::setfill('.') << sTrackID_AncestorTrackIDMap[trackID] << "] ##\n";
        std::cout << "## Ancestor PDG:       [" << std::setw(25) << std::setfill('.') << sTrackID_AncestorPDGCodeMap[trackID] << "] ##\n";
        std::cout << "## Ancestor level:     [" << std::setw(25) << std::setfill('.') << sTrackID_AncestorLevelMap[trackID] << "] ##\n";
        std::cout << "## Progeny  [.....level] [...TrackID] [.......PDG] ##\n";
        auto progeny = sTrackID_ProgenyTrackIDMap[trackID];
        auto particle_level = sTrackID_AncestorLevelMap[trackID];
        for(auto progeny_track_id : progeny)
        {
            std::cout << "##          [" << std::setw(10) << std::setfill('.') << sTrackID_AncestorLevelMap[progeny_track_id] - particle_level << "] [";
            std::cout << std::setw(10) << std::setfill('.') << progeny_track_id << "] [";
            std::cout << std::setw(10) << std::setfill('.') << sTrackID_PDGCodeMap[progeny_track_id] << "] ##\n";
        }
        std::cout << "#####################################################" << std::endl;
    }
    void SimulationWrangler::PrintEdepData(EdepID_t edepID)
    {
        auto edep = (*sMCSimEnergyDepositHandle)[edepID];
        std::cout << "## MCSimEnergyDeposit ###############################\n";
        std::cout << "## EdepID:             [" << std::setw(25) << std::setfill('.') << edepID << "] ##\n";
        std::cout << "## TrackID:            [" << std::setw(25) << std::setfill('.') << edep.TrackID() << "] ##\n";
        std::cout << "## PDG:                [" << std::setw(25) << std::setfill('.') << sTrackID_PDGCodeMap[edep.TrackID()] << "] ##\n";
        std::cout << "## Energy [MeV]:       [" << std::setw(25) << std::setfill('.') << edep.Energy() << "] ##\n";
        std::cout << "## Process:            [" << std::setw(25) << std::setfill('.') << ProcessTypeToString[sEdepID_ProcessMap[edepID]] << "] ##\n";
        std::cout << "## MidPoint [x,y,z]:   [" << std::setw(7) << std::setfill('.') << edep.MidPointX() << ", ";
        std::cout << std::setw(7) << std::setfill('.') << edep.MidPointY() << ", ";
        std::cout << std::setw(7) << std::setfill('.') << edep.MidPointZ() << "] ##\n"; 
        std::cout << "#####################################################" << std::endl;
    }
    void SimulationWrangler::PrintDetSimData(DetSimID_t detsimID)
    {

    }
    void SimulationWrangler::ResetEvent()
    {
        sMCTruthHandles.clear();
        sPrimaries.clear();
        sEnergyDepositPointCloud.clear();
        sWirePlanePointCloud.clear();
        sWirePlaneTrackTopology.clear();
        sOpDetPointCloud.clear();

        sTrackID_GeneratorLabelMap.clear();
        sTrackID_ParticleIDMap.clear();
        sTrackID_PDGCodeMap.clear();
        sTrackID_ProcessMap.clear();
        sTrackID_EndProcessMap.clear();
        sTrackID_EnergyMap.clear();
        sTrackID_StartTimeMap.clear();
        sTrackID_DaughterTrackIDMap.clear();
        sTrackID_ProgenyTrackIDMap.clear();
        sTrackID_AncestryTrackIDMap.clear();
        sTrackID_EdepIDMap.clear();
        sTrackID_EdepProcessMap.clear();
        sTrackID_DetSimIDMap.clear();
        sTrackID_RandomDetSimIDMap.clear();

        sTrackID_ParentTrackIDMap.clear();
        sTrackID_ParentPDGCodeMap.clear();
        
        sTrackID_AncestorTrackIDMap.clear();
        sTrackID_AncestorLevelMap.clear();
        sTrackID_AncestorPDGCodeMap.clear();

        sEdepID_ProcessMap.clear();
        sEdepID_DetSimIDMap.clear();

        sChannelID_TDC_DetSimIDMap.clear();
        sDetSimID_EdepIDMap.clear();

        sTrackID_LabelingFunctionMap.clear();
    }
    void SimulationWrangler::ProcessEvent(
        const Parameters& config, art::Event const& event
    )
    {
        ResetEvent();
        if(sProcessMCTruth) 
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing MCTruth"
            );
            ProcessMCTruth(event, config().labels.get_PSet());
        }
        if(sProcessMCParticles)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing MCParticles"
            );
            ProcessMCParticles(event, config().LArGeantProducerLabel());
        }
        if(sProcessSimEnergyDeposits)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing SimEnergyDeposits"
            );
            ProcessSimEnergyDeposits(event, 
                config().SimEnergyDepositProducerLabel(), config().SimEnergyDepositInstanceLabel()
            );
        }
        if(sProcessSimChannels)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing SimChannels"
            );
            ProcessSimChannels(event, 
                config().SimChannelProducerLabel(), config().SimChannelInstanceLabel()
            );
        }
        if(sProcessRawDigits)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing RawDigits"
            );
            ProcessRawDigits(event,
                config().RawDigitProducerLabel(), config().RawDigitInstanceLabel()
            );
        }
        if(sProcessHits)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing Hits"
            );
            ProcessHits(event,
                config().SimChannelProducerLabel(), config().SimChannelInstanceLabel()
            );
        }
        if(sProcessOpDetBacktrackerRecords)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing OpDetBacktrackerRecords"
            );
            ProcessOpDetBacktrackerRecords(event, 
                config().OpDetBacktrackerRecordProducerLabel()
            );
        }
        if(sProcessOpDetWaveforms)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "processing OpDetWaveforms"
            );
            ProcessOpDetWaveforms(event,
                config().OpDetWaveformProducerLabel()
            );
        }
    }
    void SimulationWrangler::ProcessMCTruth(
        art::Event const& event, fhicl::ParameterSet const& generator_labels
    )
    {
        std::map<std::string, art::InputTag> generator_tags;
        for(std::string const& name : generator_labels.get_names()) {
            generator_tags[name] = generator_labels.get<art::InputTag>(name);
        }
        for(auto const& [key, tag] : generator_tags)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "collecting simb::MCTruth from input_tag <" + 
                tag.label() + ">"
            );
            if(!event.getByLabel(tag, sMCTruthHandle))
            {
                Logger::GetInstance("SimulationWrangler")->warning(
                    "no input_tag matching " + tag.label() + 
                    " for simb::MCTruth"
                );
            }
            else 
            {
                sMCTruthHandles.emplace_back(event.getHandle<std::vector<simb::MCTruth>>(
                    tag
                ));
                sMCTruthHandleLabels.emplace_back(key);
                if(!sMCTruthHandles.back().isValid()) 
                {
                    Logger::GetInstance("SimulationWrangler")->error(
                        "data product " + tag.label() + 
                        " for simb::MCTruth is invalid!"
                    );
                    exit(0);
                }
            }
        }
    }
    void SimulationWrangler::ProcessMCParticles(
        art::Event const& event, art::InputTag input_tag
    )
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "collecting simb::MCParticle from label <" + 
            input_tag.label() + ">"
        );
        if(!event.getByLabel(input_tag, sMCParticleHandle))
        {
            Logger::GetInstance("SimulationWrangler")->error(
                "no label matching " + input_tag.label() + 
                " for simb::MCParticle!"
            );
            exit(0);
        }
        else 
        {
            sMCParticleHandle = event.getHandle<std::vector<simb::MCParticle>>(
                input_tag
            );
            if(!sMCParticleHandle.isValid()) 
            {
                Logger::GetInstance("SimulationWrangler")->error(
                    "data product " + input_tag.label() + 
                    " for simb::MCParticle is invalid!"
                );
                exit(0);
            }
        }
        Int_t particle_index = 0;
        Logger::GetInstance("SimulationWrangler")->trace(
            "creating particle track ID maps for " +
            std::to_string((*sMCParticleHandle).size()) + 
            " <simb::MCParticle>s."
        );
        for (auto particle : *sMCParticleHandle)
        {
            sTrackID_ParticleIDMap[particle.TrackId()] = particle_index;
            particle_index += 1;

            sTrackID_GeneratorLabelMap[particle.TrackId()] = GeneratorLabel::None;
            sTrackID_PhysicsMacroLabelMap[particle.TrackId()] = PhysicsMacroLabel::Undefined;
            sTrackID_UniquePhysicsMacroLabelMap[particle.TrackId()] = -1;
            sTrackID_PDGCodeMap[particle.TrackId()] = particle.PdgCode();
            sTrackID_ProcessMap[particle.TrackId()] = TrajectoryStringToProcessType[particle.Process()];
            sTrackID_EndProcessMap[particle.TrackId()] = TrajectoryStringToProcessType[particle.EndProcess()];
            sTrackID_EnergyMap[particle.TrackId()] = particle.E();
            sTrackID_StartTimeMap[particle.TrackId()] = particle.T();

            sTrackID_ParentTrackIDMap[particle.TrackId()] = particle.Mother();
            std::vector<Int_t> daughters = {};
            for(auto ii = 0; ii < particle.NumberDaughters(); ii++) {
                daughters.emplace_back(particle.Daughter(ii));
            }
            sTrackID_DaughterTrackIDMap[particle.TrackId()] = daughters;
            sTrackID_ProgenyTrackIDMap[particle.TrackId()] = {};
            sTrackID_DescendantTrackIDMap[particle.TrackId()] = daughters;

            // construct ancestry map
            std::vector<Int_t> ancestry = {};
            Int_t mother = particle.Mother();
            Int_t track_id = particle.TrackId();
            Int_t level = 0;
            while (mother != 0)
            {
                level += 1;
                track_id = mother;
                ancestry.emplace_back(mother);
                mother = sTrackID_ParentTrackIDMap[track_id];
                if(level > 1) {
                    sTrackID_ProgenyTrackIDMap[mother].emplace_back(particle.TrackId());
                    sTrackID_DescendantTrackIDMap[mother].emplace_back(particle.TrackId());
                }
            }
            sTrackID_AncestorLevelMap[particle.TrackId()] = level;
            sTrackID_AncestryTrackIDMap[particle.TrackId()] = ancestry;

            if (level == 0) {
                sPrimaries.emplace_back(particle.TrackId());
                sTrackID_ParentPDGCodeMap[particle.TrackId()] = 0;
                sTrackID_AncestorPDGCodeMap[particle.TrackId()] = 0;
                sTrackID_AncestorTrackIDMap[particle.TrackId()] = 0;
            }
            else {
                sTrackID_ParentPDGCodeMap[particle.TrackId()] = sTrackID_PDGCodeMap[particle.Mother()];
                sTrackID_AncestorPDGCodeMap[particle.TrackId()] = sTrackID_PDGCodeMap[track_id];
                sTrackID_AncestorTrackIDMap[particle.TrackId()] = track_id;
            }

            // initialize edep maps
            sTrackID_EdepIDMap[particle.TrackId()] = {};
            sTrackID_EdepProcessMap[particle.TrackId()] = {};
            sTrackID_DetSimIDMap[particle.TrackId()] = {};
        }
        for(size_t jj = 0; jj < sMCTruthHandles.size(); jj++)
        {
            for(auto truth : *sMCTruthHandles[jj])
            {
                /**
                 * MCTruth stores MCParticles starting with trackID = 0,
                 * rather than Geant4 which starts with trackID = 1.
                */
                if(truth.NParticles() == 0)
                {
                    Logger::GetInstance("SimulationWrangler")->trace(
                        "MCTruth for " + sMCTruthHandleLabels[jj] + 
                        " contains no simulated particles."
                    );
                    continue;
                }
                Logger::GetInstance("SimulationWrangler")->trace(
                    "adding labels of type " + 
                    sMCTruthHandleLabels[jj] + 
                    " for " + std::to_string(truth.NParticles()) + 
                    " particles starting with track ID = " + 
                    std::to_string(truth.GetParticle(0).TrackId()+1)
                );
            }
        }
        for(auto primary : sPrimaries)
        {
            auto particle = (*sMCParticleHandle)[sTrackID_ParticleIDMap[primary]];
            auto position = particle.Position();
            Int_t pdg_code = particle.PdgCode();
            bool found = false;
            for(size_t jj = 0; jj < sMCTruthHandles.size(); jj++)
            {
                for(auto truth : *sMCTruthHandles[jj])
                {
                    for(Int_t ii = 0; ii < truth.NParticles(); ii++)
                    {
                        /**
                         * For some reason, the cosmics generator uses a different precision than the
                         * MCParticle coming out of Geant4, so if you don't correct for this
                         * you will miss the association between Cosmics and their associated
                         * MCParticle TrackID.  We therefore round everything to six full digits to
                         * coincide with the Geant4 side of things.
                         */
                        auto truth_position = truth.GetParticle(ii).Position();
                        if(
                            round(truth_position[0]*pow(10,6))/pow(10,6) == round(position[0]*pow(10,6))/pow(10,6) &&
                            round(truth_position[1]*pow(10,6))/pow(10,6) == round(position[1]*pow(10,6))/pow(10,6) &&
                            round(truth_position[2]*pow(10,6))/pow(10,6) == round(position[2]*pow(10,6))/pow(10,6)
                        )
                        {
                            sTrackID_GeneratorLabelMap[primary] = sGeneratorMap[sMCTruthHandleLabels[jj]];
                            found = true;
                            break;
                        }
                    }
                    if(found) {
                        break;
                    }
                }
                if(found) {
                    break;
                }
            }
            if(!found)
            {
                Logger::GetInstance("SimulationWrangler")->warning(
                    "couldn't find mc truth for primary with TrackID: " +
                    std::to_string(primary) + " - PDGCode: " +
                    std::to_string(pdg_code) + " - Position: (" +
                    std::to_string(position[0]) + "," + std::to_string(position[1]) +
                    "," + std::to_string(position[2])
                ); 
            }
        }
    }
    void SimulationWrangler::ProcessSimEnergyDeposits(art::Event const& event, 
        art::InputTag producer_label, art::InputTag instance_label
    )
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "collecting sim::SimEnergyDeposit from label <" + 
            producer_label.label() + ":" + instance_label.label() + ">"
        );
        if(!event.getByLabel(
            art::InputTag(producer_label.label(), instance_label.label()), 
            sMCSimEnergyDepositHandle
        ))
        {
            Logger::GetInstance("SimulationWrangler")->error(
                "no label matching " + producer_label.label() + ":" + 
                instance_label.label() + " for sim::SimEnergyDeposit!"
            );
            exit(0);
        }
        else 
        {
            sMCSimEnergyDepositHandle = event.getHandle<std::vector<sim::SimEnergyDeposit>>(
                art::InputTag(
                    producer_label.label(), instance_label.label()
                )
            );
            if(!sMCSimEnergyDepositHandle.isValid()) 
            {
                Logger::GetInstance("SimulationWrangler")->error(
                    "data product " + producer_label.label() + ":" + 
                    instance_label.label() + " for simb::SimEnergyDeposit is invalid!"
                );
                exit(0);
            }
        }
        Int_t edep_index = 0;
        Logger::GetInstance("SimulationWrangler")->trace(
            "creating particle edep ID maps for " +
            std::to_string((*sMCSimEnergyDepositHandle).size()) + 
            " <sim::SimEnergyDeposit>s."
        );
        for(auto edep : *sMCSimEnergyDepositHandle)
        {
            auto track_id = edep.TrackID();
            auto edep_t = edep.T();
            auto edep_x = edep.X();
            auto edep_y = edep.Y();
            auto edep_z = edep.Z();
            auto energy = edep.E();
            auto num_photons = edep.NumPhotons();
            auto num_electrons = edep.NumElectrons();

            ProcessType process = DetermineEdepProcess(edep);

            sEnergyDepositPointCloud.AddPoint(
                track_id, edep_t, edep_x, edep_y, edep_z,
                energy, num_photons, num_electrons, process
            );
            sTrackID_EdepIDMap[track_id].emplace_back(edep_index);
            sTrackID_EdepProcessMap[track_id].emplace_back(process);
            sEdepID_DetSimIDMap[edep_index] = {};

            edep_index += 1;
        }
    }
    void SimulationWrangler::ProcessSimChannels(art::Event const& event,
        art::InputTag producer_label, art::InputTag instance_label
    )
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "collecting sim::SimChannel from label <" + 
            producer_label.label() + ":" + instance_label.label() + ">"
        );
        if(!event.getByLabel(
            art::InputTag(producer_label.label(), instance_label.label()),
            sMCSimChannelHandle
        ))
        {
            Logger::GetInstance("SimulationWrangler")->error(
                "no label matching " + producer_label.label() + ":" + 
                instance_label.label() + " for sim::SimChannel!"
            );
            exit(0);
        }
        else 
        {
            sMCSimChannelHandle = event.getHandle<std::vector<sim::SimChannel>>(
                art::InputTag(
                    producer_label.label(), instance_label.label()
                )
            );
            if(!sMCSimChannelHandle.isValid()) 
            {
                Logger::GetInstance("SimulationWrangler")->error(
                    "data product " + producer_label.label() + ":" + 
                    instance_label.label() + " for sim::SimChannel is invalid!"
                );
                exit(0);
            }
        }
    }
    void SimulationWrangler::ProcessRawDigits(art::Event const& event,
        art::InputTag producer_label, art::InputTag instance_label
    )
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "collecting raw::RawDigit from label <" + 
            producer_label.label() + ":" + instance_label.label() + ">"
        );
        if(!event.getByLabel(
            art::InputTag(producer_label.label(),instance_label.label()),
            sMCRawDigitHandle
        ))
        {
            Logger::GetInstance("SimulationWrangler")->error(
                "no label matching " + producer_label.label() + ":" + 
                instance_label.label() + " for raw::RawDigit!"
            );
            exit(0);
        }
        else 
        {
            sMCRawDigitHandle = event.getHandle<std::vector<raw::RawDigit>>(
                art::InputTag(
                    producer_label.label(), instance_label.label()
                )
            );
            if(!sMCRawDigitHandle.isValid()) 
            {
                Logger::GetInstance("SimulationWrangler")->error(
                    "data product " + producer_label.label() + ":" + 
                    instance_label.label() + " for raw::RawDigit is invalid!"
                );
                exit(0);
            }
        }
        detinfo::DetectorClocksData const clock_data(
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        ); 
        Int_t digit_index = 0;
        Logger::GetInstance("SimulationWrangler")->trace(
            "creating detector simulation and particle ID maps for " +
            std::to_string((*sMCRawDigitHandle).size()) + 
            " <raw::RawDigit>s."
        );
        sNumberOfTDCs = sMCRawDigitHandle->at(0).Samples();
        for(auto digit : *sMCRawDigitHandle)
        {
            // Get the channel number for this digit, number of samples,
            // and the pedestal value so that we can uncompress and
            // remove the pedestal.
            raw::ChannelID_t channel = digit.Channel();
            int num_samples = digit.Samples();
            int pedestal = (int)digit.GetPedestal();
            
            // uncompress the digits and remove the pedestal
            std::vector<short> uncompressed(num_samples);
            raw::Uncompress(
                digit.ADCs(), uncompressed, 
                pedestal, digit.Compression()
            );
            for (int ii = 0; ii < num_samples; ii++) {
                uncompressed[ii] -= pedestal;
            }
            sim::SimChannel truth_channel = (*sMCSimChannelHandle)[channel]; 

            // iterate over each tdc value
            for(int l=0; l < num_samples; l++) 
            {
                auto const& trackIDsAndEnergy = truth_channel.TrackIDsAndEnergies(l, l);
                /**
                 * This step distinguishes noise from true MC particles.
                 * If the input is noise, pass the detector output to a
                 * noise variable, otherwise, attach the output to the
                 * associated primary.
                 */
                if(trackIDsAndEnergy.size() == 0)
                {
                    if(std::abs(uncompressed[l]) >= sADCThreshold)
                    {
                        sWirePlanePointCloud.AddPoint(
                            clock_data,
                            trackIDsAndEnergy,
                            l,
                            channel,
                            (Int_t) (uncompressed[l]),
                            true
                        );
                        digit_index += 1;
                    }
                }
                else
                {
                    sWirePlanePointCloud.AddPoint(
                        clock_data,
                        trackIDsAndEnergy,
                        l,
                        channel,
                        (Int_t) (uncompressed[l]),
                        false
                    );
                    // associate this detector simulation with a particle
                    for(auto track : trackIDsAndEnergy)
                    {
                        if(track.trackID > 0) {
                            sTrackID_DetSimIDMap[track.trackID].emplace_back(digit_index);
                        }
                        else {
                            sTrackID_RandomDetSimIDMap[track.trackID] = {digit_index};
                        }
                    }
                    // determine the edeps associated with this detector simulation
                    sDetSimID_EdepIDMap[digit_index] = DetermineDetectorSimulationEdeps(
                        trackIDsAndEnergy,
                        digit_index
                    );
                    sChannelID_TDC_DetSimIDMap[std::make_pair(channel, l)] = digit_index;
                    digit_index += 1;
                }
            }
        }
    }
    void SimulationWrangler::ProcessHits(art::Event const& event,
        art::InputTag producer_label, art::InputTag instance_label
    )
    {
        detinfo::DetectorClocksData const clock_data(
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        );
        for(auto channel : *sMCSimChannelHandle) 
        {
            /** Make a list of unique track_ids on this channel and
             * group tdcs/ne accordingly.  This is not necessarily
             * the best way to do this, since there are cases where
             * several hits from the same particle could be on the
             * same channel.  
             */
            std::map<TrackID_t, std::vector<Int_t>> track_id_tdcs;
            std::map<TrackID_t, std::vector<Double_t>> track_id_nes;
            for(auto tdcide : channel.TDCIDEMap()) 
            {
                for(auto ide : tdcide.second)
                {   
                    if(track_id_tdcs.count(ide.trackID)) {
                        track_id_tdcs[ide.trackID].emplace_back(tdcide.first);
                        track_id_nes[ide.trackID].emplace_back(ide.numElectrons);
                    }
                    else {
                        track_id_tdcs[ide.trackID] = {tdcide.first};
                        track_id_nes[ide.trackID] = {ide.numElectrons};
                    }
                }
            }
            // Now calculate the mean of each distribution of tdc/ne values
            // for each track id.
            for(auto const& [key, val] : track_id_tdcs)
            {
                std::vector<Double_t> num_electrons = track_id_nes[key];
                Double_t tdc_mean = std::accumulate(val.begin(), val.end(), 0.0) / val.size();
                Double_t tdc_closest = 10e10;
                Double_t temp_tdc_rms = 0.0;
                for(auto tdc : val) {
                    temp_tdc_rms += (tdc - tdc_mean) * (tdc - tdc_mean);
                }
                for(auto tdc: val) {
                    if (abs(tdc - tdc_mean) < abs(tdc_closest - tdc_mean)) {
                        tdc_closest = tdc;
                    }
                }
                Double_t tdc_rms = std::sqrt(temp_tdc_rms / val.size());
                Double_t tdc_amplitude = *std::max_element(num_electrons.begin(), num_electrons.end());
                Double_t tdc_charge = std::accumulate(num_electrons.begin(), num_electrons.end(), 0.0);

                // Find the associated (channel, tdc) DetSimID.
                DetSimID_t detsim_id = sChannelID_TDC_DetSimIDMap[
                    std::make_pair(channel.Channel(), clock_data.TPCTDC2Tick(tdc_closest))
                ];
                sWirePlanePointCloud.AddHit(
                    detsim_id,
                    tdc_mean,
                    tdc_rms,
                    tdc_amplitude,
                    tdc_charge
                );
            }
        }
    }
    void SimulationWrangler::ProcessOpDetBacktrackerRecords(art::Event const& event,
        art::InputTag producer_label
    )
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "collecting sim::OpDetBacktrackerRecord from label <" + 
            producer_label.label() +  ">"
        );
        if(!event.getByLabel(
            art::InputTag(producer_label.label()),
            sMCOpDetBacktrackerRecordHandle
        ))
        {
            Logger::GetInstance("SimulationWrangler")->error(
                "no label matching " + producer_label.label() + 
                " for sim::OpDetBacktrackerRecord!"
            );
            exit(0);
        }
        else 
        {
            sMCOpDetBacktrackerRecordHandle = event.getHandle<std::vector<sim::OpDetBacktrackerRecord>>(
                art::InputTag(
                    producer_label.label()
                )
            );
            if(!sMCOpDetBacktrackerRecordHandle.isValid()) 
            {
                Logger::GetInstance("SimulationWrangler")->error(
                    "data product " + producer_label.label() +
                    " for sim::OpDetBacktrackerRecord is invalid!"
                );
                exit(0);
            }
        }
        Int_t opdetbacktracker_id = 0;
        for(auto opdet_record : *sMCOpDetBacktrackerRecordHandle)
        {
            sOpDetChannelID_OpDetBacktrackerIDMap[opdet_record.OpDetNum()] = opdetbacktracker_id;
            opdetbacktracker_id += 1;
        }
    }
    void SimulationWrangler::ProcessOpDetWaveforms(art::Event const& event,
        art::InputTag producer_label
    )
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "collecting raw::OpDetWaveform from label <" + 
            producer_label.label() + ">"
        );
        if(!event.getByLabel(
            art::InputTag(producer_label.label()),
            sMCOpDetWaveformHandle
        ))
        {
            Logger::GetInstance("SimulationWrangler")->error(
                "no label matching " + producer_label.label() + 
                " for raw::OpDetWaveform!"
            );
            exit(0);
        }
        else 
        {
            sMCOpDetWaveformHandle = event.getHandle<std::vector<raw::OpDetWaveform>>(
                art::InputTag(producer_label.label())
            );
            if(!sMCOpDetWaveformHandle.isValid()) 
            {
                Logger::GetInstance("SimulationWrangler")->error(
                    "data product " + producer_label.label() + 
                    " for raw::OpDetWaveform is invalid!"
                );
                exit(0);
            }
        }
        Logger::GetInstance("SimulationWrangler")->trace(
            "creating optical detector simulation and particle ID maps for " +
            std::to_string((*sMCOpDetWaveformHandle).size()) + 
            " <raw::OpDetWaveform>s."
        );
        detinfo::DetectorClocksData const clock_data(
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
        ); 
        for(auto waveform : *sMCOpDetWaveformHandle)
        {
            /**
             * (2) figure out how to add labels to these points.
            */
            auto adc = waveform.Waveform();
            auto channel = waveform.ChannelNumber();
            auto time_stamp = waveform.TimeStamp();
            auto time_tick = clock_data.Time2Tick(time_stamp);
            auto op_det_backtracker = sMCOpDetBacktrackerRecordHandle->at(GetOpDetBacktrackerID_OpDetChannelID(channel));
            auto op_det_map = op_det_backtracker.timePDclockSDPsMap();
            std::cout << "ADC: " << adc.size() << ", Map: " << op_det_map.size();
            for(size_t ii = 0; ii < adc.size(); ii++)
            {
                auto const& trackIDsAndEnergy = op_det_backtracker.TrackIDsAndEnergies(ii, ii);
                sOpDetPointCloud.AddPoint(
                    channel, 
                    time_tick + ii * (1.0/(150.0 * 1e6)), 
                    adc[ii], false
                );
            }
            // auto channel = waveform.ChannelNumber();
            // auto time_stamp = waveform.TimeStamp();
            // auto time_tick = clock_data.Time2Tick(time_stamp);
            // auto time_tdc = clock_data.TPCTick2TDC(time_tick);
            // std::cout << "channel: " << channel << ", time_stamp: " << time_stamp;
            // std::cout << ", tick: " << time_tick << ", tdc: " << time_tdc << ", num adcs: " << adc.size() << std::endl;
        }
    }
    TrackID_List SimulationWrangler::GetPrimaries_GeneratorLabel(GeneratorLabel label)
    {
        TrackID_List primaries;
        for(auto primary : sPrimaries)
        {
            if(sTrackID_GeneratorLabelMap[primary] == label) {
                primaries.emplace_back(primary);
            }
        }
        return primaries;
    }
    TrackID_List SimulationWrangler::GetPrimaries_PDGCode(Int_t pdg)
    {
        TrackID_List primaries;
        for(auto primary : sPrimaries)
        {
            if (sTrackID_PDGCodeMap[primary] == pdg) {
                primaries.emplace_back(primary);
            }
        }
        return primaries;
    }
    TrackID_List SimulationWrangler::GetPrimaries_AbsPDGCode(Int_t pdg)
    {
        TrackID_List primaries;
        for(auto primary : sPrimaries)
        {
            if (std::abs(sTrackID_PDGCodeMap[primary]) == std::abs(pdg)) {
                primaries.emplace_back(primary);
            }
        }
        return primaries;
    }
    TrackID_List SimulationWrangler::GetPrimaries_Process(ProcessType process)
    {
        TrackID_List primaries;
        for(auto primary : sPrimaries)
        {
            if (sTrackID_ProcessMap[primary] == process) {
                primaries.emplace_back(primary);
            }
        }
        return primaries;
    }
    TrackID_List SimulationWrangler::GetPrimaries_EndProcess(ProcessType process)
    {
        TrackID_List primaries;
        for(auto primary : sPrimaries)
        {
            if (sTrackID_EndProcessMap[primary] == process) {
                primaries.emplace_back(primary);
            }
        }
        return primaries;
    }
    TrackID_List SimulationWrangler::GetTrackID_GeneratorLabel(GeneratorLabel label)
    {
        TrackID_List particles;
        for(auto const& [key, value] : sTrackID_GeneratorLabelMap)
        {
            if(value == label) {
                particles.emplace_back(key);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::GetTrackID_PDGCode(Int_t pdg)
    {
        TrackID_List particles;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(value == pdg) {
                particles.emplace_back(key);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::GetTrackID_AbsPDGCode(Int_t pdg)
    {
        TrackID_List particles;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(std::abs(value) == std::abs(pdg)) {
                particles.emplace_back(key);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::GetTrackID_Process(ProcessType process)
    {
        TrackID_List particles;
        for(auto const& [key, value] : sTrackID_ProcessMap)
        {
            if(value == process) {
                particles.emplace_back(key);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::GetTrackID_EndProcess(ProcessType process)
    {
        TrackID_List particles;
        for(auto const& [key, value] : sTrackID_EndProcessMap)
        {
            if(value == process) {
                particles.emplace_back(key);
            }
        }
        return particles;
    }
    EdepID_Collection SimulationWrangler::GetEdepID_TrackID(TrackID_List trackIDList)
    {
        EdepID_Collection edep;
        for(auto track_id : trackIDList)
        {
            edep.emplace_back(sTrackID_EdepIDMap[track_id]);
        }
        return edep;
    }
    std::vector<EdepID_Collection> SimulationWrangler::GetEdepID_TrackID(TrackID_Collection trackIDCollection)
    {
        std::vector<EdepID_Collection> edep;
        for(auto track_id : trackIDCollection)
        {
            edep.emplace_back(GetEdepID_TrackID(track_id));
        }
        return edep;
    }
    DetSimID_Collection SimulationWrangler::GetDetSimID_TrackID(TrackID_List trackIDList)
    {
        DetSimID_Collection detsim;
        for(auto track_id : trackIDList)
        {
            detsim.emplace_back(sTrackID_DetSimIDMap[track_id]);
        }
        return detsim;
    }

    TrackID_List SimulationWrangler::GetStartTime_TrackID(TrackID_List trackIDList)
    {
        TrackID_List start_time;
        for(auto track_id : trackIDList)
        {
            start_time.emplace_back(sTrackID_StartTimeMap[track_id]);
        }
        return start_time;
    }

    std::vector<DetSimID_Collection> SimulationWrangler::GetDetSimID_TrackID(TrackID_Collection trackIDCollection)
    {
        std::vector<DetSimID_Collection> detsim;
        for(auto track_id : trackIDCollection)
        {
            detsim.emplace_back(GetDetSimID_TrackID(track_id));
        }
        return detsim;
    }
    TrackID_List SimulationWrangler::FilterTrackID_PDGCode(TrackID_List& trackIDList, Int_t pdg)
    {
        TrackID_List particles;
        for(auto track_id : trackIDList)
        {
            if(sTrackID_PDGCodeMap[track_id] == pdg) {
                particles.emplace_back(track_id);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::FilterTrackID_NotPDGCode(TrackID_List& trackIDList, Int_t pdg)
    {
        TrackID_List particles;
        for(auto track_id : trackIDList)
        {
            if(sTrackID_PDGCodeMap[track_id] != pdg) {
                particles.emplace_back(track_id);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::FilterTrackID_AbsPDGCode(TrackID_List& trackIDList, Int_t pdg)
    {
        TrackID_List particles;
        for(auto track_id : trackIDList)
        {
            if(std::abs(sTrackID_PDGCodeMap[track_id]) == std::abs(pdg)) {
                particles.emplace_back(track_id);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::FilterTrackID_NotAbsPDGCode(TrackID_List& trackIDList, Int_t pdg)
    {
        TrackID_List particles;
        for(auto track_id : trackIDList)
        {
            if(std::abs(sTrackID_PDGCodeMap[track_id]) != std::abs(pdg)) {
                particles.emplace_back(track_id);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::FilterTrackID_Process(TrackID_List& trackIDList, ProcessType process)
    {
        TrackID_List particles;
        for(auto track_id : trackIDList)
        {
            if(sTrackID_ProcessMap[track_id] == process) {
                particles.emplace_back(track_id);
            }
        }
        return particles;
    }
    TrackID_List SimulationWrangler::FilterTrackID_NotProcess(TrackID_List& trackIDList, ProcessType process)
    {
        TrackID_List particles;
        for(auto track_id : trackIDList)
        {
            if(sTrackID_ProcessMap[track_id] != process) {
                particles.emplace_back(track_id);
            }
        }
        return particles;
    }
    TrackID_Collection SimulationWrangler::FilterTrackID_PDGCode(TrackID_Collection& trackIDCollection, Int_t pdg)
    {
        TrackID_Collection particles;
        for(auto track_ids : trackIDCollection)
        {
            particles.emplace_back(FilterTrackID_PDGCode(track_ids, pdg));
        }
        return particles;
    }
    TrackID_Collection SimulationWrangler::FilterTrackID_NotPDGCode(TrackID_Collection& trackIDCollection, Int_t pdg)
    {
        TrackID_Collection particles;
        for(auto track_ids : trackIDCollection)
        {
            particles.emplace_back(FilterTrackID_NotPDGCode(track_ids, pdg));
        }
        return particles;
    }
    TrackID_Collection SimulationWrangler::FilterTrackID_AbsPDGCode(TrackID_Collection& trackIDCollection, Int_t pdg)
    {
        TrackID_Collection particles;
        for(auto track_ids : trackIDCollection)
        {
            particles.emplace_back(FilterTrackID_AbsPDGCode(track_ids, pdg));
        }
        return particles;
    }
    TrackID_Collection SimulationWrangler::FilterTrackID_NotAbsPDGCode(TrackID_Collection& trackIDCollection, Int_t pdg)
    {
        TrackID_Collection particles;
        for(auto track_ids : trackIDCollection)
        {
            particles.emplace_back(FilterTrackID_NotAbsPDGCode(track_ids, pdg));
        }
        return particles;
    }
    TrackID_Collection SimulationWrangler::FilterTrackID_Process(TrackID_Collection& trackIDCollection, ProcessType process)
    {
        TrackID_Collection particles;
        for(auto track_ids : trackIDCollection)
        {
            particles.emplace_back(FilterTrackID_Process(track_ids, process));
        }
        return particles;
    }
    TrackID_Collection SimulationWrangler::FilterTrackID_NotProcess(TrackID_Collection& trackIDCollection, ProcessType process)
    {
        TrackID_Collection particles;
        for(auto track_ids : trackIDCollection)
        {
            particles.emplace_back(FilterTrackID_NotProcess(track_ids, process));
        }
        return particles;
    }

    TrackID_Collection SimulationWrangler::GetDaughterTrackID_GeneratorLabel(GeneratorLabel label)
    {
        TrackID_Collection daughters;
        for(auto const& [key, value] : sTrackID_GeneratorLabelMap)
        {
            if(value == label) {
                daughters.emplace_back(sTrackID_DaughterTrackIDMap[key]);
            }
        }
        return daughters;
    }
    TrackID_Collection SimulationWrangler::GetDaughterTrackID_PDGCode(Int_t pdg)
    {
        TrackID_Collection daughters;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(value == pdg) {
                daughters.emplace_back(sTrackID_DaughterTrackIDMap[key]);
            }
        }
        return daughters;
    }
    TrackID_Collection SimulationWrangler::GetDaughterTrackID_AbsPDGCode(Int_t pdg)
    {
        TrackID_Collection daughters;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(std::abs(value) == std::abs(pdg)) {
                daughters.emplace_back(sTrackID_DaughterTrackIDMap[key]);
            }
        }
        return daughters;
    }
    TrackID_Collection SimulationWrangler::GetDaughterTrackID_Process(ProcessType process)
    {
        TrackID_Collection daughters;
        for(auto const& [key, value] : sTrackID_ProcessMap)
        {
            if(value == process) {
                daughters.emplace_back(sTrackID_DaughterTrackIDMap[key]);
            }
        }
        return daughters;
    }
    TrackID_Collection SimulationWrangler::GetDaughterTrackID_EndProcess(ProcessType process)
    {
        TrackID_Collection daughters;
        for(auto const& [key, value] : sTrackID_EndProcessMap)
        {
            if(value == process) {
                daughters.emplace_back(sTrackID_DaughterTrackIDMap[key]);
            }
        }
        return daughters;
    }
    TrackID_Collection SimulationWrangler::GetDaughterTrackID_TrackID(TrackID_List trackIDs)
    {
        TrackID_Collection daughters;
        for(auto track_id : trackIDs)
        {
            daughters.emplace_back(sTrackID_DaughterTrackIDMap[track_id]);
        }
        return daughters;
    }
    TrackID_Collection SimulationWrangler::GetProgenyTrackID_GeneratorLabel(GeneratorLabel label)
    {
        TrackID_Collection progeny;
        for(auto const& [key, value] : sTrackID_GeneratorLabelMap)
        {
            if(value == label) {
                progeny.emplace_back(sTrackID_ProgenyTrackIDMap[key]);
            }
        }
        return progeny;
    }
    TrackID_Collection SimulationWrangler::GetProgenyTrackID_PDGCode(Int_t pdg)
    {
        TrackID_Collection progeny;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(value == pdg) {
                progeny.emplace_back(sTrackID_ProgenyTrackIDMap[key]);
            }
        }
        return progeny;
    }
    TrackID_Collection SimulationWrangler::GetProgenyTrackID_AbsPDGCode(Int_t pdg)
    {
        TrackID_Collection progeny;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(std::abs(value) == std::abs(pdg)) {
                progeny.emplace_back(sTrackID_ProgenyTrackIDMap[key]);
            }
        }
        return progeny;
    }
    TrackID_Collection SimulationWrangler::GetProgenyTrackID_Process(ProcessType process)
    {
        TrackID_Collection progeny;
        for(auto const& [key, value] : sTrackID_ProcessMap)
        {
            if(value == process) {
                progeny.emplace_back(sTrackID_ProgenyTrackIDMap[key]);
            }
        }
        return progeny;
    }
    TrackID_Collection SimulationWrangler::GetProgenyTrackID_EndProcess(ProcessType process)
    {
        TrackID_Collection progeny;
        for(auto const& [key, value] : sTrackID_EndProcessMap)
        {
            if(value == process) {
                progeny.emplace_back(sTrackID_ProgenyTrackIDMap[key]);
            }
        }
        return progeny;
    }
    TrackID_Collection SimulationWrangler::GetProgenyTrackID_TrackID(TrackID_List trackIDs)
    {
        TrackID_Collection progeny;
        for(auto track_id : trackIDs)
        {
            progeny.emplace_back(sTrackID_ProgenyTrackIDMap[track_id]);
        }
        return progeny;
    }
    TrackID_Collection SimulationWrangler::GetDescendantTrackID_GeneratorLabel(GeneratorLabel label)
    {
        TrackID_Collection descendant;
        for(auto const& [key, value] : sTrackID_GeneratorLabelMap)
        {
            if(value == label) {
                descendant.emplace_back(sTrackID_DescendantTrackIDMap[key]);
            }
        }
        return descendant;
    }
    TrackID_Collection SimulationWrangler::GetDescendantTrackID_PDGCode(Int_t pdg)
    {
        TrackID_Collection descendant;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(value == pdg) {
                descendant.emplace_back(sTrackID_DescendantTrackIDMap[key]);
            }
        }
        return descendant;
    }
    TrackID_Collection SimulationWrangler::GetDescendantTrackID_AbsPDGCode(Int_t pdg)
    {
        TrackID_Collection descendant;
        for(auto const& [key, value] : sTrackID_PDGCodeMap)
        {
            if(std::abs(value) == std::abs(pdg)) {
                descendant.emplace_back(sTrackID_DescendantTrackIDMap[key]);
            }
        }
        return descendant;
    }
    TrackID_Collection SimulationWrangler::GetDescendantTrackID_Process(ProcessType process)
    {
        TrackID_Collection descendant;
        for(auto const& [key, value] : sTrackID_ProcessMap)
        {
            if(value == process) {
                descendant.emplace_back(sTrackID_DescendantTrackIDMap[key]);
            }
        }
        return descendant;
    }
    TrackID_Collection SimulationWrangler::GetDescendantTrackID_EndProcess(ProcessType process)
    {
        TrackID_Collection descendant;
        for(auto const& [key, value] : sTrackID_EndProcessMap)
        {
            if(value == process) {
                descendant.emplace_back(sTrackID_DescendantTrackIDMap[key]);
            }
        }
        return descendant;
    }
    TrackID_Collection SimulationWrangler::GetDescendantTrackID_TrackID(TrackID_List trackIDs)
    {
        TrackID_Collection descendant;
        for(auto track_id : trackIDs)
        {
            descendant.emplace_back(sTrackID_DescendantTrackIDMap[track_id]);
        }
        return descendant;
    }

    ProcessType SimulationWrangler::DetermineEdepProcess(const sim::SimEnergyDeposit& edep)
    {
        std::string process = "NotDefined";
        auto mc_particle = GetMCParticle_TrackID(edep.TrackID());
        simb::MCTrajectory trajectory = mc_particle.Trajectory();
        auto trajectory_processes = trajectory.TrajectoryProcesses();
        for(size_t ii = 0; ii < mc_particle.NumberTrajectoryPoints(); ii++)
        {
            if(
                edep.StartT() <= trajectory.T(ii) &&
                edep.EndT() >= trajectory.T(ii)
            )
            for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
            {
                if(trajectory_processes[jj].first == ii) {
                    process = trajectory_processes[jj].second;
                }
            }
        }
        return TrajectoryStringToProcessType[process];
    }
    std::vector<Int_t> SimulationWrangler::DetermineDetectorSimulationEdeps(
        const std::vector<sim::IDE>& det_ide,
        Int_t detsim_id
    )
    {
        std::vector<Int_t> edep_ids;
        for(auto track : det_ide)
        {
            if(track.trackID > 0)
            {
                std::vector<Int_t> candidate_edeps = sTrackID_EdepIDMap[track.trackID];
                for(auto edep_id : candidate_edeps)
                {
                    auto edep = GetMCSimEnergyDeposit_EdepID(edep_id);
                    if(
                        edep.MidPointX() == track.x && 
                        edep.MidPointY() == track.y &&
                        edep.MidPointZ() == track.z
                    ) {
                        edep_ids.emplace_back(edep_id);
                        sEnergyDepositPointCloud.edep_detsim_id[edep_id].emplace_back(detsim_id);
                        sEdepID_DetSimIDMap[edep_id].emplace_back(detsim_id);
                    }
                }
            }
        }
        return edep_ids;
    }
    void SimulationWrangler::FillTTree()
    {
        if(sSaveSimulationWrangler)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "saving SimulationWrangler info to root file."
            );
            sSimulationWranglerTree->Fill();
        }
        if(sSaveEnergyDepositPointCloud)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "saving energy deposit point cloud data to root file."
            );
            sEnergyDepositPointCloudTree->Fill();
        }
        if(sSaveWirePlanePointCloud)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "saving wire plane point cloud data to root file."
            );
            sWirePlanePointCloudTree->Fill();
        }
        if(sSaveWirePlaneTrackTopology)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "saving wire plane track topology info to root file."
            );
            sWirePlaneTrackTopologyTree->Fill();
        }
        if(sSaveOpDetPointCloud)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "saving optical detector point cloud data to root file."
            );
            sOpDetPointCloudTree->Fill();
        }
    }     
    DetSimID_List SimulationWrangler::GetAllDetSimID_TrackID(TrackID_t track_id)
    {
        DetSimID_List detsim = sTrackID_DetSimIDMap[track_id];
        auto track_ids = sTrackID_DescendantTrackIDMap[track_id];
        for(auto particle : track_ids)
        {
            auto detsim_ids = sTrackID_DetSimIDMap[particle];
            detsim.insert(detsim.end(), detsim_ids.begin(), detsim_ids.end());
        }
        return detsim;
    }
    EdepID_List SimulationWrangler::GetAllEdepID_TrackID(TrackID_t track_id)
    {
        EdepID_List edep = sTrackID_EdepIDMap[track_id];
        auto track_ids = sTrackID_DescendantTrackIDMap[track_id];
        for(auto particle : track_ids)
        {
            auto edep_ids = sTrackID_EdepIDMap[particle];
            edep.insert(edep.end(), edep_ids.begin(), edep_ids.end());
        }
        return edep;
    }
}