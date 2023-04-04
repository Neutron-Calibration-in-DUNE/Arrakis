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
        {ProcessType::Decay,                "Decay"},
        {ProcessType::ComptonScatter,       "ComptonScatter"},
        {ProcessType::PhotoelectricEffect,  "PhotoelectricEffect"},
        {ProcessType::ElectronBremsstrahlung, "ElectronBremsstrahlung"},
        {ProcessType::ElectronIonization,   "ElectronIonization"},
        {ProcessType::PositronAnnihilation, "PositronAnnihilation"},
        {ProcessType::MuonIonization,       "MuonIonization"},
        {ProcessType::GammaConversion,      "GammaConversion"},
        {ProcessType::IonIonization,        "IonIonization"},
        {ProcessType::MuonCaptureAtRest,    "MuonCaptureAtRest"},
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
        {"Decay",               ProcessType::Decay},
        {"ComptonScatter",      ProcessType::ComptonScatter},
        {"PhotoelectricEffect", ProcessType::PhotoelectricEffect},
        {"ElectronBremsstrahlung", ProcessType::ElectronBremsstrahlung},
        {"ElectronIonization",  ProcessType::ElectronIonization},
        {"PositronAnnihilation",ProcessType::PositronAnnihilation},
        {"MuonIonization",      ProcessType::MuonIonization},
        {"GammaConversion",     ProcessType::GammaConversion},
        {"IonIonization",       ProcessType::IonIonization},
        {"MuonCaptureAtRest",   ProcessType::MuonCaptureAtRest},
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
        {"Decay",           ProcessType::Decay},
        {"compt",           ProcessType::ComptonScatter},
        {"phot",            ProcessType::PhotoelectricEffect},
        {"eBrem",           ProcessType::ElectronBremsstrahlung},
        {"eIoni",           ProcessType::ElectronIonization},
        {"annihil",         ProcessType::PositronAnnihilation},
        {"muIoni",          ProcessType::MuonIonization},
        {"conv",            ProcessType::GammaConversion},
        {"ionIoni",         ProcessType::IonIonization},
        {"muMinusCaptureAtRest",    ProcessType::MuonCaptureAtRest},
    };
    std::map<ProcessType, std::string> TrajectoryProcessTypeToString
    {
        {ProcessType::NotDefined,           "NotDefined"},
        {ProcessType::Unknown,              "Unknown"},
        {ProcessType::Primary,              "primary"},
        {ProcessType::HadronElastic,        "hadElastic"},
        {ProcessType::PiMinusInelastic,     "pi-Inelastic"},
        {ProcessType::PiPlusInelastic,      "pi+Inelastic"},
        {ProcessType::KaonMinusInelastic,   "kaon-Inelastic"},
        {ProcessType::KaonPlusInelastic,    "kaon+Inelastic"},
        {ProcessType::ProtonInelastic,      "protonInelastic"},
        {ProcessType::NeutronInelastic,     "neutronInelastic"},
        {ProcessType::CoulombScatter,       "CoulombScatter"},
        {ProcessType::NeutronCapture,       "nCapture"},
        {ProcessType::Transportation,       "Transportation"},
        {ProcessType::Decay,                "Decay"},
        {ProcessType::ComptonScatter,       "compt"},
        {ProcessType::PhotoelectricEffect,  "phot"},
        {ProcessType::ElectronBremsstrahlung, "eBrem"},
        {ProcessType::ElectronIonization,   "eIoni"},
        {ProcessType::PositronAnnihilation, "annihil"},
        {ProcessType::MuonIonization,       "muIoni"},
        {ProcessType::GammaConversion,      "conv"},
        {ProcessType::IonIonization,        "ionIoni"},
        {ProcessType::MuonCaptureAtRest,    "muMinusCaptureAtRest"},
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
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting up EnergyDepositPointCloud tree."
        );
        sEnergyDepositPointCloudTree = sTFileService->make<TTree>(
            "edep_point_cloud", "edep_point_cloud"
        );
        sEnergyDepositPointCloudTree->Branch("x",   &sEnergyDepositPointCloud.x);
        sEnergyDepositPointCloudTree->Branch("y",   &sEnergyDepositPointCloud.y);
        sEnergyDepositPointCloudTree->Branch("z",   &sEnergyDepositPointCloud.z);
        sEnergyDepositPointCloudTree->Branch("energy",   &sEnergyDepositPointCloud.energy);
        sEnergyDepositPointCloudTree->Branch("track_id",   &sEnergyDepositPointCloud.track_id);
        sEnergyDepositPointCloudTree->Branch("detsim_id",   &sEnergyDepositPointCloud.detsim_id);

        Logger::GetInstance("SimulationWrangler")->trace(
            "setting up SimulationWrangler tree."
        );
        sSimulationWranglerTree = sTFileService->make<TTree>("SimulationWrangler", "SimulationWrangler");
        sSimulationWranglerTree->Branch("generator_map", &sGeneratorMap);

        sSimulationWranglerTree->Branch("generator_label_map",  &sTrackID_GeneratorLabelMap);
        sSimulationWranglerTree->Branch("particle_id_map",      &sTrackID_ParticleIDMap);
        sSimulationWranglerTree->Branch("pdg_code_map",         &sTrackID_PDGCodeMap);
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

        sSimulationWranglerTree->Branch("edep_process_map",         &sEdepID_ProcessMap);
        sSimulationWranglerTree->Branch("edep_detsim_map",          &sEdepID_DetSimIDMap);
        sSimulationWranglerTree->Branch("detsim_edep_map",          &sDetSimID_EdepIDMap);

        Logger::GetInstance("SimulationWrangler")->trace(
            "setting up WirePlanePointCloud tree."
        );
        sWirePlanePointCloudTree = sTFileService->make<TTree>(
            "wire_plane_point_cloud", "wire_plane_point_cloud"
        );
        sWirePlanePointCloudTree->Branch("channel", &sWirePlanePointCloud.channel);
        sWirePlanePointCloudTree->Branch("wire",    &sWirePlanePointCloud.wire);
        sWirePlanePointCloudTree->Branch("tick",    &sWirePlanePointCloud.tick);
        sWirePlanePointCloudTree->Branch("tdc",     &sWirePlanePointCloud.tdc);
        sWirePlanePointCloudTree->Branch("adc",     &sWirePlanePointCloud.adc);
        sWirePlanePointCloudTree->Branch("view",    &sWirePlanePointCloud.view);
        sWirePlanePointCloudTree->Branch("energy",  &sWirePlanePointCloud.energy);
        sWirePlanePointCloudTree->Branch("source_label",    &sWirePlanePointCloud.source_label);
        sWirePlanePointCloudTree->Branch("shape_label",     &sWirePlanePointCloud.shape_label);
        sWirePlanePointCloudTree->Branch("particle_label",  &sWirePlanePointCloud.particle_label);
        sWirePlanePointCloudTree->Branch("unique_shape",    &sWirePlanePointCloud.unique_shape);
        sWirePlanePointCloudTree->Branch("unique_particle", &sWirePlanePointCloud.unique_particle);

        sGeneratorMap["Ar39Label"] =    GeneratorLabel::Ar39;
        sGeneratorMap["Ar42Label"] =    GeneratorLabel::Ar42;
        sGeneratorMap["Kr85Label"] =    GeneratorLabel::Kr85;
        sGeneratorMap["Rn222Label"] =   GeneratorLabel::Rn222;
        sGeneratorMap["CosmicsLabel"] = GeneratorLabel::Cosmics;
        sGeneratorMap["HEPevtLabel"] =  GeneratorLabel::HEPevt;
        sGeneratorMap["PNSLabel"] =     GeneratorLabel::PNS;
    }
    void SimulationWrangler::SetConfigurationParameters(const Parameters& config)
    {
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting up configuration parameters."
        );
        sSaveSimulationWrangler = config().SaveSimulationWrangler();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveSimulationWrangler: " + std::to_string(sSaveSimulationWrangler)
        );
        sSaveWirePlanePointCloud = config().SaveWirePlanePointCloud();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting SaveWirePlanePointCloud: " + std::to_string(sSaveWirePlanePointCloud)
        );
        sADCThreshold = config().ADCThreshold();
        Logger::GetInstance("SimulationWrangler")->trace(
            "setting ADCThreshold: " + std::to_string(sADCThreshold)
        );
    }
    void SimulationWrangler::SetWirePlanePointCloudLabels(
        DetSimID_t detSimID, TrackID_t trackID,
        SourceLabelInt sourceLabel, ShapeLabelInt shapeLabel, 
        ParticleLabelInt particleLabel, Int_t uniqueShape
    )
    {
        auto track_index = sWirePlanePointCloud.GetIndex_TrackID(detSimID, trackID);
        if(track_index != -1)
        {
            sWirePlanePointCloud.source_labels[detSimID][track_index] = sourceLabel;
            sWirePlanePointCloud.shape_labels[detSimID][track_index] = shapeLabel;
            sWirePlanePointCloud.particle_labels[detSimID][track_index] = particleLabel;
            sWirePlanePointCloud.unique_shapes[detSimID][track_index] = uniqueShape;
            sWirePlanePointCloud.unique_particles[detSimID][track_index] = trackID;
        }
        sWirePlanePointCloud.source_label[detSimID] = sourceLabel;
        sWirePlanePointCloud.shape_label[detSimID] = shapeLabel;
        sWirePlanePointCloud.particle_label[detSimID] = particleLabel;
        sWirePlanePointCloud.unique_shape[detSimID] = uniqueShape;
        sWirePlanePointCloud.unique_particle[detSimID] = trackID;
    }
    void SimulationWrangler::ResetEvent()
    {
        sMCTruthHandles.clear();
        sPrimaries.clear();
        sEnergyDepositPointCloud.clear();
        sWirePlanePointCloud.clear();

        sTrackID_GeneratorLabelMap.clear();
        sTrackID_ParticleIDMap.clear();
        sTrackID_PDGCodeMap.clear();
        sTrackID_ProcessMap.clear();
        sTrackID_EndProcessMap.clear();
        sTrackID_EnergyMap.clear();
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

        sDetSimID_EdepIDMap.clear();
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
    void SimulationWrangler::ProcessEvent(
        const Parameters& config, art::Event const& event
    )
    {
        ResetEvent();
        ProcessMCTruth(event, config().labels.get_PSet());
        ProcessMCParticles(event, config().LArGeantProducerLabel());
        ProcessSimEnergyDeposits(event, 
            config().IonAndScintProducerLabel(), config().IonAndScintInstanceLabel()
        );
        ProcessSimChannels(event, 
            config().SimChannelProducerLabel(), config().SimChannelInstanceLabel()
        );
        ProcessRawDigits(event,
            config().RawDigitProducerLabel(), config().RawDigitInstanceLabel()
        );
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
            sTrackID_PDGCodeMap[particle.TrackId()] = particle.PdgCode();
            sTrackID_ProcessMap[particle.TrackId()] = TrajectoryStringToProcessType[particle.Process()];
            sTrackID_EndProcessMap[particle.TrackId()] = TrajectoryStringToProcessType[particle.EndProcess()];
            sTrackID_EnergyMap[particle.TrackId()] = particle.E();

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
            sTrackID_EdepIDMap[edep.TrackID()].emplace_back(edep_index);
            ProcessType process = DetermineEdepProcess(edep);
            sTrackID_EdepProcessMap[edep.TrackID()].emplace_back(process);
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
                            (Int_t) (std::abs(uncompressed[l])),
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
                        (Int_t) (std::abs(uncompressed[l])),
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
                    digit_index += 1;
                }
            }
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
    DetSimID_Collection SimulationWrangler::GetDetSimID_TrackID(TrackID_List trackIDList)
    {
        DetSimID_Collection detsim;
        for(auto track_id : trackIDList)
        {
            detsim.emplace_back(sTrackID_DetSimIDMap[track_id]);
        }
        return detsim;
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
                        //std::cout << "\t\tMATCH" << std::endl;
                        edep_ids.emplace_back(edep_id);
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
        if(sSaveWirePlanePointCloud)
        {
            Logger::GetInstance("SimulationWrangler")->trace(
                "saving wire plane point cloud data to root file."
            );
            sWirePlanePointCloudTree->Fill();
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
}