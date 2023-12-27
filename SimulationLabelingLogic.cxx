/**
 * @file SimulationLabelingLogic.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief
 * @version 0.1
 * @date 2023-02-22
 */
#include "SimulationLabelingLogic.h"

namespace arrakis
{
    SimulationLabelingLogic *SimulationLabelingLogic::sInstance{nullptr};
    std::mutex SimulationLabelingLogic::sMutex;

    SimulationLabelingLogic *SimulationLabelingLogic::GetInstance()
    {
        std::lock_guard<std::mutex> lock(sMutex);
        if (sInstance == nullptr)
        {
            sInstance = new SimulationLabelingLogic();
        }
        return sInstance;
    }
    SimulationLabelingLogic::SimulationLabelingLogic()
    {
    }

    void SimulationLabelingLogic::SetConfigurationParameters(const Parameters &config)
    {
        /**
         * Configuration parameters are enumerated in the Arrakis.fcl file.  
         * Each of the associated SimulationLabelingLogic parameters are assigned here.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "setting up configuration parameters.");

        sInducedChannelInfluence = config().InducedChannelInfluence();
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "setting channel influence radius to: " + std::to_string(sInducedChannelInfluence) + " channels."
        );

        sInducedTDCInfluence = config().InducedTDCInfluence();
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "setting tdc influence radius to: " + std::to_string(sInducedTDCInfluence) + " tdcs."
        );

        sShowerEnergyThreshold = config().ShowerEnergyThreshold();
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "setting shower energy threshold to: " + std::to_string(sShowerEnergyThreshold) + " MeV."
        );
    }

    void SimulationLabelingLogic::ResetEvent()
    {
        mUniqueTopologyLabel = 0;
        mUniquePhysicsMicroLabel = 0;
        mUniquePhysicsMesoLabel = 0;
        mUniquePhysicsMacroLabel = 0;
    }
    Int_t SimulationLabelingLogic::UniqueTopology()
    {
        mUniqueTopologyLabel += 1;
        return mUniqueTopologyLabel;
    }
    Int_t SimulationLabelingLogic::UniquePhysicsMicro()
    {
        mUniquePhysicsMicroLabel += 1;
        return mUniquePhysicsMicroLabel;
    }
    Int_t SimulationLabelingLogic::UniquePhysicsMeso()
    {
        mUniquePhysicsMesoLabel += 1;
        return mUniquePhysicsMesoLabel;
    }
    Int_t SimulationLabelingLogic::UniquePhysicsMacro()
    {
        mUniquePhysicsMacroLabel += 1;
        return mUniquePhysicsMacroLabel;
    }
    void SimulationLabelingLogic::FillTTree()
    {
    }

    void SimulationLabelingLogic::ProcessEvent(
        const Parameters &config, art::Event const &event)
    {
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing event."
        );
        ResetEvent();
        /**
         * Process Energy Deposit Point Clouds
        */ 
        /**
         * Process Wire Plane Point Clouds
        */ 
        PrepareInitialPointClouds(config, event);
        // Process HIPS/MIPS
        ProcessElectrons(config, event);
        ProcessPositrons(config, event);
        ProcessGammas(config, event);
        ProcessMuons(config, event);
        ProcessPion0s(config, event);
        ProcessPionss(config, event);
        ProcessKaon0s(config, event);
        ProcessKaons(config, event);
        ProcessProtons(config, event);
        // Neutron captures, NR/ER
        ProcessNeutrons(config, event);
        ProcessNuclearRecoils(config, event);
        ProcessElectronRecoils(config, event);
        // Process radiologicals
        ProcessAr39(config, event);
        ProcessAr42(config, event);
        ProcessKr85(config, event);
        ProcessRn222(config, event);
        // Process cosmics
        ProcessCosmics(config, event);
        ProcessNoise(config, event);
        CleanUpPointClouds(config, event);
        /**
         * Process Wire Plane Track Topology
        */ 
        /**
         * Process Op Det Point Clouds
        */ 
        FillTTree();
    }
    PhysicsMacroLabel SimulationLabelingLogic::DeterminePhysicsMacroLabel(TrackID_t trackID)
    {
        /**
         * Find the ancestor track id for this particle,
         * then determine how to label the source based on the 
         * generator from the ancestor.
         */
        auto mc_data = SimulationWrangler::GetInstance();
        GeneratorLabel generator_label;
        if(mc_data->GetParentTrackID_TrackID(trackID) == 0) {
            generator_label = mc_data->GetGeneratorLabel_TrackID(trackID);
        }
        else {
            generator_label = mc_data->GetGeneratorLabel_TrackID(
                mc_data->GetAncestorTrackID_TrackID(trackID)
            );
        }
        if(
            generator_label == GeneratorLabel::Ar39 ||
            generator_label == GeneratorLabel::Ar42 ||
            generator_label == GeneratorLabel::Kr85 ||
            generator_label == GeneratorLabel::Rn222
        ) {
            return PhysicsMacroLabel::Radiological;
        }
        else if(generator_label == GeneratorLabel::Cosmics) {
            return PhysicsMacroLabel::Cosmics;
        }
        else if(generator_label == GeneratorLabel::PNS) {
            return PhysicsMacroLabel::PulsedNeutronSource;
        }
        else if(generator_label == GeneratorLabel::HEPevt) {
            return PhysicsMacroLabel::HEPevt;
        }
        else {
            return PhysicsMacroLabel::Undefined;
        }
    }
    void SimulationLabelingLogic::SetLabels(
        DetSimID_List detsimID, 
        EdepID_List edepID, 
        TrackID_t trackID,
        TopologyLabel topology, 
        PhysicsMicroLabel physicsMicroLabel,
        PhysicsMesoLabel physicsMesoLabel,
        PhysicsMacroLabel physicsMacroLabel,
        Int_t uniqueTopologyLabel,
        Int_t uniquePhysicsMicroLabel,
        Int_t uniquePhysicsMesoLabel,
        Int_t uniquePhysicsMacroLabel
    )
    {
        auto mc_data = SimulationWrangler::GetInstance();
        auto source_label = DeterminePhysicsMacroLabel(trackID);
        auto particle = mc_data->GetPDGCode_TrackID(trackID);
        for(size_t ii = 0; ii < detSimIDList.size(); ii++)
        {
            mc_data->SetWirePlanePointCloudLabels(
                detSimIDList[ii], 
                trackID, 
                LabelCast(topology), 
                particle, 
                LabelCast(physicsMicroLabel),
                LabelCast(physicsMesoLabel),
                LabelCast(physicsMacroLabel),
                uniqueTopologyLabel,
                uniquePhysicsMicroLabel,
                uniquePhysicsMesoLabel,
                uniquePhysicsMacroLabel
            );
        }
        for(size_t ii = 0; ii < edepIDList.size(); ii++)
        {
            // mc_data->SetEnergyDepositPointCloudLabels(
            //     edepIDList[ii], 
            //     LabelCast(source_label), 
            //     LabelCast(topology), 
            //     particle, 
            //     LabelCast(physics),
            //     unique_topology
            // );
        }
    }
    void SimulationLabelingLogic::SetLabels(
        DetSimID_Collection detSimIDCollection,
        EdepID_Collection edepIDCollection,
        TrackID_List trackIDList,
        TopologyLabel topology, 
        PhysicsMicroLabel physicsMicroLabel,
        PhysicsMesoLabel physicsMesoLabel,
        PhysicsMacroLabel physicsMacroLabel,
        Int_t uniqueTopologyLabel,
        Int_t uniquePhysicsMicroLabel,
        Int_t uniquePhysicsMesoLabel,
        Int_t uniquePhysicsMacroLabel
    )
    {
        for(size_t ii = 0; ii < detSimIDCollection.size(); ii++)
        {
            SetLabels(
                detSimIDCollection[ii],
                edepIDCollection[ii],
                trackIDList[ii], 
                topology, 
                physicsMicroLabel,
                physicsMesoLabel,
                physicsMacroLabel,
                uniqueTopologyLabel,
                uniquePhysicsMicroLabel,
                uniquePhysicsMesoLabel,
                uniquePhysicsMacroLabel
            );
        }
    }
    void SimulationLabelingLogic::SetLabels(
        std::vector<DetSimID_Collection> detSimIDCollection,
        std::vector<EdepID_Collection> edepIDCollection,
        TrackID_Collection trackIDCollection,
        TopologyLabel topology, 
        PhysicsMicroLabel physicsMicroLabel,
        PhysicsMesoLabel physicsMesoLabel,
        PhysicsMacroLabel physicsMacroLabel,
        Int_t uniqueTopologyLabel,
        Int_t uniquePhysicsMicroLabel,
        Int_t uniquePhysicsMesoLabel,
        Int_t uniquePhysicsMacroLabel
    )
    {
        for(size_t ii = 0; ii < detSimIDCollection.size(); ii++)
        {
            SetLabels(
                detSimIDCollection[ii],
                edepIDCollection[ii],
                trackIDCollection[ii], 
                topology, 
                physicsMicroLabel,
                physicsMesoLabel,
                physicsMacroLabel,
                uniqueTopologyLabel,
                uniquePhysicsMicroLabel,
                uniquePhysicsMesoLabel,
                uniquePhysicsMacroLabel
            );
        }
    }

    void SimulationLabelingLogic::PrepareInitialPointClouds(
        const Parameters &config, art::Event const &event)
    {
        auto mc_data = SimulationWrangler::GetInstance();
        auto wire_plane_point_cloud = mc_data->GetWirePlanePointCloud();
        for(size_t detsim_id = 0; detsim_id < wire_plane_point_cloud.channel.size(); detsim_id++)
        {

        }
    }

    void SimulationLabelingLogic::CleanUpPointClouds(
        const Parameters &config, art::Event const &event)
    {
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "cleaning up point clouds."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto wire_plane_point_cloud = mc_data->GetWirePlanePointCloud();
        for(size_t detsim_id = 0; detsim_id < wire_plane_point_cloud.channel.size(); detsim_id++)
        {
            /**
             * Check to see if a point has an undefined particle label, or
             * if there are multiple labels.  If there are multiple labels,
             * then there is some logic to determine which should be
             * considered the "true" label.
             */
            if(wire_plane_point_cloud.particle_label[detsim_id] == LabelCast(ParticleLabel::Undefined))
            {
                for(auto track_id : wire_plane_point_cloud.track_ids[detsim_id])
                {
                    Logger::GetInstance("SimulationLabelingLogic")->warning(
                        "undefined point: " + std::to_string(detsim_id) + " - [trackid,pdg,process,level]: [" +
                        std::to_string(track_id) + "," + std::to_string(mc_data->GetPDGCode_TrackID(track_id)) + "," +
                        std::to_string(ProcessTypeInt(mc_data->GetProcess_TrackID(track_id))) + "," + 
                        std::to_string(mc_data->GetAncestorLevel_TrackID(track_id)) + "]"
                    );
                    auto ancestry = mc_data->GetAncestryTrackID_TrackID(track_id);
                    if(ancestry.size() > 0)
                    {
                        Logger::GetInstance("SimulationLabelingLogic")->warning(
                            "ancestry [trackid,pdg,process,level,label_function]:"
                        );
                    }
                    for(auto id : ancestry)
                    { 
                        Logger::GetInstance("SimulationLabelingLogic")->warning(
                            "\t[" + std::to_string(id) + "," + std::to_string(mc_data->GetPDGCode_TrackID(id)) +
                            "," + std::to_string(ProcessTypeInt(mc_data->GetProcess_TrackID(id))) + "," + 
                            std::to_string(mc_data->GetAncestorLevel_TrackID(id)) + "," + 
                            std::to_string(mc_data->GetLabelingFunction_TrackID(id)) + "]" 
                        );
                    }
                    mc_data->PrintParticleData(track_id);
                }
            }
            /**
             * Here we check to see if there are multiple labels for a given point,
             * and if so, then we use a rule to decide what the label should be.
             * The rule is, if any of the following points exist, then they
             * are taken with priority:
             * 
             *      track > shower > blip > neutron_capture
            */
            if(wire_plane_point_cloud.topology_labels[detsim_id].size() > 1)
            {
                auto num_tracks = std::count(
                    wire_plane_point_cloud.topology_labels[detsim_id].begin(), 
                    wire_plane_point_cloud.topology_labels[detsim_id].end(), 
                    LabelCast(TopologyLabel::Track)
                );
                auto num_showers = std::count(
                    wire_plane_point_cloud.topology_labels[detsim_id].begin(), 
                    wire_plane_point_cloud.topology_labels[detsim_id].end(), 
                    LabelCast(TopologyLabel::Shower)
                );
                auto num_blips = std::count(
                    wire_plane_point_cloud.topology_labels[detsim_id].begin(), 
                    wire_plane_point_cloud.topology_labels[detsim_id].end(), 
                    LabelCast(TopologyLabel::Blip)
                );
                if(num_tracks) {
                    wire_plane_point_cloud.topology_label[detsim_id] = LabelCast(TopologyLabel::Track);
                }
                else if(num_showers) {
                    wire_plane_point_cloud.topology_label[detsim_id] = LabelCast(TopologyLabel::Shower);
                }
                else if(num_blips) {
                    wire_plane_point_cloud.topology_label[detsim_id] = LabelCast(TopologyLabel::Blip);
                }
            }
        }
    }
    void SimulationLabelingLogic::ProcessNoise(
        const Parameters &config, art::Event const &event)
    {
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing noise."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto wire_plane_point_cloud = mc_data->GetWirePlanePointCloud();
        for(size_t detsim_id = 0; detsim_id < wire_plane_point_cloud.channel.size(); detsim_id++)
        {
            /**
             * If the point corresponds to noise, then we need to do some digging to see
             * if the point actually corresponds to some induced signal from the wirecell
             * simulation.  This can happen through various effects which are not
             * accounted for by SimChannelSink, the module which saves the drifted 
             * electrons to SimChannels.
             * 
             * The energy depositions are first drifted to 10cm from the wire planes, 
             * which is chosen as a cutoff distance for induction effects.  From this 10cm
             * distance, SimChannelSink assigns track ids to SimChannels.  
             * 
             * Other effects from the detector simulation include (a) the field response,
             * (b) detector response, where induction effects can reach out to a 10 wire distance
             * and approximately 200 tdc values (~2 microseconds).  This induction influence
             * is not backtracked to the SimChannel, which results in bands around tracks and
             * showers which do not have associated track ids.  To deal with this,
             * we search within the 10 wire and 200 tdc distance for noise labeled points
             * and weight the influence of nearby track ids by their distance and ADC value.
             */
            if(
                wire_plane_point_cloud.particle_label[detsim_id] == LabelCast(ParticleLabel::Noise)
            )
            {
                Int_t current_channel = wire_plane_point_cloud.channel[detsim_id];
                Int_t current_tdc = wire_plane_point_cloud.tdc[detsim_id];
                Int_t current_view = wire_plane_point_cloud.view[detsim_id];

                DetSimID_t largest_influence = -1;
                Double_t influence = 0.0;

                /**
                 * The logic here is meant to speed up the search time for the region of influence.
                 * We want to look out 10 wires in each direction and 200 tdc ticks in each 
                 * direction.  Since the RawDigit information is saved sequentially as [channel,tdc]:
                 * [0,0],[0,1],[0,2],...,[0,N],[1,0],[1,1],[1,2],...,[1,N],[2,0],...,[M,N],
                 * we can index to the left and right by flattening the index.
                */
                size_t start = 0;
                size_t end = wire_plane_point_cloud.channel.size();
                size_t index_distance = sInducedChannelInfluence * mc_data->GetNumberOfTDCs() + sInducedTDCInfluence;
                if(detsim_id > index_distance) {
                    start = detsim_id - index_distance;
                }
                if(detsim_id + index_distance < end) {
                    end = detsim_id + index_distance;
                }
                for(size_t other_id = start; other_id < end; other_id++)
                {
                    if(
                        (std::abs(wire_plane_point_cloud.channel[other_id] - current_channel) < sInducedChannelInfluence ||
                         std::abs(wire_plane_point_cloud.tdc[other_id] - current_tdc) < sInducedTDCInfluence) &&
                         wire_plane_point_cloud.view[other_id] == current_view &&
                         wire_plane_point_cloud.particle_label[other_id] != LabelCast(ParticleLabel::Noise)
                    )
                    {
                        Double_t temp_distance = Double_t(sqrt(
                            pow(Double_t(wire_plane_point_cloud.channel[other_id] - current_channel), 2.0) + 
                            pow(Double_t(wire_plane_point_cloud.tdc[other_id] - current_tdc), 2.0)
                        ));
                        Double_t temp_influence = Double_t(std::abs(wire_plane_point_cloud.adc[other_id])) / temp_distance;
                        if(temp_influence > influence)
                        {
                            influence = temp_influence;
                            largest_influence = other_id;
                        }
                    }
                }
                /**
                 * Now pass the region of influence and the current DetSimID to 
                 * some logic to determine how to label this point.
                 */
                if(influence > 0) {
                    mc_data->SetWirePlanePointCloudLabels(
                        detsim_id, 
                        largest_influence,
                        wire_plane_point_cloud.topology_label[largest_influence],
                        wire_plane_point_cloud.particle_label[largest_influence],
                        wire_plane_point_cloud.physics_label[largest_influence],
                        wire_plane_point_cloud.physics_meso_label[largest_influence],
                        wire_plane_point_cloud.physics_macro_label[largest_influence],
                        wire_plane_point_cloud.unique_topology[largest_influence],
                        wire_plane_point_cloud.unique_physics_micro[largest_influence],
                        wire_plane_point_cloud.unique_physics_meso[largest_influence],
                        wire_plane_point_cloud.unique_physics_macro[largest_influence],
                        1   // induction flag
                    );
                }
                else {
                    Logger::GetInstance("SimulationLabelingLogic")->warning(
                        "noise point: " + std::to_string(detsim_id) + 
                        " channel: " + std::to_string(current_channel) +
                        " tdc: " + std::to_string(current_tdc) +
                        " view: " + std::to_string(current_view) + 
                        " has no nearby track ids to influence labeling."
                    );
                }
            }
        }
    }
    void SimulationLabelingLogic::ProcessShowers(
        TrackID_t trackID, 
        Int_t uniqueTopology
    )
    {
        /**
         * Showers can come in two flavors, electromagnetic and hadronic.  Only
            dealing with the electromagnetic version for now.  The methodology
            is to assume that each trackid is coming from an electron, positron
            or gamma, but we should check this just in case.  The full shower,
            if there is one, will be associated back to whichever particle
            created the shower.

            Showers currently consist of
                (1) topology label          - TopologyLabel.Shower = 3,
                (3) physics_micro labels    - PhysicsMicroLabel::ElectronIonization = 3
                                            PhysicsMicroLabel::GammaCompton = 4
                                            PhysicsMicroLabel::GammaConversion = 5
                (3) physics_meso labels     - PhysicsMesoLabel.ElectronShower = 5
                                            PhysicsMesoLabel.PositronShower = 6
                                            PhysicsMesoLabel.PhotonShower = 7
                (4) physics_macro labels    - PhysicsMacroLabel.CCNue = 1
                                            PhysicsMacroLabel.CCNuMu = 2
                                            PhysicsMacroLabel.NC = 3
                                            PhysicsMacroLabel.Cosmics = 5

            The physics_micro labels are determined by performing a series of steps.
            First, we look to see if there are any gamma conversion events within
            the hierarchy.  This is denoted as subprocess 14 within edep-sim.
            Each set of descendants from these gammas will get the physics_micro
            label PhysicsMicroLabel::GammaConversion.  We then look for any compton
            scatters or photo-electric effect.  Comptons which have the same parent
            are linked by unique_physics_micro numbers, while photo-electric effect
            have independent ones.  The photo-electric effect gets the labels
            PhysicsMicroLabel::ElectronIonization.

            The physics_meso labels are determined simply by the type of originating
            particle for the shower.  If there is no shower, then the physics_meso
            label ...

            The physics_macro labels are determined by a separate function.
        */
        auto mc_data = SimulationWrangler::GetInstance();
        auto particle_process = mc_data->GetProcess_TrackID(trackID);
        auto particle_descendants = mc_data->GetDescendantTrackID_TrackID(trackID);

        bool shower = false;
        TopologyLabel topology = TopologyLabel::Blip;
        PhysicsMesoLabel physics_meso = PhysicsMesoLabel::Undefined;
        auto unique_physics_meso = UniquePhysicsMeso();

        // grab descendants by type
        auto ionization_descendants = mc_data->FilterTrackID_Process(
            particle_descendants, ProcessType::ElectronIonization
        );
        auto bremsstrahlung_descendants = mc_data->FilterTrackID_Process(
            particle_descendants, ProcessType::ElectronBremmstrahlung
        );
        auto annihilation_descendants = mc_data->FilterTrackID_Process(
            particle_descendants, ProcessType::ElectronPositronAnnihilation
        );
        auto photoelectric_descendants = mc_data->FilterTrackID_Process(
            particle_descendants, ProcessType::PhotoelectricEffect
        );
        auto compton_descendants = mc_data->FilterTrackID_Process(
            particle_descendants, ProcessType::ComptonScatter
        );
        auto conversion_descendants = mc_data->FilterTrackID_Process(
            particle_descendants, ProcessType::GammaConversion
        );
        
        // count number of each subprocess, including the particle itself
        if (particle_process == ProcessType::ElectronIonization) {
            ionization_descendants.emplace_back(trackID);
        }
        if (particle_process == ProcessType::ElectronBremsstrahlung) {
            bremsstrahlung_descendants.emplace_back(trackID);
        }
        if (particle_process == ProcessType::ElectronPositronAnnihilation) {
            annihilation_descendants.emplace_back(trackID);
        }
        if (particle_process == ProcessType::PhotoelectricEffect) {
            photoelectric_descendants.emplace_back(trackID);
        }
        if (particle_process == ProcessType::ComptonScatter) {
            compton_descendants.emplace_back(trackID);
        }
        if (particle_process == ProcessType::GammaConversion) {
            conversion_descendants.emplace_back(trackID);
        }

        // determine what type of shower this is
        auto num_conversion = conversion_descendatns.size();
        if (num_conversion > 0) 
        {
            shower = true;
            topology = TopologyLabel::Shower;
            auto earliest_conversion = 10e19;
            auto earliest_bremsstrahlung = 10e19;
            auto conversion_times = mc_data->GetStartTime_TrackID(conversion_descendants);
            earliest_conversion = min(conversion_times);
            if (bremsstrahlung_descendants.size() > 0) 
            {
                auto bremsstrahlung_times = mc_data->GetStartTime_TrackID(bremsstrahlung_descendants);
                earliest_bremsstrahlung = min(bremsstrahlung_times);
            }
            else
            {
                earliest_bremsstrahlung = mc_data->GetStartTime_TrackID(trackID);
            }
            if (earliest_conversion <= earliest_bremsstrahlung) {
                physics_meso = PhysicsMesoLabel::PhotonShower;
            }
            else {
                physics_meso = PhysicsMesoLabel::ElectronShower;
            }
        }
        else
        {
            if (mc_data->GetAbsPDGCode_TrackID(trackID) == 111) {
                physics_meso = PhysicsMesoLabel::Pi0Decay;
            }
            else {
                physics_meso = PhysicsMesoLabel::LowEnergyIonization;
            }
        }

        auto bremsstrahlung_detsim = mc_data->GetDetSimID_TrackID(bremsstrahlung_descendants);
        auto bremsstrahlung_edep = mc_data->GetEdepID_TrackID(bremsstrahlung_descendants);
        SetLabels(
            bremsstrahlung_detsim,
            bremsstrahlung_edep,
            bremsstrahlung_descendants,
            topology,
            PhysicsMicroLabel::Bremmstrahlung,
            physics_meso,
            PhysicsMacroLabel::Undefined,
            uniqueTopology,
            UniquePhysicsMicro(),
            unique_physics_meso,
            0,
        );

        auto photoelectric_detsim = mc_data->GetDetSimID_TrackID(photoelectric_descendants);
        auto photoelectric_edep = mc_data->GetEdepID_TrackID(photoelectric_descendants);
        SetLabels(
            photoelectric_detsim,
            photoelectric_edep,
            photoelectric_descendants,
            TopologyLabel.Blip,
            PhysicsMicroLabel::PhotoElectric,
            physics_meso,
            PhysicsMacroLabel::Undefined,
            uniqueTopology,
            UniquePhysicsMicro(),
            unique_physics_meso,
            0,
        );

        auto compton_detsim = mc_data->GetDetSimID_TrackID(compton_descendants);
        auto compton_edep = mc_data->GetEdepID_TrackID(compton_descendants);
        SetLabels(
            compton_detsim,
            compton_edep,
            compton_descendants,
            TopologyLabel.Blip,
            PhysicsMicroLabel::GammaCompton,
            physics_meso,
            PhysicsMacroLabel::Undefined,
            uniqueTopology,
            UniquePhysicsMicro(),
            unique_physics_meso,
            0,
        );

        auto annihilation_detsim = mc_data->GetDetSimID_TrackID(annihilation_descendants);
        auto annihilation_edep = mc_data->GetEdepID_TrackID(annihilation_descendants);
        SetLabels(
            annihilation_detsim,
            annihilation_edep,
            annihilation_descendants,
            topology,
            PhysicsMicroLabel::Annihilation,
            physics_meso,
            PhysicsMacroLabel::Undefined,
            uniqueTopology,
            UniquePhysicsMicro(),
            unique_physics_meso,
            0,
        );

        auto conversion_detsim = mc_data->GetDetSimID_TrackID(conversion_descendants);
        auto conversion_edep = mc_data->GetEdepID_TrackID(conversion_descendants);
        SetLabels(
            conversion_detsim,
            conversion_edep,
            conversion_descendants,
            TopologyLabel.Shower,
            PhysicsMicroLabel::GammaConversion,
            physics_meso,
            PhysicsMacroLabel::Undefined,
            uniqueTopology,
            UniquePhysicsMicro(),
            unique_physics_meso,
            0,
        );

        auto gamma_descendants = mc_data->FilterTrackID_AbsPDGCode(particle_descendants, 22);
        auto hadron_elastic = mc_data->FilterTrackID_Process(gamma_descendants, ProcessType::HadronElastic);
        auto hadron_inelastic = mc_data->FilterTrackID_Process(gamma_descendants, ProcessType::HadronInelastic);
        if (
            mc_data->GetAbsPDGCode_TrackID(trackID) == 22 &&
            mc_data->GetProcess_TrackID(trackID) == 111
        ):
            hadron_elastic.emplace_back(trackID);
        if (
            mc_data->GetAbsPDGCode_TrackID(trackID) == 22 &&
            mc_data->GetProcess_TrackID(trackID) == 121
        ):
            hadron_inelastic.emplace_back(trackID);

        auto hadron_elastic_hits = self.simulation_wrangler.get_hits_trackid(hadron_elastic);
        auto hadron_elastic_segments = self.simulation_wrangler.get_segments_trackid(hadron_elastic);
        self.simulation_wrangler.set_hit_labels_list(
            hadron_elastic_hits,
            hadron_elastic_segments,
            hadron_elastic,
            TopologyLabel.Blip,
            PhysicsMicroLabel::HadronElastic,
            PhysicsMesoLabel.NuclearRecoil,
            PhysicsMacroLabel::Undefined,
            next(self.unique_physics_micro),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
            0,
        );

        auto hadron_inelastic_hits = self.simulation_wrangler.get_hits_trackid(hadron_inelastic);
        auto hadron_inelastic_segments = self.simulation_wrangler.get_segments_trackid(hadron_inelastic);
        self.simulation_wrangler.set_hit_labels_list(
            hadron_inelastic_hits,
            hadron_inelastic_segments,
            hadron_inelastic,
            TopologyLabel.Blip,
            PhysicsMicroLabel::HadronInelastic,
            PhysicsMesoLabel.NuclearRecoil,
            PhysicsMacroLabel::Undefined,
            next(self.unique_physics_micro),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
            0,
        );
    }
    void SimulationLabelingLogic::ProcessShowers(
        TrackID_List trackIDList, 
        Int_t TopologyLabel
    )
    {
        for (auto track_id : trackIDList)
        {
            ProcessShowers(track_id, TopologyLabel);
        }
    }
    void SimulationLabelingLogic::ProcessShowers(TrackID_Collection trackIDCollection)
    {
        for (auto track_id : trackIDCollection)
        {
            ProcessShowers(track_id, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessElectrons(
        const Parameters &config, art::Event const &event)
    {
        /**
         * This function is responsible for labeling primary electrons, which
         * are for the most part considered showers, unless their total
         * kinetic energy is low.  
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing primary electrons."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto electrons = mc_data->GetPrimaries_PDGCode(11);
        for (auto electron : electrons) {
            ProcessShowers(electron, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessPositrons(
        const Parameters &config, art::Event const &event)
    {
        /**
         * This function is responsible for labeling primary positrons, which
         * are for the most part considered showers, unless their total
         * kinetic energy is low. 
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing primary positrons."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto positrons = mc_data->GetPrimaries_PDGCode(-11);
        for (auto positron : positrons) {
            ProcessShowers(positron, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessGammas(
        const Parameters &config, art::Event const &event)
    {
        /**
         * This function is responsible for labeling primary gammas, which
         * are for the most part considered showers, unless their total
         * kinetic energy is low. 
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing primary gammas."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto gammas = mc_data->GetPrimaries_AnsPDGCode(22);
        for (auto gamma : gammas) {
            ProcessShowers(gamma, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessMuons(
        const Parameters &config, art::Event const &event)
    {
        /**
         * For muons, and likewise anti-muons, there are a few different
         * special physics processes which are relevant.  For one, all muon
         * ionizations are track like objects which are MIPIonizations. Some
         * electrons generated from "muIoni", or "muPairProd" are deltas, 
         * while electrons generated from "mu+-CaptureAtRest" or "Decay"
         * are Michels.  
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing muons."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto muons = mc_data->GetTrackID_AbsPDGCode(13);
        for (auto muon : muons)
        {
            Int_t muon_topology = UniqueTopology();

            auto muon_det_sim = mc_data->GetDetSimID_TrackID(muon);
            auto muon_edep = mc_data->GetEdepID_TrackID(muon);
            SetLabels(
                muon_det_sim,
                muon_edep,
                muon,
                TopologyLabel::Track,
                PhysicsMicroLabel::MIPIonization,
                PhysicsMesoLabel::MIP,
                PhysicsMacroLabel::Undefined,
                muon_topology,
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );

            auto muon_daughters = mc_data->GetDaughterTrackID_TrackID(muon);
            auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(muon_daughters, 11);
            auto decay_daughters = mc_data->FilterTrackID_Process(elec_daughters, ProcessType::Decay);
            auto capture_daughters = mc_data->FilterTrackID_Process(elec_daughters, ProcessType::MuonCaptureAtRest);
            auto michel_decay_det_sim = mc_data->GetDetSimID_TrackID(decay_daughters);
            auto michel_capture_det_sim = mc_data->GetDetSimID_TrackID(capture_daughters);
            auto michel_decay_edep = mc_data->GetEdepID_TrackID(decay_daughters);
            auto michel_capture_edep = mc_data->GetEdepID_TrackID(capture_daughters);
            SetLabels(
                michel_decay_det_sim,
                michel_decay_edep,
                decay_daughters,
                TopologyLabel::Track,
                PhysicsMicroLabel::ElectronIonization,
                PhysicsMesoLabel::MichelElectron,
                PhysicsMacroLabel::Undefined,
                muon_topology,
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            SetLabels(
                michel_capture_det_sim,
                michel_capture_edep,
                capture_daughters,
                TopologyLabel::Track,
                PhysicsMicroLabel::ElectronIonization,
                PhysicsMesoLabel::MichelElectron,
                PhysicsMacroLabel::Undefined,
                muon_topology,
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            auto michel_decay_descendants = mc_data->GetDescendantTrackID_TrackID(decay_daughters);
            auto michel_capture_descendants = mc_data->GetDescendantTrackID_TrackID(capture_daughters);
            ProcessShowers(decay_daughters);
            ProcessShowers(capture_daughters);
            
            auto delta_daughters = mc_data->FilterTrackID_Process(elec_daughters, ProcessType::MuonIonization);
            auto delta_det_sim = mc_data->GetDetSimID_TrackID(delta_daughters);
            auto delta_edep = mc_data->GetEdepID_TrackID(delta_daughters);
            SetLabels(
                delta_det_sim,
                delta_edep,
                delta_daughters,
                TopologyLabel::Track,
                PhysicsMicroLabel::ElectronIonization,
                PhysicsMesoLabel::DeltaElectron,
                PhysicsMacroLabel::Undefined,
                muon_topology,
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            auto delta_descendants = mc_data->GetDescendantTrackID_TrackID(delta_daughters);
            ProcessShowers(delta_descendants);

            auto not_decay_elec_daughters = mc_data->FilterTrackID_NotProcess(elec_daughters, ProcessType::Decay);
            auto not_muon_capture_elec_daughters = mc_data->FilterTrackID_NotProcess(
                not_decay_elec_daughters, ProcessType::MuonCaptureAtRest
            );
            auto other_elec_daughters = mc_data->FilterTrackID_NotProcess(
                not_muon_capture_elec_daughters, ProcessType::MuonIonization
            );
            ProcessShowers(other_elec_daughters, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessPion0s(
        const Parameters &config, art::Event const &event)
    {
        /** (From Wikipedia: https://en.wikipedia.org/wiki/Pion)
         * The dominant π0 decay mode, with a branching ratio of BR2γ = 0.98823,
         * is into two photons:
         *      π0 → 2γ.
         * The decay π0 → 3γ (as well as decays into any odd number of photons) is
         * forbidden by the C-symmetry of the electromagnetic interaction: The intrinsic
         * C-parity of the π0 is +1, while the C-parity of a system of n photons is (−1)n.
         * The second largest π0 decay mode ( BRγee = 0.01174 ) is the Dalitz decay
         * (named after Richard Dalitz), which is a two-photon decay with an internal
         * photon conversion resulting a photon and an electron-positron pair in the final state:
         *      π0 → γ + e− + e+.
         * The third largest established decay mode ( BR2e2e = 3.34×10−5 ) is the double-Dalitz
         * decay, with both photons undergoing internal conversion which leads to further
         * suppression of the rate:
         *      π0 → e− + e+ + e− +	e+.
         * The fourth largest established decay mode is the loop-induced and therefore suppressed
         * (and additionally helicity-suppressed) leptonic decay mode ( BRee = 6.46×10−8 ):
         *      π0 → e− + e+.
         * The neutral pion has also been observed to decay into positronium with a branching
         * fraction on the order of 10−9. No other decay modes have been established experimentally.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing neutral pions."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto pi0s = mc_data->GetPrimaries_AnsPDGCode(111);
        for (auto pi0 : pi0s) {
            ProcessShowers(pi0, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessPions(
        const Parameters &config, art::Event const &event)
    {
        /**
         * 
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing pion+s."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto pions = mc_data->GetTrackID_AbsPDGCode(211);
        for (auto pion : pions)
        {
            auto pion_daughters = mc_data->GetDaughterTrackID_TrackID(pion);
            auto pion_det_sim = mc_data->GetDetSimID_TrackID(pion);
            auto pion_edep = mc_data->GetEdepID_TrackID(pion);
            SetLabels(
                pion_det_sim,
                pion_edep,
                pion,
                TopologyLabel::Track,
                PhysicsMicroLabel::HIPIonization,
                PhysicsMesoLabel::HIP,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            ProcessShowers(pion_daughters, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessKaon0s(
        const Parameters &config, art::Event const &event)
    {
        /**
         * 
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing neutral kaons."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto kaon0s = mc_data->GetPrimaries_AnsPDGCode(311);
        for (auto kaon0 : kaon0s) {
            ProcessShowers(kaon0, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessKaons(
        const Parameters &config, art::Event const &event)
    {
        /**
         * 
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing kaon+s."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto kaons = mc_data->GetTrackID_AbsPDGCode(321);
        for (auto kaon : kaons)
        {
            auto kaon_daughters = mc_data->GetDaughterTrackID_TrackID(kaon);
            auto kaon_det_sim = mc_data->GetDetSimID_TrackID(kaon);
            auto kaon_edep = mc_data->GetEdepID_TrackID(kaon);
            SetLabels(
                kaon_det_sim,
                kaon_edep,
                kaon,
                TopologyLabel::Track,
                PhysicsMicroLabel::HIPIonization,
                PhysicsMesoLabel::HIP,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            ProcessShowers(kaon_daughters, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessProtons(
        const Parameters &config, art::Event const &event)
    {
        /**
         * 
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing protons."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto protons = mc_data->GetTrackID_AbsPDGCode(2212);
        for (auto proton : protons)
        {
            auto proton_daughters = mc_data->GetDaughterTrackID_TrackID(proton);
            auto proton_det_sim = mc_data->GetDetSimID_TrackID(proton);
            auto proton_edep = mc_data->GetEdepID_TrackID(proton);
            SetLabels(
                proton_det_sim,
                proton_edep,
                proton,
                TopologyLabel::Track,
                PhysicsMicroLabel::HIPIonization,
                PhysicsMesoLabel::HIP,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            ProcessShowers(proton_daughters, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessNeutronCaptures(
        const Parameters &config, art::Event const &event)
    {
        /**
         * Neutron captures produce a standard candle of 6.1 MeV
         * gammas, which are generated according to a cascade: 
         *
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing neutrons."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto neutrons = mc_data->GetTrackID_PDGCode(2112);
        auto elastic_neutrons = mc_data->FilterTrackID_Process(neutrons, ProcessType::HadronElastic);
        auto inelastic_neutrons = mc_data->FilterTrackID_Process(neutrons, ProcessType::HadronInelastic);

        auto elastic_neutrons_det_sim = mc_data->GetDetSimID_TrackID(elastic_neutrons);
        auto elastic_neutrons_edep = mc_data->GetEdepID_TrackID(elastic_neutrons);
        SetLabels(
            elastic_neutrons_det_sim,
            elastic_neutrons_edep,
            elastic_neutrons,
            TopologyLabel::Blip,
            PhysicsMicroLabel::HadronElastic,
            PhysicsMesoLabel::NuclearRecoil,
            PhysicsMacroLabel::Undefined,
            UniqueTopology(),
            UniquePhysicsMicro(),
            UniquePhysicsMeso(),
            0,
        );

        auto inelastic_neutrons_det_sim = mc_data->GetDetSimID_TrackID(inelastic_neutrons);
        auto inelastic_neutrons_edep = mc_data->GetEdepID_TrackID(inelastic_neutrons);
        SetLabels(
            inelastic_neutrons_det_sim,
            inelastic_neutrons_edep,
            inelastic_neutrons,
            TopologyLabel::Blip,
            PhysicsMicroLabel::HadronInelastic,
            PhysicsMesoLabel::NuclearRecoil,
            PhysicsMacroLabel::Undefined,
            UniqueTopology(),
            UniquePhysicsMicro(),
            UniquePhysicsMeso(),
            0,
        );

        for (auto neutron : neutrons)
        {
            auto neutron_daughters = mc_data->GetDaughterTrackID_TrackID(neutron);
            auto gamma_daughters = mc_data->FilterTrackID_AbsPDGCode(neutron_daughters, 22);
            auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(neutron_daughters, 22);
            auto capture_daughters = mc_data->FilterTrackID_Process(gamma_daughters, ProcessType::NeutronCapture);
            auto other_gammas = mc_data->FilterTrackID_NotProcess(gamma_daughters, ProcessType::NeutronCapture);
            for (auto capture : capture_daughters)
            {
                //Int_t neutron_label = UniqueTopology();
                for (auto gamma : capture)
                {
                    Double_t gamma_energy = mc_data->GetEnergy_TrackID(gamma, 5);
                    auto gamma_det_sim = mc_data->GetAllDetSimID_TrackID(gamma);
                    auto gamma_edep = mc_data->GetAllEdepID_TrackID(gamma);

                    auto physics_meso_label = PhysicsLabel::NeutronCaptureGammaOther;
                    if (gamma_energy == 0.00474 || gamma_energy == 0.00475)
                    {
                        physics_meso_label = PhysicsLabel::NeutronCaptureGamma474;
                    }
                    else if (gamma_energy == 0.00336 || gamma_energy == 0.00337)
                    {
                        physics_meso_label = PhysicsLabel::NeutronCaptureGamma336;
                    }
                    else if (gamma_energy == 0.00256 || gamma_energy == 0.00257)
                    {
                        physics_meso_label = PhysicsLabel::NeutronCaptureGamma256;
                    }
                    else if (gamma_energy == 0.00118 || gamma_energy == 0.00119)
                    {
                        physics_meso_label = PhysicsLabel::NeutronCaptureGamma118;
                    }
                    else if (gamma_energy == 0.00083 || gamma_energy == 0.00084)
                    {
                        physics_meso_label = PhysicsLabel::NeutronCaptureGamma083;
                    }
                    else if (gamma_energy == 0.00051 || gamma_energy == 0.00052)
                    {
                        physics_meso_label = PhysicsLabel::NeutronCaptureGamma051;
                    }
                    else if (gamma_energy == 0.00016 || gamma_energy == 0.00017)
                    {
                        physics_meso_label = PhysicsLabel::NeutronCaptureGamma016;
                    }
                    SetLabels(
                        gamma_det_sim, 
                        gamma_edep, 
                        gamma,
                        TopologyLabel::Blip,
                        PhysicsMicroLabel::GammaCompton,
                        physics_meso_label,
                        PhysicsMacroLabel::Undefined,
                        UniqueTopology(),
                        UniquePhysicsMicro(),
                        UniquePhysicsMeso(),
                        0,
                    );
                }
            }
            ProcessShowers(other_daughters, UniqueTopology());
            ProcessShowers(other_gammas, UniqueTopology());
        }
    }
    void SimulationLabelingLogic::ProcessNuclearRecoils(
        const Parameters &config, art::Event const &event)
    {
        /**
         * Nuclear recoils can come from many things, but are essentially
         * edeps created by Ar41-Ar36 particles (with PDGCodes
         * 1000180400-1000180360), but also by fission from neutron inelastic
         * interactions, which generate an alpha and Sulfur 35 or
         * Chlorine 36.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing nuclear recoils."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto ar41 = mc_data->GetTrackID_PDGCode(1000180410);
        auto ar40 = mc_data->GetTrackID_PDGCode(1000180400);
        auto ar39 = mc_data->GetTrackID_PDGCode(1000180390);
        auto ar38 = mc_data->GetTrackID_PDGCode(1000180380);
        auto ar37 = mc_data->GetTrackID_PDGCode(1000180370);
        auto ar36 = mc_data->GetTrackID_PDGCode(1000180360);
        auto ar41_daughters = mc_data->GetDaughterTrackID_TrackID(ar41);
        auto ar40_daughters = mc_data->GetDaughterTrackID_TrackID(ar40);
        auto ar39_daughters = mc_data->GetDaughterTrackID_TrackID(ar39);
        auto ar38_daughters = mc_data->GetDaughterTrackID_TrackID(ar38);
        auto ar37_daughters = mc_data->GetDaughterTrackID_TrackID(ar37);
        auto ar36_daughters = mc_data->GetDaughterTrackID_TrackID(ar36);

        auto s33 = mc_data->GetTrackID_PDGCode(1000160330);
        auto s33_daughters = mc_data->GetDaughterTrackID_TrackID(s33);
        auto s35 = mc_data->GetTrackID_PDGCode(1000160350);
        auto s35_daughters = mc_data->GetDaughterTrackID_TrackID(s35);
        auto s36 = mc_data->GetTrackID_PDGCode(1000160360);
        auto s36_daughters = mc_data->GetDaughterTrackID_TrackID(s36);

        auto cl36 = mc_data->GetTrackID_PDGCode(1000170360);
        auto cl36_daughters = mc_data->GetDaughterTrackID_TrackID(cl36);
        auto cl37 = mc_data->GetTrackID_PDGCode(1000170370);
        auto cl37_daughters = mc_data->GetDaughterTrackID_TrackID(cl37);
        auto cl39 = mc_data->GetTrackID_PDGCode(1000170390);
        auto cl39_daughters = mc_data->GetDaughterTrackID_TrackID(cl39);
        auto cl40 = mc_data->GetTrackID_PDGCode(1000170400);
        auto cl40_daughters = mc_data->GetDaughterTrackID_TrackID(cl40);

        for (auto ar : ar41)
        {
            auto ar41_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar41_edep = mc_data->GetAllEdepID_TrackID(ar);
            SetLabels(
                ar41_det_sim,
                ar41_edep,
                ar,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto ar : ar40)
        {
            auto ar40_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar40_edep = mc_data->GetAllEdepID_TrackID(ar);
            SetLabels(
                ar40_det_sim,
                ar40_edep,
                ar,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto ar : ar39)
        {
            auto ar39_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar39_edep = mc_data->GetAllEdepID_TrackID(ar);
            SetLabels(
                ar39_det_sim,
                ar39_edep,
                ar,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto ar : ar38)
        {
            auto ar38_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar38_edep = mc_data->GetAllEdepID_TrackID(ar);
            SetLabels(
                ar38_det_sim,
                ar38_edep,
                ar,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto ar : ar37)
        {
            auto ar37_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar37_edep = mc_data->GetAllEdepID_TrackID(ar);
            SetLabels(
                ar37_det_sim,
                ar37_edep,
                ar,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto ar : ar36)
        {
            auto ar36_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar36_edep = mc_data->GetAllEdepID_TrackID(ar);
            SetLabels(
                ar36_det_sim,
                ar36_edep,
                ar,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto s : s33)
        {
            auto s33_det_sim = mc_data->GetAllDetSimID_TrackID(s);
            auto s33_edep = mc_data->GetAllEdepID_TrackID(s);
            SetLabels(
                s33_det_sim,
                s33_edep,
                s,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto s : s35)
        {
            auto s35_det_sim = mc_data->GetAllDetSimID_TrackID(s);
            auto s35_edep = mc_data->GetAllEdepID_TrackID(s);
            SetLabels(
                s35_det_sim,
                s35_edep,
                s,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto s : s36)
        {
            auto s36_det_sim = mc_data->GetAllDetSimID_TrackID(s);
            auto s36_edep = mc_data->GetAllEdepID_TrackID(s);
            SetLabels(
                s36_det_sim,
                s36_edep,
                s,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto cl : cl36)
        {
            auto cl36_det_sim = mc_data->GetAllDetSimID_TrackID(cl);
            auto cl36_edep = mc_data->GetAllEdepID_TrackID(cl);
            SetLabels(
                cl36_det_sim,
                cl36_edep,
                cl,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto cl : cl37)
        {
            auto cl37_det_sim = mc_data->GetAllDetSimID_TrackID(cl);
            auto cl37_edep = mc_data->GetAllEdepID_TrackID(cl);
            SetLabels(
                cl37_det_sim,
                cl37_edep,
                cl,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto cl : cl39)
        {
            auto cl39_det_sim = mc_data->GetAllDetSimID_TrackID(cl);
            auto cl39_edep = mc_data->GetAllEdepID_TrackID(cl);
            SetLabels(
                cl39_det_sim,
                cl39_edep,
                cl,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto cl : cl40)
        {
            auto cl40_det_sim = mc_data->GetAllDetSimID_TrackID(cl);
            auto cl40_edep = mc_data->GetAllEdepID_TrackID(cl);
            SetLabels(
                cl40_det_sim,
                cl40_edep,
                cl,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::NuclearRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
    }
    void SimulationLabelingLogic::ProcessElectronRecoils(
        const Parameters &config, art::Event const &event)
    {
        /**
         * Electron recoils can come from many things, such as
         * edeps created by deuterons/tritons/alphas coming out of
         * neutron inelastic scatters.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing electron recoils."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto deuterons = mc_data->GetTrackID_PDGCode(1000010020);
        auto tritons = mc_data->GetTrackID_PDGCode(1000010030);
        auto alphas = mc_data->GetTrackID_PDGCode(1000020040);
        auto inelastic_alphas = mc_data->FilterTrackID_Process(alphas, ProcessType::NeutronInelastic);
        auto deuteron_daughters = mc_data->GetDaughterTrackID_TrackID(deuterons);
        auto triton_daughters = mc_data->GetDaughterTrackID_TrackID(tritons);
        auto inelastic_alpha_daughters = mc_data->GetDaughterTrackID_TrackID(inelastic_alphas);
        for (auto deuteron : deuterons)
        {
            auto deuteron_det_sim = mc_data->GetAllDetSimID_TrackID(deuteron);
            auto deuteron_edep = mc_data->GetAllEdepID_TrackID(deuteron);
            SetLabels(
                deuteron_det_sim,
                deuteron_edep,
                deuteron,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::ElectronRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto triton : tritons)
        {
            auto triton_det_sim = mc_data->GetAllDetSimID_TrackID(triton);
            auto triton_edep = mc_data->GetAllEdepID_TrackID(triton);
            SetLabels(
                triton_det_sim,
                triton_edep,
                triton,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::ElectronRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
        for (auto inelastic_alpha : inelastic_alphas)
        {
            auto inelastic_alpha_det_sim = mc_data->GetAllDetSimID_TrackID(inelastic_alpha);
            auto inelastic_alpha_edep = mc_data->GetAllEdepID_TrackID(inelastic_alpha);
            SetLabels(
                inelastic_alpha_det_sim,
                inelastic_alpha_edep,
                inelastic_alpha,
                TopologyLabel::Blip,
                PhysicsMicroLabel::HadronElastic,
                PhysicsMesoLabel::ElectronRecoil,
                PhysicsMacroLabel::Undefined,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
        }
    }
    void SimulationLabelingLogic::ProcessAr39(
        const Parameters &config, art::Event const &event)
    {
        /**
         * Argon-39 decays via beta decay into Potassium-39,
         * with a Q-value of 565 keV: http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=180039.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing argon-39."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto ar39 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Ar39);
        auto ar39_daughters = mc_data->GetDaughterTrackID_TrackID(ar39);
        for(auto ar : ar39)
        {
            auto ar39_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar39_edep = mc_data->GetAllEdepID_TrackID(ar);
            SetLabels(
                ar39_det_sim,
                ar39_edep,
                ar39,
                TopologyLabel::Blip,
                PhysicsMicroLabel::ElectronIonization,
                PhysicsMesoLabel::BetaDecay,
                PhysicsMacroLabel::Ar39,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            mc_data->SetLabelingFunction_TrackID(ar, LabelCast(LabelingFunction::ProcessAr39));
        }
    }
    void SimulationLabelingLogic::ProcessAr42(
        const Parameters &config, art::Event const &event)
    {
        /**
         * Argon-42 decays via beta decay into Potassium-42,
         * with a Q-value of 599 keV: http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=180042.
         * There is a subtlety in the way this decay is simulated.  The lifetime of K42 is approximately
         * 12 hours, which beta decays to Calcium-42 with a 3.5 MeV beta.  Calcium-42 is stable.
         * See here for some details: https://indico.fnal.gov/event/50121/contributions/220205/attachments/145404/185102/20210721_Decay0_Lasorak.pdf
         * 
         * If the energy of the decay primary is 600 keV or less, we label it as an
         * Ar42 decay, otherwise it's a K42 decay.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing argon-42."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto ar42 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Ar42);
        auto ar42_betas = mc_data->FilterTrackID_PDGCode(ar42, 1000180420);
        auto ar42_daughters = mc_data->GetDaughterTrackID_TrackID(ar42_betas);
        for(auto ar : ar42)
        {
            auto ar42_det_sim = mc_data->GetAllDetSimID_TrackID(ar);
            auto ar42_edep = mc_data->GetAllEdepID_TrackID(ar);
            if (mc_data->GetEnergy_TrackID(ar) < 0.600)
            {
                SetLabels(
                ar42_det_sim,
                ar42_edep,
                ar42,
                TopologyLabel::Blip,
                PhysicsMicroLabel::ElectronIonization,
                PhysicsMesoLabel::BetaDecay,
                PhysicsMacroLabel::Ar42,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            }
            else 
            {
                SetLabels(
                ar42_det_sim,
                ar42_edep,
                ar42,
                TopologyLabel::Blip,
                PhysicsMicroLabel::ElectronIonization,
                PhysicsMesoLabel::BetaDecay,
                PhysicsMacroLabel::K42,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            }
            mc_data->SetLabelingFunction_TrackID(ar, LabelCast(LabelingFunction::ProcessAr42));
        }
    }
    void SimulationLabelingLogic::ProcessKr85(
        const Parameters &config, art::Event const &event)
    {
        /**
         * Krypton-85 decays via beta decay into Rubidium 85
         * with two prominent betas with energies of 687 keV (99.56 %) and
         * 173 keV (.43 %): http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=360085.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing krypton-85."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto kr85 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Kr85);
        auto kr85_daughters = mc_data->GetDaughterTrackID_TrackID(kr85);
        for(auto kr : kr85)
        {
            auto kr85_det_sim = mc_data->GetAllDetSimID_TrackID(kr);
            auto kr85_edep = mc_data->GetAllEdepID_TrackID(kr);
            SetLabels(
                kr85_det_sim,
                kr85_edep,
                kr85,
                TopologyLabel::Blip,
                PhysicsMicroLabel::ElectronIonization,
                PhysicsMesoLabel::BetaDecay,
                PhysicsMacroLabel::Kr85,
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            mc_data->SetLabelingFunction_TrackID(kr, LabelCast(LabelingFunction::ProcessKr85));
        }
    }
    void SimulationLabelingLogic::ProcessRn222(
        const Parameters &config, art::Event const &event)
    {
        /**
         * Radon-222 decays via alpha decay through a chain that ends in lead
         * (https://en.wikipedia.org/wiki/Radon-222).  The alpha has an energy of
         * 5.5904 MeV, which bounces around locally in Argon, but quickly thermalizes
         * due to the short scattering length.  The CSDA range of a 5.5 MeV alpha in
         * Argon is about 7.5e-3 g/cm^2.  Using a density of 1.3954 g/cm^3, the
         * scattering length is (~0.005 cm) or (~50 um).
         * 
         * The simulation of radiologicals is done through a wrap of Decay0, which does
         * not save information about what decay each simulated particle came from.
         * Instead, it saves each decay product (beta, gamma, alpha) as a primary of its
         * own.  This can be seen from the following lines of the larsimrad file
         * https://github.com/LArSoft/larsimrad/blob/develop/larsimrad/BxDecay0/Decay0Gen_module.cc,
         * which has (starting at line 162):
         * if (p.is_alpha()) {
              pdg = 1000020040;
              part = simb::MCParticle(track_id, pdg, primary_str, -1, mass / 1000, 1);
            }
            else if (p.is_gamma()) {
              pdg = 22;
              part = simb::MCParticle(track_id, pdg, primary_str);
            }
            else if (p.is_positron()) {
              pdg = -11;
              part = simb::MCParticle(track_id, pdg, primary_str);
            }
            else if (p.is_electron()) {
              pdg = 11;
              part = simb::MCParticle(track_id, pdg, primary_str);
            }
            else if (p.is_neutron()) {
              pdg = 2112;
              part = simb::MCParticle(track_id, pdg, primary_str);
            }
            else if (p.is_proton()) {
              pdg = 2212;
              part = simb::MCParticle(track_id, pdg, primary_str);
            }
            else {
              p.print(std::cout);
              throw cet::exception("Decay0Gen") << "Particle above is weird, cannot recognise it.";
            }
         * We will therefore need some way to figure out which particle came from what decay.
         * The possible decays are:
         * Rn222 -> Po218 - alpha ~ 5.590 MeV
         * Po218 -> Pb214 - alpha ~ 6.115 MeV
         * Po218 -> At218 - beta  ~ 0.294 MeV 
         * At218 -> Bi214 - alpha ~ 6.874 MeV
         * At218 -> Rn218 - beta  ~ 2.883 MeV
         * Rn218 -> Po214 - alpha ~ 7.263 MeV
         * Pb214 -> Bi214 - beta  ~ 1.024 MeV
         * Bi214 -> Tl210 - alpha ~ 5.617 MeV
         * Bi214 -> Po214 - beta  ~ 3.272 MeV
         * Po214 -> Pb210 - alpha ~ 7.833 MeV
         * Tl210 -> Pb210 - beta  ~ 5.484 MeV
         * Pb210 -> Bi210 - beta  ~ 0.064 MeV
         * Pb210 -> Hg206 - alpha ~ 3.792 MeV
         * Bi210 -> Po210 - beta  ~ 1.163 MeV
         * Bi210 -> Tl206 - alpha ~ 5.037 MeV
         * Po210 -> Pb206 - alpha ~ 5.407 MeV
         * 
         * We therefore have 9 disctint alpha energies 
         * (7.833, 7.263, 6.874, 6.115, 5.617, 5.590, 5.407, 5.037, 3.792)
         * 
         * and 7 distinct beta energies
         * (5.484, 3.272, 2.883, 1.163, 1.024, 0.294, 0.064)
         * 
         * We can assign labels then based on whatever energy is closest to each primary.
         */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing radon-222."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto rn222 = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Rn222);
        for(auto rn : rn222)
        {
            Int_t index = 0;
            Double_t energy_diff = 10e10;
            for(Int_t ii = 0; ii < mRn222Energies.size(); ii++)
            {
                if (mc_data->GetPDGCode_TrackID(rn) != mRn222PDGs[ii]) {
                    continue;
                }
                if (abs(mc_data->GetEnergy_TrackID(rn) - mRn222Energies[ii]) < energy_diff)
                {
                    index = ii;
                    energy_diff = mc_data->GetEnergy_TrackID(rn) - mRn222Energies[ii];
                }
            }
            auto rn222_det_sim = mc_data->GetAllDetSimID_TrackID(rn);
            auto rn222_edep = mc_data->GetAllEdepID_TrackID(rn);
            SetLabels(
                rn222_det_sim,
                rn222_edep,
                rn222,
                TopologyLabel::Blip,
                PhysicsMicroLabel::ElectronIonization,
                mRn222PhysicsMeso[index],
                mRn222PhysicsMacro[index],
                UniqueTopology(),
                UniquePhysicsMicro(),
                UniquePhysicsMeso(),
                0,
            );
            mc_data->SetLabelingFunction_TrackID(rn, LabelCast(LabelingFunction::ProcessRn222));
        }
    }
    void SimulationLabelingLogic::ProcessCosmics(
        const Parameters &config, art::Event const &event)
    {
        /**
         * 
        */
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "processing cosmics."
        );
        auto mc_data = SimulationWrangler::GetInstance();
        auto cosmics = mc_data->GetPrimaries_GeneratorLabel(GeneratorLabel::Cosmics);
        auto electrons = mc_data->FilterTrackID_PDGCode(cosmics, 11);
        auto positrons = mc_data->FilterTrackID_PDGCode(cosmics, -11);
        auto gammas = mc_data->FilterTrackID_AbsPDGCode(cosmics, 22);
        auto neutrons = mc_data->FilterTrackID_PDGCode(cosmics, 2112);
        auto anti_neutrons = mc_data->FilterTrackID_PDGCode(cosmics, -2112);
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "number of cosmic electrons: " + std::to_string(electrons.size())
        );
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "number of cosmic positrons: " + std::to_string(positrons.size())
        );
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "number of cosmic gammas: " + std::to_string(gammas.size())
        );
        Logger::GetInstance("SimulationLabelingLogic")->trace(
            "number of cosmic neutrons: " + std::to_string(neutrons.size())
        );
         Logger::GetInstance("SimulationLabelingLogic")->trace(
            "number of cosmic anti-neutrons: " + std::to_string(anti_neutrons.size())
        );
        // for (auto electron : electrons)
        // {
        //     auto electron_daughters = mc_data->GetDaughterTrackID_TrackID(electron);
        //     auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(electron_daughters, 11);
        //     auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(electron_daughters, 11);
        //     auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
        //     auto elec_edep = mc_data->GetEdepID_TrackID(elec_daughters);

        //     auto electron_progeny = mc_data->GetProgenyTrackID_TrackID(electron);
        //     auto electron_det_sim = mc_data->GetDetSimID_TrackID(electron);
        //     auto electron_edep = mc_data->GetEdepID_TrackID(electron);

        //     // Set electron detsim labels to Shower::ElectronShower
        //     Int_t shower_label = UniqueTopology();
        //     SetLabels(
        //         electron_det_sim, electron_edep, electron,
        //         TopologyLabel::Shower, PhysicsLabel::ElectronRecoil, 
        //         shower_label
        //     );
        //     mc_data->SetLabelingFunction_TrackID(electron, LabelCast(LabelingFunction::ProcessCosmics));

        //     SetLabels(
        //         elec_det_sim, elec_edep, elec_daughters, 
        //         TopologyLabel::Shower, PhysicsLabel::ElectronRecoil, 
        //         shower_label
        //     );
        //     for (auto daughter : elec_daughters) {
        //         mc_data->SetLabelingFunction_TrackID(daughter, LabelCast(LabelingFunction::ProcessCosmics));
        //     }

        //     ProcessShowers(electron_progeny, shower_label);
        //     ProcessShowers(other_daughters, UniqueTopology());
        // }
        // for (auto positron : positrons)
        // {
        //     auto positron_daughters = mc_data->GetDaughterTrackID_TrackID(positron);
        //     auto elec_daughters = mc_data->FilterTrackID_AbsPDGCode(positron_daughters, 11);
        //     auto other_daughters = mc_data->FilterTrackID_NotAbsPDGCode(positron_daughters, 11);
        //     auto elec_det_sim = mc_data->GetDetSimID_TrackID(elec_daughters);
        //     auto elec_edep = mc_data->GetEdepID_TrackID(elec_daughters);

        //     auto positron_progeny = mc_data->GetProgenyTrackID_TrackID(positron);
        //     auto positron_det_sim = mc_data->GetDetSimID_TrackID(positron);
        //     auto positron_edep = mc_data->GetEdepID_TrackID(positron);

        //     // Set positron detsim labels to Shower::positronShower
        //     Int_t shower_label = UniqueTopology();
        //     SetLabels(
        //         positron_det_sim, positron_edep, positron, 
        //         TopologyLabel::Shower, PhysicsLabel::PositronShower, 
        //         shower_label
        //     );
        //     mc_data->SetLabelingFunction_TrackID(positron, LabelCast(LabelingFunction::ProcessCosmics));

        //     SetLabels(
        //         elec_det_sim, elec_edep, elec_daughters,
        //         TopologyLabel::Shower, PhysicsLabel::PositronShower, 
        //         shower_label
        //     );
        //     for (auto daughter : elec_daughters) {
        //         mc_data->SetLabelingFunction_TrackID(daughter, LabelCast(LabelingFunction::ProcessCosmics));
        //     }

        //     ProcessShowers(positron_progeny, shower_label);
        //     ProcessShowers(other_daughters, UniqueTopology());
        // }
    }
}