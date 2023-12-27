/**
 * @file OpDetPointCloud.h
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
#include "lardataobj/RawData/OpDetWaveform.h"

// necessary ROOT libraries
#include <TTree.h>

#include "Core.h"
#include "DetectorGeometry.h"

namespace arrakis
{
    struct OpDetPointCloud
    {
        /**
         * A OpDetPointCloud is the basic data structure that
         * we use for data from the optical detector output.  
         * 1 tdc = 6.7 ns, whereas for wires its 500 ns. 
         * 
        */
        std::vector<Int_t> channel = {};
        std::vector<Int_t> tick = {};
        std::vector<Int_t> adc = {};
        std::vector<Double_t> energy = {};

        std::vector<std::vector<Int_t>> track_ids = {};
        std::vector<std::vector<Double_t>> energies = {};
        std::vector<std::vector<Double_t>> x = {};
        std::vector<std::vector<Double_t>> y = {};
        std::vector<std::vector<Double_t>> z = {};
        std::vector<std::vector<TopologyLabelInt>> topology_labels = {};
        std::vector<std::vector<ParticleLabelInt>> particle_labels = {};
        std::vector<std::vector<PhysicsMicroLabelInt>>  physics_micro_labels = {};
        std::vector<std::vector<PhysicsMesoLabelInt>>   physics_meso_labels = {};
        std::vector<std::vector<PhysicsMacroLabelInt>>  physics_macro_labels = {};
        std::vector<std::vector<Int_t>> unique_topology_labels = {};
        std::vector<std::vector<Int_t>> unique_particle_labels = {};
        std::vector<std::vector<Int_t>> unique_physics_micro_labels = {};
        std::vector<std::vector<Int_t>> unique_physics_meso_labels = {};
        std::vector<std::vector<Int_t>> unique_physics_macro_labels = {};

        std::vector<TopologyLabelInt> topology_label = {};
        std::vector<ParticleLabelInt> particle_label = {};
        std::vector<PhysicsMicroLabelInt> physics_micro_label = {};
        std::vector<PhysicsMesoLabelInt> physics_meso_label = {};
        std::vector<PhysicsMacroLabelInt> physics_macro_label = {};

        std::vector<Int_t> unique_topology_label = {};
        std::vector<Int_t> unique_particle_label = {};
        std::vector<Int_t> unique_physics_micro_label = {};
        std::vector<Int_t> unique_physics_meso_label = {};
        std::vector<Int_t> unique_physics_macro_label = {};



        OpDetPointCloud()
        {
        }

        void clear()
        {
            channel.clear();
            tick.clear();
            adc.clear();
            energy.clear();

            track_ids.clear();
            energies.clear();
            x.clear();
            y.clear();
            z.clear();
            topology_labels.clear();
            particle_labels.clear();
            physics_micro_labels.clear();
            physics_meso_labels.clear();
            physics_macro_labels.clear();
            unique_topology_labels.clear();
            unique_particle_labels.clear();
            unique_physics_micro_labels.clear();
            unique_physics_meso_labels.clear();
            unique_physics_macro_labels.clear();
            
            topology_label.clear();
            particle_label.clear();
            physics_micro_label.clear();
            physics_meso_label.clear();
            physics_macro_label.clear();
            unique_topology_label.clear();
            unique_particle_label.clear();
            unique_physics_micro_label.clear();
            unique_physics_meso_label.clear();
            unique_physics_macro_label.clear();
        }

        // This function is used primarily with Data,
        // which doesn't have other truth information
        // associated with the points.
        void AddPoint(
            Int_t det_tick,
            Int_t det_channel,
            Int_t det_adc
        )
        {
            channel.emplace_back(det_channel);
            tick.emplace_back(det_tick);
            adc.emplace_back(det_adc);
        }

        // This function is used for simulation, where
        // truth information about RawDigits is in the
        // sim::IDE vector.
        void AddPoint(
            Int_t det_channel,
            Int_t det_tick,
            Int_t det_adc,
            bool det_noise
        )
        {
            /**
             * 
            */
            channel.emplace_back(det_channel);
            tick.emplace_back(det_tick);
            adc.emplace_back(det_adc);
            

            std::vector<Int_t> det_track_ids;
            std::vector<Double_t> det_energies;
            std::vector<Double_t> det_x;
            std::vector<Double_t> det_y;
            std::vector<Double_t> det_z;
            std::vector<TopologyLabelInt> det_shape;
            std::vector<ParticleLabelInt> det_particle;
            std::vector<PhysicsMicroLabelInt> det_physics_micro;
            std::vector<PhysicsMesoLabelInt> det_physics_meso;
            std::vector<PhysicsMacroLabelInt> det_physics_macro;
            std::vector<Int_t> det_unique_topology;
            std::vector<Int_t> det_unique_particle;
            std::vector<Int_t> det_unique_physics_micro;
            std::vector<Int_t> det_unique_physics_meso;
            std::vector<Int_t> det_unique_physics_macro;
            Double_t det_energy = 0.0;
            for(auto ide : det_ide)
            {
                det_track_ids.emplace_back(ide.trackID);
                det_energies.emplace_back(ide.energy);
                det_x.emplace_back(ide.x);
                det_y.emplace_back(ide.y);
                det_z.emplace_back(ide.z);
                if(det_noise)
                {
                    det_shape.emplace_back(LabelCast(TopologyLabel::Noise));
                    det_particle.emplace_back(LabelCast(ParticleLabel::Noise));
                    det_physics_micro.emplace_back(LabelCast(PhysicsMicroLabel::Noise));
                    det_physics_meso.emplace_back(LabelCast(PhysicsMesoLabel::Noise));
                    det_physics_macro.emplace_back(LabelCast(PhysicsMacroLabel::Noise));
                }
                else
                {
                    det_shape.emplace_back(LabelCast(TopologyLabel::Undefined));
                    det_particle.emplace_back(LabelCast(ParticleLabel::Undefined));
                    det_physics_micro.emplace_back(LabelCast(PhysicsMicroLabel::Undefined));
                    det_physics_meso.emplace_back(LabelCast(PhysicsMesoLabel::Undefined));
                    det_physics_macro.emplace_back(LabelCast(PhysicsMacroLabel::Undefined));
                }
                det_unique_topology.emplace_back(-1);
                det_unique_particle.emplace_back(-1);
                det_unique_physics_micro.emplace_back(-1);
                det_unique_physics_meso.emplace_back(-1);
                det_unique_physics_macro.emplace_back(-1);
                det_energy += ide.energy;
            }
            
            energy.emplace_back(det_energy);
            track_ids.emplace_back(det_track_ids);
            energies.emplace_back(det_energies);
            x.emplace_back(det_x);
            y.emplace_back(det_y);
            z.emplace_back(det_z);

            topology_labels.emplace_back(det_shape);
            particle_labels.emplace_back(det_particle);
            physics_micro_labels.emplace_back(det_physics_micro);
            physics_meso_labels.emplace_back(det_physics_meso);
            physics_macro_labels.emplace_back(det_physics_macro);
            unique_topology_labels.emplace_back(det_unique_topology);
            unique_particle_labels.emplace_back(det_unique_particle);
            unique_physics_micro_labels.emplace_back(det_unique_physics_micro);
            unique_physics_meso_labels.emplace_back(det_unique_physics_meso);
            unique_physics_macro_labels.emplace_back(det_unique_physics_macro);

            if(det_noise)
            {
                topology_label.emplace_back(LabelCast(TopologyLabel::Noise));
                particle_label.emplace_back(LabelCast(ParticleLabel::Noise));
                physics_micro_label.emplace_back(LabelCast(PhysicsMicroLabel::Noise));
                physics_meso_label.emplace_back(LabelCast(PhysicsMesoLabel::Noise));
                physics_macro_label.emplace_back(LabelCast(PhysicsMacroLabel::Noise));
            }
            else
            {
                topology_label.emplace_back(LabelCast(TopologyLabel::Undefined));
                particle_label.emplace_back(LabelCast(ParticleLabel::Undefined));
                physics_micro_label.emplace_back(LabelCast(PhysicsMicroLabel::Undefined));
                physics_meso_label.emplace_back(LabelCast(PhysicsMesoLabel::Undefined));
                physics_macro_label.emplace_back(LabelCast(PhysicsMacroLabel::Undefined));
            }
            unique_topology_label.emplace_back(-1);
            unique_particle_label.emplace_back(-1);
            unique_physics_micro_label.emplace_back(-1);
            unique_physics_meso_label.emplace_back(-1);
            unique_physics_macro_label.emplace_back(-1);
        }
        Int_t GetIndex_TrackID(DetSimID_t detsim, TrackID_t track_id)
        {
            for(size_t ii = 0; ii < track_ids[detsim].size(); ii++) {
                if(track_ids[detsim][ii] == track_id) {
                    return ii;
                }
            }
            return -1;
        }
    };
}