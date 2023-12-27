/**
 * @file EnergyDepositPointCloud.h
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
    struct EnergyDepositPointCloud
    {
        std::vector<Double_t> edep_t = {};
        std::vector<Double_t> edep_x = {};
        std::vector<Double_t> edep_y = {};
        std::vector<Double_t> edep_z = {};
        std::vector<Double_t> edep_energy = {};
        std::vector<Int_t> edep_num_photons = {};
        std::vector<Int_t> edep_num_electrons = {};
        std::vector<TrackID_t>  edep_track_id = {};
        std::vector<ProcessType> edep_process = {};

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

        std::vector<std::vector<DetSimID_t>> edep_detsim_id = {};

        EnergyDepositPointCloud()
        {
        }

        void clear()
        {
            edep_t.clear();
            edep_x.clear();
            edep_y.clear();
            edep_z.clear();
            edep_energy.clear();
            edep_num_photons.clear();
            edep_num_electrons.clear();
            edep_track_id.clear();
            edep_process.clear();
            edep_detsim_id.clear();

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

        void AddPoint(
            TrackID_t track_id, Double_t t, Double_t x, 
            Double_t y, Double_t z, Double_t energy, 
            Int_t num_photons, Int_t num_electrons,
            ProcessType process
        )
        {
            edep_track_id.emplace_back(track_id);
            edep_t.emplace_back(t);
            edep_x.emplace_back(x);
            edep_y.emplace_back(y);
            edep_z.emplace_back(z);
            edep_energy.emplace_back(energy);
            edep_num_photons.emplace_back(num_photons);
            edep_num_electrons.emplace_back(num_electrons);
            edep_process.emplace_back(process);
            edep_detsim_id.emplace_back(std::vector<DetSimID_t>());

            topology_label.emplace_back(LabelCast(TopologyLabel::Undefined));
            particle_label.emplace_back(LabelCast(ParticleLabel::Undefined));
            physics_micro_label.emplace_back(LabelCast(PhysicsMicroLabel::Undefined));
            physics_meso_label.emplace_back(LabelCast(PhysicsMesoLabel::Undefined));
            physics_macro_label.emplace_back(LabelCast(PhysicsMacroLabel::Undefined));
            
            unique_topology_label.emplace_back(-1);
            unique_particle_label.emplace_back(-1);
            unique_physics_micro_label.emplace_back(-1);
            unique_physics_meso_label.emplace_back(-1);
            unique_physics_macro_label.emplace_back(-1);
        }
    };
}