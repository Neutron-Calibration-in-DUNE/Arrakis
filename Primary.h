/**
 * @file Primary.h
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
#include "ParticleMaps.h"
#include "Trajectory.h"

namespace arrakis
{
    struct Primary
    {
        Int_t track_id = {0};
        Int_t pdg = {0};

        // MC Particle info.
        std::string init_process = {""};
        Double_t init_energy = {0};
        Double_t init_t = {0};
        Double_t init_x = {0};
        Double_t init_y = {0};
        Double_t init_z = {0};
        
        std::string end_process = {""};
        Double_t end_energy = {0};
        Double_t end_t = {0};
        Double_t end_x = {0};
        Double_t end_y = {0};
        Double_t end_z = {0};

        Trajectory primary_trajectory;

        // Energy deposits of the parent
        Double_t total_edep_energy = {0};
        std::vector<Double_t> edep_energy = {};
        std::vector<std::string> edep_process = {};
        std::vector<std::string> edep_volume = {};
        std::vector<std::string> edep_material = {};
        std::vector<Double_t> edep_t = {};
        std::vector<Double_t> edep_x = {};
        std::vector<Double_t> edep_y = {};
        std::vector<Double_t> edep_z = {};

        // Daughter MC Particle info.
        std::vector<Int_t> daughter_ids = {};
        std::map<Int_t, Int_t> daughter_map = {};
        std::vector<Int_t> daughter_level = {};

        std::vector<std::string> daughter_init_process = {};
        std::vector<Double_t> daughter_init_energy = {};
        std::vector<Double_t> daughter_init_t = {};
        std::vector<Double_t> daughter_init_x = {};
        std::vector<Double_t> daughter_init_y = {};
        std::vector<Double_t> daughter_init_z = {};

        std::vector<std::string> daughter_end_process = {};
        std::vector<Double_t> daughter_end_energy = {};
        std::vector<Double_t> daughter_end_t = {};
        std::vector<Double_t> daughter_end_x = {};
        std::vector<Double_t> daughter_end_y = {};
        std::vector<Double_t> daughter_end_z = {};

        std::vector<Trajectory> daughter_trajectories = {};

        // Daughter energy deposits.
        Double_t total_daughter_edep_energy = {0};
        std::vector<Int_t> daughter_edep_ids = {};
        std::vector<Double_t> daughter_edep_energy = {};
        std::vector<std::string> daughter_edep_process = {};
        std::vector<std::string> daughter_edep_volume = {};
        std::vector<std::string> daughter_edep_material = {};
        std::vector<Double_t> daughter_edep_t = {};
        std::vector<Double_t> daughter_edep_x = {};
        std::vector<Double_t> daughter_edep_y = {};
        std::vector<Double_t> daughter_edep_z = {};

        // Raw digit information.
        std::vector<Int_t> det_track_id = {};
        std::vector<Double_t> det_energy_fraction = {};
        std::vector<Double_t> det_energy = {};
        std::vector<Int_t> det_channel = {};
        std::vector<Int_t> det_tdc = {};
        std::vector<Int_t> det_adc = {};
        std::vector<Int_t> det_edep = {};
        std::vector<std::string> det_process = {};

        Primary(){}

        Primary(simb::MCParticle particle)
        {
            track_id = particle.TrackId();
            pdg = particle.PdgCode();
            init_process = particle.Process();
            init_energy = particle.E();
            init_t = particle.T();
            init_x = particle.Vx();
            init_y = particle.Vy();
            init_z = particle.Vz();
            end_process = particle.EndProcess();
            end_energy = particle.EndE();
            end_x = particle.EndX();
            end_y = particle.EndY();
            end_z = particle.EndZ();

            simb::MCTrajectory trajectory = particle.Trajectory();
            auto trajectory_processes = trajectory.TrajectoryProcesses();

            for(size_t ii = 0; ii < particle.NumberTrajectoryPoints(); ii++)
            {
                // Get the volume information for trajectory point.
                auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                   particle.Vx(ii), particle.Vy(ii), particle.Vz(ii)
                );
                std::string process = "not_defined";
                for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
                {
                    if(trajectory_processes[jj].first == ii) {
                        process = trajectory_processes[jj].second;
                    }
                }
                primary_trajectory.AddTrajectoryPoint(
                    particle.T(ii), particle.Vx(ii), particle.Vy(ii), particle.Vz(ii),
                    particle.E(ii), process, volume.volume_name,
                    volume.material_name
                );
            }

        }

        Primary(
            Int_t _track_id, Int_t _pdg, std::string _init_process,
            Double_t _init_energy, Double_t _init_t, Double_t _init_x, Double_t _init_y,
            Double_t _init_z, std::string _end_process, Double_t _end_energy,
            Double_t _end_t, Double_t _end_x, Double_t _end_y, Double_t _end_z,
            simb::MCTrajectory _trajectory
        )
        : track_id(_track_id)
        , pdg(_pdg)
        , init_process(_init_process)
        , init_energy(_init_energy)
        , init_t(_init_t)
        , init_x(_init_x)
        , init_y(_init_y)
        , init_z(_init_z)
        , end_process(_end_process)
        , end_energy(_end_energy)
        , end_t(_end_t)
        , end_x(_end_x)
        , end_y(_end_y)
        , end_z(_end_z)
        {
            // Get the volume information for the energy deposit.
            auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                _init_x, _init_y, _init_z
            );
            
        }

        void AddEdep(
            Double_t energy, std::string process,
            Double_t t, Double_t x, Double_t y, Double_t z
        )
        {
            // Get the volume information for the energy deposit.
            auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                x, y, z
            );

            edep_energy.emplace_back(energy);
            edep_process.emplace_back(process);
            edep_volume.emplace_back(volume.volume_name);
            edep_material.emplace_back(volume.material_name);
            edep_t.emplace_back(t);
            edep_x.emplace_back(x);
            edep_y.emplace_back(y);
            edep_z.emplace_back(z);
            total_edep_energy += energy;
        }

        void AddDaughter(simb::MCParticle particle, Int_t level)
        {
            daughter_ids.emplace_back(particle.TrackId());
            daughter_level.emplace_back(level);
            daughter_init_process.emplace_back(particle.Process());
            daughter_init_energy.emplace_back(particle.E());
            daughter_init_t.emplace_back(particle.T());
            daughter_init_x.emplace_back(particle.Vx());
            daughter_init_y.emplace_back(particle.Vy());
            daughter_init_z.emplace_back(particle.Vz());
            daughter_end_process.emplace_back(particle.EndProcess());
            daughter_end_energy.emplace_back(particle.EndE());
            daughter_end_t.emplace_back(particle.EndT());
            daughter_end_x.emplace_back(particle.EndX());
            daughter_end_y.emplace_back(particle.EndY());
            daughter_end_z.emplace_back(particle.EndZ());

            simb::MCTrajectory trajectory = particle.Trajectory();
            auto trajectory_processes = trajectory.TrajectoryProcesses();
            Trajectory daughter_trajectory;

            for(size_t ii = 0; ii < particle.NumberTrajectoryPoints(); ii++)
            {
                // Get the volume information for trajectory point.
                auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                   particle.Vx(ii), particle.Vy(ii), particle.Vz(ii)
                );
                std::string process = "not_defined";
                for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
                {
                    if(trajectory_processes[jj].first == ii) {
                        process = trajectory_processes[jj].second;
                    }
                }
                daughter_trajectory.AddTrajectoryPoint(
                    particle.T(ii), particle.Vx(ii), particle.Vy(ii), particle.Vz(ii),
                    particle.E(ii), process, volume.volume_name,
                    volume.material_name
                );
            }
            daughter_trajectories.emplace_back(daughter_trajectory);
        }

        void AddDaughter(
            Int_t track_id, Int_t level,
            std::string init_process, Double_t init_energy,
            Double_t init_t, Double_t init_x, Double_t init_y, Double_t init_z,
            std::string end_process, Double_t end_energy,
            Double_t end_t, Double_t end_x, Double_t end_y, Double_t end_z,
            simb::MCTrajectory trajectory
        )
        {
            // Get the volume information for the energy deposit.
            auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                init_x, init_y, init_z
            );

            daughter_ids.emplace_back(track_id);
            daughter_level.emplace_back(level);
            daughter_init_process.emplace_back(init_process);
            daughter_init_energy.emplace_back(init_energy);
            daughter_init_t.emplace_back(init_t);
            daughter_init_x.emplace_back(init_x);
            daughter_init_y.emplace_back(init_y);
            daughter_init_z.emplace_back(init_z);
            daughter_end_process.emplace_back(end_process);
            daughter_end_energy.emplace_back(end_energy);
            daughter_end_t.emplace_back(end_t);
            daughter_end_x.emplace_back(end_x);
            daughter_end_y.emplace_back(end_y);
            daughter_end_z.emplace_back(end_z);
        }

        void AddDaughterEdep(
            Int_t track_id, Double_t energy, std::string process,
            Double_t t, Double_t x, Double_t y, Double_t z
        )
        {
            // Get the volume information for the energy deposit.
            auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                x, y, z
            );

            daughter_edep_ids.emplace_back(track_id);
            daughter_edep_energy.emplace_back(energy);
            daughter_edep_process.emplace_back(process);
            daughter_edep_volume.emplace_back(volume.volume_name);
            daughter_edep_material.emplace_back(volume.material_name);
            daughter_edep_t.emplace_back(t);
            daughter_edep_x.emplace_back(x);
            daughter_edep_y.emplace_back(y);
            daughter_edep_z.emplace_back(z);
            total_daughter_edep_energy += energy;
        }

        void AddDetectorSimulation(
            Int_t track_id,
            Double_t energy_frac,
            Double_t energy,
            Int_t channel,
            Int_t tdc,
            Int_t adc
        )
        {
            det_track_id.emplace_back(track_id);
            det_energy_fraction.emplace_back(energy_frac);
            det_energy.emplace_back(energy);
            det_channel.emplace_back(channel);
            det_tdc.emplace_back(tdc);
            det_adc.emplace_back(adc);
        }
    };
}