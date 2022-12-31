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
#include "ParticleMaps.h"
#include "Trajectory.h"

namespace arrakis
{
    struct Primary
    {
        GeneratorLabel generator_label = {kNone};
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

        // Energy deposits of the parent.
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
        std::vector<Int_t> daughter_pdgs = {};
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
        std::vector<Double_t> det_energy_fraction = {};
        std::vector<Double_t> det_energy = {};
        std::vector<Int_t> det_channel = {};
        std::vector<Int_t> det_tick = {};
        std::vector<Int_t> det_adc = {};
        std::vector<Int_t> det_edep = {};
        std::vector<Double_t> det_tdc = {};
        std::vector<std::string> det_process = {};

        // Daughter raw digit information.
        std::vector<Int_t> daughter_det_track_id = {};
        std::vector<Double_t> daughter_det_energy_fraction = {};
        std::vector<Double_t> daughter_det_energy = {};
        std::vector<Int_t> daughter_det_channel = {};
        std::vector<Int_t> daughter_det_tick = {};
        std::vector<Int_t> daughter_det_adc = {};
        std::vector<Int_t> daughter_det_edep = {};
        std::vector<Double_t> daughter_det_tdc = {};
        std::vector<std::string> daughter_det_process = {};

        /**
         * Here we are trying to match up the point at which
         * the energy deposition was created with the MC Particle
         * process that caused the energy deposition.  We do this
         * by checking if the local time of the energy deposition
         * encloses the trajectory point.
        */
        std::string FindPrimaryEnergyDepositionProcess(sim::SimEnergyDeposit& edep)
        {
            std::string process = "not_found";
            for(size_t ii = 0; ii < primary_trajectory.t.size(); ii++)
            {
                if(
                    edep.StartT() <= primary_trajectory.t[ii] &&
                    edep.EndT() >= primary_trajectory.t[ii]
                )
                {
                    process = primary_trajectory.process[ii];
                }
            }
            return process;
        }
        /**
         * Same thing as above except we loop over all daughter trajectories.
        */
        std::string FindDaughterEnergyDepositionProcess(sim::SimEnergyDeposit& edep)
        {
            std::string process = "not_found";
            if(daughter_map.find(edep.TrackID()) == daughter_map.end()) { 
                return "invalid_track_id";
            }
            Int_t daughter_index = daughter_map[edep.TrackID()];
            for(size_t ii = 0; ii < daughter_trajectories[daughter_index].t.size(); ii++)
            {
                if(
                    edep.StartT() <= daughter_trajectories[daughter_index].t[ii] &&
                    edep.EndT() >= daughter_trajectories[daughter_index].t[ii])
                {
                    process = daughter_trajectories[daughter_index].process[ii];
                }
            }
            return process;
        }

        /**
         * This looks for the closest energy deposition of the primary 
         * in (x,y,z) to the specified points.
        */
        Int_t FindPrimaryEnergyDeposition(
            Double_t x, Double_t y, Double_t z
        )
        {
            Int_t edep_index = -1;
            Double_t distance = 10e10;
            for(size_t ii = 0; ii < edep_t.size(); ii++)
            {
                Double_t temp_distance = sqrt(
                    pow(edep_x[ii] - x, 2.0) + 
                    pow(edep_y[ii] - y, 2.0) + 
                    pow(edep_z[ii] - z, 2.0)
                );
                if(temp_distance < distance) 
                {
                    edep_index = ii;
                    distance = temp_distance;
                }
            }
            return edep_index;
        }
        /**
         * Same as above except for daughters.
        */
        Int_t FindDaughterEnergyDeposition(
            Int_t track_id, Double_t x, Double_t y, Double_t z
        )
        {
            Int_t edep_index = -1;
            Double_t distance = 10e10;
            if(daughter_map.find(track_id) == daughter_map.end()) { 
                std::cout << "ERROR! Daughter: " << track_id << " not found!" << std::endl;
                std::cout << "x: " << x << ", y: " << y << ", z:" << z << std::endl;
                Int_t error;
                std::cin >> error;
                exit(0);
            }
            for(size_t ii = 0; ii < daughter_edep_x.size(); ii++)
            {
                if(daughter_ids[ii] != track_id) { 
                    continue;
                }
                Double_t temp_distance = sqrt(
                    pow(daughter_edep_x[ii] - x, 2.0) + 
                    pow(daughter_edep_y[ii] - y, 2.0) + 
                    pow(daughter_edep_z[ii] - z, 2.0)
                );
                if(temp_distance < distance) 
                {
                    edep_index = ii;
                    distance = temp_distance;
                }
            }
            return edep_index;
        }

        Primary(){}

        Primary(GeneratorLabel label, simb::MCParticle& particle)
        {
            generator_label = label;
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
                if(ii == 0) {
                    process = init_process;
                }
                else if (ii == particle.NumberTrajectoryPoints() - 1) {
                    process = end_process;
                }
                else
                {   
                    for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
                    {
                        if(trajectory_processes[jj].first == ii) {
                            process = trajectory_processes[jj].second;
                            break;
                        }
                    }
                }
                primary_trajectory.AddTrajectoryPoint(
                    particle.T(ii), particle.Vx(ii), particle.Vy(ii), particle.Vz(ii),
                    particle.E(ii), process, volume.volume_name,
                    volume.material_name
                );
            }

        }

        void AddDaughter(simb::MCParticle& particle, Int_t level)
        {
            daughter_ids.emplace_back(particle.TrackId());
            daughter_pdgs.emplace_back(particle.PdgCode());
            daughter_map[particle.TrackId()] = (daughter_ids.size() - 1);
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
                if(ii == 0) {
                    process = particle.Process();
                }
                else if (ii == particle.NumberTrajectoryPoints() - 1) {
                    process = particle.EndProcess();
                }
                else
                {   
                    for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
                    {
                        if(trajectory_processes[jj].first == ii) {
                            process = trajectory_processes[jj].second;
                            break;
                        }
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

        void AddEdep(sim::SimEnergyDeposit& edep)
        {
            // Get the volume information for the energy deposit.
            auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                edep.MidPointX(), edep.MidPointY(), edep.MidPointZ()
            );
            std::string process = FindPrimaryEnergyDepositionProcess(edep);

            edep_energy.emplace_back(edep.Energy());
            edep_process.emplace_back(process);
            edep_volume.emplace_back(volume.volume_name);
            edep_material.emplace_back(volume.material_name);
            edep_t.emplace_back(edep.Time());
            edep_x.emplace_back(edep.MidPointX());
            edep_y.emplace_back(edep.MidPointY());
            edep_z.emplace_back(edep.MidPointZ());
            total_edep_energy += edep.Energy();
        }

        void AddDaughterEdep(sim::SimEnergyDeposit& edep)
        {
            // Get the volume information for the energy deposit.
            auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
                edep.MidPointX(), edep.MidPointY(), edep.MidPointZ()
            );
            std::string process = FindDaughterEnergyDepositionProcess(edep);

            daughter_edep_ids.emplace_back(edep.TrackID());
            daughter_edep_energy.emplace_back(edep.Energy());
            daughter_edep_process.emplace_back(process);
            daughter_edep_volume.emplace_back(volume.volume_name);
            daughter_edep_material.emplace_back(volume.material_name);
            daughter_edep_t.emplace_back(edep.Time());
            daughter_edep_x.emplace_back(edep.MidPointX());
            daughter_edep_y.emplace_back(edep.MidPointY());
            daughter_edep_z.emplace_back(edep.MidPointZ());
            total_daughter_edep_energy += edep.Energy();
        }

        void AddPrimaryDetectorSimulation(
            detinfo::DetectorClocksData const& clockData,
            sim::IDE ide,
            Double_t total_energy,
            Int_t channel, 
            Int_t tick,
            Int_t adc
        )
        {
            det_energy_fraction.emplace_back(ide.energy/total_energy);
            det_energy.emplace_back(ide.energy);
            det_channel.emplace_back(channel);
            det_tick.emplace_back(tick);
            det_adc.emplace_back(adc);
            det_tdc.emplace_back(clockData.TPCTick2TDC(tick));
            Int_t edep_index = FindPrimaryEnergyDeposition(
                ide.x, ide.y, ide.z
            );
            det_edep.emplace_back(edep_index);
            if(edep_index == -1) {
                det_process.emplace_back("no_associated_edep");
            }
            else {
                det_process.emplace_back(edep_process[edep_index]);
            }            
        }

        void AddDaughterDetectorSimulation(
            detinfo::DetectorClocksData const& clockData,
            sim::IDE ide,
            Double_t total_energy,
            Int_t channel,
            Int_t tick,
            Int_t adc
        )
        {
            daughter_det_track_id.emplace_back(ide.trackID);
            daughter_det_energy_fraction.emplace_back(ide.energy/total_energy);
            daughter_det_energy.emplace_back(ide.energy);
            daughter_det_channel.emplace_back(channel);
            daughter_det_tick.emplace_back(tick);
            daughter_det_adc.emplace_back(adc);
            daughter_det_tdc.emplace_back(clockData.TPCTick2TDC(tick));
            Int_t edep_index = FindDaughterEnergyDeposition(
                ide.trackID, ide.x, ide.y, ide.z
            );
            daughter_det_edep.emplace_back(edep_index);
            if(edep_index == -1) {
                daughter_det_process.emplace_back("no_associated_edep");
            }
            else {
                daughter_det_process.emplace_back(daughter_edep_process[edep_index]);
            }
        }
    };
}