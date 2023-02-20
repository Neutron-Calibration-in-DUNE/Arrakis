/**
 * @file Primary.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-20
 */
#include "Primary.h"

namespace arrakis
{
    Primary(GeneratorLabel label, simb::MCParticle& particle)
    : Node(0, NodeType::Primary)
    {
    //     generator_label = label;
    //     track_id = particle.TrackId();
    //     pdg = particle.PdgCode();
    //     init_process = particle.Process();
    //     init_energy = particle.E();
    //     init_t = particle.T();
    //     init_x = particle.Vx();
    //     init_y = particle.Vy();
    //     init_z = particle.Vz();
    //     end_process = particle.EndProcess();
    //     end_energy = particle.EndE();
    //     end_x = particle.EndX();
    //     end_y = particle.EndY();
    //     end_z = particle.EndZ();

    //     simb::MCTrajectory trajectory = particle.Trajectory();
    //     auto trajectory_processes = trajectory.TrajectoryProcesses();

    //     for(size_t ii = 0; ii < particle.NumberTrajectoryPoints(); ii++)
    //     {
    //         // Get the volume information for trajectory point.
    //         auto volume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
    //             particle.Vx(ii), particle.Vy(ii), particle.Vz(ii)
    //         );
    //         std::string process = "not_defined";
    //         if(ii == 0) {
    //             process = init_process;
    //         }
    //         else if (ii == particle.NumberTrajectoryPoints() - 1) {
    //             process = end_process;
    //         }
    //         else
    //         {   
    //             for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
    //             {
    //                 if(trajectory_processes[jj].first == ii) {
    //                     process = trajectory_processes[jj].second;
    //                 }
    //             }
    //         }
    //         primary_trajectory.AddTrajectoryPoint(
    //             particle.T(ii), particle.Vx(ii), particle.Vy(ii), particle.Vz(ii),
    //             particle.E(ii), process, volume.volume_name,
    //             volume.material_name
    //         );
    //     }
    // }
    }
}