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
    Primary::Primary()
    : Node(0, NodeType::Primary)
    {
    }
    Primary::~Primary()
    {
    }

    Primary::Primary(GeneratorLabel label, simb::MCParticle& particle)
    : Node(0, NodeType::Primary)
    {
        mGeneratorLabel = label;
        mTrackID = particle.TrackId();
        mPdgCode = particle.PdgCode();
        mInitProcess = particle.Process();
        mInitEnergy = particle.E();
        mInitT = particle.T();
        mInitX = particle.Vx();
        mInitY = particle.Vy();
        mInitZ = particle.Vz();
        mEndProcess = particle.EndProcess();
        mEndEnergy = particle.EndE();
        mEndX = particle.EndX();
        mEndY = particle.EndY();
        mEndZ = particle.EndZ();

        // check for trajectory points
        simb::MCTrajectory trajectory = particle.Trajectory();
        auto trajectory_processes = trajectory.TrajectoryProcesses();
        if(particle.NumberTrajectoryPoints() > 0)
        {
            mTrajectoryStep = new TrajectoryStep(
                1,
                particle,
                trajectory,
                trajectory_processes
            );
        }

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