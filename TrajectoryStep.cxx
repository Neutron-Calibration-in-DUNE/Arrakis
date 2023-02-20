/**
 * @file TrajectoryStep.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-20
 */
#include "TrajectoryStep.h"

namespace arrakis
{
    TrajectoryStep::TrajectoryStep()
    : Node(0, NodeType::TrajectoryStep)
    {
    }

    TrajectoryStep::~TrajectoryStep()
    {
    }

    TrajectoryStep::TrajectoryStep(
        size_t step, simb::MCParticle& particle, 
        simb::MCTrajectory& trajectory, 
        simb::MCTrajectory::ProcessMap& trajectory_processes
    )
    : Node(step, NodeType::TrajectoryStep)
    {
        // Get the volume information for trajectory point.
        mDetectorVolume = DetectorGeometry::GetInstance("PrimaryData")->GetVolume(
            particle.Vx(step), particle.Vy(step), particle.Vz(step)
        );
        mT = particle.T(step);
        mX = particle.Vx(step);
        mY = particle.Vy(step);
        mZ = particle.Vz(step);

        for(size_t jj = 0; jj < trajectory_processes.size(); jj++)
        {
            if(trajectory_processes[jj].first == step) {
                mProcess = trajectory_processes[jj].second;
            }
        }
        std::cout << "particle: " << particle.PdgCode() << "," << particle.TrackId() << std::endl;
        std::cout << "trajectory step: " << step << "," << mProcess << "," << mT << "," << mX << "," << mY << "," << mZ << std::endl;
        if(step < particle.NumberTrajectoryPoints() -1)
        {
            mNextTrajectoryStep = new TrajectoryStep(
                step + 1,
                particle, 
                trajectory,
                trajectory_processes
            );
        }
    }
}