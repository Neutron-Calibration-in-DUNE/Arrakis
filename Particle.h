/**
 * @file Node.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-24
 */
#pragma once
#include <memory>

namespace arrakis
{
    namespace mcdata
    {
        enum class ParticleType
        {
            None = -1,
            Primary = 0,
            Daughter = 1,
        };

        enum class TrajectoryPointType
        {
            NotDefined = -1,
            Unknown = 0,
            
        };

        const std::map<TrajectoryPointType, std::string> TrajectoryPointTypeToString
        {
            {TrajectoryPointType::NotDefined, "not_defined"},
            {TrajectoryPointType::Unknown, "Unknown"},
        };
        const std::map<std::string, TrajectoryPointType> StringToTrajectoryPointType
        {
            {"not_defined", TrajectoryPointType::NotDefined},
            {"Unknown", TrajectoryPointType::Unknown},
        };

        struct Particle
        {
            Int_t track_id = {0};
            std::vector<TrajectoryPointType> trajectory_type = {};
            
        };
    }
}