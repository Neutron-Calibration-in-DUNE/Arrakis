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
        struct Particle
        {
            Int_t track_id = {0};

        };
    }
}