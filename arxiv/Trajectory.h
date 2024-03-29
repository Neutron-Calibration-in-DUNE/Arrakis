/**
 * @file Trajectory.h
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

namespace arrakis
{
    struct Trajectory
    {
        std::vector<Double_t> t = {};
        std::vector<Double_t> x = {};
        std::vector<Double_t> y = {};
        std::vector<Double_t> z = {};
        std::vector<Double_t> energy = {};
        std::vector<std::string> process = {};
        std::vector<std::string> volume = {};
        std::vector<std::string> material = {};

        Trajectory()
        {
        }

        Trajectory(
            Double_t _t, Double_t _x, Double_t _y, Double_t _z,
            Double_t _energy, std::string _process, std::string _volume, 
            std::string _material
        )
        {
            t.emplace_back(_t);
            x.emplace_back(_x);
            y.emplace_back(_y);
            z.emplace_back(_z);
            energy.emplace_back(_energy);
            process.emplace_back(_process);
            volume.emplace_back(_volume);
            material.emplace_back(_material);
        }

        void AddTrajectoryPoint(
            Double_t _t, Double_t _x, Double_t _y, Double_t _z,
            Double_t _energy, std::string _process, std::string _volume, 
            std::string _material
        )
        {
            t.emplace_back(_t);
            x.emplace_back(_x);
            y.emplace_back(_y);
            z.emplace_back(_z);
            energy.emplace_back(_energy);
            process.emplace_back(_process);
            volume.emplace_back(_volume);
            material.emplace_back(_material);
        }

        void PrintTrajectory()
        {
            std::cout << "\n ----- Trajectory ----- \n";
            std::cout << "(t,x,y,z) - energy - process - volume - material\n";
            for(size_t ii = 0; ii < t.size(); ii++)
            {
                std::cout << "(" << t[ii] << ", " << x[ii];
                std::cout << ", " << y[ii] << ", " << z[ii];
                std::cout << ") - " << energy[ii] << " - ";
                std::cout << process[ii] << " - " << volume[ii];
                std::cout << " - " << material[ii] << "\n";
            }
            std::cout << " -----    End     -----" << std::endl;
        }
    };
}