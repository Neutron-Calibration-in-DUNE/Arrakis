/**
 * @file PrimaryData.h
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
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "ParticleMaps.h"

namespace arrakis
{
    struct Primary
    {
        Int_t track_id = {0};
        Int_t pdg = {0};

        std::string init_process = {""};
        Double_t init_energy = {0};
        Double_t init_x = {0};
        Double_t init_y = {0};
        Double_t init_z = {0};
        
        std::string end_process = {""};
        Double_t end_energy = {0};
        Double_t end_x = {0};
        Double_t end_y = {0};
        Double_t end_z = {0};

        std::vector<Double_t> edep_energy = {};
        std::vector<Double_t> edep_x = {};
        std::vector<Double_t> edep_y = {};
        std::vector<Double_t> edep_z = {};

        std::vector<Int_t> daughter_ids = {};
        std::vector<Int_t> daughter_level = {};

        std::vector<std::string> daughter_init_process = {};
        std::vector<Double_t> daughter_init_energy = {};
        std::vector<Double_t> daughter_init_x = {};
        std::vector<Double_t> daughter_init_y = {};
        std::vector<Double_t> daughter_init_z = {};

        std::vector<std::string> daughter_end_process = {};
        std::vector<Double_t> daughter_end_energy = {};
        std::vector<Double_t> daughter_end_x = {};
        std::vector<Double_t> daughter_end_y = {};
        std::vector<Double_t> daughter_end_z = {};

        std::vector<Int_t> daughter_edep_ids = {};
        std::vector<Double_t> daughter_edep_energy = {};
        std::vector<Double_t> daughter_edep_x = {};
        std::vector<Double_t> daughter_edep_y = {};
        std::vector<Double_t> daughter_edep_z = {};

        Primary(){}
        Primary(
            Int_t _track_id, Int_t _pdg, std::string _init_process,
            Double_t _init_energy, Double_t _init_x, Double_t _init_y,
            Double_t _init_z, std::string _end_process, Double_t _end_energy,
            Double_t _end_x, Double_t _end_y, Double_t _end_z
        )
        : track_id(_track_id)
        , pdg(_pdg)
        , init_process(_init_process)
        , init_energy(_init_energy)
        , init_x(_init_x)
        , init_y(_init_y)
        , init_z(_init_z)
        , end_process(_end_process)
        , end_energy(_end_energy)
        , end_x(_end_x)
        , end_y(_end_y)
        , end_z(_end_z)
        {}

        void AddEdep(
            Double_t energy, 
            Double_t x, Double_t y, Double_t z
        )
        {
            edep_energy.emplace_back(energy);
            edep_x.emplace_back(x);
            edep_y.emplace_back(y);
            edep_z.emplace_back(z);
        }

        void AddDaughter(
            Int_t track_id, Int_t level,
            std::string init_process, Double_t init_energy,
            Double_t init_x, Double_t init_y, Double_t init_z,
            std::string end_process, Double_t end_energy,
            Double_t end_x, Double_t end_y, Double_t end_z
        )
        {
            daughter_ids.emplace_back(track_id);
            daughter_level.emplace_back(level);
            daughter_init_process.emplace_back(init_process);
            daughter_init_energy.emplace_back(init_energy);
            daughter_init_x.emplace_back(init_x);
            daughter_init_y.emplace_back(init_y);
            daughter_init_z.emplace_back(init_z);
            daughter_end_process.emplace_back(end_process);
            daughter_end_energy.emplace_back(end_energy);
            daughter_end_x.emplace_back(end_x);
            daughter_end_y.emplace_back(end_y);
            daughter_end_z.emplace_back(end_z);
        }

        void AddDaughterEdep(
            Int_t track_id, Double_t energy, 
            Double_t x, Double_t y, Double_t z
        )
        {
            daughter_edep_ids.emplace_back(track_id);
            daughter_edep_energy.emplace_back(energy);
            daughter_edep_x.emplace_back(x);
            daughter_edep_y.emplace_back(y);
            daughter_edep_z.emplace_back(z);
        }
    };

    class PrimaryData
    {
    public:
        PrimaryData(
            bool SavePrimaryData, bool SavePrimaryDataEdeps,
            bool SavePrimaryDataSimChannel, bool SavePrimaryDataRawTPC
        );
        ~PrimaryData();

        void ResetEvent();
        void ProcessEvent(
            ParticleMaps particle_maps,
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<raw::RawDigit>>& rawTPC
        );

        Int_t FindPrimary(Int_t track_id);

    private:
        bool mSavePrimaryData = {false};
        bool mSavePrimaryDataEdeps = {false};
        bool mSavePrimaryDataSimChannel = {false};
        bool mSavePrimaryDataRawTPC = {false};

        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mTTree;

        std::vector<Primary> mPrimaries;
        Primary mPrimary;

    };
}