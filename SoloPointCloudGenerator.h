/**
 * @file SoloPointCloudGenerator.h
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
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "Logger.h"
#include "ParticleMaps.h"
#include "PrimaryData.h"

namespace arrakis
{
    struct SoloPointCloud
    {
        Int_t event_id = {-1};
        Int_t point_cloud_id = {-1};
        std::vector<Int_t> view = {};
        std::vector<Double_t> wire = {};
        std::vector<Double_t> channel = {};
        std::vector<Double_t> tick = {};
        std::vector<Double_t> tdc = {};
        std::vector<Double_t> adc = {};
        std::vector<Double_t> energy = {};
        Double_t total_energy = {0};
        std::string group_label = {"none"};
        Int_t group_label_id = {-1};
        std::string label = {"none"};
        Int_t label_id = {-1};

        bool all_deposited = false;
        bool all_lar = false;
        bool same_apa = false;
        Double_t lar_edep_fraction = {0};

        SoloPointCloud() {}
    };


    class SoloPointCloudGenerator
    {
    public:
        SoloPointCloudGenerator();
        ~SoloPointCloudGenerator();

        void ProcessEvent(
            ParticleMaps* particle_maps, PrimaryData* primary_data
        );

        void CollectStatistics(
            Primary primary, SoloPointCloud& solo_point_cloud
        );

        void ProcessAr39(
            Primary ar39, ParticleMaps* particle_maps
        );

        void ProcessSingleNeutron(
            Primary neutron, ParticleMaps* particle_maps
        );

        void ProcessPNS(
            Primary neutron, ParticleMaps* particle_maps
        );

        void ProcessLES(
            Primary primary, ParticleMaps* particle_maps
        );

        void ProcessNeutron(
            Primary neutron, ParticleMaps* particle_maps
        );

        void ProcessGamma(
            Primary gamma, ParticleMaps* particle_maps
        );
    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mTTree;

        Int_t mPointCloudID = {0};

        Double_t mEdepEnergyThreshold = {0.01};

        std::map<GeneratorLabel, std::string> mGeneratorLabelNameMap;
        std::map<GeneratorLabel, Int_t> mGeneratorLabelIDMap;

        std::vector<SoloPointCloud> mSoloPointClouds;
        SoloPointCloud mSoloPointCloud;
    };
}