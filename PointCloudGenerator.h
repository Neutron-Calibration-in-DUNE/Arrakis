/**
 * @file PointCloudGenerator.h
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
#include "DetectorPointCloud.h"
#include "ParticleMaps.h"
#include "PrimaryData.h"
#include "SoloPointCloud.h"

namespace arrakis
{
    class PointCloudGenerator
    {
    public:
        PointCloudGenerator();
        ~PointCloudGenerator();

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

        void ProcessJunk(
            Junk junk
        );
    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mSoloPointCloudTTree;
        TTree *mDetectorPointCloudTTree;

        Int_t mPointCloudID = {0};

        Double_t mEdepEnergyThreshold = {0.01};

        std::map<GeneratorLabel, std::string> mGeneratorLabelNameMap;
        std::map<GeneratorLabel, Int_t> mGeneratorLabelIDMap;

        std::vector<SoloPointCloud> mSoloPointClouds;
        DetectorPointCloud mDetectorPointCloud;
        SoloPointCloud mSoloPointCloud;
    };
}