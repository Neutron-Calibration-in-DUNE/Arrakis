/**
 * @file DetectorGeometry.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 */
#pragma once
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

#include <string>
#include <vector>
#include <memory>

#include "Logger.h"

namespace arrakis 
{
    // list of materials in the detector
    enum MaterialList 
    {

    };

    enum VolumeType 
    {
        World,
        Cryostat,
        TPC,
    };

    // struct for detector volume information
    struct DetectorVolume
    {
        VolumeType volume_type;
        std::string volume_name;
        std::string material_name;
        double material;

        DetectorVolume() {}
        
        DetectorVolume(
            VolumeType volumeType, 
            std::string volumeName, 
            std::string materialName, 
            double material
        )
        : volume_type(volumeType)
        , volume_name(volumeName)
        , material_name(materialName)
        , material(material)
        {}
    };
    // struct for bounding boxes
    struct BoundingBox
    {
        double x_min = 0; double x_max = 0;
        double y_min = 0; double y_max = 0;
        double z_min = 0; double z_max = 0;

        double width()  { return x_max - x_min; }
        double height() { return y_max - y_min; }
        double length() { return z_max - z_min; }

        void setBox(geo::BoxBoundedGeo const& Box) {
            x_min = Box.MinX(); x_max = Box.MaxX();
            y_min = Box.MinY(); y_max = Box.MaxY();
            z_min = Box.MinZ(); z_max = Box.MaxZ();
        }
        void setBox(double xmin, double xmax,
                    double ymin, double ymax,
                    double zmin, double zmax)
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
        }
        BoundingBox() {}
        BoundingBox(double xs[2], double ys[2], double zs[2])
        {
            x_min = xs[0]; x_max = xs[1];
            y_min = ys[0]; y_max = ys[1];
            z_min = zs[0]; z_max = zs[1];
        }
        BoundingBox(double vals[6])
        {
            x_min = vals[0]; x_max = vals[1];
            y_min = vals[2]; y_max = vals[3];
            z_min = vals[4]; z_max = vals[4];
        }
        BoundingBox(double xmin, double xmax,
                    double ymin, double ymax,
                    double zmin, double zmax)
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
        }
        BoundingBox(geo::BoxBoundedGeo const& Box) {
            x_min = Box.MinX(); x_max = Box.MaxX();
            y_min = Box.MinY(); y_max = Box.MaxY();
            z_min = Box.MinZ(); z_max = Box.MaxZ();
        }
    };

    /**
     * @brief Singleton class for storing detector geometry information.
     * 
     */
    class DetectorGeometry
    {
    private:
        static DetectorGeometry * sInstance;
        static std::mutex sMutex;

    protected:
        DetectorGeometry(const std::string name);
        ~DetectorGeometry() {}
        std::string sName;

    public:
        // this singleton cannot be cloned
        DetectorGeometry(DetectorGeometry &other) = delete;
        // singleton should also not be assignable
        void operator=(const DetectorGeometry &) = delete;

        // static method that controls access to the singleton
        // instance
        static DetectorGeometry* GetInstance(const std::string& name);

        std::string Name() const {
            return sName;
        }

        unsigned int NumberOfUChannels() { return mNumberOfUChannels; }
        unsigned int NumberOfVChannels() { return mNumberOfVChannels; }
        unsigned int NumberOfZChannels() { return mNumberOfZChannels; }

        // find channel boundaries for each view
        unsigned int UChannelMin() { return mUChannelMin; }
        unsigned int UChannelMax() { return mUChannelMax; }
        unsigned int VChannelMin() { return mVChannelMin; }
        unsigned int VChannelMax() { return mVChannelMax; }
        unsigned int ZChannelMin() { return mZChannelMin; }
        unsigned int ZChannelMax() { return mZChannelMax; }

        unsigned int NumberOfAPAs() { return mNumberOfAPAs; }
        unsigned int NumberOfChannelsPerAPA() { return mNumberOfChannelsPerAPA; }

        geo::View_t View(geo::PlaneID const& pid) { return mGeometryCore->View(pid); }

        // getters
        std::string GetWorldName()      { return mWorldName; }
        BoundingBox GetWorldBox()       { return mWorldBox; }

        std::string GetDetectorName()   { return mDetectorName; }
        BoundingBox GetDetectorBox()    { return mDetectorBox; }

        std::string GetCryostatName()   { return mCryostatName; }
        BoundingBox GetCryostatBox()    { return mCryostatBox; }

        int GetNumberOfTPCs()                   { return mNumberOfTPCs; }
        std::vector<std::string> GetTPCNames()  { return mTPCNames; }
        BoundingBox GetTotalTPCBox()            { return mTotalTPCBox; }
        BoundingBox GetTotalActiveTPCBox()      { return mTotalActiveTPCBox; }
        double GetTotalTPCMass()                { return mTotalTPCMass; }
        std::vector<double> GetTPCMasses()      { return mTPCMasses; }
        std::vector<double> GetTPCDriftDistances()  { return mTPCDriftDistances; }

        std::string GetTPCName(const size_t i);
        BoundingBox GetTPCBox(const size_t i);
        BoundingBox GetActiveTPCBox(const size_t i);
        double GetTPCMass(const size_t i);
        double GetTPCDriftDistance(const size_t i);
        
        // get volume information for a point
        DetectorVolume GetVolume(std::vector<double> position);
        DetectorVolume GetVolume(double x, double y, double z);

        // function for finding total tpc volumes
        void FindTotalTPCBoxes();
        // fill the geometry ttree
        void FillTTree();
        
    private:
        art::ServiceHandle<geo::Geometry> mGeometryService;
        geo::GeometryCore const* mGeometryCore;
        art::ServiceHandle<art::TFileService> mTFileService;

        // TPC // Number of channels in each planes
        unsigned int mNumberOfUChannels;
        unsigned int mNumberOfVChannels;
        unsigned int mNumberOfZChannels;

        // find channel boundaries for each view
        unsigned int mUChannelMin;
        unsigned int mUChannelMax;
        unsigned int mVChannelMin;
        unsigned int mVChannelMax;
        unsigned int mZChannelMin;
        unsigned int mZChannelMax;

        unsigned int mNumberOfAPAs; //Number of APAs
        unsigned int mNumberOfChannelsPerAPA; //Number of channels in each APA

        TTree *mGeometryTree;
        size_t mTriggerOffset;
        
        std::map<std::string,VolumeType> mVolumeTypeMap;

        // world volume
        std::string mWorldName;
        BoundingBox mWorldBox;
        // detector volume
        std::string mDetectorName;
        BoundingBox mDetectorBox;
        // cryostat volume
        std::string mCryostatName;
        BoundingBox mCryostatBox;

        // tpc volumes
        int mNumberOfTPCs;
        std::vector<std::string> mTPCNames;
        std::vector<BoundingBox> mTPCBoxes;
        std::vector<BoundingBox> mActiveTPCBoxes;
        std::vector<double> mTPCMasses;
        std::vector<double> mTPCDriftDistances;

        // full tpc volume
        BoundingBox mTotalTPCBox;
        BoundingBox mTotalActiveTPCBox;
        double mTotalTPCMass;
        ////////////////////////////////////////////////
        // detector material variables
        ////////////////////////////////////////////////
        // we will need to ask Geant4 about material 
        // properties for the detector volume
        // at each point of interest.  This requires holding 
        // this information in a
        // TGeoMaterial object, which is part of ROOT.
        const TGeoMaterial *mMaterial;
        geo::Point_t mMaterialPOI;
    };
}