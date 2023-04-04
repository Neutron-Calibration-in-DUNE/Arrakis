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

        // setting boundaries
        void setBox(geo::BoxBoundedGeo const& Box) 
        {
            x_min = Box.MinX(); x_max = Box.MaxX();
            y_min = Box.MinY(); y_max = Box.MaxY();
            z_min = Box.MinZ(); z_max = Box.MaxZ();
        }
        void setBox(
            double xmin, double xmax,
            double ymin, double ymax,
            double zmin, double zmax
        )
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
        }

        // constructors
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
        BoundingBox(
            double xmin, double xmax,
            double ymin, double ymax,
            double zmin, double zmax
        )
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
        }
        BoundingBox(geo::BoxBoundedGeo const& Box) 
        {
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
    public:
        // this singleton cannot be cloned
        DetectorGeometry(DetectorGeometry &other) = delete;
        // singleton should also not be assignable
        void operator=(const DetectorGeometry &) = delete;

        // static method that controls access to 
        // the singleton instance
        static DetectorGeometry* GetInstance();

        unsigned int NumberOfUChannels() { return sNumberOfUChannels; }
        unsigned int NumberOfVChannels() { return sNumberOfVChannels; }
        unsigned int NumberOfZChannels() { return sNumberOfZChannels; }

        // find channel boundaries for each view
        unsigned int UChannelMin() { return sUChannelMin; }
        unsigned int UChannelMax() { return sUChannelMax; }
        unsigned int VChannelMin() { return sVChannelMin; }
        unsigned int VChannelMax() { return sVChannelMax; }
        unsigned int ZChannelMin() { return sZChannelMin; }
        unsigned int ZChannelMax() { return sZChannelMax; }

        unsigned int NumberOfAPAs() { return sNumberOfAPAs; }
        unsigned int NumberOfChannelsPerAPA() { return sNumberOfChannelsPerAPA; }

        geo::View_t View(raw::ChannelID_t const channel) { return sGeometryCore->View(channel); }
        std::vector<geo::WireID> ChannelToWire(raw::ChannelID_t const channel) { return sChannelToWireIDMap[channel]; }
        Double_t GetWirePitch(geo::View_t view) { return sWirePitchMap[view]; }

        // getters
        std::string GetWorldName()      { return sWorldName; }
        BoundingBox GetWorldBox()       { return sWorldBox; }

        std::string GetDetectorName()   { return sDetectorName; }
        BoundingBox GetDetectorBox()    { return sDetectorBox; }

        std::string GetCryostatName()   { return sCryostatName; }
        BoundingBox GetCryostatBox()    { return sCryostatBox; }

        int GetNumberOfTPCs()                   { return sNumberOfTPCs; }
        std::vector<std::string> GetTPCNames()  { return sTPCNames; }
        BoundingBox GetTotalTPCBox()            { return sTotalTPCBox; }
        BoundingBox GetTotalActiveTPCBox()      { return sTotalActiveTPCBox; }
        double GetTotalTPCMass()                { return sTotalTPCMass; }
        std::vector<double> GetTPCMasses()      { return sTPCMasses; }
        std::vector<double> GetTPCDriftDistances()  { return sTPCDriftDistances; }

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
    
    protected:
        DetectorGeometry();
        ~DetectorGeometry() {}
        
    private:
        static DetectorGeometry * sInstance;
        static std::mutex sMutex;

        art::ServiceHandle<geo::Geometry> sGeometryService;
        geo::GeometryCore const* sGeometryCore;
        art::ServiceHandle<art::TFileService> sTFileService;

        // TPC // Number of channels in each planes
        unsigned int sNumberOfUChannels;
        unsigned int sNumberOfVChannels;
        unsigned int sNumberOfZChannels;

        // find channel boundaries for each view
        unsigned int sUChannelMin;
        unsigned int sUChannelMax;
        unsigned int sVChannelMin;
        unsigned int sVChannelMax;
        unsigned int sZChannelMin;
        unsigned int sZChannelMax;

        unsigned int sNumberOfAPAs; //Number of APAs
        unsigned int sNumberOfChannelsPerAPA; //Number of channels in each APA

        std::map<raw::ChannelID_t, std::vector<geo::WireID>> sChannelToWireIDMap;
        std::map<geo::View_t, geo::Length_t> sWirePitchMap;

        TTree *sGeometryTree;
        size_t sTriggerOffset;
        
        std::map<std::string,VolumeType> sVolumeTypeMap;

        // world volume
        std::string sWorldName;
        BoundingBox sWorldBox;
        // detector volume
        std::string sDetectorName;
        BoundingBox sDetectorBox;
        // cryostat volume
        std::string sCryostatName;
        BoundingBox sCryostatBox;

        // tpc volumes
        int sNumberOfTPCs;
        std::vector<std::string> sTPCNames;
        std::vector<BoundingBox> sTPCBoxes;
        std::vector<BoundingBox> sActiveTPCBoxes;
        std::vector<double> sTPCMasses;
        std::vector<double> sTPCDriftDistances;

        // full tpc volume
        BoundingBox sTotalTPCBox;
        BoundingBox sTotalActiveTPCBox;
        double sTotalTPCMass;

        ////////////////////////////////////////////////
        // detector material variables
        ////////////////////////////////////////////////
        // we will need to ask Geant4 about material 
        // properties for the detector volume
        // at each point of interest.  This requires holding 
        // this information in a
        // TGeoMaterial object, which is part of ROOT.
        const TGeoMaterial *sMaterial;
        geo::Point_t sMaterialPOI;
    };
}