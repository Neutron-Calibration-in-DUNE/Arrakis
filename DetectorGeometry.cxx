/**
 * @file DetectorGeometry.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 */
#include "DetectorGeometry.h"

namespace arrakis 
{
    namespace geometry
    {
        DetectorGeometry* DetectorGeometry::sInstance{nullptr};
        std::mutex DetectorGeometry::sMutex;

        DetectorGeometry *DetectorGeometry::GetInstance()
        {
            std::lock_guard<std::mutex> lock(sMutex);
            if (sInstance == nullptr) {
                sInstance = new DetectorGeometry();
            }
            return sInstance;
        }

        std::string DetectorGeometry::GetTPCName(const size_t i) 
        {
            if (i < sTPCNames.size()) { 
                return sTPCNames[i]; 
            }
            else { 
                return sTPCNames[0]; 
            }
        }
        BoundingBox DetectorGeometry::GetTPCBox(const size_t i) 
        {
            if (i < sTPCBoxes.size()) { 
                return sTPCBoxes[i]; 
            }
            else { 
                return sTPCBoxes[0]; 
            }
        }
        BoundingBox DetectorGeometry::GetActiveTPCBox(const size_t i) 
        {
            if (i < sActiveTPCBoxes.size()) { 
                return sActiveTPCBoxes[i]; 
            }
            else { 
                return sActiveTPCBoxes[0]; 
            }
        }
        double DetectorGeometry::GetTPCMass(const size_t i) 
        {
            if (i < sTPCMasses.size()) { 
                return sTPCMasses[i]; 
            }
            else { 
                return sTPCMasses[0]; 
            }
        }
        double DetectorGeometry::GetTPCDriftDistance(const size_t i) 
        {
            if (i < sTPCDriftDistances.size()) { 
                return sTPCDriftDistances[i]; 
            }
            else { 
                return sTPCDriftDistances[0]; 
            }
        }

        DetectorGeometry::DetectorGeometry()
        {
            Logger::GetInstance("geometry")->trace(
                "setting up detector geometry info."
            );
            // set up the geometry interface
            sGeometryCore = lar::providerFrom<geo::Geometry>();

            // initialize TTrees
            sGeometryTree = sTFileService->make<TTree>("geometry", "geometry");
            sGeometryTree->Branch("world_name", &sWorldName);
            sGeometryTree->Branch("world_box_ranges", &(sWorldBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            sGeometryTree->Branch("detector_name", &sDetectorName);
            sGeometryTree->Branch("detector_box_ranges", &(sDetectorBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            sGeometryTree->Branch("cryostat_name", &sCryostatName);
            sGeometryTree->Branch("cryostat_box_ranges", &(sCryostatBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            sGeometryTree->Branch("number_of_tpcs", &sNumberOfTPCs);
            sGeometryTree->Branch("tpc_names", &sTPCNames);
            for (int i = 0; i < sNumberOfTPCs; i++) 
            {
                sGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_name").c_str(), &(sTPCNames[i]));
                sGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_box_ranges").c_str(), &(sTPCBoxes[i]), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
                sGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_mass").c_str(), &(sTPCMasses[i]));
                sGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_drift_distance").c_str(), &(sTPCDriftDistances[i]));
            }
            sGeometryTree->Branch("tpc_masses", &sTPCMasses);
            sGeometryTree->Branch("tpc_drift_distances", &sTPCDriftDistances);
            sGeometryTree->Branch("total_tpc_box_ranges", &(sTotalTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            sGeometryTree->Branch("total_active_tpc_box_ranges", &(sTotalActiveTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            sGeometryTree->Branch("total_tpc_mass", &sTotalTPCMass);

            // get detector clock data
            auto const clock_data = 
                art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
            sTriggerOffset = trigger_offset(clock_data);

            std::set<geo::View_t> const views = sGeometryCore->Views();
            for(auto view : views) {
                sWirePitchMap[view] = sGeometryCore->WirePitch(view);
            }

            // Accquiring geometry data
            sNumberOfAPAs=sGeometryCore->NTPC()*sGeometryCore->Ncryostats()/2; //No. of APAs
            sNumberOfChannelsPerAPA = sGeometryCore->Nchannels()/sNumberOfAPAs; //No. of channels per APA

            // taken from dune35t module a way to organise the channel mapping:
            // loop through channels in the first APA to find the channel boundaries for each view
            // will adjust for desired APA after
            sUChannelMin = 0;
            sZChannelMax = sNumberOfChannelsPerAPA - 1;
            for (unsigned int c = sUChannelMin + 1; c < sZChannelMax; c++ ){
                if ( sGeometryCore->View(c) == geo::kV && sGeometryCore->View(c-1) == geo::kU ){
                    sVChannelMin = c;
                    sUChannelMax = c - 1;
                }
                if ( sGeometryCore->View(c) == geo::kZ && sGeometryCore->View(c-1) == geo::kV ){
                    sZChannelMin = c;
                    sVChannelMax = c-1;
                }
            }

            //Number of channels in each view
            sNumberOfUChannels = sUChannelMax - sUChannelMin+1; //U
            sNumberOfVChannels = sVChannelMax - sVChannelMin+1; //V
            sNumberOfZChannels = sZChannelMax - sZChannelMin+1; //Z (collection plane)

            for(auto channel : sGeometryCore->ChannelsInTPCs())
            {
                sChannelToWireIDMap[channel] = sGeometryCore->ChannelToWire(channel);
            }

            // collect world info
            sWorldName = sGeometryCore->GetWorldVolumeName();
            sWorldBox.setBox(sGeometryCore->WorldBox());

            // create name-volumetype map for world
            sMaterialPOI.SetCoordinates(sWorldBox.x_min,sWorldBox.y_min,sWorldBox.z_min);
            std::string volumeName = sGeometryCore->VolumeName(sMaterialPOI);
            sVolumeTypeMap[volumeName] = VolumeType::World;

            // collect detector info
            sDetectorName = sGeometryCore->DetectorName();
            sDetectorBox.setBox(
                -sGeometryCore->DetHalfWidth(), sGeometryCore->DetHalfWidth(),
                -sGeometryCore->DetHalfHeight(), sGeometryCore->DetHalfHeight(),
                0, sGeometryCore->DetLength()
            );

            // collect cryostat info
            // for now, assuming analysis is done over a single cryostat
            geo::CryostatGeo const& Cryo = sGeometryCore->Cryostat();
            sCryostatName = std::string(Cryo.ID());
            sCryostatBox.setBox(Cryo.Boundaries());

            // create name-volumetype map for cryostat
            sMaterialPOI.SetCoordinates(sCryostatBox.x_min, sCryostatBox.y_min, sCryostatBox.z_min);
            volumeName = sGeometryCore->VolumeName(sMaterialPOI);
            sVolumeTypeMap[volumeName] = VolumeType::Cryostat;

            // iterate over all TPCs
            sNumberOfTPCs  = sGeometryCore->TotalNTPC();
            for (geo::TPCGeo const& TPC : sGeometryCore->IterateTPCs())
            {
                sTPCNames.emplace_back(TPC.ID());
                sTPCBoxes.emplace_back(BoundingBox(TPC.BoundingBox()));
                sActiveTPCBoxes.emplace_back(BoundingBox(TPC.ActiveBoundingBox()));
                sTPCMasses.emplace_back(TPC.ActiveMass());
                sTPCDriftDistances.emplace_back(TPC.DriftDistance());
                // create name-volumetype map for this tpc
                sVolumeTypeMap[sGeometryCore->VolumeName(TPC.GetCenter())] = VolumeType::TPC;
            }
            // find the total TPC and total Active TPC volumes
            FindTotalTPCBoxes();
            sTotalTPCMass = sGeometryCore->TotalMass();    
        }

        // get volume information for a point
        DetectorVolume DetectorGeometry::GetVolume(std::vector<double> position)
        {
            return GetVolume(position[0], position[1], position[2]);
        }

        // get volume information for a point
        DetectorVolume DetectorGeometry::GetVolume(double x, double y, double z)
        {
            sMaterialPOI.SetCoordinates(x,y,z);

            // get the volume information
            std::string volumeName = sGeometryCore->VolumeName(sMaterialPOI);
            VolumeType volumeType = sVolumeTypeMap[volumeName];

            // get the current material information
            sMaterial = sGeometryService->Material(sMaterialPOI);
            double material = sMaterial->GetZ();
            std::string materialName = sMaterial->GetName();

            // return the constructed volume 
            return DetectorVolume(volumeType, volumeName, materialName, material);
        }
        
        // get total tpc volume information
        void DetectorGeometry::FindTotalTPCBoxes()
        {
            double x_min = 0; double x_max = 0;
            double y_min = 0; double y_max = 0;
            double z_min = 0; double z_max = 0;
            for (size_t i = 0; i < sTPCBoxes.size(); i++) 
            {
                if (sTPCBoxes[i].x_min < x_min) x_min = sTPCBoxes[i].x_min;
                if (sTPCBoxes[i].x_max > x_max) x_max = sTPCBoxes[i].x_max;
                if (sTPCBoxes[i].y_min < y_min) y_min = sTPCBoxes[i].y_min;
                if (sTPCBoxes[i].y_max > y_max) y_max = sTPCBoxes[i].y_max;
                if (sTPCBoxes[i].z_min < z_min) z_min = sTPCBoxes[i].z_min;
                if (sTPCBoxes[i].z_max > z_max) z_max = sTPCBoxes[i].z_max;
            }
            sTotalTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
            x_min = 0; x_max = 0;
            y_min = 0; y_max = 0;
            z_min = 0; z_max = 0;
            for (size_t i = 0; i < sActiveTPCBoxes.size(); i++) 
            {
                if (sActiveTPCBoxes[i].x_min < x_min) x_min = sActiveTPCBoxes[i].x_min;
                if (sActiveTPCBoxes[i].x_max > x_max) x_max = sActiveTPCBoxes[i].x_max;
                if (sActiveTPCBoxes[i].y_min < y_min) y_min = sActiveTPCBoxes[i].y_min;
                if (sActiveTPCBoxes[i].y_max > y_max) y_max = sActiveTPCBoxes[i].y_max;
                if (sActiveTPCBoxes[i].z_min < z_min) z_min = sActiveTPCBoxes[i].z_min;
                if (sActiveTPCBoxes[i].z_max > z_max) z_max = sActiveTPCBoxes[i].z_max;
            }
            sTotalActiveTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
        }
        void DetectorGeometry::FillTTree()
        {
            Logger::GetInstance("geometry")->trace(
                "saving geometry info to root file."
            );
            // add geometry info
            sGeometryTree->Fill();
        }
    }
}