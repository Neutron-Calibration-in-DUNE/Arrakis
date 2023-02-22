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
            if (i < mTPCNames.size()) { 
                return mTPCNames[i]; 
            }
            else { 
                return mTPCNames[0]; 
            }
        }
        BoundingBox DetectorGeometry::GetTPCBox(const size_t i) 
        {
            if (i < mTPCBoxes.size()) { 
                return mTPCBoxes[i]; 
            }
            else { 
                return mTPCBoxes[0]; 
            }
        }
        BoundingBox DetectorGeometry::GetActiveTPCBox(const size_t i) 
        {
            if (i < mActiveTPCBoxes.size()) { 
                return mActiveTPCBoxes[i]; 
            }
            else { 
                return mActiveTPCBoxes[0]; 
            }
        }
        double DetectorGeometry::GetTPCMass(const size_t i) 
        {
            if (i < mTPCMasses.size()) { 
                return mTPCMasses[i]; 
            }
            else { 
                return mTPCMasses[0]; 
            }
        }
        double DetectorGeometry::GetTPCDriftDistance(const size_t i) 
        {
            if (i < mTPCDriftDistances.size()) { 
                return mTPCDriftDistances[i]; 
            }
            else { 
                return mTPCDriftDistances[0]; 
            }
        }

        DetectorGeometry::DetectorGeometry(const std::string name)
        : sName(name)
        {
            // set up the geometry interface
            mGeometryCore = lar::providerFrom<geo::Geometry>();
            // initialize TTrees
            mGeometryTree = mTFileService->make<TTree>("geometry", "geometry");
            // get detector clock data
            auto const clock_data = 
                art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
            mTriggerOffset = trigger_offset(clock_data);

            std::set<geo::View_t> const views = mGeometryCore->Views();
            for(auto view : views) {
                mWirePitchMap[view] = mGeometryCore->WirePitch(view);
            }

            // Accquiring geometry data
            mNumberOfAPAs=mGeometryCore->NTPC()*mGeometryCore->Ncryostats()/2; //No. of APAs
            mNumberOfChannelsPerAPA = mGeometryCore->Nchannels()/mNumberOfAPAs; //No. of channels per APA

            // taken from dune35t module a way to organise the channel mapping:
            // loop through channels in the first APA to find the channel boundaries for each view
            // will adjust for desired APA after
            mUChannelMin = 0;
            mZChannelMax = mNumberOfChannelsPerAPA - 1;
            for (unsigned int c = mUChannelMin + 1; c < mZChannelMax; c++ ){
                if ( mGeometryCore->View(c) == geo::kV && mGeometryCore->View(c-1) == geo::kU ){
                    mVChannelMin = c;
                    mUChannelMax = c - 1;
                }
                if ( mGeometryCore->View(c) == geo::kZ && mGeometryCore->View(c-1) == geo::kV ){
                    mZChannelMin = c;
                    mVChannelMax = c-1;
                }
            }

            //Number of channels in each view
            mNumberOfUChannels = mUChannelMax - mUChannelMin+1; //U
            mNumberOfVChannels = mVChannelMax - mVChannelMin+1; //V
            mNumberOfZChannels = mZChannelMax - mZChannelMin+1; //Z (collection plane)

            for(auto channel : mGeometryCore->ChannelsInTPCs())
            {
                mChannelToWireIDMap[channel] = mGeometryCore->ChannelToWire(channel);
            }

            // collect world info
            mWorldName = mGeometryCore->GetWorldVolumeName();
            mWorldBox.setBox(mGeometryCore->WorldBox());
            // create name-volumetype map mor world
            mMaterialPOI.SetCoordinates(mWorldBox.x_min,mWorldBox.y_min,mWorldBox.z_min);
            std::string volumeName = mGeometryCore->VolumeName(mMaterialPOI);
            mVolumeTypeMap[volumeName] = VolumeType::World;
            // collect detector inmo
            mDetectorName = mGeometryCore->DetectorName();
            mDetectorBox.setBox(
                -mGeometryCore->DetHalfWidth(), mGeometryCore->DetHalfWidth(),
                -mGeometryCore->DetHalfHeight(), mGeometryCore->DetHalfHeight(),
                0, mGeometryCore->DetLength()
            );
            // collect cryostat info
            // for now, assuming analysis is done over a single cryostat
            geo::CryostatGeo const& Cryo = mGeometryCore->Cryostat();
            mCryostatName = std::string(Cryo.ID());
            mCryostatBox.setBox(Cryo.Boundaries());
            // create name-volumetype map for cryostat
            mMaterialPOI.SetCoordinates(mCryostatBox.x_min, mCryostatBox.y_min, mCryostatBox.z_min);
            volumeName = mGeometryCore->VolumeName(mMaterialPOI);
            mVolumeTypeMap[volumeName] = VolumeType::Cryostat;
            // iterate over all TPCs
            mNumberOfTPCs  = mGeometryCore->TotalNTPC();
            for (geo::TPCGeo const& TPC : mGeometryCore->IterateTPCs())
            {
                mTPCNames.emplace_back(TPC.ID());
                mTPCBoxes.emplace_back(BoundingBox(TPC.BoundingBox()));
                mActiveTPCBoxes.emplace_back(BoundingBox(TPC.ActiveBoundingBox()));
                mTPCMasses.emplace_back(TPC.ActiveMass());
                mTPCDriftDistances.emplace_back(TPC.DriftDistance());
                // create name-volumetype map for this tpc
                mVolumeTypeMap[mGeometryCore->VolumeName(TPC.GetCenter())] = VolumeType::TPC;
            }
            // find the total TPC and total Active TPC volumes
            FindTotalTPCBoxes();
            mTotalTPCMass = mGeometryCore->TotalMass();    
        }

        // get volume information for a point
        DetectorVolume DetectorGeometry::GetVolume(std::vector<double> position)
        {
            return GetVolume(position[0], position[1], position[2]);
        }
        // get volume information for a point
        DetectorVolume DetectorGeometry::GetVolume(double x, double y, double z)
        {
            mMaterialPOI.SetCoordinates(x,y,z);

            // get the volume information
            std::string volumeName = mGeometryCore->VolumeName(mMaterialPOI);
            VolumeType volumeType = mVolumeTypeMap[volumeName];

            // get the current material information
            mMaterial = mGeometryService->Material(mMaterialPOI);
            double material = mMaterial->GetZ();
            std::string materialName = mMaterial->GetName();

            // return the constructed volume 
            return DetectorVolume(volumeType, volumeName, materialName, material);
        }
        
        // get total tpc volume information
        void DetectorGeometry::FindTotalTPCBoxes()
        {
            double x_min = 0; double x_max = 0;
            double y_min = 0; double y_max = 0;
            double z_min = 0; double z_max = 0;
            for (size_t i = 0; i < mTPCBoxes.size(); i++) 
            {
                if (mTPCBoxes[i].x_min < x_min) x_min = mTPCBoxes[i].x_min;
                if (mTPCBoxes[i].x_max > x_max) x_max = mTPCBoxes[i].x_max;
                if (mTPCBoxes[i].y_min < y_min) y_min = mTPCBoxes[i].y_min;
                if (mTPCBoxes[i].y_max > y_max) y_max = mTPCBoxes[i].y_max;
                if (mTPCBoxes[i].z_min < z_min) z_min = mTPCBoxes[i].z_min;
                if (mTPCBoxes[i].z_max > z_max) z_max = mTPCBoxes[i].z_max;
            }
            mTotalTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
            x_min = 0; x_max = 0;
            y_min = 0; y_max = 0;
            z_min = 0; z_max = 0;
            for (size_t i = 0; i < mActiveTPCBoxes.size(); i++) 
            {
                if (mActiveTPCBoxes[i].x_min < x_min) x_min = mActiveTPCBoxes[i].x_min;
                if (mActiveTPCBoxes[i].x_max > x_max) x_max = mActiveTPCBoxes[i].x_max;
                if (mActiveTPCBoxes[i].y_min < y_min) y_min = mActiveTPCBoxes[i].y_min;
                if (mActiveTPCBoxes[i].y_max > y_max) y_max = mActiveTPCBoxes[i].y_max;
                if (mActiveTPCBoxes[i].z_min < z_min) z_min = mActiveTPCBoxes[i].z_min;
                if (mActiveTPCBoxes[i].z_max > z_max) z_max = mActiveTPCBoxes[i].z_max;
            }
            mTotalActiveTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
        }
        void DetectorGeometry::FillTTree()
        {
            // add geometry info
            mGeometryTree->Branch("world_name", &mWorldName);
            mGeometryTree->Branch("world_box_ranges", &(mWorldBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            mGeometryTree->Branch("detector_name", &mDetectorName);
            mGeometryTree->Branch("detector_box_ranges", &(mDetectorBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            mGeometryTree->Branch("cryostat_name", &mCryostatName);
            mGeometryTree->Branch("cryostat_box_ranges", &(mCryostatBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            mGeometryTree->Branch("number_of_tpcs", &mNumberOfTPCs);
            mGeometryTree->Branch("tpc_names", &mTPCNames);
            for (int i = 0; i < mNumberOfTPCs; i++) 
            {
                mGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_name").c_str(), &(mTPCNames[i]));
                mGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_box_ranges").c_str(), &(mTPCBoxes[i]), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
                mGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_mass").c_str(), &(mTPCMasses[i]));
                mGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_drift_distance").c_str(), &(mTPCDriftDistances[i]));
            }
            mGeometryTree->Branch("tpc_masses", &mTPCMasses);
            mGeometryTree->Branch("tpc_drift_distances", &mTPCDriftDistances);
            mGeometryTree->Branch("total_tpc_box_ranges", &(mTotalTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            mGeometryTree->Branch("total_active_tpc_box_ranges", &(mTotalActiveTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            mGeometryTree->Branch("total_tpc_mass", &mTotalTPCMass);
            mGeometryTree->Fill();
        }
    }
}