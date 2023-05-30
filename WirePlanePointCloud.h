/**
 * @file WirePlanePointCloud.h
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
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
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

#include "Core.h"
#include "DetectorGeometry.h"

namespace arrakis
{
    struct WirePlanePointCloud
    {
        /**
         * A WirePlanePointCloud is the basic data structure that
         * we use for data from the detector output.  Technically,
         * since each variable (besides the energy) is discrete,
         * and contiguous, the data is automatically voxelized. Perhaps
         * then a more proper name is a *WirePlaneVoxelGrid*.
         * 
         * The main variables of interest are the signals (adc) as
         * a function of (channel, tick) or equivalently (channel, tdc).
         * To each point in the discriminating variables adc(channel, tdc),
         * we assign a set of labels which if known would allow
         * unambiguous analysis of the data.
         * 
         *  (a) source - this label identifies what the generator of the
         *               associated primary was.  For example, all particles
         *               whose ancestor was a cosmic ray would get the label
         *               'cosmic', while all radiologicals get 'radiological'.
         *               Some special classes are 'PNS', which stands for the
         *               Pulsed Neutron Source, and 'HEPevt', which comes from
         *               any HEPevt provided file.
         * 
         *  (b) shape - shape refers to a high level description of the geometry
         *               of the associated charge depositions.  Long one-dimensional
         *               tracks are aptly called 'track', while electron/position/
         *               photon showers are called 'shower'.  Other current shapes 
         *               are 'blip', which are low energy activity, and a specific
         *               blip 'neutron_capture'.
         * 
         *  (c) particle - the class particle most often refers to the actual 
         *               particle that left behind the energy deposition, however
         *               this is not always the case and can be subtle.  
         * 
         *  (d) unique_shape - this is a clustering label, meant to identify
         *               unique instances of a shape, i.e. since when the network
         *               learns to identify the shape 'track', it doesn't also learn
         *               whether two pixels belong to the same track or not.  That
         *               is what this label is meant to do.
         * 
         *  (e) unique_particle - this is also a clustering label with the same logic
         *               as unique_shape but now with respect to the 'particle' label.  
        */
        std::vector<Int_t> channel = {};
        std::vector<Int_t> wire = {};
        std::vector<Int_t> tick = {};
        std::vector<Int_t> tdc = {};
        std::vector<Int_t> adc = {};
        std::vector<Int_t> view = {};
        std::vector<Double_t> energy = {};

        std::vector<std::vector<Int_t>> track_ids = {};
        std::vector<std::vector<Double_t>> energies = {};
        std::vector<std::vector<Double_t>> x = {};
        std::vector<std::vector<Double_t>> y = {};
        std::vector<std::vector<Double_t>> z = {};
        std::vector<std::vector<SourceLabelInt>> source_labels = {};
        std::vector<std::vector<ShapeLabelInt>> shape_labels = {};
        std::vector<std::vector<ParticleLabelInt>> particle_labels = {};
        std::vector<std::vector<Int_t>> unique_sources = {};
        std::vector<std::vector<Int_t>> unique_shapes = {};
        std::vector<std::vector<Int_t>> unique_particles = {};

        std::vector<SourceLabelInt> source_label = {};
        std::vector<ShapeLabelInt> shape_label = {};
        std::vector<ParticleLabelInt> particle_label = {};

        std::vector<Int_t> unique_source = {};
        std::vector<Int_t> unique_shape = {};
        std::vector<Int_t> unique_particle = {};



        WirePlanePointCloud()
        {
        }

        void clear()
        {
            channel.clear();
            wire.clear();
            tick.clear();
            tdc.clear();
            adc.clear();
            view.clear();
            energy.clear();

            track_ids.clear();
            energies.clear();
            x.clear();
            y.clear();
            z.clear();
            source_labels.clear();
            shape_labels.clear();
            particle_labels.clear();
            unique_sources.clear();
            unique_shapes.clear();
            unique_particles.clear();
            
            source_label.clear();
            shape_label.clear();
            particle_label.clear();

            unique_source.clear();
            unique_shape.clear();
            unique_particle.clear();
        }

        // This function is used primarily with Data,
        // which doesn't have other truth information
        // associated with the points.
        void AddPoint(
            detinfo::DetectorClocksData const& clock_data,
            Int_t det_tick,
            Int_t det_channel,
            Int_t det_adc
        )
        {
            auto wires = DetectorGeometry::GetInstance()->ChannelToWire(det_channel);
            auto det_view = DetectorGeometry::GetInstance()->View(det_channel);
            Double_t wire_multiple = DetectorGeometry::GetInstance()->GetWirePitch(det_view);

            channel.emplace_back(det_channel);
            wire.emplace_back(wires[det_view].Wire * wire_multiple);
            tick.emplace_back(det_tick);
            tdc.emplace_back(clock_data.TPCTick2TDC(det_tick));
            adc.emplace_back(det_adc);
            view.emplace_back(det_view);
        }

        // This function is used for simulation, where
        // truth information about RawDigits is in the
        // sim::IDE vector.
        void AddPoint(
            detinfo::DetectorClocksData const& clock_data,
            std::vector<sim::IDE> det_ide,
            Int_t det_tick,
            Int_t det_channel,
            Int_t det_adc,
            bool det_noise
        )
        {
            /**
             * 
            */
            auto wires = DetectorGeometry::GetInstance()->ChannelToWire(det_channel);
            auto det_view = DetectorGeometry::GetInstance()->View(det_channel);
            Double_t wire_multiple = DetectorGeometry::GetInstance()->GetWirePitch(det_view);

            channel.emplace_back(det_channel);
            wire.emplace_back(wires[det_view].Wire * wire_multiple);
            tick.emplace_back(det_tick);
            tdc.emplace_back(clock_data.TPCTick2TDC(det_tick));
            adc.emplace_back(det_adc);
            view.emplace_back(det_view);

            std::vector<Int_t> det_track_ids;
            std::vector<Double_t> det_energies;
            std::vector<Double_t> det_x;
            std::vector<Double_t> det_y;
            std::vector<Double_t> det_z;
            std::vector<SourceLabelInt> det_source;
            std::vector<ShapeLabelInt> det_shape;
            std::vector<ParticleLabelInt> det_particle;
            std::vector<Int_t> det_unique_shape;
            std::vector<Int_t> det_unique_particle;
            Double_t det_energy = 0.0;
            for(auto ide : det_ide)
            {
                det_track_ids.emplace_back(ide.trackID);
                det_energies.emplace_back(ide.energy);
                det_x.emplace_back(ide.x);
                det_y.emplace_back(ide.y);
                det_z.emplace_back(ide.z);
                if(det_noise)
                {
                    det_source.emplace_back(LabelCast(SourceLabel::Noise));
                    det_shape.emplace_back(LabelCast(ShapeLabel::Noise));
                    det_particle.emplace_back(LabelCast(ParticleLabel::Noise));
                }
                else
                {
                    det_source.emplace_back(LabelCast(SourceLabel::Undefined));
                    det_shape.emplace_back(LabelCast(ShapeLabel::Undefined));
                    det_particle.emplace_back(LabelCast(ParticleLabel::Undefined));
                }
                det_unique_shape.emplace_back(-1);
                det_unique_particle.emplace_back(-1);
                det_energy += ide.energy;
            }
            energy.emplace_back(det_energy);
            track_ids.emplace_back(det_track_ids);
            energies.emplace_back(det_energies);
            x.emplace_back(det_x);
            y.emplace_back(det_y);
            z.emplace_back(det_z);
            source_labels.emplace_back(det_source);
            shape_labels.emplace_back(det_shape);
            particle_labels.emplace_back(det_particle);
            unique_shapes.emplace_back(det_unique_shape);
            unique_particles.emplace_back(det_unique_particle);

            if(det_noise)
            {
                source_label.emplace_back(LabelCast(SourceLabel::Noise));
                shape_label.emplace_back(LabelCast(ShapeLabel::Noise));
                particle_label.emplace_back(LabelCast(ParticleLabel::Noise));
            }
            else
            {
                source_label.emplace_back(LabelCast(SourceLabel::Undefined));
                shape_label.emplace_back(LabelCast(ShapeLabel::Undefined));
                particle_label.emplace_back(LabelCast(ParticleLabel::Undefined));
            }
            unique_shape.emplace_back(-1);
            unique_particle.emplace_back(-1);
        }
        Int_t GetIndex_TrackID(DetSimID_t detsim, TrackID_t track_id)
        {
            for(size_t ii = 0; ii < track_ids[detsim].size(); ii++) {
                if(track_ids[detsim][ii] == track_id) {
                    return ii;
                }
            }
            return -1;
        }
    };
}