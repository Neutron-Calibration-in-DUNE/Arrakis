/**
 * @file OpDetPointCloud.h
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
#include "lardataobj/RawData/OpDetWaveform.h"

// necessary ROOT libraries
#include <TTree.h>

#include "Core.h"
#include "DetectorGeometry.h"

namespace arrakis
{
    struct OpDetPointCloud
    {
        /**
         * A OpDetPointCloud is the basic data structure that
         * we use for data from the optical detector output.  
         * 
        */
        std::vector<Int_t> channel = {};
        std::vector<Int_t> tick = {};
        std::vector<Int_t> adc = {};
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



        OpDetPointCloud()
        {
        }

        void clear()
        {
            channel.clear();
            tick.clear();
            adc.clear();
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
            Int_t det_tick,
            Int_t det_channel,
            Int_t det_adc
        )
        {
            channel.emplace_back(det_channel);
            tick.emplace_back(det_tick);
            adc.emplace_back(det_adc);
        }

        // This function is used for simulation, where
        // truth information about RawDigits is in the
        // sim::IDE vector.
        void AddPoint(
            Int_t det_channel,
            Int_t det_tick,
            Int_t det_adc,
            bool det_noise
        )
        {
            /**
             * 
            */
            channel.emplace_back(det_channel);
            tick.emplace_back(det_tick);
            adc.emplace_back(det_adc);

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
            // for(auto ide : det_ide)
            // {
            //     det_track_ids.emplace_back(ide.trackID);
            //     det_energies.emplace_back(ide.energy);
            //     det_x.emplace_back(ide.x);
            //     det_y.emplace_back(ide.y);
            //     det_z.emplace_back(ide.z);
            //     if(det_noise)
            //     {
            //         det_source.emplace_back(LabelCast(SourceLabel::Noise));
            //         det_shape.emplace_back(LabelCast(ShapeLabel::Noise));
            //         det_particle.emplace_back(LabelCast(ParticleLabel::Noise));
            //     }
            //     else
            //     {
            //         det_source.emplace_back(LabelCast(SourceLabel::Undefined));
            //         det_shape.emplace_back(LabelCast(ShapeLabel::Undefined));
            //         det_particle.emplace_back(LabelCast(ParticleLabel::Undefined));
            //     }
            //     det_unique_shape.emplace_back(-1);
            //     det_unique_particle.emplace_back(-1);
            //     det_energy += ide.energy;
            // }
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