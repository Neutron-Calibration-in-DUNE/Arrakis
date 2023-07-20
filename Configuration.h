/**
 * @file    Configuration.h
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   A struct for holding LArSoft configuration parameters
 *          for the Arrakis module.
 * @version 0.1
 * @date 2022-02-15
 */
#pragma once

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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

namespace arrakis
{
    struct Configuration
    {
        fhicl::Atom<std::string> ProcessType
        {
            fhicl::Name("ProcessType"),
            fhicl::Comment("Whether the data to process is 'simulation' or 'data'.")
        };
        fhicl::Atom<bool> ProcessMCTruth
        {
            fhicl::Name("ProcessMCTruth"),
            fhicl::Comment("Whether to process MCTruth.")
        };
        fhicl::Atom<bool> ProcessMCParticles
        {
            fhicl::Name("ProcessMCParticles"),
            fhicl::Comment("Whether to process MCParticles.")
        };
        fhicl::Atom<bool> ProcessSimEnergyDeposits
        {
            fhicl::Name("ProcessSimEnergyDeposits"),
            fhicl::Comment("Whether to process SimEnergyDesposits.")
        };
        fhicl::Atom<bool> ProcessSimChannels
        {
            fhicl::Name("ProcessSimChannels"),
            fhicl::Comment("Whether to process SimChannels.")
        };
        fhicl::Atom<bool> ProcessRawDigits
        {
            fhicl::Name("ProcessRawDigits"),
            fhicl::Comment("Whether to process RawDigits.")
        };
        fhicl::Atom<bool> ProcessHits
        {
            fhicl::Name("ProcessHits"),
            fhicl::Comment("Whether to process Hits.")
        };
        fhicl::Atom<bool> ProcessOpDetBacktrackerRecords
        {
            fhicl::Name("ProcessOpDetBacktrackerRecords"),
            fhicl::Comment("Whether to process OpDetBacktrackerRecords.")
        };
        fhicl::Atom<bool> ProcessOpDetWaveforms
        {
            fhicl::Name("ProcessOpDetWaveforms"),
            fhicl::Comment("Whether to process OpDetWaveforms.")
        };

        /**
         * @brief Producer and Instance labels.
         * 
         */
        fhicl::Atom<art::InputTag> LArGeantProducerLabel
        {
            fhicl::Name("LArGeantProducerLabel"),
            fhicl::Comment("Tag of the input data product for the largeant side of the simulation.")
        };
        
        fhicl::Atom<art::InputTag> SimEnergyDepositProducerLabel
        {
            fhicl::Name("SimEnergyDepositProducerLabel"),
            fhicl::Comment("Tag of the input data product for the SimEnergyDeposit side of the simulation.")
        };
        fhicl::Atom<art::InputTag> SimEnergyDepositInstanceLabel
        {
            fhicl::Name("SimEnergyDepositInstanceLabel"),
            fhicl::Comment("Tag of the input data product for the SimEnergyDeposit side of the simulation.")
        };
        fhicl::Atom<art::InputTag> SimChannelProducerLabel
        {
            fhicl::Name("SimChannelProducerLabel"),
            fhicl::Comment("Tag of the input data product for the SimChannelProducerLabel.")
        };

        fhicl::Atom<art::InputTag> SimChannelInstanceLabel
        {
            fhicl::Name("SimChannelInstanceLabel"),
            fhicl::Comment("Tag of the input data product for the SimChannelInstanceLabel.")
        };

        fhicl::Atom<art::InputTag> RawDigitProducerLabel
        {
            fhicl::Name("RawDigitProducerLabel"),
            fhicl::Comment("Tag of the input data product for the RawDigitProducerLabel.")
        };

        fhicl::Atom<art::InputTag> RawDigitInstanceLabel
        {
            fhicl::Name("RawDigitInstanceLabel"),
            fhicl::Comment("Tag of the input data product for the RawDigitInstanceLabel.")
        };

        fhicl::Atom<art::InputTag> OpDetBacktrackerRecordProducerLabel
        {
            fhicl::Name("OpDetBacktrackerRecordProducerLabel"),
            fhicl::Comment("Tag of the input data product for the OpDetBacktrackerRecordProducerLabel.")
        };

        fhicl::Atom<art::InputTag> OpDetWaveformProducerLabel
        {
            fhicl::Name("OpDetWaveformProducerLabel"),
            fhicl::Comment("Tag of the input data product for the OpDetWaveformProducerLabel.")
        };

        /**
         * @brief Generator labels for various 
         * particle generators used in the simulation.
         */
        struct GeneratorLabels
        {
            fhicl::Atom<art::InputTag> Ar39Label
            {
                fhicl::Name("Ar39Label"),
                fhicl::Comment("Tag of the input data product for the Ar39Label.")
            };
            fhicl::Atom<art::InputTag> Ar42Label
            {
                fhicl::Name("Ar42Label"),
                fhicl::Comment("Tag of the input data product for the Ar42Label.")
            };
            fhicl::Atom<art::InputTag> Kr85Label
            {
                fhicl::Name("Kr85Label"),
                fhicl::Comment("Tag of the input data product for the Kr85Label.")
            };
            fhicl::Atom<art::InputTag> Rn222Label
            {
                fhicl::Name("Rn222Label"),
                fhicl::Comment("Tag of the input data product for the Rn222Label.")
            };
            fhicl::Atom<art::InputTag> BeamLabel
            {
                fhicl::Name("BeamLabel"),
                fhicl::Comment("Tag of the input data product for the BeamLabel.")
            };
            fhicl::Atom<art::InputTag> CosmicsLabel
            {
                fhicl::Name("CosmicsLabel"),
                fhicl::Comment("Tag of the input data product for the CosmicsLabel.")
            };
            fhicl::Atom<art::InputTag> HEPevtLabel
            {
                fhicl::Name("HEPevtLabel"),
                fhicl::Comment("Tag of the input data product for the HEPevt.")
            };
            fhicl::Atom<art::InputTag> PNSLabel
            {
                fhicl::Name("PNSLabel"),
                fhicl::Comment("Tag of the input data product for the PNSLabel.")
            };
        };
        fhicl::Table<GeneratorLabels> labels
        {
            fhicl::Name("GeneratorLabels"),
            fhicl::Comment("Input Tag Table for accessing GeneratorLabels for various generator labels.")
        };
        
        fhicl::Atom<bool> SaveMeta
        {
            fhicl::Name("SaveMeta"),
            fhicl::Comment("Whether to save meta info arrays.")
        };
        fhicl::Atom<bool> SaveGeometry
        {
            fhicl::Name("SaveGeometry"),
            fhicl::Comment("Whether to save geometry arrays.")
        };
        fhicl::Atom<bool> SaveEnergyDepositPointCloud
        {
            fhicl::Name("SaveEnergyDepositPointCloud"),
            fhicl::Comment("Whether to save the energy deposit point clouds.")
        };
        fhicl::Atom<bool> SaveSimulationWrangler
        {
            fhicl::Name("SaveSimulationWrangler"),
            fhicl::Comment("Whether to save MC Maps to the output file.")
        };
        fhicl::Atom<bool> SaveWirePlaneHits
        {
            fhicl::Name("SaveWirePlaneHits"),
            fhicl::Comment("Whether to save the wire plane hits.")
        };
        fhicl::Atom<bool> SaveWirePlanePointCloud
        {
            fhicl::Name("SaveWirePlanePointCloud"),
            fhicl::Comment("Whether to save the wire plane point clouds.")
        };
        fhicl::Atom<bool> SaveWirePlaneTrackTopology
        {
            fhicl::Name("SaveWirePlaneTrackTopology"),
            fhicl::Comment("Whether to save the wire plane track topology.")
        };
        fhicl::Atom<bool> SaveOpDetPointCloud
        {
            fhicl::Name("SaveOpDetPointCloud"),
            fhicl::Comment("Whether to save the optical detector point clouds.")
        };

        /**
         * Parameters for the simulation labeling logic class.
         */
        fhicl::Atom<bool> SaveWirePlaneInductionFlag
        {
            fhicl::Name("SaveWirePlaneInductionFlag"),
            fhicl::Comment("Whether to save the flag that labels induction signals.")
        };
        fhicl::Atom<bool> SaveWirePlaneUnderThreshold
        {
            fhicl::Name("SaveWirePlaneUnderThreshold"),
            fhicl::Comment("Whether to save the true points which are under threshold.")
        };
        fhicl::Atom<std::string> NeutronCaptureGammaDetail
        {
            fhicl::Name("NeutronCaptureGammaDetail"),
            fhicl::Comment("The level of detail for labeling neutron capture gammas, can be either 'simple', 'medium' of 'full'.")
        };
        fhicl::Atom<Int_t> ADCThreshold
        {
            fhicl::Name("ADCThreshold"),
            fhicl::Comment("ADC threshold value to use for filling 2D arrays.")
        };    
        fhicl::Atom<Int_t> InducedChannelInfluence
        {
            fhicl::Name("InducedChannelInfluence"),
            fhicl::Comment("The amount of channels a nearby wire effects through induction.  This is to account for the wirecell simulation which is not back tracked to sim::IDE's in SimChannel.")
        };
        fhicl::Atom<Int_t> InducedTDCInfluence
        {
            fhicl::Name("InducedTDCInfluence"),
            fhicl::Comment("The amount of tdc values influenced by a surrounding signal.  This is to account for the wirecell simulation which is not back tracked to sim::IDE's in SimChannel.")
        };
    };

    using Parameters = art::EDAnalyzer::Table<Configuration>;
}