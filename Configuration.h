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

        /**
         * @brief Producer and Instance labels.
         * 
         */
        fhicl::Atom<art::InputTag> LArGeantProducerLabel
        {
            fhicl::Name("LArGeantProducerLabel"),
            fhicl::Comment("Tag of the input data product for the largeant side of the simulation.")
        };
        
        fhicl::Atom<art::InputTag> IonAndScintProducerLabel
        {
            fhicl::Name("IonAndScintProducerLabel"),
            fhicl::Comment("Tag of the input data product for the IonAndScint side of the simulation.")
        };
        fhicl::Atom<art::InputTag> IonAndScintInstanceLabel
        {
            fhicl::Name("IonAndScintInstanceLabel"),
            fhicl::Comment("Tag of the input data product for the IonAndScint side of the simulation.")
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
            fhicl::Atom<art::InputTag> CosmicsLabel
            {
                fhicl::Name("CosmicsLabel"),
                fhicl::Comment("Tag of the input data product for the CosmicsLabel.")
            };
            fhicl::Atom<art::InputTag> SingleNeutronLabel
            {
                fhicl::Name("SingleNeutronLabel"),
                fhicl::Comment("Tag of the input data product for the SingleNeutronLabel.")
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
        struct MelangeParameters
        {
            fhicl::Atom<art::InputTag> NeutronCaptureGammaDetail
            {
                fhicl::Name("NeutronCaptureGammaDetail"),
                fhicl::Comment("The level of detail for labeling neutron capture gammas, can be either 'simple', 'medium' of 'full'.")
            };
        };
        fhicl::Table<MelangeParameters> melange_parameters
        {
            fhicl::Name("MelangeParameters"),
            fhicl::Comment("Parameters for the melange algorithms.")
        };

        fhicl::Atom<double> ADCThreshold
        {
            fhicl::Name("ADCThreshold"),
            fhicl::Comment("ADC threshold value to use for filling 2D arrays.")
        };    

    };

    using Parameters = art::EDAnalyzer::Table<Configuration>;
}