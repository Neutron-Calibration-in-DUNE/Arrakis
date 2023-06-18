# Arrakis

Arrakis is a LArSoft module for creating training data for ML tasks.  

### Table of Contents

1. [ Installing Arrakis ](#install)
2. [ Using Arrakis ](#usage)
	* [ LArSoft Configuration ](#config)

<a name="install"></a>
## Installing Arrakis
Currently, Arrakis can only be installed by cloning this repository and compiling it with LArSoft.  This can be done easily by using the associated [LArSoftArrakis](https://github.com/Neutron-Calibration-in-DUNE/LArSoftArrakis) repository.  To include in your own local LArSoft install, clone this repository into one of the *larsoft* or *dune* namespace products,
```{bash}
/local_larsoft_dir/srcs/<namespace>/:$ git clone https://github.com/Neutron-Calibration-in-DUNE/Arrakis
```
Then, you'll have to edit the CMakeLists.txt in the *namespace* directory to include the lines,
```{cmake}
add_subdirectory(Arrakis)
```

We are currently working to make Arrakis part of a LArSoft namespace.

<a name="usage"></a>
## Using Arrakis

<a name="config"></a>

### LArSoft Configuration
There are several FHiCL parameters that must be specified at run time.  These include the *Producer* and *Instance* labels for the various products that Arrakis uses.  Other required parameters are the names of the various generators that are used in the simulation.  Both of these are shown below from the [Arrakis.fcl](https://github.com/Neutron-Calibration-in-DUNE/Arrakis/blob/main/Arrakis.fcl) file,
```{yaml}
LArGeantProducerLabel:      "largeant"
SimEnergyDepositProducerLabel:   "IonAndScint"
SimEnergyDepositInstanceLabel:   "priorSCE"
SimChannelProducerLabel:    "tpcrawdecoder"
SimChannelInstanceLabel:    "simpleSC"
RawDigitProducerLabel:      "tpcrawdecoder"
RawDigitInstanceLabel:      "daq"
OpDetWaveformProducerLabel: "opdigi"

GeneratorLabels: 
{
    Ar39Label:          "Ar39" 
    Ar42Label:          "Ar42"
    Kr85Label:          "Kr85"
    Rn222Label:         "Rn222"
    BeamLabel:          "Beam"
    CosmicsLabel:       "Cosmics"
    HEPevtLabel:        "HEPevt"
    PNSLabel:           "PNS"
}
```

