# Surface Acoustic Wave Device SPICE Generator

Script for generating models for surface acoustic wave devices in SPICE. The output is divided in two files, one for the hierarchical blocks (out.cir) and other for a testbench with a sample circuit and simulation directives (tb.cir). 

For graphical usability of the model, I included a sample LTSPICE file (tb.asc) that imports the circuit generated from the script and performs common analysis.

## Dependencies

* Python3;
* Numpy.

## Usage

This software is spread between two python files. While SAW_FUNCTIONS.py carries all the necessary functions for the script, SAW_SPICE_GEN.py is responsible for the model generation. Before running the program, make sure that both files are in the **same folder**. Then, execute it with:

```
python SAW_SPICE.GEN.py
```

If you have more than one python version in your computer, you may need to specify it when running, such as

```
python3 SAW_SPICE.GEN.py
```

Then, the program will ask you to input a few parameters necessary for evaluation of the models. If you are not sure about any of them, you can just press _Enter_ in your keyboard and the software will use default values. The most important parameters, however, are the center frequency (f0) and the number of finger pairs in one IDT, and they must be known if you desire that the model minimally replicates your device's behavior.

When generating the model, the program produces many different hierarchical blocks, according to what is better described in my Master's degree dissertation "MODELAGEM DE DISPOSITIVOS DE ONDAS ACÚSTICAS DE SUPERFÍCIE COM FOCO EM AMBIENTES DE SIMULAÇÃO" or, in english, "SIMULATION-ORIENTED MODELLING OF SURFACE ACOUSTIC WAVE DEVICES". Depending on when you read this documentation, it may be easier or harder for this manuscript to be found online, so you can always contact me directly via email.

When the program is done, you can test the model in any SPICE simulator (or compatible). I would recommend LTSPICE because it is free. If you use ADS or something different, you can either look how to import SPICE models or you could read the text files (out.cir and tb.cir) and implement the circuit in your graphic environment. I have included a graphical LTSPICE example, where the model is implemented as text but the rest of the circuit is a schematic as you are probably used, so check out tb.asc if you don't know where to start from.

In the following section, I give a general description of the model created by this script.

## Description of the model

The sample circuit, present in the tb.cir file, is shown in the figure.

<p align="center">
  <img src="images/fullcircuit.png" width="900" />
</p>

The final model for the device is contained in the DLINE block, and it can be seen in the figure. Inside of it, there are two IDTs, which themselves are composed of smaller FING models, connected in series, determined by the amount of finger pairs chosen by the user, during the script execution. There are two alternatives for the FING model, the first is the distributed components model, which has the best behavior in general, with good frequency and time domain responses and it is the recommended for most uses.

However, you can also use the lumped component model, which approximates the A and B impedances (dependent on the tangent and cosecant of the frequency) with LC pairs through the Foster or Mittag-Leffler theorems. Since this is an approximation of periodic functions, it is necessary to choose how many periods will be approximated (n=1 to 5) for Foster or any number n>1 for Mittag-Leffler. The lumped component approach results in an unstable model in the time domain simulations, but can still be used in the frequency domain. The value of n is responsible for how many periods appear in the frequency response, as seen below.

<p align="center">
  <img src="images/broad.png" width="300" />
</p>

The SENSE block, between IDTs, allows for good customization of the response by the user. In that hierarchical block, the values of Ap and tp are defined as parameters, that are controlled by the tb.cir file. By changing these values, the user can freely attenuate or delay the response of the device, which is interesting for simulating perturbations in a sensor. The limitation of this circuit is that Ap can only model **attenuations**, so if you want to increase the gain level, you can change the electroacoustic coupling parameter (K^2) during the execution of the script. 

The capacitor Cf can also be modified through a parameter in tb.cir. This capacitance models the electromagnetic coupling between IDTs and is responsible for:
* Electromagnetic pulse at the beginning of the time-domain response;
* Imperfections in the frequency domain response, as seen in the figures below.

<p align="center">
  <img src="images/capdiff.png" width="600" />
</p>

PS: these graphs are in different scales.

If you want more information, you can look at our published papers.

## Published papers

* Evaluation of an Equivalent Circuit Model for Simulation of Surface Acoustic Wave Sensors (INSCIT 2021) - https://ieeexplore.ieee.org/document/9557258
* Equivalent Circuit Models for SAW Delay Line Sensors (IEEE Sensors 2022) - https://ieeexplore.ieee.org/abstract/document/9768799

## Developed at

<p align="center">
  <img src="images/ufs_horizontal_positiva.png" width="300" />
</p>

by Raphael Cardoso (cardosodeoliveir@gmail.com)

Feel free to contact me via email if you have any questions or there is any bugs in this program.
