# Surface Acoustic Wave Device SPICE Generator

Script for generating models for surface acoustic wave devices in SPICE. The output is divided in two files, one for the hierarchical blocks (default out.cir) and other for a testbench with a sample circuit and simulation directives (tb.cir). 

## Dependencies

* Python3
* Numpy

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

When generating the model, the program produces many different hierarchical blocks, according to what is better described in my Master's degree dissertation "MODELAGEM DE DISPOSITIVOS DE ONDAS ACÚSTICAS DE SUPERFÍCIE COM FOCO EM AMBIENTES DE SIMULAÇÃO" or, in english "SIMULATION-ORIENTED MODELLING OF SURFACE ACOUSTIC WAVE DEVICES". Depending on when you read this documentation, it may be easier or harder for this manuscript to be found online, so you can always contact me directly via email.

In the following section, I give a general description of the model created by this script.

## Description of the model

# Developed at

<p align="center">
  <img src="imagens/ufs_horizontal_positiva.png" width="300" />
</p>

by Raphael Cardoso (cardosodeoliveir@gmail.com)
