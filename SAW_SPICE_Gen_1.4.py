#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 16:06:49 2021
UNIVERSIDADE FEDERAL DE SERGIPE

SAW SPICE Generator 1.4
Author: Raphael Cardoso

CHANGELOG:
    1.X -> 2.0: GUI (TODO)
    1.3 -> 1.4: BETTER USER INPUT METHOD
    1.2 -> 1.3: SENSE BLOCK
    1.1 -> 1.2: TRANSMISSION LINE MODEL
    1.0 -> 1.1: MITTAG LEFFLER'S THEOREM
    0.0 -> 1.0: FOSTER THEOREM, LUMPED IMPEDANCE MODEL
"""

import numpy as np
from SAW_FUNCTIONS import *

# ======================== START OF THE CODE ==================================

# =============================================================================
print("Welcome to the SAW SPICE Generator 1.4!")
print("This program was written by Raphael Cardoso")
print("From Universidade Federal de Sergipe (UFS)")
print("------------------------------------")
f0 = input("Please enter the center frequency in MHz (empty for 100 MHz): ")
if (f0==''):
    f0 = 100e6
else:
    f0 = float(f0)*1e6

Ns = input("Please enter number of finger pairs (empty for 10): ")
if (Ns==''):
    Ns=10
else:
    Ns = int(Ns)

tl = input("Please enter the time of flight of the acoustic signal in microseconds (empty for 1 us): ")
if (tl==''):
    tl=1e-6
else:
    tl=float(tl)*1e-6

Ct = input("Please enter the IDT capacitance in pF (empty for 1 pF): ")
if (Ct==''):
    Ct = 1e-12
else:
    Ct=float(Ct)*1e-12

Rt = input("Please enter the IDT resistance in Ohms (empty for 10 Ohms): ")
if (Rt==''):
    Rt=10
else:
    Rt = float(Rt)

kk = input("Please enter the K^2 parameter in percentage (empty for 1%): ")
if (kk==''):
    kk=1/100
else:
    kk= float(kk)/100

modelType = input("Please choose between lumped (L) or distributed (D) component model for finger pairs (empty for D): ")
if (modelType == ''):
    modelType = 'D'

if (modelType == 'L'):
    n = input("Please enter the number of LC pairs you want to use from 1 to 5 (empty for 3): ")
    if (n==''):
        n=3
    else:
        n = int(n)
print("------------------------------------")
print("Calculating . . .")
# =============================================================================

w0 = 2*np.pi*f0
Cs = Ct/Ns
R0 = 2*np.pi/(w0*Cs*kk)
Rs = Rt*Ns

print("Characteristic impedance in MOhm: ", R0/1e6)

blocks_name = "out.cir"
print("File in which the hierarchical blocks are stored: ", blocks_name)
f = open(blocks_name, "w") # Open file in write mode - will erase all
f.write(".title SPICE Blocks of SAW Device")
f.write("\n")
f.close()

testbench_name = "tb.cir"
print("File in which the test circuit is stored: ", testbench_name)
f = open(testbench_name, "w") # Open file in write mode - will erase all
f.write(".title SPICE Testbench of SAW Device")
f.write("\n")
f.close()


# Initializing SUBCKT names
dline_subckt_name = "DLINE"
idt_subckt_name = "IDT"
fing_subckt_name = "FING"
tan_subckt_name = "TAN"
csc_subckt_name = "CSC"
sense_subckt_name = "SENSE"

# ============================ Generating testbench ===========================
testbench_gen(testbench_name, blocks_name, dline_subckt_name, f0, tl)

# ============================ Generating DLINE block =========--==============

append_DLINE(blocks_name, dline_subckt_name, idt_subckt_name, sense_subckt_name, R0)

# ============================ Generating IDT block ===========================

append_IDT(blocks_name, idt_subckt_name, fing_subckt_name, Ns)

# ============================ Generating FING block ==========================

if (modelType=='L'):
    append_fing_bhata(blocks_name, fing_subckt_name, tan_subckt_name, csc_subckt_name, Rs/2, Cs/2)
else:
    append_fing_ltra(blocks_name, fing_subckt_name, R0, f0, Rs/2, Cs/2)

# =========================== Generating SENSE block ==========================

append_SENSE(blocks_name, sense_subckt_name, R0, tl)

# ======================== LUMPED MODEL ONLY - LC PAIRS =======================
if (modelType=='L'):
# ============================ Generating TAN block ===========================
    L_tan = np.zeros(n)
    C_tan = np.zeros(n)
    
    L_tan, C_tan = foster_tan_gen(n, R0, w0)
    #L_tan, C_tan = mittag_tan_gen(n, R0, w0)
    
    append_LC_series(blocks_name, tan_subckt_name, C_tan, L_tan)

# ========================== Generating CSC block =============================
    
    L_csc = np.zeros(2*n+1)
    C_csc = np.zeros(2*n+1)
    
    L_csc, C_csc = foster_csc_gen(n,R0,w0)
    #L_csc, C_csc = mittag_csc_gen(n,R0,w0)
    
    append_LC_series(blocks_name, csc_subckt_name, C_csc, L_csc)

# ========================== END ==============================================
print("Finished with no errors!")
