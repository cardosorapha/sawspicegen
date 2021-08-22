#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 11:29:02 2021

UNIVERSIDADE FEDERAL DE SERGIPE

SAW SPICE Generator 1.4
Author: Raphael Cardoso

This file contains the functions instantiated
in SAW_SPICE_GEN.py
"""

import numpy as np

# ================== FUNCTION DECLARATIONS ====================

# =============================================================
# Function: append_LC_series                                  =
# Description: adds a new subcircuit of LC pairs to the file  =
# No inductor: L=np.inf -- No capacitor: C=0                  =
#                                                             =
# OBS: The LC pairs will be in SERIES                         =
# PIN ORDER: IN OUT                                           =
# =============================================================
def append_LC_series(ofile_name, subckt_name, C, L):
    n = C.size    
    f = open(ofile_name, "a") # Open file in append mode - will not erase all
    f.write("\n") # New line
    f.write(".SUBCKT " + subckt_name + " 1 "+ str(n+1)) # First line of subcircuit
    f.write("\n") # New line
    m = 0 # Counter
    while(m<n):
        if (abs(C[m])>0): #Zero C has infinite Z (open)
            block = "C" + str(m+1) + " " + str(m+1) + " " + str(m+2) + " " + str(C[m])
            # ====== Name ============== IN ============= OUT ========= Value in sci notation ====
            f.write(block)
            f.write("\n") # New line
        
        if (abs(L[m])<np.inf): #Infinite L has infinite Z (open)
            block = "L" + str(m+1) + " " + str(m+1) + " " + str(m+2) + " " + str(L[m])
            # ====== Name ============== IN ============= OUT ========= Value in sci notation ====
            f.write(block)
            f.write("\n") # New line
        
        m = m + 1
    f.write(".ENDS\n")
    f.close()    
    
# =============================================================
# Function: append_IDT                                        =
# Description: adds a new subcircuit of finger pairs          =
# The pairs will be acoustically in series and electrically   =
# in paralel.                                                 =
#                                                             =
# PIN ORDER: V+ V- Ain Aout                                   =
# =============================================================
def append_IDT(ofile_name, subckt_name, fing_subckt_name, Ns):
    f = open(ofile_name, "a")
    f.write("\n") # New line
    f.write(".SUBCKT " + subckt_name + " Vp" + " Vm" + " 1"+ " " + str(Ns+1)) 
    # Pin order: V+ V- Ain Aout 
    f.write("\n") # New line
    m = 1 # Counter
    while(m<=Ns):
        block = "X" + str(m) + " Vp" + " Vm" + " " + str(m) + " " + str(m+1) + " " + fing_subckt_name
        f.write(block)
        f.write("\n")
        m = m + 1
    f.write(".ENDS\n")
    f.close()
    
# =============================================================
# Function: append_fing_bhata                                 =
# Description: adds a new subcircuit composing tan and csc    =
# This model is according to Bhatacharyya's implementation    =
# and explained in our INSCIT paper (2021)                    =
#                                                             =
# PIN ORDER: V+ V- Ain Aout                                   =
# =============================================================
def append_fing_bhata(ofile_name, subckt_name, tan_subckt_name, csc_subckt_name, Rs, Cs):
    f = open(ofile_name, "a")
    f.write("\n") # New line
    f.write(".SUBCKT " + subckt_name + " Vp" + " Vm" + " 1"+ " 7")
    f.write("\n") # New line
    f.write("X1 1 2 " + tan_subckt_name + "\n")
    f.write("X2 2 3 " + csc_subckt_name + "\n")
    f.write("X3 2 4 " + tan_subckt_name + "\n") # 4 is the center node
    
    f.write("C1 3 0 " + str(Cs/2) + "\n")
    f.write("R1 3 Vp " + str(Rs) + "\n")
    
    f.write("X4 4 5 " + tan_subckt_name + "\n")
    f.write("X5 5 6 " + csc_subckt_name + "\n")
    f.write("X6 5 7 " + tan_subckt_name + "\n")
    
    f.write("C2 6 0 " + str(Cs/2) + "\n")
    f.write("R2 6 Vm " + str(Rs) + "\n")
    
    f.write(".ENDS\n")
    f.close()
    
# =============================================================
# Function: append_fing_ltra                                  =
# Description: adds a new subcircuit composing lossless       =
# transmission lines, and is capable of transient simulations =
#                                                             =
#                                                             =
# PIN ORDER: V+ V- Ain Aout                                   =
# =============================================================
def append_fing_ltra(ofile_name, subckt_name, R0, f0, Rs, Cs):
    td = 0.5/f0; # Time between fingers of different polarity
    f = open(ofile_name, "a")
    f.write("\n") # New line
    f.write(".SUBCKT " + subckt_name + " Vp" + " Vm" + " 1"+ " 5")
    f.write("\n") # New line
    f.write("T1 1 2 3 2 " + "Z0=" +str(R0) + " td=" + str(td)+ "\n")
    f.write("T2 3 4 5 4 " + "Z0=" +str(R0) + " td=" + str(td)+ "\n")
    
    f.write("C1 2 0 " + str(Cs/2) + "\n")
    f.write("R1 2 Vp " + str(Rs) + "\n")
    
    f.write("C2 4 0 " + str(Cs/2) + "\n")
    f.write("R2 4 Vm " + str(Rs) + "\n")
    
    f.write(".ENDS\n")
    f.close()

# =============================================================
# Function: append_sense                                      =
# Description: adds a new subcircuit to facilitate adjustments=
# of gain and phase by the user, through the parameters       =
#  Ap, tp                                                     =
#                                                             =
# PIN ORDER: 1 5                                              =
# =============================================================
def append_SENSE(ofile_name, subckt_name, R0, tl):
    f = open(ofile_name, "a")
    f.write("\n") # New line
    f.write(".param k1={10**(-Ap/20)} \n")
    f.write(".param k2={(5/k1+sqrt(25/(k1**2)-16))/2} \n")
    f.write(".SUBCKT " + subckt_name + " 1 5" + "\n")
    f.write("T1 1 0 2 0 " + "Z0=" + str(R0) + " td={" + str(tl) + "+tp}" + "\n")
    f.write("R1 2 0 " + str(R0) + "\n")
    f.write("F1 0 2 V1 {k2}" + "\n")
    f.write("V1 3 0 0" + "\n")
    f.write("E1 4 3 2 0 {k2}" + "\n")
    f.write("R2 4 5 " + str(R0) + "\n")
    f.write(".ENDS \n")
    f.close()

# =============================================================
# Function: append_DLINE                                      =
# Description: adds a complete subcircuit for the delay line  =
# Currently without feedforward capacitor or SENSE subcircuit =
#                                                             =
#                                                             =
# PIN ORDER: Vi+ Vi- Vo+ Vo-                                  =
# =============================================================
def append_DLINE(ofile_name, subckt_name, idt_subckt_name, sense_subckt_name, R0):
    f = open(ofile_name, "a")
    f.write("\n") # New line
    f.write(".SUBCKT " + subckt_name + " Vip Vim Vop Vom")
    f.write("\n") # New line
    f.write("R1 1 0 " + str(R0) + "\n")
    f.write("X1 Vip Vim 1 2 " + idt_subckt_name + "\n")
    f.write("X2 2 3 " + sense_subckt_name + "\n")
    f.write("X3 Vop Vom 3 4 " + idt_subckt_name + "\n")
    f.write("R2 4 0 " + str(R0) + "\n")
    f.write("C1 Vip Vop {Cf}" + "\n")
    f.write(".ENDS \n")
    f.close()
    
    

# =============================================================
# Function: foster_tan_gen                                    =
# Description: Generates the LC pairs in series according to  =
# the foster theorem given R0 and w0 of the tangent function  =
#                                                             =
#                                                             =
# Returns: tuple of (L,C)                                     =
# =============================================================
def foster_tan_gen(n, R0, w0):
    C = np.zeros(n)
    L = np.zeros(n)
    if (n==1):
        C[0] = (2/3)/(R0*w0)
        L[0] = (3/2)*R0/w0
    
    if (n==2):
        C[0] = ((2*8)/(3*7))/(R0*w0)
        L[0] = ((3*7)/(2*8))*R0/w0
        
        C[1] = ((2*8)/(7*5))/(R0*w0)
        L[1] = ((7*5)/(9*2*8))*R0/w0
    
    if (n==3):
        C[0] = ((8*24*2)/(3*15*11))/(R0*w0)
        L[0] = ((3*15*11)/(8*24*2))*R0/w0
        
        C[1] = ((8*16*2)/(5*7*11))/(R0*w0)
        L[1] = ((5*7*11)/(8*16*2*9))*R0/w0
        
        C[2] = ((24*16*2)/(21*9*11))/(R0*w0)
        L[2] = ((21*9*11)/(24*16*2*25))*R0/w0
        
    if (n==4):
        C[0] = ((8*24*48*2)/(3*15*35*15))/(R0*w0)
        L[0] = ((3*15*35*15)/(8*24*48*2))*R0/w0
        
        C[1] = ((8*16*40*2)/(5*7*27*15))/(R0*w0)
        L[1] = ((5*7*27*15)/(8*16*40*2*9))*R0/w0
        
        C[2] = ((24*16*24*2)/(21*9*11*15))/(R0*w0)
        L[2] = ((21*9*11*15)/(24*16*24*2*25))*R0/w0
        
        C[3] = ((48*40*24*2)/(45*33*13*15))/(R0*w0)
        L[3] = ((45*33*13*15)/(48*40*24*2*49))*R0/w0
    
    if (n==5):
        C[0] = ((8*24*48*80*2)/(3*15*35*63*19))/(R0*w0)
        L[0] = ((3*15*35*63*19)/(8*24*48*80*2))*R0/w0
        
        C[1] = ((8*16*40*72*2)/(5*7*27*55*19))/(R0*w0)
        L[1] = ((5*7*27*55*19)/(8*16*40*72*2*9))*R0/w0
        
        C[2] = ((24*16*24*56*2)/(21*9*11*39*19))/(R0*w0)
        L[2] = ((21*9*11*39*19)/(24*16*24*56*2*25))*R0/w0
        
        C[3] = ((48*40*24*32*2)/(45*33*13*15*19))/(R0*w0)
        L[3] = ((45*33*13*15*19)/(48*40*24*32*2*49))*R0/w0
        
        C[4] = ((80*72*56*32*2)/(77*65*45*17*19))/(R0*w0)
        L[4] = ((77*65*45*17*19)/(80*72*56*32*2*81))*R0/w0
        
    return (L,C)

# =============================================================
# Function: mittag_tan_gen                                    =
# Description: Generates the LC pairs in series according to  =
# the Mittag-Leffler's theorem given the                      =
# R0 and w0 of the tangent function                           =
#                                                             =
# Returns: tuple of (L,C)                                     =
# =============================================================
def mittag_tan_gen(n,R0,w0):
    C = np.zeros(n)
    L = np.zeros(n)    
    m = 0 # Counter
    while(m<n):
        mt = 2*m+1
        L[m] = 4*R0/(np.pi*w0*mt**2)
        C[m] = np.pi/(4*R0*w0)
        m = m+1
    return (L,C)

# =============================================================
# Function: foster_csc_gen                                    =
# Description: Generates the LC pairs in series according to  =
# the foster theorem given R0 and w0 of the cosecant function =
#                                                             =
#                                                             =
# Returns: tuple of (L,C)                                     =
# =============================================================
def foster_csc_gen(n,R0,w0):
    C = np.zeros(2*n+1)
    L = np.zeros(2*n+1)
    
    if (n==1):
        C[0] = -(4/3)/(R0*w0)
        L[0] = -(3/4)*R0/w0
        
        C[1] = (16/5)/(R0*w0)
        L[1] = np.inf # The capacitor is alone
        
        C[2] = ((4*4)/(5*3))/(R0*w0)
        L[2] = ((5*3)/(4*4*4))*R0/w0
        
    if (n==2):
        C[0] = -((4*8)/(7*3))/(R0*w0)
        L[0] = -((7*3)/(4*8))*R0/w0
        
        C[1] = -((4*8)/(7*5))/(R0*w0)
        L[1] = -((7*5)/(4*8*9))*R0/w0
        
        C[2] = ((4*4*16)/(9*1*9))/(R0*w0)
        L[2] = np.inf # The capacitor is alone
        
        C[3] = ((4*12*4)/(9*3*5))/(R0*w0)
        L[3] = ((9*3*5)/(4*4*12*4))*R0/w0
        
        C[4] = ((4*16*12)/(9*15*7))/(R0*w0)
        L[4] = ((9*15*7)/(4*16*12*16))*R0/w0
    
    if (n==3):
        C[0] = -((4*8*24)/(11*3*15))/(R0*w0)
        L[0] = -((11*3*15)/(4*8*24))*R0/w0
        
        C[1] = -((4*8*16)/(5*7*11))/(R0*w0)
        L[1] = -((5*7*11)/(4*8*16*9))*R0/w0
        
        C[2] = -((4*24*16)/(21*9*11))/(R0*w0)
        L[2] = -((21*9*11)/(4*24*16*25))*R0/w0
        
        C[3] = ((4*4*16*36)/(13*1*9*25))/(R0*w0)
        L[3] = np.inf # The capacitor is alone
        
        C[4] = ((4*4*12*32)/(13*3*5*21))/(R0*w0)
        L[4] = ((13*3*5*21)/(4*4*4*12*32))*R0/w0
        
        C[5] = ((4*16*12*20)/(13*15*7*9))/(R0*w0)
        L[5] = ((13*15*7*9)/(16*4*16*12*20))*R0/w0
        
        C[6] = ((4*36*32*20)/(13*35*27*11))/(R0*w0)
        L[6] = ((13*35*27*11)/(36*4*36*32*20))*R0/w0
    
    if (n==4):
        C[0] = -((8*24*48*4)/(3*15*35*15))/(R0*w0)
        L[0] = -((3*15*35*15)/(8*24*48*4))*R0/w0
        
        C[1] = -((8*16*40*4)/(5*7*27*15))/(R0*w0)
        L[1] = -((5*7*27*15)/(8*16*40*4*9))*R0/w0
        
        C[2] = -((24*16*24*4)/(21*9*11*15))/(R0*w0)
        L[2] = -((21*9*11*15)/(24*16*24*4*25))*R0/w0
        
        C[3] = -((48*40*24*4)/(45*33*13*15))/(R0*w0)
        L[3] = -((45*33*13*15)/(48*40*24*4*49))*R0/w0
        
        C[4] = ((4*4*16*36*64)/(17*1*9*25*49))/(R0*w0)
        L[4] = np.inf # The capacitor is alone
        
        C[5] = ((4*12*32*60*4)/(3*5*21*45*17))/(R0*w0)
        L[5] = ((3*5*21*45*17)/(4*12*32*60*4*4))*R0/w0
        
        C[6] = ((16*12*20*48*4)/(15*7*9*33*17))/(R0*w0)
        L[6] = ((15*7*9*33*17)/(16*12*20*48*4*16))*R0/w0
        
        C[7] = ((36*32*20*28*4)/(35*27*11*13*17))/(R0*w0)
        L[7] = ((35*27*11*13*17)/(36*32*20*28*4*36))*R0/w0
        
        C[8] = ((64*60*48*28*4)/(63*55*39*15*17))/(R0*w0)
        L[8] = ((63*55*39*15*17)/(64*60*48*28*4*64))*R0/w0
    
    if (n==5):
        C[0] = -((8*24*48*80*4)/(3*15*35*63*19))/(R0*w0)
        L[0] = -((3*15*35*63*19)/(8*24*48*80*4))*R0/w0
        
        C[1] = -((8*16*40*72*4)/(5*7*27*55*19))/(R0*w0)
        L[1] = -((5*7*27*55*19)/(8*16*40*72*4*9))*R0/w0
        
        C[2] = -((24*16*24*56*4)/(21*9*11*39*19))/(R0*w0)
        L[2] = -((21*9*11*39*19)/(24*16*24*56*4*25))*R0/w0
        
        C[3] = -((48*40*24*32*4)/(45*33*13*15*19))/(R0*w0)
        L[3] = -((45*33*13*15*19)/(48*40*24*32*4*49))*R0/w0
        
        C[4] = -((80*72*56*32*4)/(77*65*45*17*19))/(R0*w0)
        L[4] = -((77*65*45*17*19)/(80*72*56*32*4*81))*R0/w0
        
        
        C[5] = ((4*4*16*36*64*100)/(21*1*9*25*49*81))/(R0*w0)
        L[5] = np.inf # The capacitor is alone
        
        C[6] = ((4*12*32*60*96*4)/(3*5*21*45*77*21))/(R0*w0)
        L[6] = ((3*5*21*45*77*21)/(4*12*32*60*96*4*4))*R0/w0
        
        C[7] = ((16*12*20*48*84*4)/(15*7*9*33*65*21))/(R0*w0)
        L[7] = ((15*7*9*33*65*21)/(16*12*20*48*84*4*16))*R0/w0
        
        C[8] = ((36*32*20*28*64*4)/(35*27*11*13*45*21))/(R0*w0)
        L[8] = ((35*27*11*13*45*21)/(36*32*20*28*64*4*36))*R0/w0
        
        C[9] = ((64*60*48*28*36*4)/(63*55*39*15*17*21))/(R0*w0)
        L[9] = ((63*55*39*15*17*21)/(64*60*48*28*36*4*64))*R0/w0
        
        C[10] = ((100*96*84*64*36*4)/(99*91*75*51*19*21))/(R0*w0)
        L[10] = ((99*91*75*51*19*21)/(100*96*84*64*36*4*100))*R0/w0
    
    return (L,C)

# =============================================================
# Function: mittag_csc_gen                                    =
# Description: Generates the LC pairs in series according to  =
# the Mittag-Leffler's theorem given the                      =
# R0 and w0 of the cosecant function                          =
#                                                             =
# Returns: tuple of (L,C)                                     =
# =============================================================
def mittag_csc_gen(n,R0,w0):
    C = np.zeros(2*n+1)
    L = np.zeros(2*n+1)    
    
    C[0] = np.pi/(w0*R0)
    L[0] = np.inf # The capacitor is alone
    
    m = 1 # Counter
    while(m<=n):
        L[2*m] = R0/(2*np.pi*w0*m**2)
        C[2*m] = np.pi/(2*R0*w0)
        
        L[2*m-1] = -2*R0/(np.pi*w0*(4*m**2-4*m+1))
        C[2*m-1] = -np.pi/(2*R0*w0);
        m = m+1
    return (L,C)

# =============================================================
# Function: impedance_LC                                      =
# Description: Evaluates the impedance of the given LC pair   =
# it only works with SCALAR L and C values                    =
#                                                             =
#                                                             =
# Returns: vector Z of impedance over w                       =
# =============================================================
def impedance_LC(w, L, C):
    Z = (1/C)*1j*w/((1j*w)**2+1/(L*C))
    return Z

# =============================================================
# Function: testbench_gen                                     =
# Description: Generates the testbench file with a sample     =
# circuit and simulation directives in the chosen file        =
#                                                             =
#                                                             =
# =============================================================
def testbench_gen(ofile_name, blocks_name, dline_subckt_name, f0, tl):
    f = open(ofile_name, "a")
    f.write(".include " + blocks_name + "\n")
    f.write(".param Ap=0 \n")
    f.write(".param tp=0 \n")
    f.write(".param Cf=10f \n")
    f.write("V1 Vi 0 AC 1 PULSE(0 5 0 1n 1n 10n) DC 0 \n")
    f.write("R1 Vi Di 50 \n")
    f.write("X1 Di 0 Dop Dom " + dline_subckt_name + "\n")
    f.write("R2 Dop Dom 50 \n")
    f.write("\n")
    f.write(".AC LIN 10000 " + str(f0/100) + " " + str(2*f0) + "\n")
    f.write("*.TRAN 1N " + str(5*tl) + " uic")
    f.close()
    
    