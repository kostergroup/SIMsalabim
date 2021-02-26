# SIMsalabim

SIMsalabim: A 1D drift-diffusion simulator for semiconductor devices (LEDs, solar cells, diodes, organics, perovskites). 

## Table of Contents
1. [Introduction](#introduction)
2. [Copyright and license](#copyright-and-license)
3. [Instructional videos](#instructional-videos)
4. [How to cite](#how-to-cite)
5. [How to compile](#how-to-compile)
6. [Scientific publications based on SIMsalabim](#scientific-publications-based-on-simsalabim)

## Introduction

SIMsalabim can be used to simulate current-voltage (JV) characteristics of semiconductor devices. It includes the effects of generation, recombination and trapping of electrons and holes, the effect of ions and dopants, and self-consistently solves the electric field that results from all charged species.

## Copyright and license

SIMsalabim is licensed under the GNU Lesser General Public Licence version. The details of this licence can be found in the files COPYING.txt and COPYING.LESSER.txt. 
Several authors (all from the University of Groningen) have contributed to the code: 
- Dr T.S. (Tejas) Sherkar
- V.M. (Vincent) Le Corre
- M. (Marten) Koopmans
- F. (Friso) Wobben
- Prof. Dr. L.J.A. (Jan Anton) Koster

## Instructional videos

In order to help new users to get the most out of SIMsalabim, we have recorded a few instructional videos:

1 [Installation and basic use](https://youtu.be/0mtpGJMnbFE "Installation and basic use")


## How to cite

The essentials of this code have been outlined in several publications. 
The original paper showing how this code can be used for organic solar cells is: 
- L.J.A. Koster, E. C. P. Smits, V. D. Mihailetchi, and P. W. M. Blom, Device model for the operation of polymer/fullerene bulk heterojunction solar cells, Phys. Rev. B **72**, 085205 (2005).

For the use of perovskite solar cells, the following paper offers a detailed description of the modelling assumptions:
- T.S. Sherkar, C. Momblona, L. Gil-Escrig, J. Ávila, M. Sessolo, H. Bolink, and L.J.A. Koster, Recombination in Perovskite Solar Cells: Significance of Grain Boundaries, Interface Traps and Defect Ions, ACS Energy Lett. **2**, 1214 (2017).

## Pre-compiled binaries

SIMsalabim comes as a pre-compiled binary for WIN or Linux (64 bits), see folder 'Binaries'. This avoids having to download and install the FPC compiler and compiling the code which is described below. There is no guarantee that these will work, but it is worth a try. Simply copy the binary to the main SIMsalabim folder. The Linux one might require a change of permissions to run, something like:
<pre><code>
    chmod +x SIMsalabim
 </code></pre>


## How to compile

1) SIMsalabim uses fpc (free pascal compiler)
- download and install fpc from https://www.freepascal.org/
- You'll probably want to use (install?) some form of integrated development environment (IDE). Popular IDEs include Lazarus and Geany.

2) Compile SIMsalabim (note: on Linux this is case-sensitve):
- open a terminal window and go to the SIMsalabim/Units folder
- within Units, compile the units in the following order:
<pre><code>
    fpc TypesAndConstants.pp
    fpc NumericalUtils.pp
    fpc InputOutputUtils.pp
 </code></pre>
- now go back to the SIMsalabim parent folder and compile SIMsalabim:
 <pre><code>
 fpc SIMsalabim
  </code></pre>

You may see this warning message:
/usr/bin/ld.bfd: warning: link.res contains output sections; did you forget -T?
Just ignore this.

3) Running SIMsalabim:
- all parameters are specified in device_parameters.txt
- on Linux run SIMsalabim by entering in the terminal: 
<pre><code>
./SIMsalabim
 </code></pre>

on windows:
<pre><code>
SIMsalabim.exe
 </code></pre>

- all parameters listed in device_parameters.txt can also be changed via the command line. This will override the respective parameter value in device_parameters.txt. Simply add the parameter name and value to the command line after a dash (-).
Example: change of thickness (L) and JV output file (JV_file):
<pre><code>
./SIMsalabim -L 345E-9 -JV_file anotherJV.dat
</code></pre>
- multiple output files will be generated (see device_parameters) and log file.

4) Other remarks
- in general, input files and the log file have extension .txt, whereas output files that can be used for plotting (J-V curves, for example) have extension .dat.
- all changes are shown and commented on in the file Docs/Change_log.txt. Here we document not just the changes, but also explain and motivate some of the choices in naming, physical models, etc. 

## Scientific publications based on SIMsalabim

List of publications (not complete):

- V.M. Le Corre, T.S. Sherkar, M. Koopmans, and L.J.A. Koster, Identification of the Dominant Recombination Process for Perovskite Solar Cells Based on Machine Learning, Cell Rep. Phys. Sci. 2, 100346 (2021).

- Y. Firdaus, V.M. Le Corre, S. Karuthedath, W. Liu, A. Markina, W. Huang, S. Chattopadhyay, M.M. Nahid, M.I. Nugraha, A. Seitkhan, A. Basu, Y. Lin, I. McCulloch, H. Ade, J. Labram, F. Laquai, D. Andrienko, L.J.A. Koster, and T.D. Anthopoulos, Long-range exciton diffusion in molecular non-fullerene acceptors, Nature Comm. 11, 5220 (2020).

- D. Hu, Q. Yang, H. Chen, F. Wobben, V.M. Le Corre, R. Singh, L.J.A. Koster, Z. Kan, Z. Xiao, and S. Lu, 15.34% Efficiency All-Small-Molecule Organic Solar Cells with Improved Fill Factor Enabled by a Fullerene Additive, Energy Environ. Sci. 13, 2134 (2020).

- L. Hou, J. Lv, F. Wobben, V.M. Le Corre, H. Tang, R. Singh, M. Kim, F. Wang, H. Sun, W. Chen, Z. Xiao, M. Kumar, T. Xu, W. Zhang, I. McCulloch, T. Duan, H. Xie, L.J.A. Koster, S. Lu, and Z. Kan, Effects of Fluorination on Fused Ring Electron Acceptor for Active Layer Morphology, Exciton Dissociation, and Charge Recombination in Organic Solar Cells, ACS Appl. Mater. Interfaces 12, 56231 (2020).

- E.A. Duijnstee, V.M. Le Corre, M.B. Johnston, L.J.A. Koster, J. Lim, and H.J. Snaith, Understanding dark current-voltage characteristics in metal-halide perovskite single crystals, Phys. Rev. Appl. 15, 014006 (2021).

- D. Neher, J. Kniepert , A. Elimelech, and L.J.A. Koster, A new Figure of Merit for Organic Solar Cells with Transport-limited Photocurrents, Sci. Rep. 6, 24861 (2016).

- S. Shao, M. Abdu-Aguye, T.S. Sherkar, H.-H. Fang, G. ten Brink, B.J. Kooi, L.J.A. Koster, and M.A. Loi, The effect of the microstructure on trap-assisted recombination and light soaking phenomenon in hybrid perovskite solar cells, Adv. Funct. Mater. 26, 8094 (2016).

- T.S. Sherkar, C. Momblona, L. Gil-Escrig, H.J. Bolink, and L.J.A. Koster, Improving perovskite solar cells: Insights from a validated device model, Adv. Energy Mater. 1602432 (2017).
 
- T.S. Sherkar, C. Momblona, L. Gil-Escrig, J. Ávila, M. Sessolo, H. Bolink, and L.J.A. Koster, Recombination in Perovskite Solar Cells: Significance of Grain Boundaries, Interface Traps and Defect Ions, ACS Energy Lett. 2, 1214 (2017).

- V.M. Le Corre, A. Rahimi Chatri, N.Y. Doumon, and L.J.A. Koster, Charge carrier extraction in organic solar cells governed by steady-state mobilities, Adv. Energy Mater. 1701138 (2017).

- F.J.M. Colberts, M.M. Wienk, R. Heuvel, W. Li, V.M. Le Corre, L.J.A. Koster, and R.A.J. Janssen, Bilayer-Ternary Polymer Solar Cells Fabricated Using Spontaneous Spreading on Water, Adv. Energy Mater. 8, 1802197 (2018).

- N.Y. Doumon, M.V. Dryzhov, F.V. Houard, V.M. Le Corre, A. Rahimi Chatri, P. Christodoulis, and L.J.A. Koster, Photostability of Fullerene and Non-Fullerene Polymer Solar Cells: The Role of the Acceptor, ACS Appl. Mater. Interfaces 11, 8310 (2019).

- V.M. Le Corre,  M. Stolterfoht, L. Perdigón Toro, M. Feuerstein, C. Wolff, L. Gil-Escrig, H.J. Bolink, D. Neher, and L.J.A. Koster, Charge transport layers limiting the efficiency of perovskite solar cells: how to optimize conductivity, doping and thickness, ACS Appl. Energy Mater. 2, 6280 (2019). 

- E. A. Duijnstee, J.M. Ball, V.M. Le Corre, L.J.A. Koster, H.J. Snaith, and J. Lim, Towards Understanding Space-charge Limited Current Measurements on Metal Halide Perovskites, ACS Energy Lett. 5, 376 (2020).

- D. Bartesaghi and L. J. A. Koster, The effect of large compositional inhomogeneities on the performance of organic solar cells: A numerical study, Adv. Funct. Mater. 25, 2013 (2015).

- J. Kniepert, I. Lange, J. Heidbrink, J. Kurpiers, T. Brenner, L.J.A. Koster, and D. Neher, Effect of Solvent Additive on Generation, Recombination and Extraction in PTB7:PCBM Solar Cells: A conclusive Experimental and Numerical Simulation Study, J. Phys. Chem. C 119, 8310 (2015).

- D. Bartesaghi, I. del Carmen Pérez, J. Kniepert, S. Roland, M. Turbiez, D. Neher, and L.J.A. Koster, Competition between recombination and extraction of free carriers determines the fill-factor of organic solar cells, Nature Comm. 6, 7083 (2015).

- N.J. van der Kaap, I. Katsouras, K. Asadi, P.W.M. Blom, L.J.A. Koster, and D.M. de Leeuw, Charge transport in disordered semiconducting polymers driven by nuclear tunneling, Phys. Rev. B 93, 140206(R) (2016).

- D. Bartesaghi, M. Turbiez, and L. J. A. Koster, Charge Transport and Recombination in PDPP5T:[70]PCBM Organic Solar Cells: the Influence of Morphology, Org. Elec. 15, 3191 (2014).

- G. A. H. Wetzelaer, N. J. van der Kaap, L. J. A. Koster, and P. W. M. Blom, Quantifying bimolecular recombination in organic solar cells in steady-state, Adv. Energy Mater. 3, 1130 (2013).

- J. Kniepert, I. Lange, N. J. van der Kaap, L. J. A. Koster, and D. Neher, A conclusive view on charge generation, recombination and extraction in as-prepared and annealed P3HT:PCBM blends: a combined experimental-simulation work, Adv. Energy Mater. 4, 1301401 (2014).

- L. J. A. Koster, M. Kemerink, M. M. Wienk, K. Maturová, and R. A. J. Janssen, Quantifying bimolecular recombination losses in bulk heterojunction solar cells, Adv. Mater. 63, 1670 (2011).

- G. A. H. Wetzelaer, L. J. A. Koster, and P. W. M. Blom, Validity of the Einstein relation in disordered organic semiconductors, Phys. Rev. Lett. 107, 066605 (2011).

- L.J.A. Koster, S. E. Shaheen, and J. C. Hummelen, Pathways to a new efficiency regime for organic solar cells, Adv. Energy Mater. 2, 1246 (2012). Listed in Adv. Energy Mater. top-5 most accessed articles in May 2012
- D.J. Wehenkel, L.J.A. Koster, M. M. Wienk, and R. A. J. Janssen, Influence of injected charge carriers on photocurrents in polymer solar cells, Phys. Rev. B 85, 125203 (2012).
- M. Kuik, L. J. A. Koster, A. G. Dijkstra, G. A. H. Wetzelaer, and P. W. M. Blom, Non-radiative recombination losses in polymer light-emitting diodes, Org. Elec. 13, 969 (2012).

- M. Kuik, L. J. A. Koster, G. A. H. Wetzelaer, and P. W. M. Blom, Trap-assisted recombination in disordered organic semiconductors, Phys. Rev. Lett. 107, 256805 (2011).

- L. J. A. Koster, Charge carrier mobility in disordered organic blends for photovoltaics, Phys. Rev. B 81, 205318 (2010).

- M. M. Mandoc, W. Veurman, L. J. A. Koster, M. M. Koetse, J. Sweelssen, B. de Boer, and P. W. M. Blom, Charge transport in MDMO-PPV : PCNEPV all-polymer solar cells, J. Appl. Phys. 101, 104512 (2007).

- V. D. Mihailetchi, H. Xie, B. de Boer, L. J. A. Koster, and P. W. M. Blom, Charge Transport and Photocurrent Generation in Poly(3-hexylthiophene):Methanofullerene Bulk-Heterojunction Solar Cells, Adv. Funct. Mater. 16, 699 (2006).L. J. A. Koster, W. J. van Strien, W. J. E. Beek, and P. W. M. Blom, Device operation of conjugated polymer/zinc oxide bulk heterojunction solar cells, Adv. Funct. Mater. 17, 1297 (2007). 

- M. M. Mandoc, L. J. A. Koster, and P. W. M. Blom, Optimum charge carrier mobility in organic solar cells, Appl. Phys. Lett. 90, 133504 (2007).

- M. M. Mandoc, W. Veurman, L. J. A. Koster, B. de Boer, and P. W. M. Blom, Origin of the Reduced Fill Factor and Photocurrent in MDMO-PPV:PCNEPV All-Polymer Solar Cells, Adv. Funct. Mater. 17, 2167 (2007).J. D. Kotlarski, M. Lenes, L. J. A. Koster, L. H. Slooff, and P. W. M. Blom, Combined optical and electrical modeling of polymer:fullerene bulk heterojunction solar cells, J. Appl. Phys. 103, 084502 (2008).

- L. J. A. Koster, V. D. Mihailetchi, and P. W. M. Blom, Ultimate efficiency of polymer/fullerene bulk heterojunction solar cells, Appl. Phys. Lett. 88, 093511 (2006).

- M. Lenes, L. J. A. Koster, V. D. Mihailetchi, and P. W. M. Blom, Thickness dependence of the efficiency of polymer:fullerene bulk heterojunction solar cells, Appl. Phys. Lett. 88, 243502 (2006).

- L. J. A. Koster, V. D. Mihailetchi, and P. W. M. Blom, Bimolecular recombination in polymer/fullerene bulk heterojunction solar cells, Appl. Phys. Lett. 88, 052104 (2006).
 
 - V. D. Mihailetchi, H. Xie, B. de Boer, L. J. A. Koster, L. M. Popescu, J. C. Hummelen, and P. W. M. Blom, Origin of the enhanced performance in poly(3-exylthiophene):methanofullerene solar cells using slow drying, Appl. Phys. Lett. 89, 012107 (2006).
  
  - L. J. A. Koster, E. C. P. Smits, V. D. Mihailetchi, and P. W. M. Blom, Device model for the operation of polymer/fullerene bulk heterojunction solar cells, Phys. Rev. B 72, 085205 (2005).
  
  - L. J. A. Koster, V. D. Mihailetchi, R. Ramaker, and P. W. M. Blom, Light intensity dependence of open-circuit voltage of polymer:fullerene solar cells, Appl. Phys. Lett. 86, 123509 (2005).
   
   - L. J. A. Koster, V. D. Mihailetchi, H. Xie, and P. W. M. Blom, Origin of the light intensity dependence of the short-circuit current of polymer/fullerene solar cells, Appl. Phys. Lett. 87, 203502 (2005).






