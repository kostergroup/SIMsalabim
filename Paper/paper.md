---
title: 'SIMsalabim: An open-source drift-diffusion simulator for semiconductor devices'
tags:
    - semiconductor
    - solar cell
    - drift-diffusion
    - diode
authors:
    - name: Marten Koopmans
      affiliation: 1
    - name: Vincent M. Le Corre
      affiliation: 1
    - name: L. Jan Anton Koster^[corresponding author]
      affiliation: 1
affiliations:
    - name: Zernike Institute for Advanced Materials, University of Groningen, Nijenborgh 4, 9747 AG Groningen, The Netherlands
      index: 1
date: 24 June 2021
bibliography: paper.bib
---

# Summary
The drift-diffusion technique is a well-established way of describing charge transport in semiconductors. While in some special cases it is possible to analytically solve the problem, typically numerical solutions are the best achievable result, especially in more complicated semiconductor devices, such as solar cells.
Numerical solutions can be found, however, and often prove to closely align with experimental observations.

SIMsalabim is a  software that allows researchers to simulate solar cells or other semiconductor devices to access parameters that might be hard to observe in an experiment. SIMsalabim can do both steady-state and time-dependent simulations including effects of illumination, mobile ions, recombination, trapping, and dielectric mismatch. While SIMsalabim is not restricted to the simulation of solar cells, some relevant solar cell experiments that can be simulated include transient photovoltage and current-voltage sweeps (including hysteresis). SIMsalabim can be used to simulate many different device types, but the user should make sure that the physical model implemented in SIMsalabim fits the device and operating regime that is simulated. Because the physical model is constantly evolving, we refer the reader to the manual for what is currently implemented.

# Statement of need
SIMsalabim is a numerical 1D drift-diffusion simulation software implemented in Pascal. Pascal enables a combination of speed offered by low-level languages (e.g. C) while enforcing code readability and simplicity. The code is therefore relatively accessible to the typical material scientists who only have experience in higher-level languages (e.g. Python, MATLAB). SIMsalabim allows setting of parameters either in a file or the command line, meaning that it is both scriptable and user-friendly. The scriptable nature allows for automated fitting of experiments, or scanning of extremely large parameter spaces [@LeCorre2021Feb].

Multiple purposes can be fulfilled using SIMsalabim. It can be used to help in designing experiments by allowing physics-based estimates of experimental parameters. 
Another use is checking the validity of a simplified physical model or formula for a wide range of different input parameters, representing a broad class of devices [@VincentTLs2019]. 
Due to the fast pace of research in fields such as that of perovskite solar cells, it is important that the user can perform adaptations to the codebase to include novel effects, without waiting on another entity to implement those.

SIMsalabim is currently used by multiple research groups in the solar cell community. Researchers interested in the modelling of semiconductors can use SIMsalabim for a combination of features: flexible device architecture, fully transient ions/charge-carriers/trapping-detrapping, open-circuit voltage tracking, rigorous recombination treatment, high numerical stability in extreme cases, and scriptable high-performance-computing friendly design.

As drift-diffusion modelling is a well established technique, there are more pieces of software that achieve similar modelling each with their own unique strengths and weaknesses. A few alternatives are: SCAPS [@BURGELMAN2000] (freeware), IonMonger [@Courtier2019] (open-source), SETFOS [@SETFOS] (proprietary), Driftfusion [@driftfusion] (open-source), and gpvdm [@gpvdm] (open-source).

# Recent projects involving SIMsalabim (incomplete)
- Identification of the dominant recombination process for perovskite solar cells based on machine learning [@LeCorre2021Feb]
- Understanding Dark Current-Voltage Characteristics in Metal-Halide Perovskite Single Crystals [@Elisabeth2021]
- Revealing Charge Carrier Mobility and Defect Densities in Metal Halide Perovskites via Space-Charge-Limited Current Measurements [@SpaceCharge2021]
- Toward Understanding Space-Charge Limited Current Measurements on Metal Halide Perovskites [@Elisabeth2020]
- 15.34% efficiency all-small-molecule organic solar cells with an improved fill factor enabled by a fullerene additive [@Friso2020]
- Effects of Fluorination on Fused Ring Electron Acceptor for Active Layer Morphology, Exciton Dissociation, and Charge Recombination in Organic Solar Cells [@Frisotwo2020]
- Photostability of Fullerene and Non-Fullerene Polymer Solar Cells: The Role of the Acceptor [@Nutifafa2019]
- Charge Transport Layers Limiting the Efficiency of Perovskite Solar Cells: How To Optimize Conductivity, Doping, and Thickness [@VincentTLs2019]

# References
