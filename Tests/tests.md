
# Introduction

The literature on semiconductor physics provides many analytical expressions that hold in certain limits that can be reproduced by SIMsalabim/ZimT. These expressions therefore serve as highly suitable validation tests for the codebase. We chose the tests in such a way that most of the functionality is covered by the test cases. In this document, we describe the tests, while we provide the test results in graph form and the used simulation parameters in the test-specific folders. In all tests, we specify in what limits SIMsalabim/ZimT should produce the same results as the analytical expression it is tested against.

# How to run the tests
Each folder contains the input files need to reproduce the test. Copy the input file(s) to either the SIMsalabim or ZimT folder. Note that ZimT needs an additional input file (tVG*.txt). Running the program generates output that can be compared with the sample output in the test folder. A graphical comparison of the latest version with the analytical results is also provided. 

# SIMsalabim tests

## Test 1: Photoconductivity of an insulator

This test is based on the work by Sokel and Hughes (R. Sokel and R.C. Hughes, J. Appl. Phys. **53**, 7414 (1982)). They provide a formula (their Eq. 49) for the current density in an insulator that is able to absorb light. It assumes that the electric field in the device is constant, i.e. there is a neglible density of electrons and holes. The file JV.dat contains the simulated and analytical result. The agreement is excellent for all voltages.

## Test 2: Space-charge-limited-current with field-dependent mobility
This test assesses the accuracy of the Poisson solver. It simulates a so-called electron-only diode (a diode where the injection of holes is blocked) with an electron mobility that depends on the electric field. The analytical result by Murgatroyd is used as a reference (see  P. N. Murgatroyd, J. Phys. D **3**, 151 (1970)). The analytical result by Murgatroyd neglects the effect of diffusion and, hence, the simulated current density should be lower than the analytical result, especially at low voltages. This is indeed the case: at high voltages, the simulation and analytical result overlap, while at low voltages the simulated current is significantly smaller.

## Test 3: Double-injection with very weak recombination

Rosenberg and Lampert (L.M. Rosenberg and M.A. Lampert, J. Appl. Phys. **41**, 508 (1970)) derive two limiting cases for double injection into an insulator, corresponding to injection of electrons and holes from Ohmic contacts. The semiconductor/insulator does not contain any trapping and recombination is of the direct type. The plasma limit is obtained when the recombination is very weak. The simulated and analytical results are quite close, especially for high voltages. Note, the analytical derivation does not consider any built-in voltage (and no diffusion) which is present in the simulations. As a result, the simulation should be below the analytical result, specifically at low voltages where diffusion and built-in voltage are relevant. 

## Test 4: Space-charge-limited diode with diffusion regime
Test 4 is another test of the space-charge-limited regime. This time, the mobility is taken as a constant (i.e. not dependent on the electric field) and we compare the high- and low-field behaviour with the analytical results as descirbed in J.A. Röhr, T. Kirchartz, and J. Nelson, J. Phys.: Condens. Matter **29**, 205901 (2017).
At high voltages, the simulation should approach the Mott-Gruney law: This is indeed the case, despite the fact that the simulation includes a built-in voltage. At low voltages, the current is linear as a consequence of background charge density that stems from the contacts. This regime is described by the moving-electrode equation. The agreement in both limits is excellent.

## Test 5: Space-charge-limited-diode with trapping
In this test, we assess the ability of SIMsalabim to simulate a diode with a single trap level (situated mid-gap). As a reference, we also show the simulated current-voltage curve without traps. To obtain the reference, either put Bulk_tr to zero in the parameter file, or run SIMsalabim with -Bulk_tr 0. The inclusion of traps, as expected, strongly reduces the current at low voltages. The trap-filled-limit (V<sub>TFL</sub>) is obtained as outlined in V.M. Le Corre, E.A. Duijnstee, O.  El Tambouli, J.M. Ball, H.J. Snaith, J. Lim, and L.J.A. Koster, ACS Energy Lett. **6**,  1087 (2021): two tangents are used (dashed lines in graph), one in the steepest part of the curve and one in the quadratic regime. They cross at 1.2 V, which agrees very nicely with the anticipated value for V<sub>TFL</sub> of 1 V which is based on the input parameters.


# ZimT tests

## Test 6: Open-circuit voltage in steady-state
ZimT can be used to solve for the open-circuit voltage of a solar cell, in transient cases as well as steady-state. The open-circuit in steady-state is compared with an analytical solution given in L.J.A. Koster, V.D. Mihailetchi, R. Ramaker, and P.W.M. Blom, Appl. Phys. Lett. **86**, 123509 (2005). The agreement is excellent.


## Test 7: RC-time
This tests ZimT by simulating an RC-circuit. The parameters define a capacitor with a series resistance. The voltage is changed from 0V (initial condition) to 1V. The RC-time, based on the parameters, is 5 μs. A fit to the simulated data yields an RC time of 5.36 μs, which is satisfactory.


## Test 8: Transient Photo-voltage (TPV)
This test ZimT's ability to simulate the transient decay of the open-circuit voltage (Voc) upon a small change in the light intensity. We start at a generation rate of electron-hole pairs of Gehp = 1.2E26 m<sup>-3</sup>s<sup>-1</sup> at time zero, followed by a quick reduction to Gehp = 1.0E26 m<sup>-3</sup>s<sup>-1</sup>. The Voc will thus decay a little bit. The decay rate can be estimated based on A. Rahimi Chatri, S. Torabi, V.M. Le Corre, and L.J.A. Koster, ACS Appl. Mater. Interfaces **10**, 12013 (2018). Their formula (Eq. 9) predicts a decay time of 7.09 μs for the parameters used in this test. Fitting to the simulated data yield a decay time of 6.96 μs.









