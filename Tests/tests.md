# SIMsalabim Tests

The literature on semiconductor physics provides many analytical expressions that hold in certain limits that can be reproduced by SimSS/ZimT. These expressions therefore serve as highly suitable validation tests for the codebase. We chose the tests in such a way that most of the functionality is covered by the test cases. In this document, we describe the tests, while we provide the test results in graph form and the used simulation parameters in the test-specific folders. In all tests, we specify in what limits SimSS/ZimT should produce the same results as the analytical expression it is tested against.

# How to run the tests
There is a script to automatically asses the SIMsalabim project that is located in the 'Tests' folder. These tests test the code against physics that apply in certain regimes of operation. For every test case there is a folder that contains the input parameters required for the test to run. In this folder, a plot showing the result of the test will be generated during the test.

A Python 3 interpretor is installed by default on most operating systems, however on Windows it should first be installed (for example from the Microsoft Store).

The test script also requires a few non-default python packages. To install these simply run the command `pip install numpy scipy matplotlib pandas` in the terminal on your computer. 

The test script also requires the free Pascal compiler if not installed already. For details on how to install the free Pascal compiler see README.txt.

To run the automated tests, run `python run_tests.py` in a terminal from the 'Tests' folder.


# SimSS tests

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
In this test, we assess the ability of SimSS to simulate a diode with a single trap level (situated mid-gap). As a reference, we also show the simulated current-voltage curve without traps. To obtain the reference, either put N_t_bulk to zero in the layer parameter file, or run SimSS with -l1.N_t_bulk 0. The inclusion of traps, as expected, strongly reduces the current at low voltages. The trap-filled-limit (V<sub>TFL</sub>) is obtained as outlined in V.M. Le Corre, E.A. Duijnstee, O.  El Tambouli, J.M. Ball, H.J. Snaith, J. Lim, and L.J.A. Koster, ACS Energy Lett. **6**,  1087 (2021): two tangents are used (dashed lines in graph), one in the steepest part of the curve and one in the quadratic regime. They cross at 1.2 V, which agrees very nicely with the anticipated value for V<sub>TFL</sub> of 1 V which is based on the input parameters.


# ZimT tests

## Test 6: Open-circuit voltage in steady-state
ZimT can be used to solve for the open-circuit voltage of a solar cell, in transient cases as well as steady-state. The open-circuit in steady-state is compared with an analytical solution given in L.J.A. Koster, V.D. Mihailetchi, R. Ramaker, and P.W.M. Blom, Appl. Phys. Lett. **86**, 123509 (2005). The agreement is excellent.


## Test 7: RC-time
This tests ZimT by simulating an RC-circuit. The parameters define a capacitor with a series resistance. The voltage is changed from 0V (initial condition) to 1V. The RC-time, based on the parameters, is 5 μs. A fit to the simulated data yields an RC time of 5.36 μs, which is satisfactory.


## Test 8: Transient Photo-voltage (TPV)
This test ZimT's ability to simulate the transient decay of the open-circuit voltage (Voc) upon a small change in the light intensity. We start at a generation rate of electron-hole pairs of G_ehp = 1.2E26 m<sup>-3</sup>s<sup>-1</sup> at time zero, followed by a quick reduction to Gehp = 1.0E26 m<sup>-3</sup>s<sup>-1</sup>. The Voc will thus decay a little bit. The decay rate can be estimated based on A. Rahimi Chatri, S. Torabi, V.M. Le Corre, and L.J.A. Koster, ACS Appl. Mater. Interfaces **10**, 12013 (2018). Their formula (Eq. 9) predicts a decay time of 7.09 μs for the parameters used in this test. Fitting to the simulated data yields a very similar decay time.

## Test 9: Transient Photo-current with bulk traps (TPC)
This test assesses ZimT's ability to simulate the emptying of bulk traps: The simulation starts with a solar cell at 0 V under illumination. Quickly, the light is switched off and free carriers are swept out of the device (the mobilities are very high). Trapped electrons will be slowly released from the traps and this is what generates the current in the ms-regime. The parameters are chosen such that the detrapping time should be 0.050 s. Fitting to ZimT's results indeed yields a decay time sufficiently close to this.

## Test 10: Transient Photo-current with interface traps (TPC)
This test assesses ZimT's ability to simulate the emptying of interface traps between adjacent layers. The input parameters are chosen such that de-trapping should have a lifetime of 25 ms. Fitting to ZimT's result yields a lifetime of 25.8 ms, i.e. sufficiently close to the input value.

## Test 11: Transient Photo-current for estimating Urbach energy
This tests ZimT's ability to simulate detrapping of charges from bulk traps that are distributed exponentially from the conduction band defined by an Urbach energy of 100 meV. The Urbach energy can be fitted by transforming time and current output as described in MacKenzie, et al., J. Phys. Chem. C **117 (24)**, 12407-12414 (2013). First we run the steady state solver under illumination, effectively filling the traps present. Then we switch off the light and track detrapping charges through the current. We then manipulate the current and time values to estimate the DOS(E) and fit the Urbach energy, which yields 95 meV, sufficiently close to the input value.

## Test 12: Transient Photo-current for estimating Urbach energy
This tests ZimT's ability to simulate detrapping of charges from interface traps that are distributed exponentially from the valence band, defined by an Urbach of 70 meV. We defined 10 layers for this purpose. The Urbach energy can be fitted by transforming time and current output as described in MacKenzie, et al., J. Phys. Chem. C **117 (24)**, 12407-12414 (2013). First we run the steady state solver under illumination, effectively filling the traps present. Then we switch off the light and track detrapping charges through the current. We then manipulate the current and time values to estimate the DOS(E) and fit the Urbach energy, which yields 66 meV, i.e. sufficiently close to the input value.

## Test 13: Generation profile using the TransferMatrix method
This tests the creation of a generation profile using the TransferMatrix method based on the layers in the device structure, which are characterized by width and n,k properties of the layer material. The resulting generation profile is compared to a generation profile for the same device calculated using the script from the work of G.F. Burkhard, E.T. Hoke and M.D. McGehee, Adv. Mater. **22**, 3293 (2010).












