unit DDTypesAndConstants;
{provides drift-diffusion-specific types and constants}

{
SIMsalabim:a 1D drift-diffusion simulator 
Copyright (c) 2021, 2022, 2023, 2024, S. Heester, Dr T.S. Sherkar, Dr V.M. Le Corre, Dr M. Koopmans,
F. Wobben, and Prof. Dr. L.J.A. Koster, University of Groningen
This source file is part of the SIMsalabim project.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License and 
the GNU Lesser General Public License along with this program.


The SIMsalabim project can be found on Github at https://github.com/kostergroup/SIMsalabim 
email:l.j.a.koster@rug.nl
surface mail: 
L.J.A. Koster
Zernike Institute for Advanced Materials
Nijenborgh 4, 9747 AG Groningen, the Netherlands
}

{$MODE OBJFPC} {force OBJFPC mode}

INTERFACE

USES TypesAndConstants; {provides a couple of types}

CONST
    DDTypesAndConstantsVersion = '5.13'; {version of this unit}
    defaultParameterFile = 'simulation_setup.txt'; {name of file with parameters}
    q = 1.6022e-19;  	{C} {elementary charge}
    k = 1.3807e-23;     {J/K} {Boltzmann's constant}
    h = 6.62606957e-34; {Js}  {Planck's constant}
    c = 2.99792458e8;	{m/s} {speed of light}
    eps_0 = 8.8542e-12; {F/m} {dielectric constant of vacuum}
    tolReal = 1e-6;     {tolerance: when are 2 floats equal?}

{Magic numbers, used in the code:}
    {Name                    Where                        What does it do?}
    threshold_err = 0.1; {IN: Find_Solar_Cell_Parameters, defines max relative error when displaying solar cell parameters}
    MaxInterpolationOrder = 3; {IN: Find_Solar_Cell_Parameters, maximum interpolation order in estimating Jsc, Voc}
	minExpData = 4; {IN: SimSS /  Read_Experimental_JV, min number of experimental JV points}
	MinGridPointsPerLayer = 5; {IN: Make_Grid, mininum number of grid points per layer}
    myDelims = [#0..' ', ';', '''', '"', '`', '*', '=']; {IN Tidy_up_parameterFile_Exit, note: not a decimal point or comma!}
    minErr = 1E-4; {IN: Calc_and_Output_Soalr_Cell_Parameters, lower bound to error in parameters}
    tab1 = 6; {IN: Calc_and_Output_Solar_Cell_Parameters, 1st position of tab in table with solar cell parameters}
    tab2 = 28; {IN: Calc_and_Output_Solar_Cell_Parameters, 2nd position of tab in table with solar cell parameters}
    tab3 = 50; {IN: Calc_and_Output_Solar_Cell_Parameters, 3rd position of tab in table with solar cell parameters}
    MinCountStatic = 5; {number of time steps we consider to be static}
    MinCountChangeSmall = 3; {IN: determine_convergence, used if loop hardly changes.}
    MinCountJSmall = 3; {IN: determine_convergence, used if we're simulating really small currents.}
    TolRomb = 0.01; {IN Calc_Dissociation, sets relative tolerance of Romberg integration}
    MaxRombIt = 15; {IN Calc_Dissociation, max. # of iterations in Romberg integration}
    LowerLimBraun = 0.01; {IN Calc_Dissociation, lower limit of integration in Braun model, should be non-zero}
    UpperLimBraun = 20; {IN Calc_Dissociation, upper limit of integration in Braun model }
    minDeltaLambda = 0.5E-9; {IN: Read_AM_From_File, Read_nk_Material_From_File, minimum spacing between wavelengths in supplied spectrum/n,k files}


TYPE 
    TProgram = (ZimT, SimSS); {which program are we using?}
    
    TJVList = ARRAY OF RECORD {for storing the JV curve in SimSS}
                        Vint, Vext : myReal;
                        Jint, Jext : myReal;
                        UpdateIons, Store, Use : BOOLEAN
                      END;

	TTrapArray = ARRAY[0..Max_NP + 1] OF ARRAY OF myReal;

	
	TTrapDistLayer = RECORD {stores the energy and distribution of trap levels in a layer}
						NLevels, cwe : INTEGER; {cwe: charge of trap when empty}
						en, Nt : Row; {energy and trap density of a level}
						nt0, pt0 : Row; {nt0 and pt0 detrapping for bulk}
						nt0_L, nt0_R, pt0_L, pt0_R : Row {for interfaces, we have a left and right version}
					 END;


    TRec = RECORD {for storing the linearization of f_ti and f_tb at a grid point}
			direct, bulk, int, {direct (band to band), bulk SRH, interface SRH recombination respectively}
            {the rest are auxiliary variables:}
			dir_cont_rhs, dir_cont_m, bulk_cont_rhs, bulk_cont_m,
			int_cont_lo, int_cont_up, int_cont_m, int_cont_rhs : vector;        
           END;	

    TState = RECORD {stores all variables that either define a state or change in time}
		    G_frac, tijd : myReal; {set by input}
		    dti, Vint, Vext, Jint, Jext, errJ : myReal; {derivative variables or result of simulation}
		    SimType, convIndex : INTEGER; 
		    V, Vgn, Vgp, n, p, nion, pion, Gm, 
		    Jn, Jp, Jnion, Jpion, JD, mun, mup,
		    gen, Lang, SRH, diss_prob, 
		    Ntb_charge, Nti_charge : vector; {charge density of bulk/interface trap}
		    f_tb, f_ti : TTrapArray; {f_tb/i: occupancy of bulk/interface trap}
		    f_ti_numer, f_ti_inv_denom : TTrapArray; {numerator and inverse of denominator of f_ti}
		    Rn, Rp : TRec;
		    UpdateIons : BOOLEAN;
		END;

    TFitMode = (linear, logarithmic);

	TIonicRegions = ARRAY OF RECORD 
						istart, ifinish : INTEGER;
						AvC : myReal;
					END;

    TStaticVars = RECORD {stores all parameters that are calculated from input, but don't change during the simulation}
		    NcLoc, ni, eps, h, x, nid, pid, E_CB, E_VB, orgGm, mu_n_ion, mu_p_ion : vector;
		    lid : ShortIntVector;
		    i0, i1 : intArray;
		    Ntb, Nti : ARRAY OF TTrapDistLayer; {array over Layers}
		    Ltot, Lgen, epsi, V0, VL, Vt, Vti : myReal;
		    NLayers, NJV : INTEGER;
		    Traps_int_poisson : BOOLEAN;
		    NegIonRegion, PosIonRegion : TIonicRegions
		END;												  
																

    TLayerParameters =  RECORD {all input parameters and ones that are directly derived from the input (e.g. booleans), per layer}
			    {order of variables: same as in input file, but grouped by type}
			    {Floating point:}
				L, eps_r, E_c, E_v, N_c, N_D, N_A, 
				mu_n, mu_p, gamma_n, gamma_p,
				nu_int_n, nu_int_p, N_t_int, E_t_int, C_n_int, C_p_int,
				N_anion, N_cation, mu_anion, mu_cation,
				G_ehp, P0, a, 
				k_f, preLangevin, k_direct,
				N_t_bulk, C_n_bulk, C_p_bulk, E_t_bulk : myReal; 
			    {integers:}		
				mobnDep, mobpDep, intTrapType,
				Use_gen_profile, thermLengDist, bulkTrapType : INTEGER; 
			    {derived booleans:}
				layerGen, fieldDepG, negIons, posIons, 
				ionsMayEnter, useLangevin,
				bulkTrapFromFile, intTrapFromFile : BOOLEAN;
			    {filenames and other strings:}
				layerFile, nkLayer, bulkTrapFile, intTrapFile : ANSISTRING
			END;																				


    TInputParameters =  RECORD {all input parameters and ones that are directly derived from the input (e.g. booleans)}
			    {order of variables: same as in input file, but grouped by type}
				lyr : ARRAY OF TLayerParameters; {all parameters for a layer}
			    {Floating point:}
				T, W_L, W_R, L_TCO, L_BE, lambda_min, lambda_max, 
				S_n_R, S_n_L, S_p_L, S_p_R, R_series, R_shunt, G_frac, tolPois, maxDelV,
				tolDens, couplePC, minAcc, maxAcc, grad, tolVint,
				Vpre, Vmin, Vmax, Vstep, Vacc, timeout, fitThreshold : myReal;
			    {integers:}
				Use_gen_profile, maxItPois, 
				NP, currDiffInt, maxItSS, maxItTrans, 
				failureMode, outputRatio, Vdist, Vscan : INTEGER; 
			    {derived booleans:}
				fixIons, ignoreNegDens, autoTidy, StoreVarFile, limitDigits, 
				preCond, untilVoc, useExpData, autoStop, pauseAtEnd : BOOLEAN;
			    {filenames and other strings:}
				nkSubstrate, nkTCO, nkBE, spectrum, 
				genProfile, logFile, tVGFile, 
				varFile, tJFile, expJV, JVFile, scParsFile	: ANSISTRING;
			    {Misc: }
				fitMode : TFitMode;				
			END;																				

	
	TSCPar = RECORD {for storing solar cell parameters}
			Jsc, Vmpp, MPP, FF, Voc, ErrJsc, ErrVmpp, ErrMPP, ErrFF, ErrVoc : myReal;
			calcSC, calcMPP, calcFF, calcOC : BOOLEAN
		END;

    TLinFt = RECORD {for storing the linearization of f_ti and f_tb at a grid point}
				f_ti_lo, f_ti_up, f_ti_m, f_ti_rhs : vector;
                f_tb_m : vector;
             END;							


IMPLEMENTATION

BEGIN 

END.
