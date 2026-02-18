unit HelperUnit;


{
This unit supplies some routines to help the transition from v4.5x to 5.00
 
Copyright (c) 2023, 2024, 2025, 2026, S. Heester, Dr T.S. Sherkar, Dr V.M. Le Corre, 
Dr M. Koopmans, F. Wobben, and Prof. Dr. L.J.A. Koster, University of Groningen
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
Nijenborgh 3, 9747 AG Groningen, the Netherlands
}

{$MODE OBJFPC} {force OBJFPC mode}
  
INTERFACE

{$UNITPATH ../../Units/} {first tell compiler where our own units are located}

USES TypesAndConstants, 
     InputOutputUtils, 
     StrUtils,
     SysUtils,
     DDTypesAndConstants,
     Math;

TYPE

	{this is version 4.56 of TInputParameters, with the addition of NJV and mob_ion_spec:}
   TInputParameters456 =  RECORD {all input parameters and ones that are directly derived from the input (e.g. booleans)}
			    {order of variables: same as in input file, but grouped by type}
			    {Floating point:}
				T, L, eps_r, CB, VB, Nc, n_0, p_0, 
				L_TCO, L_BE, lambda_min, lambda_max,
				mun_0, mup_0, gamma_n, gamma_p,
				W_L, W_R,
				Sn_R, Sn_L, Sp_L, Sp_R, Rseries, Rshunt, L_LTL, L_RTL, Nc_LTL, 
				Nc_RTL, doping_LTL, doping_RTL, mob_LTL, mob_RTL,  
				nu_int_LTL, nu_int_RTL, eps_r_LTL, eps_r_RTL, 
				CB_LTL, VB_LTL, CB_RTL, VB_RTL, CNI, CPI,
				mobnion, mobpion,
				Gehp, Gfrac, P0, a, 
				kf, Lang_pre, kdirect,
				Bulk_tr, St_L, St_r, GB_tr, 
				Cn, Cp, ETrapSingle, tolPois, maxDelV,
				tolDens, couplePC, minAcc, maxAcc, grad, TolVint,
				Vpre, Vmin, Vmax, Vstep, Vacc, timeout, rms_threshold : myReal; 
			    {integers:}		
				MaxItPois, 
				mob_n_dep, mob_p_dep,
				Use_gen_profile, ThermLengDist, mob_ion_spec, ion_red_rate, num_GBs, Tr_type_L, Tr_type_R, Tr_type_B, 
				NP, CurrDiffInt, MaxItSS, MaxItTrans, 
				FailureMode, OutputRatio, Vdistribution, Vscan, NJV, {Note: NJV is defined here. In 456, it's in TStaticVars}
				Pause_at_end : INTEGER; 
			    {derived booleans:}
				Field_dep_G, negIonsMove, posIonsMove, TLsGen,
				IonsInTLs, TLsTrap, UseLangevin,
				BulkTrapFromFile, IntTrapFromFile, 
				IgnoreNegDens, AutoTidy, StoreVarFile, LimitDigits, 
				PreCond, until_Voc, UseExpData, AutoStop : BOOLEAN;
			    {filenames and other strings:}
				nk_substrate, nk_TCO, nk_active, nk_BE, spectrum, 
				nk_LTL, nk_RTL, Gen_profile, BulkTrapFile, IntTrapFile, 
				log_file, tVG_file, 
				Var_file, tj_file, ExpJV, JV_file, scPars_file	: ANSISTRING;
			    {Misc: }
				rms_mode : TFitMode;
			END;

    {this is version 5.00 of TLayerParameters}
    TLayerParameters500 =  RECORD {all input parameters and ones that are directly derived from the input (e.g. booleans), per layer}
			    {order of variables: same as in input file, but grouped by type}
			    {Floating point:}
				L, eps_r, CB, VB, Nc, n_0, p_0, 
				mun_0, mup_0, gamma_n, gamma_p,
				nu_int, St, ETrapInt, CnInt, CpInt,
				CNI, CPI, mobnion, mobpion,
				Gehp, P0, a, 
				kf, Lang_pre, kdirect,
				Bulk_tr, CnBulk, CpBulk, ETrapBulk : myReal; 
			    {integers:}		
				MaxItPois, 
				mob_n_dep, mob_p_dep, Tr_type_Int,
				Use_gen_profile, ThermLengDist, Tr_type_B : INTEGER; 
			    {derived booleans:}
				layerGen, Field_dep_G, negIonsMove, posIonsMove, 
				IonsMayEnter, UseLangevin,
				BulkTrapFromFile, IntTrapFromFile : BOOLEAN;
			    {filenames and other strings:}
				layerFile, nk_layer, BulkTrapFile, IntTrapFile	: ANSISTRING
			END;	
    	
    {version 5.00 of TInputParameters, except for the addition of NJV and NLayers:}		
    TInputParameters500 =  RECORD {all input parameters and ones that are directly derived from the input (e.g. booleans)}
			    {order of variables: same as in input file, but grouped by type}
				lyr : ARRAY OF TLayerParameters500; {all parameters for a layer}
			    {Floating point:}
				T, W_L, W_R, L_TCO, L_BE, lambda_min, lambda_max, 
				Sn_R, Sn_L, Sp_L, Sp_R, Rseries, Rshunt, Gehp, Gfrac, tolPois, maxDelV,
				tolDens, couplePC, minAcc, maxAcc, grad, TolVint,
				Vpre, Vmin, Vmax, Vstep, Vacc, timeout, rms_threshold : myReal;
			    {integers:}
				NLayers, Use_gen_profile, MaxItPois, 
				NP, CurrDiffInt, MaxItSS, MaxItTrans, 
				FailureMode, OutputRatio, Vdistribution, NJV, Vscan, 
				Pause_at_end : INTEGER; 
			    {derived booleans:}
				FixIons, IgnoreNegDens, AutoTidy, StoreVarFile, LimitDigits, 
				PreCond, until_Voc, UseExpData, AutoStop : BOOLEAN;
			    {filenames and other strings:}
				nk_substrate, nk_TCO, nk_active, nk_BE, spectrum, 
				Gen_profile, log_file, tVG_file, 
				Var_file, tj_file, ExpJV, JV_file, scPars_file	: ANSISTRING;
			    {Misc: }
				rms_mode : TFitMode;				
			END;				
			
PROCEDURE Read_Parameters_456(parameterFile : STRING; VAR msg : ANSISTRING; VAR par : TInputParameters456; ProgName : String);
{Routine directly from 456 with slight modifications: no stv, etc.}			

procedure write_simulation_setup_500(parnew : TInputParameters500; fn, path, ProgName : string); 
{produces general setup file}

procedure write_layer_500(lyr : TLayerParameters500; path : STRING); 
{produces layer file for main absober}

PROCEDURE Tidy_Up_Parameter_Files(parameterFile, path : STRING; QuitWhenDone : BOOLEAN; CONSTREF par : TInputParameters500);
{This procedure cleans up the parameter files: the main file and the files with the parameters of each layer}


IMPLEMENTATION


PROCEDURE Read_Parameters_456(parameterFile : STRING; VAR msg : ANSISTRING; VAR par : TInputParameters456; ProgName : String);
{Routine directly from 456 with slight modifications: no stv, etc.}
VAR dumint : INTEGER; {a dummy integer variable}
    dumstr : STRING; {again, a dummy variable}
    inv : TEXT;
    ZimT, SimSS	: BOOLEAN; 
BEGIN
    {use 2 booleans to check if we're using ZimT or SimSS}
    ZimT:= (ProgName='ZimT');
    SimSS:= (ProgName='SimSS');
    
    IF NOT FileExists(parameterFile) {the file with input par. is not found}
		THEN Stop_Prog('Could not find file '+parameterFile+'.', EC_FileNotFound);
    ASSIGN(inv, parameterFile);
    RESET(inv);

	WITH par DO BEGIN
{**General**************************************************************************}
		Get_Float(inv, msg,'T',T);  {abs. temperature, K}
		Get_Float(inv, msg, 'L', L); {device thickness, m}
		Get_Float(inv, msg, 'eps_r', eps_r);  {relative dielectric constant}
		Get_Float(inv, msg, 'CB', CB);  {eV, conduction band edge}
		Get_Float(inv, msg, 'VB', VB);  {eV, valence band edge}
		Get_Float(inv, msg, 'Nc',Nc);  {effective DOS, m^-3}
		Get_Float(inv, msg, 'n_0', n_0);  {ionised n-doping density, m^-3}
		Get_Float(inv, msg, 'p_0', p_0);  {ionised p-doping density, m^-3}

{**Optics****************************************************************************}
		Get_Float(inv, msg, 'L_TCO', L_TCO); {m, thickness of the TCO. Set to 0 if layer is not used}
		Get_Float(inv, msg, 'L_BE', L_BE); {m, thickness of back electrode, must be >0}
		Get_String(inv, msg, 'nk_substrate', nk_substrate); {name of file with n,k values of substrate}
		Get_String(inv, msg, 'nk_TCO', nk_TCO); {name of file with n,k values of TCO}
		Get_String(inv, msg, 'nk_active', nk_active); {name of file with n,k values of active layer}
		Get_String(inv, msg, 'nk_BE', nk_BE); {name of file with n,k values of back electrode}
		Get_String(inv, msg, 'spectrum', spectrum); {name of file that contains the spectrum}
		Get_Float(inv, msg, 'lambda_min', lambda_min); {m, lower bound wavelength}
		Get_Float(inv, msg, 'lambda_max', lambda_max); {m, upper bound wavelength}

{**Mobilities************************************************************************}
		Get_Float(inv, msg, 'mun_0', mun_0); {electron zero-field mobility, m^2/Vs}
		Get_Float(inv, msg, 'mup_0', mup_0); {hole zere-field mobility, m^2/Vs}
		Get_Integer(inv, msg, 'mob_n_dep', mob_n_dep);{dependence of elec mobility, 0 : const. mob, 1 : field-dep}
		Get_Integer(inv, msg, 'mob_p_dep', mob_p_dep);  {dependence of hole mobility, 0 : const. mob, 1 : field-dep}
		Get_Float(inv, msg,'gamma_n', gamma_n); {field depedence of mobility, eV/(V/m)^0.5}
		Get_Float(inv, msg, 'gamma_p', gamma_p); {field depedence of mobility, eV/(V/m)^0.5}

{**Contacts**************************************************************************}
		Get_Float(inv, msg, 'W_L', W_L); {eV, work function left electrode (= cathode)}
		Get_Float(inv, msg, 'W_R', W_R); {eV, work function right electrode (= anode)}
		Get_Float(inv, msg, 'Sn_L', Sn_L); {m/s, surface recombination of electrons at the left electrode.}
		Get_Float(inv, msg, 'Sp_L', Sp_L); {m/s, surface recombination of holes at the left electrode.}
		Get_Float(inv, msg, 'Sn_R', Sn_R); {m/s, surface recombination of electrons at the right electrode.}
		Get_Float(inv, msg, 'Sp_R', Sp_R); {m/s, surface recombination of holes at the right electrode.}
		Get_Float(inv, msg, 'Rshunt', Rshunt); {Ohms m2, shunt resistance. Use negative value for infinite Rshunt}
		Get_Float(inv, msg, 'Rseries', Rseries); {Ohms m2, series resistance}

{**Transport Layers************************************************************************}
		Get_Float(inv, msg, 'L_LTL', L_LTL); {m, thickness left TL}
		Get_Float(inv, msg, 'L_RTL', L_RTL); {m, thickness right TL}
		Get_String(inv, msg, 'nk_LTL', nk_LTL); {name of file with n,k values of the left TL}
		Get_String(inv, msg, 'nk_RTL', nk_RTL); {name of file with n,k values of the right TL}
		Get_Float(inv, msg, 'Nc_LTL', Nc_LTL); {m^-3, DOS of left TL}
		Get_Float(inv, msg, 'Nc_RTL', Nc_RTL); {m^-3, DOS of right TL}
		Get_Float(inv, msg, 'doping_LTL', doping_LTL);  {m^-3, doping in left TL if >0 p-type doping if <0 n-type doping}
		Get_Float(inv, msg, 'doping_RTL', doping_RTL);  {m^-3, doping in right TL if >0 p-type doping if <0 n-type doping}
		Get_Float(inv, msg, 'mob_LTL', mob_LTL); {m2/Vs, mobility of left TL}
		Get_Float(inv, msg, 'mob_RTL', mob_RTL); {m2/Vs, mobility of right TL}
		Get_Float(inv, msg, 'nu_int_LTL', nu_int_LTL); {m/s, interface transfer velocity, left TL}
		Get_Float(inv, msg, 'nu_int_RTL', nu_int_RTL); {m/s, interface transfer velocity, right TL}	
		Get_Float(inv, msg, 'eps_r_LTL', eps_r_LTL); {relative dielectric constant left TL}
		Get_Float(inv, msg, 'eps_r_RTL', eps_r_RTL); {relative dielectric constant right TL}
		Get_Float(inv, msg, 'CB_LTL', CB_LTL); {eV, conduction band left TL}
		Get_Float(inv, msg, 'CB_RTL', CB_RTL); {eV, conduction band right TL}
		Get_Float(inv, msg, 'VB_LTL', VB_LTL); {eV, valence left TL}
		Get_Float(inv, msg, 'VB_RTL', VB_RTL); {eV, valence right TL}
		Get_Integer(inv, msg, 'TLsGen', dumint); {TLsGen, TLs generate elctrons and holes? yes(1)/no(0), overrides the profile}
		TLsGen:=ROUND(dumint)=1;
		Get_Integer(inv, msg, 'TLsTrap', dumint); {TLsTrap, traps in TLs yes(1)/no(0)}
		TLsTrap:=ROUND(dumint)=1;
		Get_Integer(inv, msg, 'IonsInTLs', dumint); {can ions, if any, move into the TLs? yes(1)/no(<>1)}
		IonsInTLs:=ROUND(dumint)=1;

{**Ions*******************************************************************}
		Get_Float(inv, msg, 'CNI', CNI); {m^-3, concentration of negative ions}
		Get_Float(inv, msg, 'CPI', CPI); {m^-3, concentration of positive ions}
		IF ZimT THEN Get_Float(inv, msg, 'mobnion', mobnion);{mobility of negative ions}
		IF ZimT THEN Get_Float(inv, msg, 'mobpion', mobpion);{mobility of negative ions}
		IF SimSS THEN Get_Integer(inv, msg, 'mob_ion_spec', mob_ion_spec);{mobile ion species: -1: negative, 0: both, 1: positive ions}
		IF ZimT THEN negIonsMove:=(CNI>0) AND (mobnion>0); {are there any moving negative ions?}
		IF ZimT THEN posIonsMove:=(CPI>0) AND (mobpion>0); {are there any moving positive ions?}
		IF SimSS THEN negIonsMove:=(CNI>0) AND (dumint <= 0); {do the negative ions move?}
		IF SimSS THEN posIonsMove:=(CPI>0) AND (dumint >= 0); {do the positive ions move?}
		IF SimSS THEN Get_Integer(inv, msg, 'ion_red_rate', ion_red_rate);{number of voltage steps after which ions redistribute, }

{**Generation and recombination******************************************************}
		IF SimSS 
		THEN BEGIN
			Get_Float(inv, msg, 'Gehp', Gehp);  {generation rate of electron-hole pairs, m^-3/s}
			Get_Float(inv, msg, 'Gfrac', Gfrac);
		END;
		Get_String(inv, msg, 'Gen_profile', Gen_profile); {name of file generation profile (or 'none')}
		CASE lowercase(Trim(Gen_profile)) OF
			'none' : Use_gen_profile := 0; {Uniform generation}
			'calc' : Use_gen_profile := 1; {Calculate generation profile using the Transfermatrix method}
		ELSE
			Use_gen_profile := 2; {Use an user-defined generation profile}
		END;
		Get_Integer(inv, msg, 'Field_dep_G', dumint);  {field-dependent G, true or false}
		Field_dep_G:=(ROUND(dumint) = 1);
    	Get_Float(inv, msg, 'P0', P0); {0<=P0<1, fraction of quenched excitons that directly yield free carriers}
		Get_Float(inv, msg, 'a', a); {thermalization length, Braun model used, m}
		Get_Integer(inv, msg, 'ThermLengDist', ThermLengDist);
		Get_Float(inv, msg, 'kf', kf); {decay rate of CT state, 1/s}
		Get_Float(inv, msg, 'kdirect', kdirect); {m3/s, direct (band-to-band, bimolecular) recombination rate}
		Get_Float(inv, msg, 'Lang_pre', Lang_pre); {Langevin prefactor}
		Get_Integer(inv, msg, 'UseLangevin', dumint);
		UseLangevin:=(ROUND(dumint) = 1); {Calculate recombination using Langevin equation (1) or direct input (<>1, kdirect is used))}

{**Trapping**************************************************************************}
		{** Bulk traps}
		Get_Float(inv, msg,'Bulk_tr', Bulk_tr); {m^-3, trap density (in bulk)}
		{** Interface traps}
		Get_Float(inv, msg, 'St_L', St_L); {m^-2, left interface trap density}
		Get_Float(inv, msg, 'St_R', St_R); {m^-2, right interface trap density}
		{** Grain boundaries}
		Get_Integer(inv, msg,'num_GBs',num_GBs); {number of grain boundaries}
		Get_Float(inv, msg, 'GB_tr', GB_tr); {m^-2, grain boundary trap density per grain boundary}
		{** traps coefficients}
		Get_Float(inv, msg, 'Cn', Cn); {m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)}
		Get_Float(inv, msg, 'Cp', Cp); {m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)}
		Get_Float(inv, msg, 'ETrapSingle', ETrapSingle); {eV, energy level of all traps}
		Get_String(inv, msg, 'BulkTrapFile', BulkTrapFile); {name of file with bulk trap energy profile (or 'none'). If specified, overrides ETrapSingle}
		BulkTrapFromFile:= lowercase(Trim(BulkTrapFile))<>'none'; {use the profile if BulkTrapFile isn't 'none'}	
		Get_String(inv, msg, 'IntTrapFile', IntTrapFile); {name of file with interface trap energy profile (or 'none'). If specified, overrides ETrapSingle}
		IntTrapFromFile:= lowercase(Trim(IntTrapFile))<>'none'; {use the profile if IntTrapFile isn't 'none'}	
		Get_Integer(inv, msg, 'Tr_type_L', Tr_type_L); {Trap type for the left interface: -1: acceptor, 0: neutral, 1: donor}	
		Get_Integer(inv, msg, 'Tr_type_R', Tr_type_R); {Trap type for the right interface: -1: acceptor, 0: neutral, 1: donor}	
		Get_Integer(inv, msg, 'Tr_type_B', Tr_type_B); {Trap type of bulk and grain boundary traps: -1: acceptor, 0: neutral, 1: donor}	

{**Numerical Parameters**************************************************************}
		Get_Integer(inv, msg, 'NP', NP); {number of grid points}
		Get_Float(inv, msg, 'tolPois', tolPois); {abs. tolerance of Poisson solver}
		Get_Float(inv, msg, 'maxDelV', maxDelV); {maximum change (in Vt) of the potential per loop}
		Get_Integer(inv, msg, 'MaxItPois', MaxItPois); {Max. number of loops Poisson solver}
		Get_Integer(inv, msg, 'MaxItSS', MaxItSS); {max. number it. steady-state loops}
		IF ZimT THEN Get_Integer(inv, msg, 'MaxItTrans', MaxItTrans); {max. number it. transient solver}
		Get_Integer(inv, msg, 'CurrDiffInt', CurrDiffInt); {Calc. current from differential (1) or integral (2) expression}
		Get_Float(inv, msg, 'tolDens', tolDens); {relative tolerance of density solver}
		Get_Float(inv, msg, 'couplePC', couplePC); {>= 0, coupling between Poisson equation and continuity equations}
		Get_Float(inv, msg, 'minAcc', minAcc); {>0, min. acceleration parameter}
		Get_Float(inv, msg, 'maxAcc', maxAcc); {<2, max. acceleration parameter}
		Get_Integer(inv, msg, 'IgnoreNegDens', dumint);
		IgnoreNegDens:= dumint=1; {whether(1) or not(<>1) to ignore negative densities}
		Get_Integer(inv, msg, 'FailureMode', FailureMode); {how treat failed (t,V,G) points: 0: stop, 1: ignore, 2: skip}
		Get_Float(inv, msg, 'grad', grad); {gradient of grid, increase grad for smaller h[1]}
		IF ZimT THEN Get_Float(inv, msg, 'TolVint', TolVint); {V, tolerance in internal voltage (Vint)}

{**Voltage range of simulation*******************************************************}
		IF SimSS {this entire block is only relevant to SimSS}
		THEN BEGIN
			Get_Integer(inv, msg, 'Vdistribution', Vdistribution); {type of V distribution, 1=linear, 2=logarithmic}
			Get_Integer(inv, msg, 'PreCond', dumint); {Pre-condition in light(1)/dark(0)}
			PreCond:=ROUND(dumint)=1;
			Get_Float(inv, msg, 'Vpre', Vpre); {V, pre-conditioned voltage}
			Get_Integer(inv, msg, 'Vscan', Vscan); {integer, direction of voltage scan: up = 1, down = -1}
			Get_Float(inv, msg, 'Vmin', Vmin); {V, minimum voltage in JV characteristic}
			Get_Float(inv, msg, 'Vmax', Vmax); {V, max. voltage in JV}
			Get_Float(inv, msg, 'Vstep', Vstep); {V, voltage step}
			Get_Float(inv, msg, 'Vacc', Vacc); {accumulation voltage for logarithmic JV, should be outside [Vmin, Vmax]}
			Get_Integer(inv, msg, 'NJV', NJV); {Number of JV points, for logarithmic JV}
			Get_Integer(inv, msg, 'until_Voc', dumint); {if 1 then SimSS stops at Voc}
			until_Voc:=(dumint=1);
		END;

{**User interface********************************************************************}
		Get_Float(inv, msg, 'timeout', timeout); {s, max run time. use negative value for unlimited run time.}
		Get_Integer(inv, msg, 'Pause_at_end', Pause_at_end);  {pause at the end of the simulation yes(1) or no (0)}
		Get_Integer(inv, msg, 'AutoTidy', dumint);
		AutoTidy:=dumint = 1;	{if 1, then we will always tidy up the device_parameter file}
		IF SimSS 
		THEN BEGIN
			Get_Integer(inv, msg, 'UseExpData', dumint);
			UseExpData:=dumint = 1; {if 1 then  SimSS will try to read ExpJV and use it}
			Get_String(inv, msg, 'ExpJV', ExpJV); {name of file with experimental JV points}
			Get_String(inv, msg, 'rms_mode', dumstr); {lin or log: use J or log(J) in calc. of rms error}
			dumstr:=lowercase(dumstr);
			IF NOT ((dumstr='lin') OR (dumstr='log')) THEN Stop_Prog('rms_mode has to be either lin or log.', EC_InvalidInput);
			IF dumstr='lin' THEN rms_mode:=linear ELSE rms_mode:=logarithmic;
			Get_Float(inv, msg, 'rms_threshold', rms_threshold); {threshold of fraction converged points in calc. rms error}
		END;		
		IF ZimT THEN
		BEGIN
			Get_Integer(inv, msg, 'AutoStop', dumint);
			AutoStop:= dumint=1; {stop ZimT if change of system stops chaning, yes(1) or no (<>1).	}
			Get_String(inv, msg, 'tVG_file', tVG_file); {name of file that specifies time t, voltage V and gen. rate G}
			Get_String(inv, msg, 'tj_file', tj_file); {name of file with (t, J, V, G)}
		END;
		IF SimSS THEN Get_String(inv, msg, 'JV_file', JV_file); {name of file with simulated JV points}
		Get_String(inv, msg, 'Var_file', Var_file); {name of file with internal variables}
		Get_Integer(inv, msg, 'LimitDigits', dumint); {if 1, then number of digits in output is limited}
		LimitDigits:=dumint = 1;
		Get_Integer(inv, msg, 'OutputRatio', OutputRatio); {output (ZimT: J to screen and) variables to var_file every OutputRatio timesteps/voltages}	
		IF SimSS THEN StoreVarFile:=OutputRatio>0;
		IF ZimT THEN StoreVarFile:=lowercase(Trim(Var_file))<>'none'; {only store var_file if Var_file isn't 'none'}    
		IF SimSS THEN Get_String(inv, msg, 'scPars_file', scPars_file); {name of file with solar cell parameters}
		Get_String(inv, msg, 'log_file', log_file); { name of log file}
    END; {WITH par statement}

    CLOSE(inv);
    WRITELN('Read parameters from ',parameterFile);

END;

function formatFlt(val : myReal) : string;
{formats a float (val) into a pretty string.}
begin
    formatFlt:=FloatToStrF(val, ffGeneral, 3, 0)
end;

procedure write_simulation_setup_500(parnew : TInputParameters500; fn, path, ProgName : string); 
{produces general setup file}
var outp : text;
    i : integer;
begin
    writeln('Will create simulation setup file for ',ProgName);    
    assign(outp, path + fn);
    rewrite(outp);    
    writeln(outp,'** ',ProgName,' Simulation Setup:');
    writeln(outp,'** Don''t change the order of the parameters, comments can be added anywhere,');
    writeln(outp,'** but only after an ''*''. Use ''**'' if you want your comment to be left-justified.');
    writeln(outp,'** version: 5.00');
    writeln(outp);

    with parnew do begin

        writeln(outp, '**General***************************************************************************');
        writeln(outp, 'T = ',formatFlt(T), '                   * K, absolute temperature');
        writeln(outp);
     
        writeln(outp, '**Layers****************************************************************************');
        writeln(outp,'l1 = ',lyr[1].layerFile,' * parameter file for layer 1, mandatory');        
        for i:=2 to NLayers do
            writeln(outp,'l',i,' = ',lyr[i].layerFile,'    * parameter file for layer ',i);
        writeln(outp);

        writeln(outp, '**Contacts**************************************************************************');
        writeln(outp, 'W_L = ',formatFlt(W_L), '        * eV, work function left electrode');
        writeln(outp, 'W_R = ',formatFlt(W_R), '         * eV, work function right electrode');
        writeln(outp, 'Sn_L = ',formatFlt(Sn_L), '        * m/s, surface recombination of electrons at the left electrode');
        writeln(outp, 'Sp_L = ',formatFlt(Sp_L), '        * m/s, surface recombination of holes at the left electrode');
        writeln(outp, 'Sn_R = ',formatFlt(Sn_R), '        * m/s, surface recombination of electrons at the right electrode');
        writeln(outp, 'Sp_R = ',formatFlt(Sp_R), '        * m/s, surface recombination of holes at the right electrode');
        writeln(outp, '                                        * nb: use negative values if Sn/pR/L should be infinite');
        writeln(outp, 'Rshunt =', formatFlt(Rshunt), '                        * Ohms m2, shunt resistance. Use negative value for infinite Rshunt');
        writeln(outp, 'Rseries =', formatFlt(Rseries), '                        * Ohms m2, series resistance.');
        writeln(outp);

        writeln(outp, '**Optics****************************************************************************');

		if ProgName='SimSS' then begin
            writeln(outp, 'Gehp = ',formatFlt(Gehp), '        * m^-3 s^-1, average generation rate of electron-hole pairs');			
            writeln(outp, 'Gfrac = ',formatFlt(Gfrac), '        * fraction of Gmax used in solar cell');
		end;
        writeln(outp, 'Gen_profile = ',Gen_profile,'             * name of file generation profile (or ''none'' or ''calc'') ');   
        writeln(outp, 'L_TCO = ',formatFlt(L_TCO), '       * m, thickness of the TCO. Set to 0 if layer is not used');
        writeln(outp, 'L_BE = ',formatFlt(L_BE), '   * m, thickness of back electrode, must be >0');
        writeln(outp, 'nk_substrate = ',nk_substrate,'       * name of file with n,k values of substrate');    
        writeln(outp, 'nk_TCO = ',nk_TCO,'       * name of file with n,k values of TCO'); 
        writeln(outp, 'nk_BE = ',nk_BE,'       * name of file with n,k values of back electrode');   
        writeln(outp, 'spectrum = ', spectrum, '         * name of file that contains the spectrum');
        writeln(outp, 'lambda_min = ',formatFlt(lambda_min),'                  * m, lower bound wavelength');   
        writeln(outp, 'lambda_max = ',formatFlt(lambda_max),'                  * m, upper bound wavelength');        
        writeln(outp);


        writeln(outp, '**Numerical Parameters**************************************************************');
        writeln(outp, 'NP = ',NP,'                           * integer, number of grid points, must be at least 5 per layer.');    
        writeln(outp, 'tolPois = ',formatFlt(tolPois), '           * V, abs. tolerance of iterative Poisson solver');
        writeln(outp, 'maxDelV = ',formatFlt(maxDelV),'                   * maximum change (in Vt) of the potential per loop');
        writeln(outp, 'MaxItPois = ', MaxItPois, '                  * max. number it. Poisson loop');
        writeln(outp, 'MaxItSS = ', MaxItSS, '                       * max. number it. main loop');
        if ProgName='ZimT' then writeln(outp, 'MaxItTrans = ', MaxItTrans, '* max. number it. transient solver');
        writeln(outp, 'CurrDiffInt = ', CurrDiffInt, '                  * Calc. current from differential (1) or integral (2) expression');
        writeln(outp, 'tolDens = ',formatFlt(tolDens),'                    * relative tolerance of density solver');
        writeln(outp, 'couplePC = ',formatFlt(couplePC),'              * >= 0, coupling between Poisson equation and continuity equations');
        writeln(outp, 'minAcc = ', formatFlt(minAcc),'                   * >0, min. acceleration parameter');
        writeln(outp, 'maxAcc = ',formatFlt(maxAcc),'         * <2, max. acceleration parameter');
        writeln(outp, 'IgnoreNegDens = ',ord(IgnoreNegDens),'            * whether(1) or not(<>1) to ignore negative densities');
        writeln(outp, 'FailureMode = ',FailureMode,'                   * how treat failed (t,V,G) points: 0: stop, 1: ignore, 2: skip');
        writeln(outp, 'grad = ',formatFlt(grad),'          * determines shape of exp. grid, increase grad for smaller h[1]');
        if ProgName='ZimT' then writeln(outp, 'TolVint = ',formatFlt(TolVint), '         * V, tolerance internal voltage (Vint)');
        writeln(outp);        

        if ProgName='SimSS' then begin
            writeln(outp, '**Voltage range of simulation*******************************************************');
            writeln(outp,'Vdistribution = ',Vdistribution, '                   * 1 for uniform (specified by Vstep), 2 for logarithmic (specified by Vacc and NJV)');
            writeln(outp,'PreCond = ',ord(PreCond), '              * pre-conditioning, yes(1)/no(0)');
            writeln(outp,'Vpre = ',formatFlt(Vpre),'             * V, pre-conditioned voltage');
            writeln(outp,'FixIons = ',ord(FixIons),' * fix ions at first applied voltage? yes(1) or no (0).');  
            writeln(outp,'Vscan = ',Vscan, '                   * integer, 1 for forward sweep direction, -1 for reverse sweep');
            writeln(outp,'Vmin = ', formatFlt(Vmin), '         * V');
            writeln(outp,'Vmax = ', formatFlt(Vmax),'          * V');
            writeln(outp,'Vstep = ', formatFlt(Vstep),'         * V');
            writeln(outp,'Vacc = ', formatFlt(Vacc), '  * V, point of accumulation of row of V''s, note: Vacc should be');
            writeln(outp,'                           * slightly larger than Vmax or slightly lower than Vmin');
            writeln(outp,'NJV = ',NJV, '         * number of JV points in logarithmic distribution');
            writeln(outp,'until_Voc = ',ord(until_Voc), '              * if 1 then SimSS will stop at Voc');
            writeln(outp)
      end;        


        writeln(outp, '**User interface********************************************************************');
        writeln(outp,'timeout = ',formatFlt(timeout),'                   * s, max run time, use negative value for unlimited run time.');
        writeln(outp,'Pause_at_end = ',ord(Pause_at_end),'        * pause at the end of the simulation yes(1) or no (0)');
        writeln(outp,'AutoTidy = ',ord(AutoTidy),'        * if 1, then the program will always tidy up this file');
      
        if ProgName='SimSS' then begin
            writeln(outp,'UseExpData = ', ord(UseExpData), '             * if 1, SimSS will try to read JV_Exp and use it');
            writeln(outp,'ExpJV = ',ExpJV,'       * name of file with experimental JV characteristics');
            if rms_mode = linear then
                writeln(outp,'rms_mode = lin          * lin or log: use J or log(J) in calc. of rms error')
            else 
                writeln(outp,'rms_mode = log          * lin or log: use J or log(J) in calc. of rms error');
            writeln(outp,'rms_threshold = ',formatFlt(rms_threshold),'       * threshold of fraction converged points in calc. rms error');
            writeln(outp,'JV_file = ',JV_file,'             * name of the file with simulated JV characteristics')    
        end 
        else begin
            writeln(outp,'AutoStop = ', ord(AutoStop), '                    * stop ZimT if change of system stops changing, yes(1) or no (<>1).');
            writeln(outp,'tVG_file = ', tVG_file, '                    * name of file that specifies time t, voltage V and gen. rate G');
            writeln(outp,'tj_file = ', tj_file, '                       * name of file with (t, V, G, J, range)')
        end;
        
        writeln(outp,'Var_file = ',Var_file, '                       * name of the file with (x,V,n,p,Jn,etc) or none for no file.');
        writeln(outp,'LimitDigits = ', ord(LimitDigits),'               * if 1, then number of digits in output is limited');
        if ProgName='ZimT' then           
            writeln(outp,'OutputRatio = ',OutputRatio, '           * Output J to screen and variables to var_file every OutputRatio timesteps')
        else begin
            writeln(outp,'OutputRatio = ',OutputRatio, '            * Output to Var_file every OutputRatio voltages');
            writeln(outp,'scPars_file = ', scPars_file,'           * name of file with solar cell parameters')
        end;

        writeln(outp,'log_file = ', log_file,'                * name of log file')   

    end;
    close(outp)
end;

procedure write_layer_500(lyr : TLayerParameters500; path : string); 
{produces layer file for main absober}
var outp : text;
begin
    with lyr do begin
        assign(outp, path + layerFile);
        rewrite(outp);    
        writeln(outp,'** SIMsalabim Layer parameters:');
        writeln(outp,'** Don''t change the order of the parameters, comments can be added anywhere,');
        writeln(outp,'** but only after an ''*''. Use ''**'' if you want your comment to be left-justified.');
        writeln(outp,'** version: 5.00');
        writeln(outp);
   
        writeln(outp,'**General**************************************************************************');
        writeln(outp,'L = ',formatFlt(L),'* m, device length/thickness');
        writeln(outp,'eps_r = ',formatFlt(eps_r),'* relative dielectric constant');
        writeln(outp,'CB = ',formatFlt(CB),' * eV, conduction band edge ');                         
        writeln(outp,'VB = ',formatFlt(VB),' * eV, valence band edge '); 
        writeln(outp,'Nc = ',formatFlt(Nc),' * m^-3, DOS of conduction and valence bands '); 
        writeln(outp,'n_0 = ',formatFlt(n_0),'  * m^-3, ionised n-doping'); 
        writeln(outp,'p_0 = ',formatFlt(p_0),' * m^-3, ionised p-doping '); 
        writeln(outp);
        
        writeln(outp,'**Mobilities************************************************************************');
        writeln(outp,'mun_0 = ',formatFlt(mun_0),' * m^2/Vs, zero field mobility '); 
        writeln(outp,'mup_0 = ',formatFlt(mup_0),' * m^2/Vs, zero field mobility '); 
        writeln(outp,'mob_n_dep = ',mob_n_dep, '* 0 : const. mob, 1 : field-dependent');
        writeln(outp,'mob_p_dep = ',mob_p_dep, '* 0 : const. mob, 1 : field-dependent');
        writeln(outp,'gamma_n = ',formatFlt(gamma_n),'   * (m/V)^0.5, field dependence of mob, Poole-Frenkel form'); 
        writeln(outp,'gamma_p = ',formatFlt(gamma_p),'   * (m/V)^0.5, field dependence of mob, Poole-Frenkel form'); 
        writeln(outp);
        
        writeln(outp,'**Interface-layer-to-right**********************************************************');
        writeln(outp,'nu_int = ',formatFlt(nu_int),'  * m/s, interface transfer velocity, to layer to the right'); 
        writeln(outp,'St = ',formatFlt(St),'  * m^-2, trap density at interface with layer to the right'); 
        writeln(outp,'ETrapInt = ',formatFlt(ETrapInt),' * eV, energy level of traps at interface '); 
        writeln(outp,'IntTrapFile = ', IntTrapFile, '* name of file with interface trap energy profile (or ''none''). If specified, overrides EtrapSingle');
        writeln(outp,'Tr_type_Int = ', Tr_type_Int,' * Trap type for the right interface: -1: acceptor, 0: neutral, 1: donor');
        writeln(outp,'CnInt =',formatFlt(CnInt),'  * m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)');    
        writeln(outp,'CpInt = ',formatFlt(CpInt),'  * m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)');           
        writeln(outp);
                    
        writeln(outp,'**Ions******************************************************************************');
        writeln(outp,'CNI = ',formatFlt(CNI),'* m^-3, concentration of negative ions');    
        writeln(outp,'CPI = ',formatFlt(CPI),' * m^-3, concentration of positive ions '); 
        writeln(outp,'mobnion = ',formatFlt(mobnion),'  * m^2/Vs, mobility of negative ions (take 0 if they don''t move)'); 
        writeln(outp,'mobpion = ',formatFlt(mobpion),' * m^2/Vs, mobility of positive ions (take 0 if they don''t move) ');
        writeln(outp,'IonsMayEnter = ',ord(IonsMayEnter),'* may ions enter from other layers? yes(1) or no(<>1)');      
        writeln(outp);                   
                      
        writeln(outp,'**Generation and recombination******************************************************');
        writeln(outp,'layerGen = ',ord(layerGen), '* does this layer generate electron/hole pairs? yes(1) or no (0)');
        writeln(outp,'nk_layer = ',nk_layer,' * name of file with n,k values of this layer');
        writeln(outp,'Field_dep_G = ',ord(Field_dep_G), '  * field dependent generation yes (1) or no (0)'); 
        writeln(outp,'P0 = ',formatFlt(P0),' * 0<=P0<1, fraction of quenched excitons that direcltly yield free carriers '); 
        writeln(outp,'a = ',formatFlt(a),'  * m, charge separation distance, Braun model used '); 
        writeln(outp,'ThermLengDist = ', ThermLengDist,'* distribution of a, 1 for delta function, 2 for Gaussian');    
        writeln(outp,'                                       * 3 for exponential and 4 for r^2 exponential 5 for r^4 Gaussian');    
        writeln(outp,'kf = ',formatFlt(kf),' * 1/s, decay rate '); 
        writeln(outp,'kdirect = ',formatFlt(kdirect),' * m3/s, direct (band-to-band, bimolecular) recombination rate '); 
        writeln(outp,'Lang_pre = ',formatFlt(Lang_pre),'  * Langevin recombination prefactor ');
        writeln(outp,'UseLangevin = ',ord(UseLangevin), '    * (1) use Langevin to calc. recombination or not (<>1, kdirect is used)'); 
        writeln(outp);
                             
        writeln(outp,'**Bulk trapping**************************************************************************');
        writeln(outp,'Bulk_tr = ',formatFlt(Bulk_tr),'  * m^-3, trap density (in bulk) '); 
        writeln(outp,'CnBulk = ',formatFlt(CnBulk),' * m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)'); 
        writeln(outp,'CpBulk = ',formatFlt(CpBulk),'  * m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band) '); 
        writeln(outp,'ETrapBulk = ',formatFlt(ETrapBulk),'   * eV, energy level of all traps'); 
        writeln(outp,'BulkTrapFile = ',BulkTrapFile,'   * name of file with bulk trap energy profile (or ''none''). If specified, overrides EtrapBulk'); 
        writeln(outp,'Tr_type_B = ',Tr_type_B, ' * Trap type of bulk traps: -1: acceptor, 0: neutral, 1: donor')

    end; {with lyr statement}
    close(outp)

end;


PROCEDURE Tidy_Up_File(FileName, path : STRING);
{This procedure does the actual work of tidying the parameter files.}
VAR inp, outp : TEXT;
	temp : ANSISTRING;
    line, dumstr, outline : STRING;
    max_pos, linecount, pos_asterix, i : INTEGER;
BEGIN
    {open the original parameterFile}
    ASSIGN(inp, FileName);
    RESET(inp);
	temp:=''; {this is an ansistring (so no length limit) where we store the input file and manipulate it later on}

    max_pos:=0;
    linecount:=0;
    
    WHILE NOT EOF(inp) DO
    BEGIN
		READLN(inp, line);
		INC(linecount);
		line:=TRIM(line); {strips blank characters at the beginning and end of line}
		{now we have read a line, we are going to process it to form outline.}
		{outline is what we write in the temp output file.}
	
		{if the line is empty, or it is a comment (*) then simply copy:}
		IF (LeftStr(line, 1)='*') OR (LeftStr(line, 1)='') 
			THEN outline:=line;
		
		IF (LeftStr(line, 1)<>'*') AND (LENGTH(line)<>0) THEN {we should have a line that contains a parameter}
		BEGIN
			{first check if there is an =}
			dumstr:='';
			dumstr:=ExtractWord(1, line, myDelims); {Extracts the first word from a string, should contain the parameter name}
			dumstr:=dumstr + ' = '+ ExtractWord(2, line, myDelims); {this puts the = in the right spot, including white space}
			pos_asterix:=POS('*', line);
			{a line with a parameter should contain some units, description, etc.}
			{so if there is no asterix, or there is no comment after the asterix, we have a problem:}
			IF (pos_asterix=0) OR (TRIM(RightStr(line, LENGTH(line)-pos_asterix))='') THEN {if this is not the case, we will point this out:}
			BEGIN
				WRITELN('A line that contains a parameter should also contain some unit,');
				WRITELN('description, etc, and always after a *.');
				WRITELN('In file ',FileName,', the line that starts with "',LeftStr(line, MIN(LENGTH(line),20)),'"');
				WRITELN('does not obey this rule.');				
				WRITELN('This is something you need to fix yourself.');
				Stop_Prog('See line number '+IntToStr(linecount)+'.', EC_DevParCorrupt);
			END;
			{now we are sure there is an *}
			dumstr:=dumstr + ' * ' + TrimLeft(RightStr(line, LENGTH(line)-pos_asterix));
			pos_asterix:=POS('*', dumstr);
			outline:=dumstr;
			max_pos:=MAX(max_pos, pos_asterix);
		END;
		
		temp:=temp + outline + LineEnding; {add line to our temp version of the parameter file}
		
    END; {while loop reading input file}
    
    CLOSE(inp);
    {now we know where to put the * (max_pos)}
    
    {start writing into the real parameterFile:}
    ASSIGN(outp, FileName);
    REWRITE(outp);
    
    WHILE LENGTH(temp) <> 0 DO
    BEGIN
		line:=Copy2StringDel(temp, LineEnding); {Deletes and returns all characters in a string till a given string (not included).}
		IF (LeftStr(line, 2)<>'**') AND (LENGTH(line)<>0) THEN {we should have a line that contains a parameter}
		BEGIN
			pos_asterix:=POS('*', line); {now we know there must be an asterix}
			dumstr:=TrimRight(LeftStr(line, pos_asterix-1)); {this is the part that contains the name of the parameter}
			WRITE(outp, dumstr);
			FOR i:=1 TO max_pos - 1 - LENGTH(dumstr) DO WRITE(outp, ' '); {write spaces}
			dumstr:=TrimLeft(RightStr(line, LENGTH(line)-pos_asterix)); 
			WRITELN(outp, '* ',dumstr); {now add * unit, etc.}
		END
		ELSE WRITELN(outp, line);
	
    END;
    
    CLOSE(outp);

	WRITELN('Tidied-up ',FileName,'.');
END;

PROCEDURE Tidy_Up_Parameter_Files(parameterFile, path : STRING; QuitWhenDone : BOOLEAN; CONSTREF par : TInputParameters500);
{This procedure cleans up the parameter files: the main file and the files with the parameters of each layer}
{It makes sure that all the * are aligned, that every parameter has a 
unit or other description/comment and left-aligns any line that starts with **. It does this by reading the original device
parameter file line by line and writing corrected lines to a temp file. Once the correct position of the descriptions 
(starting with *) has been found, it uses the temp file to create the new, tidy parameter file. Lastly, the temp file is
removed and the program exits.}
VAR i : INTEGER;	
BEGIN
	{first, do the main paramter file:}
	Tidy_Up_File(parameterFile, path);
	
	{next, loop over files that contain the parameters of the layers:}
	FOR i:=1 TO par.NLayers DO
		Tidy_Up_File(par.lyr[i].layerFile, path);
		
	IF QuitWhenDone THEN Stop_Prog('Done cleaning parameter files.', EC_Warning)
END;




BEGIN 

END.
