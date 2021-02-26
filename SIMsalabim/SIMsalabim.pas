PROGRAM SIMsalabim;

{
SIMsalabim:a 1D drift-diffusion simulator 
Copyright (c) 2020 Dr T.S. Sherkar, V.M. Le Corre, M. Koopmans,
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


{If we're using windows we need to make it a console application:}
{$IFDEF WINDOWS}
{$APPTYPE CONSOLE}
{$ENDIF}

{$MODE DELPHI} {force Delphi mode}
{Note: we use Delphi mode for the passing of functions in RombergIntegration.
The default FPC mode and also ObjFPC mode requires you to prefix
function given as argument with @. If you don't like this behavior,
you can use either ($mode tp) or ($mode delphi) directive to change
to Turbo Pascal or Delphi mode, respectively.}

USES TypesAndConstants IN 'Units/',
InputOutputUtils IN 'Units/', {for DelWhite and DelWhite1}
NumericalUtils IN 'Units/',
Math,  {for power, min and max functions}
StrUtils, {for DelSpace}
SysUtils; {for TRIM function in Init_Elec_Mob_Table}


CONST
    version = '3.86';   {version, 1.00 = 10-03-2004}
    parameter_file = 'device_parameters.txt'; {name of file with parameters}

{Magic numbers, used in the code:}
    {Name                    Where                        What does it do?}
    threshold_err = 0.2; {IN: Find_Solar_Cell_Parameters, defines max relative error when displaying solar cell parameters}
    InterpolationOrder = 4; {IN: Find_Solar_Cell_Parameters, interpolation order in estimating Jsc, Voc}
    temp_file = '.temp.txt';{IN: Tidy_up_parameter_file_Exit, a temporary file used to store the device parameters}
    myDelims = [#0..' ', ';', '/', '\', '''', '"', '`', '*', '=']; {IN Tidy_up_parameter_file_Exit, note: not a decimal point or comma!}
    nd = 25; {IN: main, it limits the number of digits to this value to prevent unnecessarily large output files}
    minErr = 1E-4; {IN: Calc_and_Output_Soalr_Cell_Parameters, lower bound to error in parameters}
    tab1 = 6; {IN: Calc_and_Output_Soalr_Cell_Parameters, 1st position of tab in table with solar cell parameters}
    tab2 = 28; {IN: Calc_and_Output_Soalr_Cell_Parameters, 2nd position of tab in table with solar cell parameters}
    tab3 = 50; {IN: Calc_and_Output_Soalr_Cell_Parameters, 3rd position of tab in table with solar cell parameters}


VAR i, VCount, MainIt, i1, i2 : INTEGER;

	V, Vgn, Vgp, E_VB, E_CB, n, p, nid, pid, empty, oldn, oldp, Jn, Jp, rec, h, x,
    mun, mup, Dn, Dp, Gm, gen, Lang, SRHB, diss_prob, eps, Ntb,  Nti, Nti_charge, cwe_i, charged_traps, nion, pion, 
    NcLoc, nt0, pt0, ni, n_contn, p_contp : vector;
    
{for tabulated electron mobility:}
    n_vals, Fn_vals, p_vals, Fp_vals : row; {values of n and F for tabulated electron mobility}
    n_mob_tab, p_mob_tab : Table; {the table of mobilities}
    n_points, p_points, Fn_points, Fp_points, flip_count_deln, flip_count_delp : INTEGER; {number of points in F and n direction}

	check_Poisson, Conv_Main, extra_F, extra_ln_mob, Field_dep_G, Use_gen_profile,
    UseLangevin, ImageLowering, Traps, Traps_int, until_Voc, quit_Voc,
	TLsAbsorb, TLsTrap, IonsInTLs, PreCond, resetNegDens, AutoTidy, UseExpData : BOOLEAN;

	Egap, Jtot, oldJtot, totn, totp, Vt, Vti, F, Braun_rec, V0, VL, Av_Diss,
	recLan, recSRH, phi_left, phi_right, phi_left_eff, phi_right_eff, epsi : myReal;

	T, L, eps_r, CB, VB, Nc, n_0, p_0, mun_0, mup_0, beta_n, beta_p,
	gamma_n, gamma_p, W_L, W_R, L_LTL, L_RTL, Nc_LTL, Nc_RTL, 
	doping_LTL , doping_RTL, mob_LTL , mob_RTL, nu_int_LTL, nu_int_RTL,
	eps_r_LTL, eps_r_RTL, CB_LTL, VB_LTL, CB_RTL, VB_RTL,
	RoughLeft, RoughRight, Sn_L, Sp_L, Sn_R, Sp_R, Rshunt, Rseries,
	eta_n, eta_p, Gmax, Gfrac, P0, a, kf, kdirect, Lang_pre, Bulk_tr, St_L,
	St_R, Etrap, cwe_b, GB_tr, Cn, Cp, nu_ni, nu_pi,
	tolPois, maxDelV, accPois, accDens,tolMain, grad, TolRomb, Vpre, 
	Vmin, Vmax, Vstep, Vacc, LowerLimBraun, UpperLimBraun,
	CIM, rms_threshold, del_totn, del_totp, old_totn, old_totp, old_del_totn,
	old_del_totp, red_fac, max_del, max_found_deln, max_found_delp, test_p : myReal; {paramaters from input file}

	Jbimo, JSRH_bulk, JSRH_LI, JSRH_RI, Jph, Jn_l, Jn_r, Jp_l, Jp_r : myReal; {recombination current}

	MaxItPois, MaxItMain, MaxRombIt, mob_n_dep, mob_p_dep, Dn_dep, Dp_dep,
	ThermLengDist, Tr_type_L, Tr_type_R, Tr_type_B, Conv_Var, Vdistribution,
	NJV, mob_ion_spec, ion_red_rate, num_GBs, NP, Vscan,
	Pause_at_end, Exit_after_fail : INTEGER; { also form input file}

	rms_mode : Tfitmode; {parameter from input file}

	inv, uitv , log : TEXT; {the input and output files}

	delDens, Va, Vaold : myReal;

	Gen_profile, ExpJV, log_file, JV_file, Var_file, n_file, p_file, scPars_file : STRING;

	JVSim, JVExp : JVList; {stores the current-voltage characteristics}
    
	MsgStr : ANSISTRING = ''; {Ansistrings have no length limit, init string to ''}

PROCEDURE DisplayHelpExit;
{displays a short help message and exits}
BEGIN
    WRITELN('SIMsalabim can be used with the following options:');
    WRITELN;
    WRITELN('All parameters can be set via ',parameter_file);
    WRITELN('or via the command line, which overrides the values in the file.');
    WRITELN('Example: ./SIMsalabim -T 400 -Var_file Var.dat');
    WRITELN;
    WRITELN('To tidy-up the parameter file, use ''-tidy''');
    WRITELN;
    Stop_Prog('''-h''     : displays this help message');
END;

FUNCTION Norm(vec : vector; istart, ifinish : integer) : myReal;
VAR i : INTEGER;
    ans : myReal;
    {computes the infinity norm of vector vec}
BEGIN
    ans:=0;
    FOR i:=istart TO ifinish DO
	IF ABS(vec[i]) > ans THEN ans:=ABS(vec[i]);
    Norm:=ans
END;


FUNCTION Correct_version_parameter_file : BOOLEAN;
{check for version info:
we do this by checking if there is a line that contains both the
string 'version' and the string that contains the version number of the program (version).
Note, this is not very strict. Example: if the version in the
parameter_file is 'version: 13.57' and
SIMsalabim is version '3.57', then this function will return true even
though it is not correct. However,
I don't think that is a problem as it is much more like to run
SIMsalabim with a parameter_file that is just slightly
older or newer ('3.57' versus '3.52' for example}
VAR inp : TEXT;
    found_it : BOOLEAN;
    line : STRING;
BEGIN
    IF NOT FileExists(parameter_file) {the file with input par. is not found}
	THEN Stop_Prog('Could not find file '+parameter_file+'.');
    ASSIGN(inp, parameter_file);
    RESET(inp);
    found_it:=FALSE;
    WHILE NOT EOF(inp) AND NOT found_it DO
    BEGIN
		READLN(inp, line); {read a line from the file}
		line:=LOWERCASE(line); {convert all characters to lower case to make search easier}
		IF (POS('version', line) > 0) AND (POS(version, line) > 0) THEN
			found_it:=TRUE; {we have found the version number and it is correct!}
		{POS returns the index of Substr in S, if S contains Substr. In case Substr isn't found, 0 is returned.}
    END;
    CLOSE(inp);
    Correct_version_parameter_file:=found_it;
END;

PROCEDURE Init_Generation_Profile(Use_gen_profile : BOOLEAN; gen_file : STRING; VAR genrate : vector);
VAR inv : TEXT; {file with generation profile}
    a, gr : ARRAY[0..5000] OF myReal; {a : x-coordinate, gr: generation rate}
    counter, i, j : INTEGER;
    Gav : myReal;
BEGIN
    IF Gmax=0
	THEN FOR i:=0 TO NP+1 DO Gm[i]:=0 {if Gmax=0 then profile is irrelevent}
	ELSE BEGIN {Gmax <>0}
	    IF (NOT Use_gen_profile)
		THEN FOR i:=0 TO NP+1 DO Gm[i]:=Gmax {constant Gmax=generation rate of e-h pairs}

		ELSE 
		BEGIN {optical interference effects are taken into account}
		    IF NOT FileExists(gen_file) {file with gen. profile is not found}
			THEN Stop_Prog('Could not find file '+gen_file);
		    ASSIGN(inv, gen_file); {once we get here, we're sure the file exists}
		    RESET(inv);

		    {the input file contains the generation profile, however, the grid points
		    in this file may not correspond to our grid here}
		    counter:=-1; {count the number of points in the input file}
		    WHILE NOT(EOF(inv)) DO {read all the generation rates from file}
		    BEGIN counter:=counter+1; READLN(inv, a[counter], gr[counter]) END;
		    {a: x-coordinate, gr: corresponding generation rate}
		    CLOSE(inv);
		    IF counter=-1 THEN Stop_Prog('The file generation_profile.txt is empty.');

		    { rescale a to ensure that a[counter]=L: }
		    FOR i:=0 TO counter DO a[i]:=a[i] * L/a[counter];

		    {Now interpolate gr to get genrate on x[i] grid:}
		    i:=0; {counter for x[i] grid}
		    FOR j:=0 TO counter-1 DO
				WHILE (x[i] < a[j+1]) AND (i<NP+1) DO
				BEGIN
					genrate[i]:=gr[j] + (x[i]-a[j]) * (gr[j+1]-gr[j])/(a[j+1]-a[j]);
					{we're using linear interpolation here}
					i:=i+1
				END;
		END; {reading generation profile from file}

	    {if the transport layers (if there are any) don't absorb we need set the generation
	    rate to zero. This overrides the generation profile (if any).}
	    IF NOT TLsAbsorb THEN
	    BEGIN
			{left transport layer, offset generalised potentials:}
			IF i1>0 THEN { I would have thought that one shouldn't need this if statement}
				FOR i:=0 TO i1 DO
					genrate[i]:=0;
			{right transport layer, offset generalised potentials:}
			IF i2<NP+1 THEN
				FOR i:=i2 TO NP+1 DO
					genrate[i]:=0;
	    END;

	    {Rescale such that Gmax equals the average of the profile over the
	    length of the absorber. The latter equals L if TLsAbsorb or there
	    are no transport layers, and L-L_LTL-L_RTL otherwise. }
	    Gav:=Average(genrate, h, 0, NP+1);
	    IF NOT TLsAbsorb THEN Gav:=Gav * L / (L-L_LTL-L_RTL);
	    FOR i:=0 TO NP+1 DO
			genrate[i]:=genrate[i] * (Gmax/Gav)
	END; {Gmax<>0}
END;

PROCEDURE Prepare_Log_File(VAR log : TEXT; MsgStr : ANSISTRING);
BEGIN
    ASSIGN(log, log_file); 

    REWRITE(log);
    WRITELN(log,'Version ', version);
    WRITELN(log,'Size of reals used in simulation: ',SizeOf(myReal),' bytes');
    
    {Read_Parameters may have added something to MsgStr (values from the command line):}
    IF Length(MsgStr) > 0 THEN
    BEGIN
		WRITELN(log, 'Values from command line:');
		WRITE(log, MsgStr) {MsgStr contains any messages from Read_Parameter. No need for writeln as it either empty, or end with LineEnding}
    END;
    FLUSH(log);
END;

PROCEDURE Tidy_up_parameter_file(QuitWhenDone : BOOLEAN);
{This procedure cleans up the parameter file: it makes sure that all the * are aligned, that every parameter has a 
unit or other description/comment and left-aligns any line that starts with **. It does this by reading the original device
parameter file line by line and writing corrected lines to a temp file. Once the correct position of the descriptions 
(starting with *) has been found, it uses the temp file to create the new, tidy parameter file. Lastly, the temp file is
removed and the program exits.}
VAR inp, outp : TEXT;
    line, dumstr, outline : STRING;
    max_pos, linecount, pos_asterix, i : INTEGER;
BEGIN
    {open the original parameter_file}
    ASSIGN(inp, parameter_file);
    RESET(inp);
    {and create a temp output file:}
    ASSIGN(outp, temp_file);
    REWRITE(outp);
    
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
				WRITELN('The line that starts with "',LeftStr(line, MIN(LENGTH(line),20)),'" does not obey this rule.');				
				WRITELN('This is something you need to fix yourself.');
				Stop_Prog('See line number '+IntToStr(linecount)+'.');
			END;
			{now we are sure there is an *}
			dumstr:=dumstr + ' * ' + TrimLeft(RightStr(line, LENGTH(line)-pos_asterix));
			pos_asterix:=POS('*', dumstr);
			outline:=dumstr;
			max_pos:=MAX(max_pos, pos_asterix);
		END;
		WRITELN(outp, outline);
    END; {while loop reading input file}
    
    CLOSE(inp);
    CLOSE(outp);
    {now we know where to put the * (max_pos)}
    
    {now start reading from the temp file:}
    ASSIGN(inp, temp_file);
    RESET(inp);
    {and start writing into the real parameter_file:}
    ASSIGN(outp, parameter_file);
    REWRITE(outp);
    
    WHILE NOT EOF(inp) DO
    BEGIN
		READLN(inp, line);
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
	
    END; {while loop reading input file}
    
    CLOSE(inp);
    CLOSE(outp);
    
    {now delete old parameter_file and replace with temp_file}
    ERASE(inp); {the file should be assigned with Assign, but not opened with Reset or Rewrite}
    RENAME(outp, parameter_file); {file must be assigned, but not opened}
    IF QuitWhenDone 
		THEN Stop_Prog('Replaced old '+parameter_file+' with a tidy version.')
		ELSE WRITELN('Tidied-up ',parameter_file,'.');
END;

PROCEDURE Read_Parameters(VAR msg : ANSISTRING);
VAR dummy : myReal; {dummy is used for reading integer values, these are read as
		    myReal and then converted to integers}
    dumint : INTEGER; {a dummy integer variable}
    dumstr : STRING; {again, a dummy variable}
BEGIN
    IF NOT FileExists(parameter_file) {the file with input par. is not found}
	THEN Stop_Prog('Could not find file '+parameter_file+'.');
    ASSIGN(inv, parameter_file);
    RESET(inv);
    
    {**General**************************************************************************}
    Get_Float(inv, msg,'T',T);  {abs. temperature, K}
    Get_Float(inv, msg, 'L', L); {Read_Number(inv, L);  {device thickness, m}
    Get_Float(inv, msg, 'eps_r', eps_r);  {relative dielectric constant}
    Get_Float(inv, msg, 'CB', CB);  {eV, conduction band edge}
    Get_Float(inv, msg, 'VB', VB);  {eV, valence band edge}
    Get_Float(inv, msg, 'Nc',Nc);  {effective DOS, m^-3}
    Get_Float(inv, msg, 'n_0', n_0);  {ionised n-doping density, m^-3}
    Get_Float(inv, msg, 'p_0', p_0);  {ionised p-doping density, m^-3}
    Get_Float(inv, msg, 'beta_n', beta_n); {(m/V)^0.5, field  activation factor of n_0}
    Get_Float(inv, msg, 'beta_p', beta_p); {(m/V)^0.5, field  activation factor of p_0}

    {**Mobilities************************************************************************}
    Get_Float(inv, msg, 'mun_0', mun_0); {electron zero-field mobility, m^2/Vs}
    Get_Float(inv, msg, 'mup_0', mup_0); {hole zere-field mobility, m^2/Vs}
    Get_Integer(inv, msg, 'mob_n_dep', mob_n_dep);{dependence of elec mobility, 0 : const. mob, 1 : field-dep., 2 : get mobility from table specified in n_file}
    Get_Integer(inv, msg, 'mob_p_dep', mob_p_dep);  {dependence of hole mobility, 0 : const. mob, 1 : field-dep., 2 : get mobility from table specified in p_file}
    Get_Float(inv, msg,'gamma_n', gamma_n); {field depedence of mobility, eV/(V/m)^0.5}
    Get_Float(inv, msg, 'gamma_p', gamma_p); {Poole-Frenkel form}
    Get_String(inv, msg, 'n_file', n_file); {read name of file with tabulated electron mobility}
    Get_String(inv, msg, 'p_file', p_file); {read name of file with tabulated hole mobility}
    Get_Integer(inv, msg, 'extra_F', dumint); {whether (1) or not (0) to extrapolate F beyond its max. value in n_file or p_file}
    extra_F:=ROUND(dumint) = 1;
    Get_Integer(inv, msg, 'extra_ln_mob', dumint); {whether (1) or not (0) to inter-/extrapolate ln(mob) instead of mob in n_file or p_file}
    extra_ln_mob:=ROUND(dumint) = 1;
    Get_Integer(inv, msg, 'Dn_dep', dumint); {read Dn_dep}
    Dn_dep:=ROUND(dumint); {0 : Einstein, 1 : value (eta_n) times Einstein}
    Get_Integer(inv, msg, 'Dp_dep', dumint); {read Dp_dep}
    Dp_dep:=ROUND(dumint); {0 : Einstein, 1 : value (eta_p) times Einstein}
    Get_Float(inv, msg, 'eta_n', eta_n); {multiplication factor Dn; Dn = eta_n Vt mun}
    Get_Float(inv, msg, 'eta_p', eta_p); {multiplication factor Dp; Dp = eta_p Vt mup}

    {**Contacts**************************************************************************}
    Get_Float(inv, msg, 'W_L', W_L); {eV, work function left electrode (= cathode)}
    Get_Float(inv, msg, 'W_R', W_R); {eV, work function right electrode (= anode)}
    
    phi_left:=W_L-CB;
    phi_right:=VB-W_R;

    V0:=0.5 * (VB+CB) - W_L;
    VL:=0.5 * (VB+CB) - W_R;
    
    Get_Integer(inv, msg, 'ImageLowering', dumint);
    ImageLowering:= (dumint = 1); {whether (1) or not (<>1) to use image-force lowering of barriers.
				  Note: only works for electrons at left and holes at right electrode
				  and only if surface recombination is infinite. }
    Get_Float(inv, msg, 'RoughLeft', RoughLeft); {extra barrier lowering due to roughness of cathode}
    Get_Float(inv, msg, 'RoughRight', RoughRight); {extra barrier lowering due to roughness of cathode}
    Get_Float(inv, msg, 'Sn_R', Sn_R); {m/s, surface recombination of electrons at the right electrode.}
    Get_Float(inv, msg, 'Sp_R', Sp_R); {m/s, surface recombination of holes at the right electrode.}
    Get_Float(inv, msg, 'Sn_L', Sn_L); {m/s, surface recombination of electrons at the left electrode.}
    Get_Float(inv, msg, 'Sp_L', Sp_L); {m/s, surface recombination of holes at the left electrode.}
    Get_Float(inv, msg, 'Rshunt', Rshunt); {Ohms m2, shunt resistance. Use negative value for infinite Rshunt}
    Get_Float(inv, msg, 'Rseries', Rseries); {Ohms m2, series resistance}

    {**Transport Layers******************************************************************}
    {**syntax: [name variable]_LTL/RTL for left, resp., right transport layer}
    Get_Float(inv, msg, 'L_LTL', L_LTL); {m, thickness left TL}
    Get_Float(inv, msg, 'L_RTL', L_RTL); {m, thickness right TL}
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
    Get_Integer(inv, msg, 'TLsAbsorb', dumint); {TLsAbsorb, TLs absorb yes(1)/no(0), overrides the profile}
    TLsAbsorb:=ROUND(dumint)=1;
    Get_Integer(inv, msg, 'TLsTrap', dumint); {TLsTrap, traps in TLs yes(1)/no(0)}
    TLsTrap:=ROUND(dumint)=1;
    Get_Integer(inv, log, 'IonsInTLs', dumint); {can ions, if any, move into the TLs? yes(1)/no(<>1)}
    IonsInTLs:=ROUND(dumint)=1;
   
    {**Ions*******************************************************************}
    Get_Float(inv, msg, 'CIM', CIM); {m^-3, concentration of ions}
    Get_Float(inv, msg, 'mob_ion_spec', dummy);{mobile ion species: -1: negative, 0: both, 1: positive ions}
    mob_ion_spec := ROUND(dummy);
    Get_Integer(inv, msg, 'ion_red_rate', ion_red_rate);{number of voltage steps after which ions redistribute, }
    
    {**Generation and recombination******************************************************}
    Get_Float(inv, msg, 'Gmax', Gmax);  {maximum generation rate, m^-3/s}
    Get_Float(inv, msg, 'Gfrac', Gfrac);
    Gmax:=Gmax*Gfrac;
    Get_String(inv, msg, 'Gen_profile', Gen_profile); {name of file generation profile (or 'none')}
    Use_gen_profile:= lowercase(Trim(Gen_profile))<>'none'; {use the profile if Gen_profile isn't 'none'}
    Get_Integer(inv, msg, 'Field_dep_G', dumint);  {field-dependent G, true or false}
    Field_dep_G:=(ROUND(dumint) = 1);
    Get_Float(inv, msg, 'P0', P0); {0<=P0<1, fraction of quenched excitons that direcltly yield free carriers}
    Get_Float(inv, msg, 'a', a); {thermalization length, Braun model used, m}
    Get_Integer(inv, msg, 'ThermLengDist', dumint);
    ThermLengDist:=ROUND(dumint);
    Get_Float(inv, msg, 'kf', kf); {decay rate of CT state, 1/s}
    Get_Float(inv, msg, 'kdirect', kdirect); {m3/s, direct (band-to-band, bimolecular) recombination rate}
    Get_Float(inv, msg, 'Lang_pre', Lang_pre); {Langevin prefactor}
    Get_Integer(inv, msg, 'UseLangevin', dumint);
    UseLangevin:=(ROUND(dumint) = 1); {Calculate recombination using Langevin equation (1) or direct input (<>1, kdirect is used))}

    {**Trapping**************************************************************************}
    {** Bulk traps}
    Get_Float(inv, msg,'Bulk_tr', Bulk_tr); {m^-3, trap density (in bulk)}
    {** Surface traps}
    Get_Float(inv, msg, 'St_L', St_L); {m^-2, left interface trap density}
    Get_Float(inv, msg, 'St_R', St_R); {m^-2, right interface trap density}
    {** Grain boundaries}
    Get_Integer(inv, msg,'num_GBs',num_GBs); {number of grain boundaries}
    Get_Float(inv, msg, 'GB_tr', GB_tr); {m^-2, grain boundary trap density per grain boundary}
    {** traps coefficients}
    Get_Float(inv, msg, 'Cn', Cn); {m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)}
    Get_Float(inv, msg, 'Cp', Cp); {m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)}
    nu_ni := Cn; {m/s, interface recombination velocity for electrons at interfaces and grain boundaries}
    nu_pi := Cp; {m/s, interface recombination velocity for holes at interfaces and grain boundaries}
    Get_Float(inv, msg, 'Etrap', Etrap); {eV, energy level of all traps}
    Get_Integer(inv, msg, 'Tr_type_L', dumint);
    Tr_type_L:=ROUND(dumint); {Trap type for the left interface: -1: acceptor, 0: neutral, 1: donor}	
    Get_Integer(inv, msg, 'Tr_type_R', dumint);
    Tr_type_R:=ROUND(dumint); {Trap type for the right interface: -1: acceptor, 0: neutral, 1: donor}	
    Get_Integer(inv, msg, 'Tr_type_B', dumint);
    Tr_type_B:=ROUND(dumint); {Trap type of bulk and grain boundary traps: -1: acceptor, 0: neutral, 1: donor}	
    Traps:=(Bulk_tr<>0) OR (St_L<>0) OR (St_R<>0) OR ((GB_tr<>0) AND (num_GBs>0)); {are there any traps?}
    Traps_int:=(St_L<>0) OR (St_R<>0) OR ((GB_tr<>0) AND (num_GBs>0)); {are there any interface traps?}
    
    {**Numerical Parameters**************************************************************}
    Get_Integer(inv, msg, 'NP', NP); {number of grid points}
    Get_Float(inv, msg, 'tolPois', tolPois); {abs. tolerance of Poisson solver}
    Get_Float(inv, msg, 'maxDelV', maxDelV); {maximum change (in Vt) of the potential per loop}
    Get_Float(inv, msg, 'accPois', accPois); {SOR/SUR acceleration parameter for Poisson, (0<accPois<2)}
    Get_Float(inv, msg, 'accDens', accDens); {SOR/SUR acceleration parameter for densities, (0<accDens<2)}
    Get_Integer(inv, msg, 'resetNegDens', dumint);
    resetNegDens:=(dumint = 1); {whether(1) or not(<>1) to reset points with a negative density}
    Get_Float(inv, msg, 'tolMain', tolMain); {rel. tolerance of main loop}
    Get_Float(inv, msg, 'Conv_VAR', dummy); {1 selects current to be monitored for convergence in main loop}
    Conv_Var:=ROUND(dummy); {<> selects the Slotboom variables (cf. densities)}
    Get_Float(inv, msg, 'MaxItPois', dummy); {Max. number of loops Poisson solver}
    MaxItPois:=ROUND(dummy);
    Get_Float(inv, msg, 'MaxitMain', dummy); {Max. number of main loops}
    MaxItMain:=ROUND(dummy);
    Get_Float(inv, msg, 'grad', grad); {gradient of grid, increase grad for smaller h[1]}
    Get_Float(inv, msg, 'TolRomb', TolRomb); {rel. tolerance of Romberg integration}
    Get_Float(inv, msg, 'MaxRombIt', dummy);
    MaxRombIt:=ROUND(dummy); {max. # Romberg iterations}
    Get_Float(inv, msg, 'LowerLimBraun', LowerLimBraun); {lower limit of integration over distribution Braun}
    Get_Float(inv, msg, 'UpperLimBraun', UpperLimBraun); {upper limit}

    {**Voltage range of simulation*******************************************************}
    Get_Float(inv, msg, 'Vdistribution', dummy);
    Vdistribution:=ROUND(dummy); {type of V distribution, 1=linear, 2=logarithmic}
    Get_Integer(inv, msg, 'PreCond', dumint); {Pre-condition in light(1)/dark(0)}
    PreCond:=ROUND(dumint)=1;
    Get_Float(inv, msg, 'Vpre', Vpre); {V, pre-conditioned voltage}
    Get_Integer(inv, msg, 'Vscan', Vscan); {integer, direction of voltage scan: up > 0, down < 0}
    Get_Float(inv, msg, 'Vmin', Vmin); {V, minimum voltage in JV characteristic}
    Get_Float(inv, msg, 'Vmax', Vmax); {V, max. voltage in JV}
    Get_Float(inv, msg, 'Vstep', Vstep); {V, voltage step}
    Get_Float(inv, msg, 'Vacc', Vacc); {accumulation voltage for logarithmic JV, should be outside [Vmin, Vmax]}
    Get_Integer(inv, msg, 'NJV', NJV); {Number of JV points, for logarithmic JV}
    IF (Vstep<>0) AND (Vdistribution=1) THEN
	NJV:=TRUNC((Vmax - Vmin)/Vstep + 1e-10) + 1; {needed for setting length of Jdat and Vdat}
    {1e-10 is needed to get right value}
    Get_Integer(inv, msg, 'until_Voc', dumint); {if 1 then SIMsalabim stops at Voc}
    until_Voc:=(dumint=1);

    {**User interface********************************************************************}
    Get_Integer(inv, msg, 'Pause_at_end', Pause_at_end); {pause at the end of the simulation yes(1) or no (0)}
    Get_Integer(inv, msg, 'AutoTidy', dumint);
    AutoTidy:=dumint = 1;	{if 1, then SIMsalabim will always tidy up the device_parameter file}
    Get_Integer(inv, msg, 'UseExpData', dumint);
    UseExpData:=dumint = 1; {if 1 then  SIMsalabim will try to read ExpJV and use it}
    Get_String(inv, msg, 'ExpJV', ExpJV); {name of file with experimental JV points}
    Get_String(inv, msg, 'rms_mode', dumstr); {lin or log: use J or log(J) in calc. of rms error}
    dumstr:=lowercase(dumstr);
    IF NOT ((dumstr='lin') OR (dumstr='log')) THEN Stop_Prog('rms_mode has to be either lin or log.');
    IF dumstr='lin' THEN rms_mode:=linear ELSE rms_mode:=logarithmic;
    Get_Float(inv, msg, 'rms_threshold', rms_threshold); {threshold of fraction converged points in calc. rms error}
    Get_Integer(inv, msg, 'Exit_after_fail', Exit_after_fail); {Exit simulation when a point does not converge yes(1) or no (0)}
    Get_String(inv, msg, 'log_file', log_file); { name of log file}
    Get_String(inv, msg, 'JV_file', JV_file); {name of file with simulated JV points}
    Get_String(inv, msg, 'Var_file', Var_file);{name of file with various variables at V=Vmax}
    Get_String(inv, msg, 'scPars_file', scPars_file);{solar cell parameter file}

    CLOSE(inv);
END;

PROCEDURE Check_Parameters;
{performs a number of checks on the parameters. Just to ensure that they are 
valid, consistent, make sense}
BEGIN
    {when adding new check, please keep the order in line with the device_parameter file}
    {Check first if Vt has been initialised. This should have happened in the main code}
    IF (Vt=0) THEN Stop_Prog('Vt needs to be initialised before calling Check_Parameters.');
    {checks on general parameters:}
    IF CB >= VB THEN Stop_Prog('CB should be smaller than VB.');
    IF (CB<0) OR (VB<0) THEN Stop_Prog('CB and VB should be positive.');

    {checks on mobilities:}
    IF NOT (mob_n_dep IN [0, 1, 2]) THEN Stop_Prog('Invalid mob_dep_n selected.');
    IF NOT (mob_p_dep IN [0, 1, 2]) THEN Stop_Prog('Invalid mob_dep_p selected.');
    {checks on contacts:}
    {part of this can only be done after checking the TLs!}
    IF Rshunt=0 THEN Stop_Prog('Rshunt cannot be zero, use positive (negative) value for finite (infinite) shunt resistance.');
    {checks on transport layers:}
    IF L_LTL+L_RTL>=L THEN Stop_Prog('Sum of L_LTL and L_RTL (transport layers) should be less than total thickness L.');
    {more checks on the energy levels of the TLs:}
    {first check left:}
    IF L_LTL>0 THEN
	BEGIN
		IF W_L<CB_LTL THEN Stop_Prog('W_L cannot be smaller than CB_LTL.');
		IF W_L>VB_LTL THEN Stop_Prog('W_L cannot be larger than VB_LTL.');
		IF CB_LTL>= VB_LTL THEN Stop_Prog('CB_LTL must be smaller than VB_LTL');
    END
    ELSE BEGIN {no TLs:}
		IF W_L<CB THEN Stop_Prog('W_L cannot be smaller than CB.');
		IF W_L>VB THEN Stop_Prog('W_L cannot be larger than VB.');
    END;
    {check the right:}
    IF L_RTL > 0 THEN
    BEGIN
		IF W_R<CB_RTL THEN Stop_Prog('W_R cannot be smaller than CB_RTL.');
		IF W_R>VB_RTL THEN Stop_Prog('W_R cannot be larger than VB_RTL.');
		IF CB_RTL>= VB_RTL THEN Stop_Prog('CB_RTL must be smaller than VB_RTL');
    END
    ELSE
    BEGIN
        IF W_R<CB THEN Stop_Prog('W_R cannot be smaller than CB.');
        IF W_R>VB THEN Stop_Prog('W_R cannot be larger than VB.');
    END;

    {check whether there are a possible number of traps (negative not allowed)}
    IF (St_L < 0) OR (St_R < 0) OR (Bulk_tr < 0) OR (GB_tr < 0) THEN Stop_Prog('Negative trap density not allowed.');

    {check trap types}
    IF NOT ((Tr_type_L <= 1) AND (Tr_type_L >= -1)) THEN Stop_Prog('Invalid left interface trap type.');
    IF NOT ((Tr_type_R <= 1) AND (Tr_type_R >= -1)) THEN Stop_Prog('Invalid right interface trap type.');
    IF NOT ((Tr_type_B <= 1) AND (Tr_type_B >= -1)) THEN Stop_Prog('Invalid bulk trap type.');        

    {check whether the level at which Etrap sits makes sense (i.e. falls within the bands everywhere in
    the device}
    IF Etrap <= CB THEN Stop_Prog('Trap energy level should be lower than the conduction band.');
    IF Etrap >= VB THEN Stop_Prog('Trap energy level should be higher than the valence band.');  

    IF Etrap <= CB_LTL THEN Stop_Prog('Trap energy level should be lower than the LTL conduction band.');
    IF Etrap <= CB_RTL THEN Stop_Prog('Trap energy level should be lower than the RTL conduction band.');    
    IF Etrap >= VB_LTL THEN Stop_Prog('Trap energy level should be higher than the LTL valence band.');
    IF Etrap >= VB_RTL THEN Stop_Prog('Trap energy level should be higher than the RTL valence band.');     
    
    {the following check is needed as long as we don't treat the barriers correctly:}
    IF ImageLowering AND ((L_LTL>0) AND (CB<>CB_LTL) OR (L_RTL>0) AND (VB<>VB_RTL)) 
		THEN Stop_Prog('Cannot use ImageLowering and TLs with CB/VB different from main layer at the same time.');
    {checks on ions:}
    IF NOT(ABS(mob_ion_spec) IN [0,1]) THEN Stop_Prog('Invalid mob_ion_spec selected.');
    {note: pascal set cannot contain negative numbers, hence the ABS()}	
    IF (ion_red_rate<0) THEN Stop_Prog('Scan rate cannot be lower than zero, should be: 0 <= ion_red_rate < NJV');
    {checks on generation and recombination parameters}
    IF NOT (ThermLengDist IN [1,2,3,4,5]) THEN Stop_Prog('Invalid ThermLengDist selected.');
    IF (P0>=1) OR (P0<0) THEN Stop_Prog('Invalid value of P0, should be: 0<=P0<1');
    IF (P0<>0) AND (Field_dep_G = FALSE) THEN Stop_Prog('P0 should be zero if not using field dependent generation');
    {checks on numerical parameters:}
    IF (NP<=5) OR (NP>Max_NP) THEN Stop_Prog('Invalid number of grid points (NP) selected, must be >=5 and <'+IntToStr(Max_NP)+'.');
    IF maxDelV<=0 THEN Stop_Prog('maxDelV should be positive.');
    {check if value of accPois makes any sense:}
    IF (accPois>=2) or (accPois<=0) THEN Stop_Prog('Invalid value of accPois selected.');
    {check if value of accDens makes any sense:}
    IF (accDens>=2) or (accDens<=0) THEN Stop_Prog('Invalid value of accDens selected.');  
    {checks on voltages:}
    IF NOT (Vdistribution IN [1,2]) THEN Stop_Prog('Invalid voltage distribution selected.');
    IF Vscan=0 THEN Stop_Prog('Vscan must be positive or negative, cannot be zero.');
    IF Vmin/Vt < -22657 THEN Stop_Prog('Vmin is too small.');
    IF Vmax/Vt > 22657 THEN Stop_Prog('Vmax is too big');
    {Note: upper and lower limits of voltage depend on G, these values seem reasonable though}
    {of course they also depend on the size of myReal. The numbers here are OK for doubles}
    IF Vmin > Vmax THEN Stop_Prog('Vmin should not be greater than Vmax.');
    {now check for redundancy of pre-bias:}
    IF PreCond THEN
    BEGIN
		IF Rseries>0 THEN WarnUser('Pre-bias voltage does not take Rseries into account, so Vpre=Vint.');
		IF CIM=0 THEN Stop_Prog('Do not use a pre-bias without any ions, makes no sense.');
		IF (Vscan=1) AND (Vpre=Vmin) THEN Stop_Prog('Pre-bias voltage is equal to Vmin, makes no sense.');
		IF (Vscan=-1) AND (Vpre=Vmax) THEN Stop_Prog('Pre-bias voltage is equal to Vmax, makes no sense.');
    END;
    IF (CIM<>0) AND (NOT (ion_red_rate in [0,1])) AND (Vdistribution =2) 
		THEN Stop_Prog('Do not use Vdistribution=2 with ion_red_rate other than 0 or 1.');
	IF (Vacc >= Vmin) AND (Vacc <= Vmax) AND (Vdistribution = 2)
		THEN {Vacc is not valid} Stop_Prog('Invalid Vacc selected, must be outside [Vmin, Vmax].');
    IF (Vdistribution=1) AND (Vstep <= 0) THEN Stop_Prog('Vstep should be positive.');	
    IF (ABS(Vmin-Vmax) < 1e-10) AND (Vdistribution=2) {to avoid infinite loop of Va}
		THEN Stop_Prog('Do not use Vdistribution=2 when Vmin = Vmax.');	
    {checks on user-interface:}
    IF (Gmax * Gfrac <> 0) AND (V0 <> VL) AND UseExpData AND (rms_mode=logarithmic) {this is a weird combination, warn user}
		THEN WarnUser('You are fitting a solar cell with rms_mode=log.');
    IF UseExpData AND until_Voc THEN
        Stop_Prog('You cannot use until_Voc = 1 and UseExpData = 1 at the same time.');
    IF ((rms_threshold<=0) OR (rms_threshold>1)) AND UseExpData THEN
        Stop_Prog('rms_threshold has to be larger than 0 but not larger than 1.');
END;

PROCEDURE DefineLayers(VAR NcLoc, ni, nid, pid, eps : vector; VAR i1, i2 : INTEGER);
VAR i : INTEGER;
BEGIN
    {define insulating layers:}
    {need i1, i2 to update the potential}
    i1:=0;
    WHILE x[i1+1]<L_LTL DO INC(i1);
	i2:=NP+1;
    WHILE x[i2-1]>L-L_RTL DO DEC(i2);
    {i1 is the last point in the left transport layer (or 0 if there isn't any)
    i2 is the first point in the right transport layer (or NP+1 if there is none)}

    {define different properties of the layers:}
    {first define bulk}
    FOR i:=0 TO NP+1 DO 
    BEGIN
		NcLoc[i]:=Nc; {bulk value of Nc}
		ni[i]:=Nc*EXP(-Egap/(2*Vt)); {equilibrium concentration}
		E_CB[i]:= CB;
		E_VB[i]:= VB;
		nid[i]:=n_0;
		pid[i]:=p_0;
		eps[i]:=eps_r * eps_0;
    END;    
    
    {now the transport layers}
    
    {left TL:}
    IF i1>0 THEN
	FOR i:=0 TO i1 DO
	BEGIN
	    NcLoc[i]:=Nc_LTL;
	    E_CB[i]:= CB_LTL;
	    E_VB[i]:= VB_LTL;
	    ni[i]:=Nc_LTL*EXP(-(VB_LTL-CB_LTL)/(2*Vt));
	    nid[i]:=0;
	    pid[i]:=0;	    
	    IF doping_LTL<=0 THEN nid[i]:=-doping_LTL ELSE pid[i]:=doping_LTL;
	    eps[i]:=eps_r_LTL * eps_0
	END;
    
    {right TL:}	
    IF i2<NP+1 THEN
	FOR i:=i2 TO NP+1 DO
	BEGIN
	    NcLoc[i]:=Nc_RTL;
	    E_CB[i]:= CB_RTL;
	    E_VB[i]:= VB_RTL;
	    ni[i]:=Nc_RTL*EXP(-(VB_RTL-CB_RTL)/(2*Vt));
	    nid[i]:=0;
	    pid[i]:=0;	    	    
	    IF doping_RTL<=0 THEN nid[i]:=-doping_RTL ELSE pid[i]:=doping_RTL;
	    eps[i]:=eps_r_RTL * eps_0;
	END;		

    
END;

PROCEDURE Make_Grid(VAR k, x : vector);
VAR i, ending : INTEGER;
    norm : myReal;
    {Makes a exponential symmetric grid, normalized to unity}
    {k[i] = (x[i+1] - x[i])/L}
    {and initialises the array with x-positions}
BEGIN
    FOR i:=0 TO NP+1 DO k[i]:=1;
    {note: we assign a value to k[np+1], but it doesn't have any meaning!}
    ending:=ROUND(NP/4); {defines the exponential part of the grid}
    FOR i:=0 TO ending DO
    BEGIN
		k[i]:=EXP(grad*(1/ending - 1/(i+1)));
		k[NP-i]:=k[i]
    END;
    norm:=0;
    FOR i:=0 TO NP DO norm:=norm + k[i];
    FOR i:=0 TO NP+1 DO k[i]:=k[i]/norm;
    x[0]:=0;    {Now calculate the positions}
    FOR i:=1 TO NP+1 DO x[i]:=x[i-1] + L*h[i-1];
END;

FUNCTION Applied_voltage(Vcount : INTEGER) : myReal;
{computes the applied voltage based on Vcount (=i)}
VAR i : INTEGER;

    FUNCTION LogarithmicV(i : INTEGER) : myReal;
    {computes a logarithmic V distribution, stepsize becomes zero at Vacc, so
    Vacc should lie outside [Vmin, Vmax] }
    VAR d : myReal;
    BEGIN
		d:=Vacc-Vmax;
		LogarithmicV:=Vacc - d*EXP((1-i/(NJV-1))*LN((Vacc-Vmin)/d))
    END;

BEGIN

    CASE Sign(Vscan) OF
      -1 : i:=NJV+1-Vcount; {scan down}
      1 : i:=Vcount {scan up}
    END;

    CASE Vdistribution of
      1 : Applied_voltage:=Vmin + Vstep*(i-1);
      2 : Applied_voltage:=LogarithmicV(i-1);
    END;
END;

PROCEDURE Read_Experimental_JV(VAR JVExp : JVList; VAR NJV : INTEGER; VAR log : TEXT);
{reads the experiment JV curve and determines Vmin, Vmax, Vstep, and Vscan}
{note: This overrides the Vmin,max,step,scan,distribution in the parameter file AND the command line.}
{we could have used a linked list instead...}
VAR i : INTEGER;
    inp : TEXT;
    dumstr : STRING;
    V, J : ARRAY[1..maxExpData] OF myReal;
BEGIN
    IF NOT FileExists(ExpJV) THEN Stop_Prog('Could not find file '+ExpJV, FALSE);
    ASSIGN(inp, ExpJV);
    RESET(inp);
    
    {first try to read a header}
    TRY 
		READLN(inp, dumstr);
    EXCEPT
		Stop_Prog('Cannot read a header from '+ExpJV+', is it empty?', FALSE);
    END;

    i:=0;
    WHILE (NOT EOF(inp)) AND (i < maxExpData) DO
    BEGIN
        INC(i);
        TRY
            READLN(inp, V[i], J[i]);
        EXCEPT
            Stop_Prog('The experimental JV curve in '+ExpJV+' can only contain voltage and current density, etc.', FALSE);
        END
    END;

    CLOSE(inp);

    {now check a number of things:}
    NJV:=i; {this overrides NJV calculated from Vmin,Vmax,Vstep, or from direct input}
    IF NJV < minExpData THEN Stop_Prog('Not enough experimental data points.', FALSE);
    IF NJV = maxExpData THEN Stop_Prog('It looks like the experimental JV file contains more points than I can handle.', FALSE);

    {now we can copy the arrays V and J into ExpJV}
    SetLength(JVExp, NJV+1); {the regular voltages start at 1}
    FOR i:=1 TO NJV DO {note, we only use part of the JVList record}
    BEGIN
        JVExp[i].Vint:=V[i];
        JVExp[i].Vext:=V[i];
        JVExp[i].Jint:=J[i];
        JVExp[i].Jext:=J[i];
        JVExp[i].Use:=TRUE
    END;

    {Vmin, Vmax: we will need these later on}
    Vmin:=MIN(JVexp[1].Vint, JVexp[NJV].Vint); {we use min, max function as Vscan might be -1}
    Vmax:=MAX(JVexp[1].Vint, JVexp[NJV].Vint);

    {now document what happened and write to log file:}
    WRITELN(log);
    WRITELN(log, 'Read experiment JV curve from ', ExpJV,'.');
    WRITELN(log, 'This overrides the voltage distribution that is in ',parameter_file);
    WRITELN(log, 'and any voltage parameters passed via the command line.');
    WRITELN(log, 'Vmin: ',Vmin:6:4,' Vmax: ',Vmax:6:4);
    WRITELN(log);
    WRITELN('Read experimental JV curve from ', ExpJV,'.');
END;



PROCEDURE Init_Voltages_and_Tasks(VAR JVSim, JVExp : JVList; VAR log : TEXT);
VAR i, k : INTEGER;
BEGIN
    IF UseExpData THEN Read_Experimental_JV(JVExp, NJV, log); {calling this overrides NJV calculated from Vmin,Vmax,Vstep, or from direct input}
    
    SetLength(JVSim, NJV+1); {if there is a pre-bias, we need an extra data point, so we make the array one longer than NJV}
    IF PreCond THEN {the pre-bias will be stored in point 0}
	WITH JVSim[0] DO {pre-bias point}
	BEGIN
	    Vint:=Vpre;
	    UpdateIons:=TRUE; {yes, ions are moving. We checked in Read_Parameters that if PreCond then CIM <>0}
	    Store:=FALSE; {however, don't store the pre-bias point}
	END;
    
    k:=ORD(NOT PreCond); {so if PreCond=true then k=0, we will use this below}
    {now for the rest of the voltages:}
    FOR i:=1 TO NJV DO {the regular voltages start at 1}
    BEGIN
		WITH JVSim[i] DO
		BEGIN
			IF UseExpData
				THEN Vint:=JVExp[i].Vint
				ELSE Vint:=Applied_voltage(i);
			IF ion_red_rate=0 THEN {means ions are fixed}
            BEGIN
				UpdateIons:=(i=1) AND (NOT PreCond);
				{so if the ions are fixed, but we're doing the first point (i=1)
				then update ions unless we already did a preconditioning}
				Store:=TRUE {by default, store this JV point. If this point does not converge, we will set Store to FALSE}
			END
			ELSE BEGIN {so ions are not always fixed}
				UpdateIons:=(CIM<>0) AND ( (i-k) MOD ion_red_rate = 0 );
				{note the k! whether we update the ions also depends on the pre-conditioning}
				Store:=(CIM=0) OR ( (i-k+1) MOD ion_red_rate = 0)
			END;
		END;
    END;
    
    {now write to log to show what we are going to do:}
    WRITELN(log);
    WRITELN(log,'The following voltages will be simulated:');
    WRITELN(log,'  i       V    UpdateIons   Store');
    FOR i:=ORD(NOT PreCond) TO NJV DO {note: if not PreCond we start at i=1}
	WITH JVSim[i] DO
	    WRITELN(log,i:3,'  ',Vint:8:3,'    ',UpdateIons:5,'      ',Store);
    WRITELN(log);
    FLUSH(log);	
END;

PROCEDURE Init_Elec_Mob_Table(VAR n_vals, Fn_vals : row; VAR n_points, Fn_points : INTEGER; VAR mob_tab : Table);
{reads the elec. mobility as a function of F and n from input file n_file}
VAR i, j : INTEGER;
    dum : myReal;
    input : TEXT;
    dumstr : STRING[3];
    ascending : boolean;
BEGIN
    {first check if file exists:}
    IF (mob_n_dep = 2) AND (NOT FileExists(n_file)) THEN Stop_Prog('Could not find file '+n_file);
    {open file:}
    ASSIGN(input, n_file);
    RESET(input);
    {get the parameters:}
    Read_Number(input, dum);  {read number of points in n direction}
    n_points:=ROUND(dum);
    Read_Number(input, dum);  {read number of points in F direction}
    Fn_points:=ROUND(dum);

    SetLength(n_vals, n_points+1); {SetLength makes the array start at 0, but we'll start at 1}
    SetLength(Fn_vals, Fn_points+1);
    SetLength(mob_tab, n_points+1, Fn_points+1);

    READ(input, dumstr);
    IF TRIM(dumstr) <> 'n\F' THEN Stop_Prog('Not sure about table in file '+n_file+' as it does not have format n\F');

    {read the values of F}
    FOR i:=1 TO Fn_points DO
		READ(input, Fn_vals[i]);
    READLN(input);

    {now check wheter F values are ascending}
    ascending:=TRUE;
    FOR i:=1 TO Fn_points-1 DO
		ascending:=ascending AND (Fn_vals[i] < Fn_vals[i+1]);
    IF NOT ascending THEN stop_prog('F_values not ascending');

    {now read table}
    FOR i:=1 TO n_points DO
    BEGIN
		{1st number is the density:}
		READ(input, dum);
		n_vals[i]:=Nc * dum; {note: here we use Nc bulk value as the mobility table is only valid for the main absorber layer, not the TLs}
		{now read the corresponding mobilities:}
		FOR j:=1 TO Fn_points DO
		BEGIN
			READ(input, mob_tab[i,j]);
			IF extra_ln_mob THEN mob_tab[i,j]:=ln(mob_tab[i,j]) {if extra_ln_mob then take the ln of the mob}
		END;
		READLN(input);  {and go to next line}
    END;

    {now check wheter n values are ascending}
    ascending:=TRUE;
    FOR i:=1 TO n_points-1 DO
		ascending:=ascending AND (n_vals[i] < n_vals[i+1]);
    IF NOT ascending THEN stop_prog('n_values not ascending');

    CLOSE(input);
    WRITELN('Read MOB_TAB');
    IF extra_ln_mob THEN WRITELN('USING LN(MOB) FOR INTER-/EXTRAPOLATION');
END;

PROCEDURE Init_Hole_Mob_Table(VAR p_vals, Fp_vals : row; VAR p_points, Fp_points : INTEGER; VAR mob_tab : Table);
{reads the hole mobility as a function of F and p from input file p_file}
VAR i, j : INTEGER;
    dum : myReal;
    input : TEXT;
    dumstr : STRING[3];
    ascending : boolean;
BEGIN
    {first check if file exists:}
    IF (mob_p_dep = 2) AND (NOT FileExists(p_file)) THEN Stop_Prog('Could not find file '+p_file);
    {open file:}
    ASSIGN(input, p_file);
    RESET(input);
    {get the parameters:}
    Read_Number(input, dum);  {read number of points in p direction}
    p_points:=ROUND(dum);
    Read_Number(input, dum);  {read number of points in F direction}
    Fp_points:=ROUND(dum);

    SetLength(p_vals, p_points+1); {SetLength makes the array start at 0, but we'll start at 1}
    SetLength(Fp_vals, Fp_points+1);
    SetLength(mob_tab, p_points+1, Fp_points+1);

    READ(input, dumstr);
    IF TRIM(dumstr) <> 'p\F' THEN Stop_Prog('Not sure about table in file '+p_file+' as it does not have format p\F');

    {read the values of F}
    FOR i:=1 TO Fp_points DO
		READ(input, Fp_vals[i]);
    READLN(input);

    {now check wheter F values are ascending}
    ascending:=TRUE;
    FOR i:=1 TO Fp_points-1 DO
		ascending:=ascending AND (Fp_vals[i] < Fp_vals[i+1]);
    IF NOT ascending THEN stop_prog('F_values not ascending');

    {now read table}
    FOR i:=1 TO p_points DO
    BEGIN
		{1st number is the density:}
		READ(input, dum);
		p_vals[i]:=Nc * dum; {note: here we use Nc bulk value as the mobility table is only valid for the main absorber layer, not the TLs}
		{now read the corresponding mobilities:}
		FOR j:=1 TO Fp_points DO
		BEGIN
			READ(input, mob_tab[i,j]);
			IF extra_ln_mob THEN mob_tab[i,j]:=LN(mob_tab[i,j]) {if extra_ln_mob then take the ln of the mob}
		END;
		READLN(input);  {and go to next line}
    END;

    {now check wheter p values are ascending}
    ascending:=TRUE;
    FOR i:=1 TO p_points-1 DO
		ascending:=ascending AND (p_vals[i] < p_vals[i+1]);
    IF NOT ascending THEN stop_prog('n_values not ascending');

    CLOSE(input);
    WRITELN('Read p_mob_tab');
    IF extra_ln_mob THEN WRITELN('USING LN(MOB) FOR INTER-/EXTRAPOLATION');
END;


PROCEDURE Trap_Distribution();
{Trap distribution along the thickness of the device}
VAR i, gb_nr, cur_idx : INTEGER;
    gb_centre_pos, inv_int_trap_width, cwe_dummy : myReal;
    new_gb_found : BOOLEAN;

{grain_size,dum : myReal;}
BEGIN    
    {we concider a trap filled when it contains an electron}
    {a donor type trap is changed when empty, and an acceptor trap is charged when filled.}
    CASE Tr_type_B OF
      -1 : cwe_b := 0;
      0	 : cwe_b := 0; {we will not use this variable in this case, but it should have a value regardless}
      1	 : cwe_b := 1;
    END;
      
    FOR i:=0 TO NP+1 DO Ntb[i] := Bulk_tr;
    FOR i :=0 TO NP+1 DO charged_traps[i]:= 1;    
    
    IF (St_L > 0) THEN {if there is traps at the left interface}
    BEGIN
		inv_int_trap_width := 2 / (L * (0.5 * h[i1-1] + h[i1] + 0.5 * h[i1+1]));
	
		Nti[i1] := St_L * inv_int_trap_width;
		Nti[i1+1] := St_L * inv_int_trap_width;

        CASE Tr_type_L OF
			-1 : cwe_dummy := 0;
			0  : BEGIN cwe_dummy := 0; charged_traps[i1]:= 0; charged_traps[i1+1]:= 0 END; {we will not use this variable in this case, but it should have a value regardless}
			1  : cwe_dummy := 1;
		END;
		cwe_i[i1] := cwe_dummy;
		cwe_i[i1+1] := cwe_dummy;
    END;
    
    IF (St_R > 0) THEN {if there is traps at the right interface}
    BEGIN	
		inv_int_trap_width := 2 / (L * (0.5*h[i2] + h[i2-1] + 0.5*h[i2-2]));
	
		Nti[i2-1] := St_R * inv_int_trap_width;
		Nti[i2] := St_R * inv_int_trap_width;

        CASE Tr_type_R OF
			-1 : cwe_dummy := 0;
			0  : BEGIN cwe_dummy := 0; charged_traps[i2]:= 0; charged_traps[i2-1]:= 0 END;{we will not use this variable in this case, but it should have a value regardless}
			1  : cwe_dummy := 1;
		END;
	
		cwe_i[i2-1] := cwe_dummy;
		cwe_i[i2] := cwe_dummy;
    END;
    
    IF (num_GBs > 0) AND (GB_tr <> 0) THEN {include grain boundary recombination}
    BEGIN
		cur_idx := i1;
        CASE Tr_type_B OF
			-1 : cwe_dummy := 0;
			0  : BEGIN cwe_dummy := 0; FOR i:=i1+2 TO i2-2 DO charged_traps[i]:= 0 END; {we will not use this variable in this case, but it should have a value regardless}
			1  : cwe_dummy := 1;	
		END;
		
		FOR gb_nr := 1 TO num_GBs DO
		BEGIN 
			new_gb_found := FALSE;
			gb_centre_pos := (L - L_RTL - L_LTL) / (num_GBs + 1) * gb_nr + L_LTL;

			{find the start and end idx of the grain boundary and set trap density}
			WHILE new_gb_found = FALSE DO 
			BEGIN 
				IF (gb_centre_pos >= x[cur_idx]) AND (gb_centre_pos < x[cur_idx+1]) THEN 
				BEGIN 
					new_gb_found := TRUE;
					inv_int_trap_width := 2 / (L * (0.5 * h[cur_idx-1] + h[cur_idx] + 0.5 * h[cur_idx+1]));
					Nti[cur_idx] := abs(GB_tr) * inv_int_trap_width;
					Nti[cur_idx+1] := abs(GB_tr) * inv_int_trap_width;
		    
					cwe_i[cur_idx] := cwe_dummy;
					cwe_i[cur_idx+1] := cwe_dummy;
				END;
				inc(cur_idx);
			END
		END   
    END;

    {There may not be two touching interfaces, as this is conflicting with how the equations are derived.}
    FOR i := 1 TO NP DO
	IF (Nti[i-1] <> 0) AND (Nti[i] <> 0) AND (Nti[i+1] <> 0) THEN
	    Stop_Prog('There are multiple consecutive interfaces, decrease interface density.');
    
    IF NOT TLsTrap THEN {whether traps in transport layers:}
    BEGIN
		{left transport layer}
		IF i1 > 0 THEN
			FOR i := 0 TO i1 DO
				Ntb[i] := 0;
		{right transport layer}
		IF i2 < NP+1 THEN
			FOR i:=i2 TO NP+1 DO
				Ntb[i] := 0;
    END;
END;

PROCEDURE Calc_Ion_Distribution(VAR nion, pion, V : vector);
{calculates the ion distribution in steady state. This is either constant 
(if the ionic species doesn't move), or it follows from requiring that
the associated ionic particle current be zero.} 
VAR nrmn, nrmp, fac, sum_nion, sum_pion : myReal;
    istart, ifinish : INTEGER;
BEGIN
    sum_nion:=0; {sums of ion concentration}
    sum_pion:=0;
    
    {ions are limited to the middle layer (which can take up the whole device i=0...NP+1.}
    IF (i1>0) AND (NOT IonsInTLs) THEN istart:=i1+1 ELSE istart:=0;
    IF (i2<NP+1) AND (NOT IonsInTLs) THEN ifinish:=i2-1 ELSE ifinish:=NP+1;
    
    nion[istart]:=1;
    pion[istart]:=1;

    FOR i:=istart+1 TO ifinish DO
    BEGIN
		fac:=B(Vti*(V[i-1]-V[i]))/B(Vti*(V[i]-V[i-1]));
		{if an ionic species moves, then we use the expressions for the 
		ion currents to make sure that they are zero. If they don't move,
		then we simply set them to 1 and normalise later. This yields the profile of neg/pos ions.}
		CASE mob_ion_spec OF
			-1 : BEGIN nion[i]:=nion[i-1]*fac; pion[i]:=1 END;
			0 : BEGIN nion[i]:=nion[i-1]*fac; pion[i]:=pion[i-1]/fac; END;
			1 : BEGIN nion[i]:=1; pion[i]:=pion[i-1]/fac; END;
		END;
	END;
    
    {nomalize concentrations}
    FOR i:=istart+1 TO ifinish DO {note: we start at istart + 1 as we access ion[i-1]}
    BEGIN
		{if the grid is non-uniform, we need to take this into account}
		sum_nion:=sum_nion + 0.5*(nion[i]+nion[i-1])*h[i-1]*L; 
		sum_pion:=sum_pion + 0.5*(pion[i]+pion[i-1])*h[i-1]*L; 
    END;
    
    IF IonsInTLs THEN
    BEGIN
		nrmn:=CIM*L/sum_nion;
        nrmp:=CIM*L/sum_pion;
    END
    ELSE
    BEGIN
		nrmn:=CIM*(L-L_LTL-L_RTL)/sum_nion;
        nrmp:=CIM*(L-L_LTL-L_RTL)/sum_pion;
    END;
	
    FOR i:=istart TO ifinish DO
    BEGIN
		nion[i]:=nion[i]*nrmn;
		pion[i]:=pion[i]*nrmp
    END;
    {now the total number of ions should be equal to the concentrion CMI times (L-L_LTL-L_RTL)}
END;


FUNCTION Calc_f_tb() : vector;
VAR i : INTEGER;
BEGIN
    FOR i := 0 TO NP+1 DO
    BEGIN
		IF (Ntb[i] <> 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
			Calc_f_tb[i] := (Cn*n[i] + Cp*pt0[i]) / (Cn*(n[i]+nt0[i]) + Cp*(p[i]+pt0[i]))
		ELSE
			Calc_f_tb[i] := 0
    END;

END;



PROCEDURE Calc_SRH_bulk_prefactor(VAR srhp : vector; V, n, p : vector);
{Calculates the denominator of the occupational probabilities of bulk traps
as in equation 4.2 of Shockley and Read Physical Review 1952, divided by Cn*Cp
Previously called Shockley read hall prefactor. In can be used to calculate 
bulk trap SRH recombination by R_SRH = srhp * (np - n_int^2), where R_SRH 
is the SRH recombination rate.}
VAR i : INTEGER;
BEGIN
    IF (Cn= 0) OR (Cp = 0) OR (NOT Traps)
	THEN FOR i:=0 TO NP+1 DO srhp[i]:=0 {no SRH recombination}
	ELSE
	BEGIN 
	    FOR i:=1 TO NP DO srhp[i]:=(Ntb[i])/ ( (n[i] + nt0[i]) /Cp + (p[i] + pt0[i])/Cn);
	END;
END;


FUNCTION Calc_inv_f_ti_denom(): vector;
VAR i, j	     : INTEGER;
    denom_trap_occup : myReal;
BEGIN
    Calc_inv_f_ti_denom[0] := 0;
    Calc_inv_f_ti_denom[NP+1] := 0;    
    FOR i := 1 TO NP DO
    BEGIN
		IF (Nti[i] <> 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
		BEGIN
			denom_trap_occup:=0;
			FOR j:=-1 TO 1 DO denom_trap_occup := denom_trap_occup + nu_ni * (Nti[i+j] * n[i+j] + Nti[i+j] * nt0[i+j]) + nu_pi * (Nti[i+j] * p[i+j] + Nti[i+j] * pt0[i+j]);
			calc_inv_f_ti_denom[i] := 1/denom_trap_occup;
		END
		ELSE Calc_inv_f_ti_denom[i] := 0;
    END;
END;


FUNCTION Calc_f_ti(inv_denom_SRH_int : vector): vector;
VAR i, j : INTEGER;
    f_t	 : myReal;
BEGIN
    Calc_f_ti[0] := 0;
    Calc_f_ti[NP+1] := 0;
    
    FOR i := 1 TO NP DO
    BEGIN
		IF Nti[i] <> 0 THEN
		BEGIN
			f_t := 0;
			FOR j:=-1 TO 1 DO f_t := f_t + nu_ni * (Nti[i+j] * n[i+j]) + nu_pi * (Nti[i+j] * pt0[i+j]);
			Calc_f_ti[i] := f_t * inv_denom_SRH_int[i];
		END
		ELSE Calc_f_ti[i] := 0;
    END;   
END;


FUNCTION Calc_f_ti_numer(): vector;
VAR i, j : INTEGER;
    f_t	 : myReal;
BEGIN
    Calc_f_ti_numer[0] := 0;
    Calc_f_ti_numer[NP+1] := 0;
    
    FOR i := 1 TO NP DO
    BEGIN
		IF Nti[i] <> 0 THEN
		BEGIN
			f_t := 0;
			FOR j:=-1 TO 1 DO f_t := f_t + nu_ni * (Nti[i+j] * n[i+j]) + nu_pi * (Nti[i+j] * pt0[i+j]);
			Calc_f_ti_numer[i] := f_t;
		END
		ELSE Calc_f_ti_numer[i] := 0;
    END;   
END;


PROCEDURE Solve_Poisson(VAR V, n, p, nid, pid : vector; UpdateIons : BOOLEAN; VAR conv : BOOLEAN);
VAR it, i											: INTEGER;
    delV, rhs, lower, upper, main, Ntb_nr_charges, f_ti, f_tb, inv_denom_trap_occup, f_ti_numer	: vector;
    fac, f_ti_first_ord_lo, f_ti_first_ord_up, f_ti_first_ord_main, f_tb_numer, f_tb_inv_denom,
    f_tb_first_ord_main                                                                         : myReal;
BEGIN
    FOR i:=1 TO NP DO delV[i]:=1; {init delV}
    {    FOR i:=1 TO NP DO Vnew[i]:=1; {init delV}
    delV[0]:=0; {delV=0 at contacts, by definition, since we know V(0,L)}
    delV[NP+1]:=0;
    it:=0; 

    FOR i:=0 TO NP+1 DO {init vectors}
    BEGIN
		f_tb[i] := 0;
		f_ti[i] := 0;
		Ntb_nr_charges[i] := 0;
		Nti_charge[i] := 0;
    END;

    
    WHILE (Norm(delV, 0, NP+1) > tolPois) AND (it < MaxItPois) DO
    BEGIN
		IF (bulk_tr <> 0) AND NOT (Tr_type_B = 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN f_tb := Calc_f_tb;
		IF (Traps_int) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
		BEGIN
			inv_denom_trap_occup := Calc_inv_f_ti_denom();
			f_ti_numer := Calc_f_ti_numer();
			f_ti := Calc_f_ti(inv_denom_trap_occup);
		END;
	
		{init traps stuff (f_tb and inv denom)}
		FOR i:=1 TO NP DO
		BEGIN
			{non-interface traps, so the ones that do not capture from neighbouring points.}
			{implement using eq 4.2 and 4.3 from shockley and read paper}
			IF (Ntb[i] <> 0) AND NOT (Tr_type_B = 0) THEN Ntb_nr_charges[i] := (cwe_b - f_tb[i]) * Ntb[i]; {note that we have f_tb minus charges or 1-f_tb positive charges} 	    
			IF (Traps_int) THEN Nti_charge[i] := (cwe_i[i] - f_ti[i]) * 0.5 * Nti[i]; {each side of the interface hosts half the traps}

			{linearize f_ti in delV}
			IF (Nti[i-1] = 0) THEN
				f_ti_first_ord_lo := 0
			ELSE
			BEGIN
				f_ti_first_ord_lo := Cn * n[i-1] * Nti[i-1] * Vti * inv_denom_trap_occup[i];
				f_ti_first_ord_lo := f_ti_first_ord_lo - (Cn * n[i-1] - Cp * p[i-1]) * Nti[i-1] * Vti * f_ti_numer[i] * SQR(inv_denom_trap_occup[i]);
				f_ti_first_ord_lo := f_ti_first_ord_lo * 0.5 * Nti[i-1];
			END;

			IF (Nti[i] = 0) THEN
				f_ti_first_ord_main := 0
			ELSE
			BEGIN
				{the factor Vti is multiplied in the poisson equation with for this main diagonal}
				f_ti_first_ord_main := Cn * n[i] * Nti[i] * inv_denom_trap_occup[i];
				f_ti_first_ord_main := f_ti_first_ord_main - (Cn * n[i] - Cp * p[i]) * Nti[i] * f_ti_numer[i] * SQR(inv_denom_trap_occup[i]);
				f_ti_first_ord_main := f_ti_first_ord_main * 0.5 * Nti[i];
			END;

			IF (Nti[i+1] = 0) THEN
				f_ti_first_ord_up := 0
			ELSE
			BEGIN
				f_ti_first_ord_up := Cn * n[i+1] * Nti[i+1] * Vti * inv_denom_trap_occup[i];
				f_ti_first_ord_up := f_ti_first_ord_up - (Cn * n[i+1] - Cp * p[i+1]) * Nti[i+1] * Vti * f_ti_numer[i] * SQR(inv_denom_trap_occup[i]);
				f_ti_first_ord_up := f_ti_first_ord_up * 0.5 * Nti[i+1];
			END;

    	    {linearize f_tb in delV}
			f_tb_numer := Cn*n[i] + Cp*pt0[i];
			IF NOT (Tr_type_B = 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
			BEGIN
				f_tb_inv_denom := 1 / (Cn*(n[i]+nt0[i]) + Cp*(p[i]+pt0[i]));
				f_tb_first_ord_main := Cn * n[i+1] * f_tb_inv_denom;
				f_tb_first_ord_main := f_tb_first_ord_main - (Cn * n[i] - Cp * p[i]) * f_tb_numer * SQR(f_tb_inv_denom);
				f_tb_first_ord_main := (cwe_b - f_tb_first_ord_main) * Ntb[i];
			END
			ELSE
				f_tb_first_ord_main := 0;

			{ this factor repeats often in the equations}
			fac :=  h[i]*h[i-1]*L*(h[i] + h[i-1]) * L * 0.5 * q;
	    
			rhs[i] := fac * (n[i] + nion[i] + pid[i] -p[i] - pion[i] - nid[i] - Nti_charge[i] * charged_traps[i] - ABS(Tr_type_B) * Ntb_nr_charges[i])
					- eps[i-1]*h[i]*V[i-1] - eps[i]*h[i-1]*V[i+1] + (eps[i]*h[i-1] + eps[i-1]*h[i])*V[i];
			lower[i]:= eps[i-1]*h[i] - fac * f_ti_first_ord_lo * charged_traps[i];
			upper[i]:= eps[i]*h[i-1] - fac * f_ti_first_ord_up * charged_traps[i];
			{ while unintuitive, the main diagonal adds all charges regardless of sign. This comes from the difference in partial derivatives with respect to delV.}
			main[i]:=-(h[i-1]*eps[i]+h[i]*eps[i-1]) - fac * (n[i]+p[i] + nion[i]+pion[i] + f_ti_first_ord_main * charged_traps[i] + ABS(Tr_type_B) * f_tb_first_ord_main)*Vti;
		END;

		Tridiag(delV, lower, main, upper, rhs, 1, NP); {solve for delV}
		FOR i:=1 TO NP DO 
		BEGIN
			delV[i]:=accPois*delV[i]; {use SOR/SUR}	    
			IF ABS(delV[i])>maxDelV*Vt THEN delV[i]:=SIGN(delV[i])*maxDelV*Vt; 
			V[i]:=V[i]+delV[i];  {and add delV to V}
			n[i]:=n[i]*EXP(delV[i]*Vti); {now update the densities: we have to do this}
			p[i]:=p[i]*EXP(-delV[i]*Vti); {in order to conserve Gummel iteration}
			IF UpdateIons THEN
			CASE mob_ion_spec OF
				-1: nion[i]:=nion[i]*EXP(delV[i]*Vti);
				0 : BEGIN
						nion[i]:=nion[i]*EXP(delV[i]*Vti);
						pion[i]:=pion[i]*EXP(-delV[i]*Vti);
					END;
				1 : pion[i]:=pion[i]*EXP(-delV[i]*Vti);
			END		
		END;
		it:=it+1;
    END;
    conv:= Norm(delV, 0, NP+1) <= tolPois  {finally, check for convergence}
END;


PROCEDURE UpdateGenPot(V, NcLoc : vector; VAR Vgn, Vgp : vector);
VAR i : INTEGER;
    facDOS : myReal;
BEGIN
    Vgn:=V;
    Vgp:=V;
    {left transport layer, offset generalised potentials:}
    IF i1>0 THEN
	FOR i:=0 TO i1 DO
	BEGIN
	    facDOS:=Vt*LN(NcLoc[i]/Nc);
	    Vgn[i]:=V[i] + CB_LTL - CB + facDOS;
	    Vgp[i]:=V[i] - VB + VB_LTL - facDOS;
	END;
    {right transport layer, offset generalised potentials:}
    IF i2<NP+1 THEN
	FOR i:=i2 TO NP+1 DO
	BEGIN
	    facDOS:=Vt*LN(NcLoc[i]/Nc);
	    Vgn[i]:=V[i] + CB_RTL - CB + facDOS;
	    Vgp[i]:=V[i] - VB + VB_RTL - facDOS;
	END;
END;

PROCEDURE Init_Pot_Dens(VAR V, Vgn, Vgp, n, p, empty : vector; Va : myReal);
{init. for V, Vgn,p, n, p at bias voltage Va}
VAR i : INTEGER;
    facDOS : myReal;
BEGIN
    FOR i:=0 TO NP+1 DO {guess new V in interior of device} {linear interpolation of V}
		V[i]:=V0 - Va/2 + x[i]*(VL-V0+Va)/L;
    UpdateGenPot(V, NcLoc, Vgn, Vgp); {update generalised potentials}

    FOR i:=0 TO NP+1 DO {guess n, p in interior of device: not critical!}
    BEGIN
		n[i]:=ni[i]*EXP( (Va * (0.5 -x[i]/L) + Vgn[i])*Vti ); {note: we use Vgn,p, NOT V}
		p[i]:=ni[i]*EXP( (Va * (x[i]/L - 0.5) - Vgp[i])*Vti );
		empty[i]:=0;  {Empty vector, used to separately calculate Lan & SRH Recombination}
    END;
    
    {now we set the boundary conditions, these are the ones that matter:} 
    {left, i=0:}
    facDOS:=Vt*LN(NcLoc[0]/Nc);
    n[0]:=NcLoc[0]*EXP((Vgn[0]-V[0]+CB-facDOS-W_L)*Vti);
    p[0]:=SQR(ni[0])/n[0]; {use law of mass action for p}
    
    {right, i=NP+1:}
    facDOS:=Vt*LN(NcLoc[NP+1]/Nc);
    p[NP+1]:=NcLoc[NP+1]*EXP((W_R+V[NP+1]-Vgp[NP+1]-VB-facDOS)*Vti);
    n[NP+1]:=SQR(ni[NP+1])/p[NP+1]; 
    
END;

PROCEDURE Calc_elec_mob(VAR mu : vector; V, n : vector);
VAR i : INTEGER;
    conc, field: myReal; {conc: concentation on x=i+1/2}
    {calculates the elec. mob. on the interleaved mesh, mun[i]=mun at x=i+1/2}
    {therefore the field or concentration on x=i+1/2 is needed}
BEGIN
    CASE mob_n_dep OF
      0 : FOR i:=0 TO NP DO mu[i]:=mun_0; { mob. is constant }
      1 : FOR i:=0 TO NP DO {field-dep. mob}
	      mu[i]:=mun_0 * EXP(gamma_n*SQRT(ABS(V[i+1]-V[i])/(L*h[i])));
      2 : BEGIN {electron mobility from table in file n_file}
			FOR i:=0 TO NP DO
			BEGIN
				conc:=0.5*(n[i] + n[i+1]);
				{note: n, V are defined on grid, but mu on interleaved mesh!}
				Field:=ABS(V[i+1]-V[i])/(L*h[i]);
				IF extra_ln_mob {decide whether to extrapolate ln(mob) or just mob}
					THEN mu[i]:=EXP(BilinearInterpolation(conc, Field, n_vals, Fn_vals, n_points, Fn_points, n_mob_tab, extra_F))
					ELSE mu[i]:=BilinearInterpolation(conc, Field, n_vals, Fn_vals, n_points, Fn_points, n_mob_tab, extra_F);
				END; {for loop}
			END; {case 2}
    END; {case selector}

    mu[NP+1]:=mu[NP]; {does not have a meaning, should not be used}
    IF i1>0 THEN
    BEGIN
		FOR i:=0 TO i1-1 DO mu[i]:=mob_LTL;
		{remember: i1 is the last point (highest index) in the right TL}
		mu[i1]:=(L*h[i1])*nu_int_LTL*Vti; {mobility AT the TL/main absorber interface}
    END;
    IF i2<NP+1 THEN
    BEGIN
		mu[i2-1]:=(L*h[i2-1])*nu_int_RTL*Vti; {mobility AT the TL/main absorber interface}
		{remember: i2 is the first point (lowest index) in the right TL}
		FOR i:=i2 TO NP+1 DO mu[i]:=mob_RTL;
    END
END;

PROCEDURE Calc_hole_mob(VAR mu : vector; V, p : vector);
VAR i : INTEGER;
    conc, field: myReal; {conc: concentration on x=i+1/2}
    {calculates the hole mob. on the interleaved mesh, mup[i]=mup at x=i+1/2}
    {therefore the field or concentration on x=i+1/2 is needed}
BEGIN
    CASE mob_p_dep OF
      0 : FOR i:=0 TO NP DO mu[i]:=mup_0; { mob. is constant }
      1 : FOR i:=0 TO NP DO {field-dep. mob}
	      mu[i]:=mup_0 * EXP(gamma_p*SQRT(ABS(V[i+1]-V[i])/(L*h[i])));
      2 : BEGIN {hole mobility from table in file p_file}
			FOR i:=0 TO NP DO
			BEGIN
				conc:=0.5*(p[i] + p[i+1]);
				{note: p, V are defined on grid, but mu on interleaved mesh!}
				Field:=ABS(V[i+1]-V[i])/(L*h[i]);
				IF extra_ln_mob {decide whether to extrapolate ln(mob) or just mob}
					THEN mu[i]:=EXP(BilinearInterpolation(conc, Field, p_vals, Fp_vals, p_points, Fp_points, p_mob_tab, extra_F))
					ELSE mu[i]:=BilinearInterpolation(conc, Field, p_vals, Fp_vals, p_points, Fp_points, p_mob_tab, extra_F);
			END; {for loop}
		END; {case 2}
    END; {case selector}
    mu[NP+1]:=mu[NP]; {does not have a meaning, should not be used}

    IF i1>0 THEN
    BEGIN
		FOR i:=0 TO i1-1 DO mu[i]:=mob_LTL;
		{remember: i1 is the last point (highest index) in the right TL}
		mu[i1]:=(L*h[i1])*nu_int_LTL*Vti; {mobility AT the TL/main absorber interface}
    END;
    IF i2<NP+1 THEN
    BEGIN
		mu[i2-1]:=(L*h[i2-1])*nu_int_RTL*Vti; {mobility AT the TL/main absorber interface}
		{remember: i2 is the first point (lowest index) in the right TL}
		FOR i:=i2 TO NP+1 DO mu[i]:=mob_RTL;
    END
END;

PROCEDURE Calc_Dn(VAR D : vector; mu, n : vector);
{Calculates the diffusion constant, just Einstein for now, Dn[i]= Dn at x=i+1/2 }
VAR i : INTEGER;
BEGIN
    CASE Dn_dep OF
      0 : FOR i:=0 TO NP DO D[i]:=mu[i]*Vt;
      1 : FOR i:=0 TO NP DO D[i]:=eta_n*mu[i]*Vt
	  ELSE Stop_Prog('Invalid Dn_dep selected.');
    END;
END;

PROCEDURE Calc_Dp(VAR D : vector; mu, p : vector);
{Calculates the diffusion constant, just Einstein for now, Dn[i]= Dn at x=i+1/2 }
VAR i : INTEGER;
BEGIN
    CASE Dp_dep OF
      0 : FOR i:=0 TO NP DO D[i]:=mu[i]*Vt;
      1 : FOR i:=0 TO NP DO D[i]:=eta_p*mu[i]*Vt
	  ELSE Stop_Prog('Invalid Dp_dep selected.');
    END;
END;

PROCEDURE Calc_Langevin_factor(VAR Lan : vector; mob_n, mob_p : vector);
{Calculates the Langevin recombination strength. Lan[i] is defined on the
regular grid, i.e., Lan[i]=Lan at x=xi. Note that mobilities are defined on the
interleaved mesh!}
VAR i : INTEGER;
    rec_mob : myReal;
BEGIN
    IF UseLangevin THEN {use Langevin formula to calculate bimolecular, direct, band-to-band recombination rate}
    FOR i:=1 TO NP DO
    BEGIN
		rec_mob:=(mob_n[i-1] + mob_n[i]+ mob_p[i-1] + mob_p[i] )/2;
		{we take mob(x=xi)=(mob(x=xi-1/2)+mob(x=xi+1/2))/2}
		Lan[i]:=Lang_pre * q * rec_mob/eps[i];
    END
    ELSE
        FOR i:=1 TO NP DO {use input value for direct recombination}
            Lan[i]:=kdirect;

    Lan[0]:=0; {no recombination (or generation) at the contacts}
    Lan[NP+1]:=0;
END;


FUNCTION DissProb_Delta(r : myReal) : myReal;
{Calculates the dissociation probability as a function of distance r}
VAR b, kdF, delE : myReal;
BEGIN
    delE:=q/(4*PI*epsi*r); {binding energy in eV}
    b:=q*ABS(F)/(8*PI*epsi*SQR(Vt)); {scaled ABOLUTE field strength}
    kdF:=3*Braun_rec/(4*PI*r*SQR(r))*EXP(-delE*Vti)*Bessel(b);
    {note: Braun_rec appears here!}
    DissProb_Delta:=kdF/(kdF + kf)
END;


FUNCTION DissProb_Gauss(r : myReal) : myReal;
{calculates the diss. prob. as a function of distance r with a Gaussian distribution
of a}
BEGIN
    DissProb_Gauss:=(4/(a*a*a*SQRT(PI)))*SQR(r)*EXP(-SQR(r/a))*DissProb_Delta(r)
END;


FUNCTION DissProb_Exp(r : myReal) : myReal;
{calculates the diss. prob. as a function of distance r with an exponential
distribution of a}
BEGIN
    DissProb_Exp:=EXP(-r/a)/a * DissProb_Delta(r)
END;


FUNCTION DissProb_SQRrExp(r : myReal) : myReal;
{calculates the diss. prob. as a function of distance r with an r^2 exponential
distribution of a}
BEGIN
    DissProb_SQRrExp:=SQR(r)*EXP(-r/a)/(2*a*a*a) * DissProb_Delta(r)
END;


FUNCTION DissProb_r4Gauss(r : myReal) : myReal;
{calculates the diss. prob. as a function of distance r with an r^4 * Gaussian
distribution of a}
BEGIN
    DissProb_r4Gauss:=8*SQR(r)*SQR(r)*EXP(-SQR(r/a))/(3*Power(a,5)*SQRT(PI))
    * DissProb_Delta(r)
END;


PROCEDURE Calc_Dissociation(VAR dp, g : vector; Lan, V, mob_n, mob_p : vector);
{calculates the field-dependent dissociation rate dp and the generation rate,
g[i]=g at x=i}
{see C. Braun, J. Chem. Phys. 80, p. 4157 (1984)}
VAR i : INTEGER;
BEGIN
    IF Field_dep_G
		THEN { calc. field dep. G }
		FOR i:=1 TO NP DO
		BEGIN
			F:=(V[i+1]-V[i-1])/(L*(h[i]+h[i-1]));
			{F= field on mesh point i, centered difference, global variable}
			Braun_rec:=Lan[i];
			{Braun recombination strength at x=xi, global variable!}
			epsi:=eps[i]; {Global variable to be used in DissProb_Delta}
			CASE ThermLengDist OF  {a = thermalization length}
				1 : dp[i]:=DissProb_Delta(a); {delta-function distribution}
				2 : dp[i]:=RombergIntegration(DissProb_Gauss, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
				3 : dp[i]:=RombergIntegration(DissProb_Exp, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
				4 : dp[i]:=RombergIntegration(DissProb_SQRrExp, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
				5 : dp[i]:=RombergIntegration(DissProb_r4Gauss, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE)
			END; {case selector}
			{total free-carrier yield is sum of direct generation (P0) and the field dependent part (1-P0)*(dp)}
			g[i]:=(P0 + (1-P0)*dp[i]) * Gm[i]
		END {for loop}
	ELSE FOR i:=1 TO NP DO BEGIN dp[i]:=0; g[i]:=Gm[i] END; {G is constant}
    dp[0]:=0;
    dp[NP+1]:=0;
    g[0]:=0;   {no generation at the contacts}
    g[NP+1]:=0
END;


PROCEDURE Init_nt0_and_pt0();
VAR i : INTEGER;
BEGIN
    nt0[0] := 0;
    nt0[NP+1] := 0;
    pt0[0] := 0;
    pt0[NP+1] := 0;
    
    FOR i := 1 TO NP DO
    BEGIN
		nt0[i]:=NcLoc[i]*EXP((E_CB[i]-Etrap)*Vti);
		pt0[i]:=NcLoc[i]*EXP((Etrap-E_VB[i])*Vti);
    END;
END;


PROCEDURE Recombi(VAR r : vector; V, n, p, g, Lan, dp, srhp : vector);
VAR i, j		       : INTEGER;
    rec_SRH_int		       : myreal;
    inv_denom_trap_occup, f_ti : vector;
BEGIN {Langevin and Shockley-Read-Hall recombination with generation term G}
    inv_denom_trap_occup := Calc_inv_f_ti_denom();
    f_ti := Calc_f_ti(inv_denom_trap_occup);
    
    FOR i:=0 TO NP+1 DO 
    BEGIN
		rec_SRH_int := 0;
		FOR j := MAX(0,i-1) TO MIN(NP+1,i+1) DO rec_SRH_int := rec_SRH_int + (Cn * n[j] * Nti[j] * (1-f_ti[j]) - Cn * nt0[j] * Nti[j] * f_ti[j]);

		{note that srhp, g, and, lan are empty depending on calculated recombination}
		r[i]:=( (1-dp[i])*Lan[i]+srhp[i] ) * (n[i]*p[i]-SQR(ni[i])) + rec_SRH_int / 2 - g[i]
    END;
END;


PROCEDURE Calc_Recombi_current(VAR Jbimo, JSRH_bulk, JSRH_LI, JSRH_RI, Jph, Jn_l, Jn_r, Jp_l, Jp_r : myReal);
VAR i				     : INTEGER;
    inv_denom_trap_occup, f_tb, f_ti : vector;
BEGIN {Langevin and Shockley-Read-Hall recombination with generation term G}
    inv_denom_trap_occup := Calc_inv_f_ti_denom();
    f_tb := Calc_f_tb();
    f_ti := Calc_f_ti(inv_denom_trap_occup);
    
    Jbimo:=0;
    JSRH_bulk:=0;
    JSRH_LI:=0;
    JSRH_RI:=0;
    Jph:=0;
    FOR i:=0 TO NP+1 DO Jbimo := Jbimo + q*L*h[i]*( (1-dp[i])*Lang[i]) * (n[i]*p[i]-SQR(ni[i]));

    {the recombination current for traps}
    IF NOT ((nu_pi = 0) OR (nu_ni = 0)) THEN
    BEGIN
		FOR i:=i1 TO i1+1 DO JSRH_LI := JSRH_LI + q*L*h[i]*(Cn * n[i] * Nti[i] * (1-f_ti[i]) - Cn * nt0[i] * Nti[i] * f_ti[i]);
		FOR i:=i2-1 TO i2 DO JSRH_RI := JSRH_RI + q*L*h[i]*(Cn * n[i] * Nti[i] * (1-f_ti[i]) - Cn * nt0[i] * Nti[i] * f_ti[i]);
		FOR i:=i1+2 TO i2-2 DO JSRH_bulk := JSRH_bulk + q*L*h[i]*(Cn * n[i] * Nti[i] * (1-f_ti[i]) - Cn * nt0[i] * Nti[i] * f_ti[i]);
		FOR i:=1 TO NP DO JSRH_bulk := JSRH_bulk + q*L*h[i]*(Cn * n[i] * Ntb[i] * (1-f_tb[i]) - Cn * nt0[i] * Ntb[i] * f_tb[i]);    
    END;
    
    FOR i:=0 TO NP+1 DO Jph := Jph-q*L*h[i]*gen[i];
    Jn_l := Jn[0];
    Jn_r := Jn[NP];
    Jp_l := Jp[0];
    Jp_r := Jp[NP];

END;


PROCEDURE SetDensitiesImageLowering(VAR V, n, p : vector; VAR phi_left_eff, phi_right_eff : myReal);
{calculates the carrier density if there is image-force-lowering.
See Scott & Malliaras, Chem. Phys. Lett. 299, 115, equation 11}
{we also account for the effect of roughness, see Palasantzas et al, PRB 60, 9157}
VAR F0,f, psi : myReal;
BEGIN
    F0:=(V[1]-V[0])/(L*h[0]); {electric field at cathode}
    IF F0>0 THEN {only if field is positive do we have lowering:}
    BEGIN
		f:=q*F0/(4*pi*eps[0]*SQR(Vt));
		psi:=1/f * (1-SQRT(1+2*SQRT(f))) + 1/SQRT(f);
		n[0]:=4*SQR(psi)*NcLoc[0]*EXP(SQRT(f)-phi_left*Vti);
		n[0]:=n[0]*EXP( RoughLeft * SQRT(q*F0/(4*pi*eps[0]))/(2*Vt)); {now account for influence of roughness}
		n[0]:=MIN(n[0],NcLoc[0]); {n[0] cannot exceed Nc}
		p[0]:=SQR(ni[i])/n[0];
    END;
    phi_left_eff:=-Vt*LN(n[0]/NcLoc[0]); {calculate effective barrier height}

    {now do the same at the anode:}
    F0:=(V[NP+1]-V[NP])/(L*h[NP]); {electric field at anode}
    IF F0>0 THEN {only if field is positive do we have lowering:}
    BEGIN
		f:=q*F0/(4*pi*eps[NP+1]*SQR(Vt));
		psi:=1/f * (1-SQRT(1+2*SQRT(f))) + 1/SQRT(f);
		p[NP+1]:=4*SQR(psi)*NcLoc[NP+1]*EXP(SQRT(f)-phi_right*Vti);
		p[NP+1]:=p[NP+1]*EXP( RoughRight * SQRT(q*F0/(4*pi*eps[NP+1]))*(0.5*Vti)); {now account for influence of roughness}
		p[NP+1]:=MIN(p[NP+1],NcLoc[NP+1]); {p[NP+1] cannot exceed Nc}
		n[NP+1]:=SQR(ni[i])/p[NP+1];
    END;
    phi_right_eff:=-Vt*LN(p[NP+1]/NcLoc[NP+1]); {calculate effective barrier height}
END;


PROCEDURE Contn(VAR n : vector; V, p, mu, D, g, Lan, dp, srhp : vector);
VAR i, j : INTEGER;
    m_srhi, rhs_srhi, denom_trap_curr, u_trap_curr_num, l_trap_curr_num, m_trap_curr_num, u_trap_curr, l_trap_curr, m_trap_curr, fac, sum_aj, sum_bj, sum_cj, sum_dj, ai, ai_min, ai_plus, ci, rhs_zero_order, lo_part_der, up_part_der, main_part_der, lo_srhi, up_srhi  : myReal;
    lo, m, u, rhs, f_ti, inv_denom_trap_occup : vector;
BEGIN

    {Set the boundary conditions for finite surface recombination on the left side of the device}
    IF (Sn_L < 0) OR (Jn[0]>0) THEN	
    BEGIN
		{Infinite surface recombination}
		lo[0]:=0;
		m[0]:= 1;
		u[0]:=0;
		rhs[0]:=n[0];
    END
    ELSE
    BEGIN
		{Finite surface recombination}
		lo[0]:= 0;
		m[0]:=-D[0]*B(mu[0]*(V[0]-V[1])/D[0])/(L*h[0]) - Sn_L;
		u[0]:=D[0]*B(mu[0]*(V[1]-V[0])/D[0])/(L*h[0]);
		rhs[0]:= - Sn_L * NcLoc[0] * EXP((E_CB[0] - W_L) *Vti);
    END;

    IF (Traps_int) THEN
    BEGIN
		inv_denom_trap_occup := Calc_inv_f_ti_denom();
		f_ti := Calc_f_ti(inv_denom_trap_occup);
    END;
    
    FOR i:=1 TO NP DO  { set up the matrix form of the continuity equation. }
    BEGIN
		IF (Nti[i] > 0) AND NOT ((nu_pi = 0) OR (nu_ni = 0)) THEN
		BEGIN
			sum_aj := 0;
			sum_bj := 0;
			sum_cj := 0;
			sum_dj := 0;
			FOR j:=-1 TO 1 DO
			BEGIN
				sum_aj := sum_aj + nu_ni * Nti[i+j] * n[i+j];
				sum_bj := sum_bj + nu_pi * Nti[i+j] * p[i+j];
				sum_cj := sum_cj + nu_ni * Nti[i+j] * nt0[i+j];
				sum_dj := sum_dj + nu_pi * Nti[i+j] * pt0[i+j];
			END;
			ci := nu_ni * Nti[i] * nt0[i];
			ai_min := nu_ni * Nti[i-1] ;
			ai := nu_ni * Nti[i] ;
			ai_plus := nu_ni * Nti[i+1];
			{calculate linearization of the recombination}
			{zero'th order in the taylor expansion, so all to rhs}
			rhs_zero_order := (-ci * (sum_aj + sum_dj) + ai*n[i] * (sum_bj + sum_cj)) * inv_denom_trap_occup[i];

			{ first order of taylor expansion}
			{derivative to ni-1, same as to ni+1}
			lo_part_der := ci * ai_min * (inv_denom_trap_occup[i] - ai_min * n[i-1]);
			lo_part_der := lo_part_der - ai_min * (ci * ((sum_aj-ai_min*n[i-1]) + sum_dj) - ai*n[i]* (sum_bj + sum_cj));
			lo_part_der := lo_part_der * SQR(inv_denom_trap_occup[i]);

			up_part_der := ci * ai_plus * (inv_denom_trap_occup[i] - ai_plus * n[i+1]);
			up_part_der := up_part_der - ai_plus * (ci * ((sum_aj-ai_plus*n[i+1]) + sum_dj) - ai*n[i]* (sum_bj + sum_cj));
			up_part_der := up_part_der * SQR(inv_denom_trap_occup[i]);

			{derivative to ni}
			main_part_der := (ci - sum_bj - sum_cj) * ai  * (inv_denom_trap_occup[i] - ai * n[i]);
			main_part_der := main_part_der - ai * ci * (sum_dj + sum_aj - ai * n[i]);
			main_part_der := main_part_der * SQR(inv_denom_trap_occup[i]);

			m_srhi := -main_part_der;
			lo_srhi :=-lo_part_der;
			up_srhi :=-up_part_der;
			rhs_srhi := rhs_zero_order - n[i-1] * lo_srhi - n[i+1] * up_srhi - n[i] * main_part_der;

		END
		ELSE BEGIN
			lo_srhi := 0;
			m_srhi := 0;
			up_srhi := 0;
			rhs_srhi := 0;
		END;
		{ calculate the terms governing conduction through traps. }
		IF (Nti[i] > 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
		BEGIN
			denom_trap_curr:=0;
			FOR j:=-1 TO 1 DO
				denom_trap_curr := denom_trap_curr + nu_pi * Nti[i+j] * p[i+j] + nu_ni * Nti[i+j] * nt0[i+j];
			u_trap_curr_num := SQR(nu_ni) * (1-f_ti[i])  * Nti[i+1] * Nti[i] * nt0[i];
			l_trap_curr_num := SQR(nu_ni) * (1-f_ti[i])  * Nti[i-1] * Nti[i] * nt0[i];
			m_trap_curr_num := SQR(nu_ni) * (1-f_ti[i]) * -Nti[i] * (Nti[i+1]*nt0[i+1] + Nti[i-1]*nt0[i-1]);
			m_trap_curr := m_trap_curr_num / denom_trap_curr;
			u_trap_curr := u_trap_curr_num / denom_trap_curr;
			l_trap_curr := l_trap_curr_num / denom_trap_curr;
		END
		ELSE BEGIN
			m_trap_curr := 0;
			u_trap_curr := 0;
			l_trap_curr := 0;
		END;
		{the continuity equation in matrix form}
		fac := 0.5*SQR(L)*h[i]*h[i-1]*(h[i]+h[i-1]);
		rhs[i]:=-fac* ((1-dp[i])*Lan[i]*SQR(ni[i])+g[i]+srhp[i]*SQR(ni[i])-rhs_srhi);
		lo[i]:=h[i]*D[i-1]*B(mu[i-1]*(V[i-1]-V[i])/D[i-1])
				+ fac* lo_srhi {off-diagonal srh interface term}
				+ fac*  l_trap_curr; {trap current term}
		m[i]:=-(h[i-1]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i]) +
				h[i]*D[i-1]*B(mu[i-1]*(V[i]-V[i-1])/D[i-1]))
				-fac*((1-dp[i])*Lan[i]*p[i]+ p[i]*srhp[i] - m_srhi)
				+ fac * m_trap_curr; {trap current term}
		u[i]:=h[i-1]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i])
				+fac* up_srhi {off-diagonal srh interface term}
				+  fac * u_trap_curr; {trap current term}
    END;

    {Set the boundary conditions on the right side of the device}
    IF (Sn_R < 0)  OR (Jn[NP]>0) THEN
    BEGIN
		{Infinite surface recombination}
		lo[NP+1]:= 0;
		m[NP+1]:= 1;
		u[NP+1]:= 0;
		rhs[NP+1]:= n[NP+1];
    END
    ELSE BEGIN
		{Finite surface recombination}
		lo[NP+1]:=-D[NP] * B(mu[NP]*(V[NP]-V[NP+1])/D[NP])/(L*h[NP]);
		m[NP+1]:=D[NP] * B(mu[NP]*(V[NP+1]-V[NP])/D[NP])/(L*h[NP]) - Sn_R;
		u[NP+1]:=0;
		rhs[NP+1]:= - Sn_R * NcLoc[NP+1] * EXP((E_CB[NP+1] - W_R)*Vti);
    END;

    Tridiag(n, lo, m, u, rhs, 0, NP+1); {Solve for the new electron densities}
    FOR i:=1 TO NP DO {and check whether n>=0}
	IF n[i] < 0 THEN
        BEGIN
	    IF resetNegDens THEN n[i]:=1 ELSE Stop_Prog('Negative electron concentration encountered!')
	END;
END;


PROCEDURE Contp(VAR p : vector; V, n,  mu, D, g, Lan, dp, srhp : vector);
VAR i, j : INTEGER;
    m_srhi, rhs_srhi, denom_trap_curr, u_trap_curr_num, l_trap_curr_num, m_trap_curr_num, u_trap_curr, l_trap_curr, m_trap_curr, fac, sum_aj, sum_bj, sum_cj, sum_dj, bi, bi_min, bi_plus, di, rhs_zero_order, lo_part_der, up_part_der, main_part_der, lo_srhi, up_srhi : myReal;
    lo, m, u, rhs, f_ti, inv_denom_trap_occup : vector;
BEGIN

    {Set the boundary conditions for finite surface recombination on the left side of the device}
    IF (Sp_L < 0) OR (Jp[0]>0) THEN
    BEGIN
		{Infinite surface recombination}
		lo[0]:=0;
		m[0]:=1;
		u[0]:=0;
		rhs[0]:=p[0];
    END
    ELSE BEGIN
		{Finite surface recombination}
		lo[0]:= 0;
		m[0]:=D[0] * B(mu[0]*(V[1]-V[0])/D[0])/(L*h[0]) - Sp_L;
		u[0]:=-D[0] * B(mu[0]*(V[0]-V[1])/D[0])/(L*h[0]);
		rhs[0]:=-Sp_L*NcLoc[0] * EXP((W_L - E_VB[0])*Vti);
    END;

    IF (Traps_int) THEN
    BEGIN
		inv_denom_trap_occup := Calc_inv_f_ti_denom();
		f_ti := Calc_f_ti(inv_denom_trap_occup);
    END;
    
    FOR i:=1 TO NP DO { set up the matrix form of the continuity equation. }
    BEGIN
		IF (Nti[i] > 0) AND NOT ((nu_pi = 0) OR (nu_ni = 0)) THEN
		BEGIN
			sum_aj := 0;
			sum_bj := 0;
			sum_cj := 0;
			sum_dj := 0;
			FOR j:=-1 TO 1 DO
			BEGIN
				sum_aj := sum_aj + nu_ni * Nti[i+j] * n[i+j];
				sum_bj := sum_bj + nu_pi * Nti[i+j] * p[i+j];
				sum_cj := sum_cj + nu_ni * Nti[i+j] * nt0[i+j];
				sum_dj := sum_dj + nu_pi * Nti[i+j] * pt0[i+j];
			END;
			di := nu_pi * Nti[i] * pt0[i];
			bi_min := nu_pi * Nti[i-1];
			bi := nu_pi * Nti[i];
			bi_plus := nu_pi * Nti[i+1] ;
			{calculate linearization of the recombination}
			{zero'th order in the taylor expansion, so all to rhs}
			rhs_zero_order := (-di * (sum_bj + sum_cj) + bi * p[i] * (sum_aj + sum_dj)) / (sum_aj+sum_bj+sum_cj+sum_dj);

			{ first order of taylor expansion}
			{derivative to ni-1, same as to ni+1}
			lo_part_der := di * bi_min * (inv_denom_trap_occup[i] - bi_min * p[i-1]);
			lo_part_der := lo_part_der - bi_min * (di * ((sum_bj-bi_min*p[i-1]) + sum_cj) - bi*n[i]* (sum_aj + sum_dj));
			lo_part_der := lo_part_der * SQR(inv_denom_trap_occup[i]);

			up_part_der := di * bi_plus * (inv_denom_trap_occup[i] - bi_plus * p[i+1]);
			up_part_der := up_part_der - bi_plus * (di * ((sum_bj-bi_plus*p[i+1]) + sum_cj) - bi*n[i]* (sum_aj + sum_dj));
			up_part_der := up_part_der * SQR(inv_denom_trap_occup[i]);

			{derivative to ni}
			main_part_der := (di - sum_aj - sum_dj) * bi  * (inv_denom_trap_occup[i] - bi * p[i]);
			main_part_der := main_part_der - bi * di * (sum_cj + sum_bj - bi * p[i]);
			main_part_der := main_part_der * SQR(inv_denom_trap_occup[i]);
	    
			m_srhi := -main_part_der;
			lo_srhi :=  -lo_part_der;
			up_srhi := -up_part_der;
			rhs_srhi := rhs_zero_order - p[i-1] * lo_srhi - p[i+1] * up_srhi - p[i] * main_part_der;	    

		END
		ELSE BEGIN
			lo_srhi := 0;
			m_srhi := 0;
			up_srhi := 0;
			rhs_srhi := 0;
		END;
		{ set up the terms governing conduction through traps. }
		IF (Nti[i] > 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
		BEGIN
			denom_trap_curr:=0;
			FOR j:=-1 TO 1 DO
				denom_trap_curr := denom_trap_curr + nu_ni * Nti[i+j] * n[i+j] + nu_pi * Nti[i+j] * pt0[i+j]; {changed}
			u_trap_curr_num := SQR(nu_pi) * f_ti[i] * Nti[i+1] * Nti[i] * pt0[i]; {changed} 
			l_trap_curr_num := SQR(nu_pi) * f_ti[i] * Nti[i-1] * Nti[i] * pt0[i]; {changed}
			m_trap_curr_num := SQR(nu_pi) * f_ti[i] * -Nti[i] * (Nti[i+1]*pt0[i+1] + Nti[i-1]*pt0[i-1]); {changed}
			m_trap_curr := m_trap_curr_num / denom_trap_curr;
			u_trap_curr := u_trap_curr_num / denom_trap_curr;
			l_trap_curr := l_trap_curr_num / denom_trap_curr;
		END
		ELSE BEGIN
			m_trap_curr := 0;
			u_trap_curr := 0;
			l_trap_curr := 0;	    
		END;
		{ The continuity equation in matrix form }
		fac := 0.5*SQR(L)*h[i]*h[i-1]*(h[i]+h[i-1]);
		rhs[i]:=-fac * ((1-dp[i])*Lan[i]*SQR(ni[i])+g[i]+srhp[i]*SQR(ni[i])-rhs_srhi);
		lo[i]:=h[i]*D[i-1]*B(mu[i-1]*(V[i]-V[i-1])/D[i-1])
				+ fac* lo_srhi {off-diagonal srh interface term}
				+ fac * l_trap_curr; {trap current term}
		m[i]:=-(h[i-1]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i]) +
				h[i]*D[i-1]*B(mu[i-1]*(V[i-1]-V[i])/D[i-1]))
				-fac*((1-dp[i])*Lan[i]*n[i]+ n[i]*srhp[i] - m_srhi)
				+ fac * m_trap_curr; {trap current term}
		u[i]:=h[i-1]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i])
				+ fac* up_srhi {off-diagonal srh interface term}
			+ fac* u_trap_curr; {trap current term}
    END;

    {Set the boundary conditions on the right side of the device}
    IF (Sp_R < 0) OR (Jp[NP]>0) THEN
    BEGIN
		{Infinite surface recombination}
		lo[NP+1]:=0;
		m[NP+1]:=1;
		u[NP+1]:=0;
		rhs[NP+1]:=p[NP+1];
    END
    ELSE BEGIN
		{Finite surface recombination}
		lo[NP+1]:=D[NP]*B(mu[NP]*(V[NP+1]-V[NP])/D[NP])/(L*h[NP]);
		m[NP+1]:=-D[NP]*B(mu[NP]*(V[NP]-V[NP+1])/D[NP])/(L*h[NP]) - Sp_R;
		u[NP+1]:=0;
		rhs[NP+1]:=-Sp_R * NcLoc[NP+1] * EXP((W_R - E_VB[NP+1])*Vti);
    END;

    Tridiag(p, lo, m, u, rhs, 0, NP+1); {Solve for the new hole densities}
    FOR i:=1 TO NP DO {and check whether p>=0}
	IF p[i] < 0 THEN BEGIN
	    IF resetNegDens THEN p[i]:=1 ELSE Stop_Prog('Negative hole concentration encountered!')
	END;
END;


PROCEDURE Calc_elec_curr(VAR Curr, V, n, D, mu : vector);
{This procedure calculates the electron current density from Selberherr eq. 6.1-39 modified for generalised Einstein relation}
VAR i, j												      : INTEGER;
    denom_trap_curr, u_trap_curr_num, m_trap_curr_num, m_trap_curr, u_trap_curr : myReal;
    f_ti												      : vector;
BEGIN
    f_ti:= calc_f_ti(calc_inv_f_ti_denom);
    FOR i:=0 TO NP DO
    BEGIN
		Curr[i]:=-q*(n[i+1]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i]) - n[i]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i]))/(L*h[i]);

		IF (Nti[i] <> 0) AND (Nti[i+1] <> 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
		BEGIN
			denom_trap_curr:=0;
			FOR j:=-1 TO 1 DO
				denom_trap_curr := denom_trap_curr + nu_pi * Nti[i+j] * p[i+j] + nu_ni * Nti[i+j] * nt0[i+j];
			u_trap_curr_num := n[i+1] * SQR(nu_ni) * (1-f_ti[i])  * Nti[i+1] * Nti[i] * nt0[i];
			m_trap_curr_num := n[i] * SQR(nu_ni) * (1-f_ti[i]) * - Nti[i] * (Nti[i+1]*nt0[i+1] + Nti[i-1]*nt0[i-1]);
			m_trap_curr := m_trap_curr_num / denom_trap_curr;
			u_trap_curr := u_trap_curr_num / denom_trap_curr;
			{The factor two is because in contn, the trap current is applied at both sides of the interface}
			Curr[i]:= Curr[i] - 2*q*L*h[i]*(m_trap_curr + u_trap_curr);
		END;
    END;
END;


PROCEDURE Calc_hole_curr(VAR Curr, V, p, D, mu : vector);
{This procedure calculates the electron current density from Selberherr eq. 6.1-41 modified for generalised Einstein relation}
VAR i,j : INTEGER;
    denom_trap_curr, u_trap_curr_num, m_trap_curr_num, m_trap_curr, u_trap_curr : myReal;
    f_ti												      : vector;
BEGIN
    f_ti:= calc_f_ti(calc_inv_f_ti_denom);
    FOR i:=0 TO NP DO
    BEGIN
		Curr[i]:=-q*(p[i]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i]) - p[i+1]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i]))/(L*h[i]);
		IF (Nti[i] <> 0) AND (Nti[i+1] <> 0) AND NOT ((nu_pi = 0) AND (nu_ni = 0)) THEN
		BEGIN
			denom_trap_curr:=0;
			FOR j:=-1 TO 1 DO
				denom_trap_curr := denom_trap_curr + nu_ni * Nti[i+j] * n[i+j] + nu_pi * Nti[i+j] * pt0[i+j]; {changed}

			u_trap_curr_num := p[i+1] *SQR(nu_pi) * f_ti[i] * Nti[i+1] * Nti[i] * pt0[i]; {changed} 
			m_trap_curr_num := p[i] * SQR(nu_pi) * f_ti[i] * - Nti[i] * (Nti[i+1]*pt0[i+1] + Nti[i-1]*pt0[i-1]); {changed}
			m_trap_curr := m_trap_curr_num / denom_trap_curr;
			u_trap_curr := u_trap_curr_num / denom_trap_curr;
			{The factor two is because in contp, the trap current is applied at both sides of the interface}
			Curr[i]:= Curr[i] + 2*q*L*h[i]*(m_trap_curr + u_trap_curr);
		END;
    END;
    
END;


FUNCTION Calc_trap_charge_dens() : vector;
VAR i				     : INTEGER;
    f_tb, f_ti, inv_denom_trap_occup : vector;
    C_trapped_bulk, C_trapped_int    : myReal;
BEGIN
    inv_denom_trap_occup := Calc_inv_f_ti_denom();
    f_tb := Calc_f_tb();
    f_ti := Calc_f_ti(inv_denom_trap_occup);    

    FOR i := 1 TO NP DO
    BEGIN
		{trap charge, positive if positive charge, negative if negative charge}
		IF (Ntb[i] <> 0) AND NOT (Tr_type_B = 0) THEN C_trapped_bulk := (cwe_b - f_tb[i]) * Ntb[i] ELSE C_trapped_bulk := 0;
		IF (Traps_int) THEN C_trapped_int := (cwe_i[i] - f_ti[i]) * 0.5* Nti[i] * charged_traps[i] ELSE C_trapped_int := 0;
		Calc_trap_charge_dens[i] := C_trapped_bulk + C_trapped_int;
    END;
END;


PROCEDURE Calc_Doping(VAR V, nid, pid : vector); {calc. ionised doping densities}
VAR field : myReal;
    i : INTEGER;
BEGIN
    IF (i2<NP+1) OR (i1>0) THEN {0 to i1: left transport layer, i2 to NP: right transport layer}
	FOR i:=i1+1 TO i2-1 DO
	BEGIN
	    {note: we approximate the field at i by the field between i and i+1,
	    just like we do in the mobility}
	    field:=(V[i+1]-V[i])/(L*h[i]);
	    IF beta_n = 0 THEN nid[i]:=n_0 ELSE nid[i]:=n_0 * EXP(beta_n*SQRT(ABS(field)));
	    IF beta_p = 0 THEN pid[i]:=p_0 ELSE pid[i]:=p_0 * EXP(beta_p*SQRT(ABS(field)));
	END
    ELSE
    BEGIN
        FOR i:=0 TO NP DO
        BEGIN
			field:=(V[i+1]-V[i])/(L*h[i]);
			IF beta_n = 0 THEN nid[i]:=n_0 ELSE nid[i]:=n_0 * EXP(beta_n*SQRT(ABS(field)));
			IF beta_p = 0 THEN pid[i]:=p_0 ELSE pid[i]:=p_0 * EXP(beta_p*SQRT(ABS(field)));
        END;

        {at NP+1 we use the field between NP and NP+1, which is the last one we calculated:}
        IF beta_n = 0 THEN nid[NP+1]:=n_0 ELSE nid[NP+1]:=n_0 * EXP(beta_n*SQRT(ABS(field)));
        IF beta_p = 0 THEN pid[NP+1]:=p_0 ELSE pid[NP+1]:=p_0 * EXP(beta_p*SQRT(ABS(field)));
    END;
END;


FUNCTION Average_Dissociation(dp : vector) : myReal;
{This function computes the average dissociation rate throughout the device}
VAR i : INTEGER;
    ans : myReal;
BEGIN
    IF Field_dep_G THEN
    BEGIN
		ans:=0;
		FOR i:=0 TO NP DO ans:=ans + (dp[i] + dp[i+1]) * h[i]/2;
		{Trapezoidal integeration is used}
		Average_Dissociation:=ans
    END
	ELSE Average_Dissociation:=1 {when Field_dep_G is false we take the av_diss to be 1}
END;

PROCEDURE Find_Solar_Cell_Parameters(JVChar : JVList; VAR SCPar : TSCPar);
{Finds Jsc, Voc, MPP (and thus FF) by interpolating the V-J data points}
VAR Nusable, i, k, i_start : INTEGER;
    c1, c2, c3 : myReal;
    x_dat, y_dat : Row;
BEGIN
    {first count how many JV points are usable}
    Nusable:=0;
    
    FOR i:=1 TO LENGTH(JVChar)-1 DO 
		IF JVChar[i].Use THEN INC(Nusable);
    {we start at i=1 as the 0th point (if any) corresponds to a pre-bias}
    {note that the length of JVChar includes index 0, so we have to stop at i=LENGHT()-1}

    
    {now create 2 arrays for the J and V points:}
    SetLength(x_dat, Nusable); {note: such arrays always start at index 0}
    SetLength(y_dat, Nusable);
    
    {now copy the usable points into x_dat and y_dat:}
    k:=0; {index couter for x_dat, y_dat arrays}
    FOR i:=1 TO LENGTH(JVChar)-1 DO {we don't need i=0 as it corresponds to Vpre (if any)}
		WITH JVChar[i] DO
			IF Use THEN 
			BEGIN
				x_dat[k]:=Vext;
				y_dat[k]:=Jext;
				INC(k)
			END;
    
    IF Nusable >= InterpolationOrder + 1 THEN {we need at 1 extra point to do the right order of interpolation}
    BEGIN
		{first calculate at short-circuit}
		Interpolation(x_dat, y_dat, 0, SCPar.Jsc, SCPar.ErrJsc, InterpolationOrder);
		SCPar.ErrJsc:=ABS(SCPar.ErrJsc);
		SCPar.calcSC:=SCPar.ErrJsc < threshold_err*ABS(SCPar.Jsc);
	
		{now calculate Voc by swapping x_dat and y_dat:}
		Interpolation(y_dat, x_dat, 0, SCPar.Voc, SCPar.ErrVoc, InterpolationOrder);
		SCPar.ErrVoc:=ABS(SCPar.ErrVoc);
		SCPar.calcOC:=SCPar.ErrVoc < threshold_err*ABS(SCPar.Voc);

		{Calculate the Maximum Power Point (MPP) and the fill-factor (FF)}
		{copy the power into y_dat:}
		FOR i:=0 TO Nusable-1 DO y_dat[i]:=-x_dat[i]*y_dat[i];
		i_start:=0;
		FOR i:=1 TO Nusable-3 DO
			IF y_dat[i] > y_dat[i_start] THEN i_start:=i;
		{i_start = i where power is maximum}
		SCPar.calcMPP:=FALSE;
		IF i_start >= 1 THEN {if i_start = 0 we didn't bracket the maximum power}
		BEGIN
			Fit_Parabola(x_dat[i_start-1], x_dat[i_start], x_dat[i_start+1], y_dat[i_start-1],
						y_dat[i_start], y_dat[i_start+1], c1, c2, c3);
			{Fit a parabola through the points around the max, store the coefficients
			in c1, c2, c3: Power[V]= c1*V^2 + c2*V + c3.}
			SCPar.Vmpp:=-c2/(2*c1); {This is the extremum of the parabola = V of MPP}
			SCPar.ErrVmpp:=0.25 * ABS(x_dat[i_start+1] - x_dat[i_start-1]); {error on Vmpp from voltage steps}
			SCPar.MPP:= c1*SQR(SCPar.Vmpp) + c2*SCPar.Vmpp + c3; {The MPP}
			SCPar.ErrMPP:=0.25 * (ABS(y_dat[i_start+1] - y_dat[i_start]) + ABS(y_dat[i_start]-y_dat[i_start-1]));
			SCPar.calcMPP:=TRUE;
	    
			{now move on to the FF:}
			SCPar.calcFF:=FALSE;
			IF (SCPar.Voc<>0) AND (SCPar.Jsc<>0)  AND SCPar.calcOC AND SCPar.calcSC THEN
			BEGIN
				SCPar.calcFF:=TRUE;
				SCPar.FF:=-SCPar.MPP/(SCPar.Voc*SCPar.Jsc); { the fill-factor}
				{now calc. the error in FF based on the propagation of errors:}
				SCPar.ErrFF:=SQRT( SQR(SCPar.ErrMPP/(SCPar.Jsc*SCPar.Voc)) + SQR(SCPar.ErrJsc*SCPar.FF/SCPar.Jsc) + SQR(SCPar.ErrVoc*SCPar.FF/SCPar.Voc));
			END;
		END;
    END;
END;


PROCEDURE Calc_and_Output_Solar_Cell_Parameters(JVExp, JVSim : JVList);
VAR SCParExp, SCParSim : TSCPar; {Experimental and simulated solar cell parameters}
    digits : INTEGER;
    uitv : TEXT;
    str : STRING;
BEGIN
    WRITELN;
    str:='Simulated';
    str:=AddChar(' ',str,tab1+Length(str));
    {add white space: AddCharR adds character ' ' to str until the length is tabpos:}
    str:=AddCharR(' ',str,tab2);
    IF UseExpData THEN
    BEGIN
		str:=str + 'Experimental';
		str:=AddCharR(' ',str,tab3); 
		str:=str + 'Deviation';
    END;
    WRITELN(str);
    
    {note: we will use Vext, Jext in calculating the Voc, Jsc, etc. as this is what is
    in the external circuit}
    Find_Solar_Cell_Parameters(JVSim, SCParSim);
    IF UseExpData 
	THEN Find_Solar_Cell_Parameters(JVExp, SCParExp)
    ELSE WITH SCParExp DO 
		BEGIN {if we are not using exp JV, then we need to set all SCParExp booleans = FALSE}
			calcSC:=FALSE; calcMPP:=FALSE; calcFF:=FALSE; calcOC:=FALSE;
		END;

    {now go by the solar parameters, one by one:}
    {each time we will calculate the number of digits (capped by minerr)}
    {if there is an experimental value we can use, we will calculate the deviation}

    {we will build up a string to write to screen. We use a string so that we can 
    easily control the layout of our table. A note on ffFixed: Fixed point format. 
    The value is converted to a string using fixed point notation. The result is 
    composed of all digits of the integer part of Value, preceded by a
    minus sign if Value is negative. Following the integer part is 
    DecimalSeparator and then the fractional part of Value, rounded off to Digits numbers. 
    If the number is too large then the result will be in scientific notation.}

    IF SCParSim.calcSC THEN 
	WITH SCParSim DO
	BEGIN	
	    str:='Jsc:  ';
	    digits:=CEIL(-LOG10(MAX(ErrJsc, minErr)));	
	    str:=str + FloatToStrF(Jsc, ffFixed, 0, digits);
	    str:=str + ' +/- ';
	    str:=str + FloatToStrF(ErrJsc, ffFixed, 0, digits);
	    str:=AddCharR(' ',str,tab2);
	    IF SCParExp.calcSC THEN
	    BEGIN
			digits:=CEIL(-LOG10(MAX(SCParExp.ErrJsc, minErr)));
			str:=str + FloatToStrF(SCParExp.Jsc, ffFixed, 0, digits);
			str:=str + ' +/- ';
			str:=str + FloatToStrF(SCParExp.ErrJsc, ffFixed, 0, digits);
			str:=AddCharR(' ',str,tab3);
			str:=str + FloatToStrF(Jsc-SCParExp.Jsc, ffFixed, 0, digits);
	    END;
	    WRITELN(str, ' A/m2');
	END;
    
    IF SCParSim.calcMPP THEN 
	WITH SCParSim DO
	BEGIN
	    str:='Vmpp: ';
	    digits:=CEIL(-LOG10(MAX(ErrVmpp, minErr)));
	    str:=str + FloatToStrF(Vmpp, ffFixed, 0, digits);
	    str:=str + ' +/- ';
	    str:=str + FloatToStrF(ErrVmpp, ffFixed, 0, digits);
	    str:=AddCharR(' ',str,tab2);
	    IF SCParExp.calcMPP THEN
	    BEGIN
			digits:=CEIL(-LOG10(MAX(SCParExp.ErrVmpp, minErr)));
			str:=str + FloatToStrF(SCParExp.Vmpp, ffFixed, 0, digits);
			str:=str + ' +/- ';
			str:=str + FloatToStrF(SCParExp.ErrVmpp, ffFixed, 0, digits);
			str:=AddCharR(' ',str,tab3);
			str:=str + FloatToStrF(Vmpp-SCParExp.Vmpp, ffFixed, 0, digits);
	    END;
	    WRITELN(str, ' V');

	    str:='MPP:  ';
	    digits:=CEIL(-LOG10(MAX(ErrMPP, minErr)));
	    str:=str + FloatToStrF(MPP, ffFixed, 0, digits);
	    str:=str + ' +/- ';
	    str:=str + FloatToStrF(ErrMPP, ffFixed, 0, digits);
	    str:=AddCharR(' ',str,tab2);
	    IF SCParExp.calcMPP THEN
	    BEGIN
			digits:=CEIL(-LOG10(MAX(SCParExp.ErrMPP, minErr)));
			str:=str + FloatToStrF(SCParExp.MPP, ffFixed, 0, digits);
			str:=str + ' +/- ';
			str:=str + FloatToStrF(SCParExp.ErrMPP, ffFixed, 0, digits);
			str:=AddCharR(' ',str,tab3);
			str:=str + FloatToStrF(MPP-SCParExp.MPP, ffFixed, 0, digits);
	    END;
	    WRITELN(str, ' W/m2');	
	END;

    IF SCParSim.calcOC THEN 
	WITH SCParSim DO
	BEGIN
	    str:='Voc:  ';
	    digits:=CEIL(-LOG10(MAX(ErrVoc, minErr)));
	    str:=str + FloatToStrF(Voc, ffFixed, 0, digits);
	    str:=str + ' +/- ';
	    str:=str + FloatToStrF(ErrVoc, ffFixed, 0, digits);
	    str:=AddCharR(' ',str,tab2);
	    IF SCParExp.calcOC THEN
	    BEGIN
			digits:=CEIL(-LOG10(MAX(SCParExp.ErrVoc, minErr)));
			str:=str + FloatToStrF(SCParExp.Voc, ffFixed, 0, digits);
			str:=str + ' +/- ';
			str:=str + FloatToStrF(SCParExp.ErrVoc, ffFixed, 0, digits);
			str:=AddCharR(' ',str,tab3);
			str:=str + FloatToStrF(Voc-SCParExp.Voc, ffFixed, 0, digits);
	    END;
	    WRITELN(str, ' V');
	END;
    
    IF SCParSim.calcFF THEN 
	WITH SCParSim DO
	BEGIN
	    str:='FF:   ';
	    digits:=CEIL(-LOG10(MAX(ErrFF, minErr)));
	    str:=str + FloatToStrF(FF, ffFixed, 0, digits);
	    str:=str + ' +/- ';
	    str:=str + FloatToStrF(ErrFF, ffFixed, 0, digits);
	    str:=AddCharR(' ',str,tab2);
	    IF SCParExp.calcFF THEN
	    BEGIN
			digits:=CEIL(-LOG10(MAX(SCParExp.ErrFF, minErr)));
			str:=str + FloatToStrF(SCParExp.FF, ffFixed, 0, digits);
			str:=str + ' +/- ';
			str:=str + FloatToStrF(SCParExp.ErrFF, ffFixed, 0, digits);
			str:=AddCharR(' ',str,tab3);
			str:=str + FloatToStrF(FF-SCParExp.FF, ffFixed, 0, digits);
	    END;
	    WRITELN(str);
	END;
    
    WRITELN; 

    
    {Store the key parameters in the scPars_file}
    WITH SCParSim DO
	IF (Voc<>0) AND (Jsc<>0) AND calcOC AND calcSC THEN
	BEGIN
	    ASSIGN(uitv, scPars_file);  {create the scPars_file}
	    REWRITE(uitv);
	    WRITELN(uitv,'Jsc ErrJsc Voc ErrVoc Vmpp ErrVmpp MPP ErrMPP FF ErrFF');
	    WRITELN(uitv,Jsc:6:4,' ',ErrJsc:6:4,' ',Voc:6:4,' ',ErrVoc:6:4,' ',Vmpp:6:4,' ',ErrVmpp:6:4,' ',MPP:6:4,' ',ErrMPP:6:4,' ',FF:6:4,' ',ErrFF:6:4);
	    CLOSE(uitv);
	END
    ELSE WRITELN('Could not determine (all) solar cell parameters.');
END;


PROCEDURE Compare_Exp_Sim_JV(JVExp, JVSim : JVList);
VAR rms, Jmin, Jmax : myReal;
    i, j, count : INTEGER;
    disgardedPoints, bracketed : BOOLEAN;
BEGIN
    {If there is series resistance we need to interpolate the experimental data to match the voltages of the simulation.}
    IF (Rseries > 0) THEN
    BEGIN
		j:=1;
		FOR i:=1 TO NJV DO
		BEGIN
			{Some simuation voltages are outside the range of the experiment data, so we can not use these}
			IF (JVSim[i].Vext < Vmin) OR (JVSim[i].Vext > Vmax) THEN JVExp[i].Use:= FALSE
			ELSE
			BEGIN
				{Set the voltage and interpolate the corresponding experimental current density.}
				JVExp[i].Vext:=JVSim[i].Vext;
				bracketed:=(JVExp[j].Vint < JVExp[i].Vext) AND (JVExp[j+1].Vint > JVExp[i].Vext);
				WHILE (NOT bracketed) AND (j < NJV-1) DO
				BEGIN
					INC(j);
					bracketed:=(JVExp[j].Vint < JVExp[i].Vext) AND (JVExp[j+1].Vint > JVExp[i].Vext);
				END;
		
				IF bracketed THEN 
					JVExp[i].Jext:=JVExp[j].Jint + (JVExp[j+1].Jint - JVExp[j].Jint) * (JVExp[i].Vext - JVExp[j].Vint) / (JVExp[j+1].Vint - JVExp[j].Vint)
				ELSE
					JVExp[i].Use:= FALSE
			END;
		END;
    END;
    
    rms:=0;
    count:=0;
    Jmin:=MaxSingle;
    Jmax:=-MaxSingle;
    disgardedPoints:=FALSE; {keep track of points that we don't use in log-fitting}
    
    FOR i:=1 TO NJV DO
		{If the simulation did not converge at we can not use it. Similarly if we could not interpolate a point in 
		the experimental data we can not use this voltage}
		IF JVSim[i].Use AND JVExp[i].Use
		THEN 
		BEGIN
			{Calculate the sum of squared residuals}
			CASE rms_mode OF
				linear : BEGIN
							rms:=rms + SQR(JVExp[i].Jext-JVSim[i].Jext);
							INC(count);
						END; {case linear}
				logarithmic : BEGIN
								IF JVExp[i].Jext*JVSim[i].Jext>=0 {we can only calc rms if both are <>0 and they have the same sign}
								THEN BEGIN
									INC(count);
									rms:=rms + SQR(LN(JVExp[i].Jext/JVSim[i].Jext))
								END
								ELSE disgardedPoints:=TRUE
							END; {case logarithmic}
			END; {case}
			{Look for the interval [Jmin,Jmax] in both the simulated and experiment data}
			IF JVExp[i].Jext < Jmin THEN Jmin:=JVExp[i].Jext;
			IF JVSim[i].Jext < Jmin THEN Jmin:=JVSim[i].Jext;
			IF JVExp[i].Jext > Jmax THEN Jmax:=JVExp[i].Jext;
			IF JVSim[i].Jext > Jmax THEN Jmax:=JVSim[i].Jext;
		END;
    

    IF (rms_threshold <= count/NJV) AND (Jmax-Jmin>tolReal) THEN {we need at least 1 data point and a measurable interval [Jmin, Jmax]}
    BEGIN
		{Here we calculate the root mean square error and normalise with respect to the interval [Jmin,Jmax].}
		CASE rms_mode OF
			linear : rms:=SQRT(rms/count)/(Jmax-Jmin);
			logarithmic : rms:=SQRT(rms/count)/ABS(LN(ABS(Jmax/Jmin))); {note: Jmax > Jmin, of course, but ABS(Jmin) can be larger than ABS(Jmax) so we need to use ABS(LN(...)) to get a positive rms}
		END;
		WRITELN('Comparing simulated and experimental JV curves:');
		WRITELN('norm_rms_error: ',rms:6:5);
		IF disgardedPoints THEN WarnUser('Not all JV points were used in computing the rms-error.');
		WRITELN;
    END
    ELSE
    BEGIN
        WRITELN('Could not compute a meaningful rms-error.');
        WRITELN('Possible reasons:');
        WRITELN('-not enough simulated points converged, check rms_threshold');
        WRITELN('-the range in experimental currents is not large enough.');
        IF rms_mode = logarithmic
	    THEN WRITELN('-too many voltages were sim. and exp. currents have different signs. Try using rms_mode = lin.');
        WRITELN
    END
END;


PROCEDURE Write_Vars_To_File(VAR V, n, p, Vgn, Vgp, rec, Lang, diss_prob, SRHB, nid, pid, mun, mup, Gm, Jn, Jp, 
				 nion, pion : vector; Va : myReal);
{writes all relevant parameters to file. } 
VAR i			     : INTEGER;
    Evac, Ec, Ev, phin, phip : myReal;
    charge_traps : vector;
    uitv : TEXT;
BEGIN
    {Create file:}
    ASSIGN(uitv, Var_file);
    REWRITE(uitv);
    FOR i:= 0 TO NP + 1 DO charge_traps[i] := 0;
    charge_traps := Calc_trap_charge_dens();
    Recombi(rec,V, n, p, empty, Lang, diss_prob, SRHB); {calc. net recombination}

    {note: we limit the length to 25 (nd) characters, otherwise, if myReal=EXTENDED, Origin cannot read the file!} 
    WRITELN(uitv, 'x V n p Evac Ec Ev phin phip ctraps nid pid nion pion mun mup rec dp Gm Jn Jp');
    FOR i:=0 TO NP+1 DO
    BEGIN
		Evac:=V[0] - V[i]; {the vacuum level. We take it zero at x=0}
		{use the generalised potentials to take care of the band-offsets}
		{but we need to correct for the effect of the DOS on Vgn/p:}
		Ec:=V[0] - Vgn[i] - CB + Vt*LN(NcLoc[i]/Nc);
		Ev:=V[0] - Vgp[i] - VB - Vt*LN(NcLoc[i]/Nc);
		{electron and hole quasi-Fermi levels:}
		phin:=Ec + Vt*LN(n[i]/NcLoc[i]);
		phip:=Ev - Vt*LN(p[i]/NcLoc[i]);
		WRITELN(uitv,x[i]:nd,' ',V[i]:nd,' ',n[i]:nd,' ',p[i]:nd,' ',Evac:nd,' ',Ec:nd,' ',Ev:nd,' ',
			phin:nd,' ',phip:nd,' ',charge_traps[i]:nd,' ',nid[i]:nd,' ',pid[i]:nd,' ',
			nion[i]:nd,' ', pion[i]:nd,' ',mun[i]:nd,' ',mup[i]:nd,' ',rec[i]:nd,' ',
			diss_prob[i]:nd,' ',Gm[i]:nd,' ',Jn[i]:nd,' ',Jp[i]:nd);
    END;
    CLOSE(uitv);

    WRITELN('The values of (V, n, p) for the last voltage are written in ', Var_file,'.');
END;


BEGIN {main program}
    WRITELN('Welcome to SIMsalabim version ',version,'.');
    WRITELN('Copyright (C) 2020 Dr T.S. Sherkar, V.M. Le Corre, M. Koopmans');
    WRITELN('F. Wobben, and Prof L.J.A. Koster, University of Groningen.');
    WRITELN;

    {if '-h' or '-H' option is given then display some help and exit:}
    IF hasCLoption('-h') THEN DisplayHelpExit;
    IF hasCLoption('-tidy') THEN Tidy_up_parameter_file(TRUE); {tidy up file and exit}
    IF NOT Correct_version_parameter_file THEN Stop_Prog('Version of SIMsalabim and '+parameter_file+' do not match.');

    {Initialisation:}
    Read_Parameters(MsgStr); {Read parameters from input file}
    Prepare_Log_File(log, MsgStr); {open log file}
    IF AutoTidy THEN Tidy_up_parameter_file(FALSE); {clean up file but don't exit!}
    Egap:=VB-CB; {define band gap}
    Vt:=k*T/q;  {thermal voltage}
    Vti:=1/Vt; {inverse of Vt, used a lot}
    Check_Parameters; {perform a number of chekcs on the paramters. Note: we need Vt and ni}
    Make_Grid(h, x); {Initialize the grid}
    DefineLayers(NcLoc, ni, nid, pid, eps, i1, i2); {define layers}
    Init_nt0_and_pt0();
    Trap_Distribution();
    Init_Voltages_and_Tasks(JVSim, JVExp, log);

    IF PreCond
	THEN VCount:=0 {Counter of the number of voltages which have been computed}
    ELSE VCount:=1; {index VCount=0 is reserved for the pre-bias voltage (if any)}

    Va:=JVSim[VCount].Vint; {get applied voltage from the list}
    Vaold:=Va; {Vaold: old value of Va, in 1st loop equal to Va}

    Init_Generation_Profile(Use_gen_profile, Gen_profile, Gm); {init. the Gm array}

    IF mob_n_dep = 2 THEN Init_Elec_Mob_Table(n_vals, Fn_vals, n_points, Fn_points, n_mob_tab); {init. the table with values of mob_n vs. F and n}
    IF mob_p_dep = 2 THEN Init_Hole_Mob_Table(p_vals, Fp_vals, p_points, Fp_points, p_mob_tab); {init. the table with values of mob_p vs. F and p}

    Init_Pot_Dens(V, Vgn, Vgp, n, p, empty, Va); {init. (generalised) potentials and densities}

    ASSIGN(uitv, JV_file);  {create the JV-file}
    REWRITE(uitv);
    WRITELN(uitv,'Vext  Vint  Jext  Jint  P  recLan  recSRH  Jbimo  JSRH_bulk  JSRH_LI  JSRH_RI  Jph  Jn_l  Jp_l  Jn_r  Jp_r');
    WRITELN('The calculation has started, please wait.');
    
    quit_Voc:=FALSE;

    WHILE (Vcount<=NJV) AND (NOT quit_Voc) DO  {loop over voltages}
    BEGIN
		Vaold:=Va; {keep old applied voltage}
		Va:=JVSim[VCount].Vint; {get new applied voltage from the list}
		FOR i:=0 TO NP+1 DO {guess new V, n, p in interior of device}
			V[i]:=V[i] + (Va-Vaold) * (x[i]-L/2)/L; {new V guessed from old V and linear term proportional to Va-Vaold}
	
		Solve_Poisson(V, n, p, nid, pid, JVSim[VCount].UpdateIons, check_Poisson);
		Calc_elec_mob(mun, V, n); {calc. mobilities with the new V}
		Calc_hole_mob(mup, V, p);
		Calc_Dn(Dn, mun, n); {calc. the diffusivities}
		Calc_Dp(Dp, mup, p);
		Calc_Langevin_factor(Lang, mun, mup);{Calc. the Langevin recombination strength}
		Calc_Dissociation(diss_prob, gen, Lang, V, mun, mup); {calculate the field-dependent generation}
		Calc_SRH_bulk_prefactor(SRHB, V, n, p); {calc. the Shockley-Read-Hall recombi. prefactor}
		Calc_Doping(V, nid, pid); {calc. ionised doping densities}

		MainIt:=0;
		delDens:=TolMain + 1;  {to ensure that delDens>TolMain in 1st loop}
		Jtot:=0;
		totn:=0.5;
		totp:=0.5; {this will ensure delDens/(totn + totp)>TolMain in 1st loop}
		old_totn:= 0;
		old_totp:= 0;
		old_del_totn:= 0;
		old_del_totp:= 0;
		FOR i:=0 TO NP+1 DO
		BEGIN
			totn:=totn + n[i];
			totp:=totp + p[i];
		END;	
		max_del:= MAX(totn,totp);
		max_found_deln := 0;
		max_found_delp := 0;
		flip_count_deln := 0;
		flip_count_delp := 0;
	
		Conv_Main:=FALSE;
	
		WHILE (NOT Conv_Main) AND (MainIt<MaxItMain) DO    {main loop}
		BEGIN
			oldn:=n; {keep the old densities to monitor the progress}
			oldp:=p;
			totn:=0; {total sum of densities throughout device}
			totp:=0;

			IF ImageLowering THEN
			BEGIN
				SetDensitiesImageLowering(V,n,p, phi_left_eff, phi_right_eff);
				V0:=0.5 * Egap - phi_left_eff; {modify V0 to account for change in barrier height}
				VL:=-0.5 * Egap + phi_right_eff; {do the same for VL}
			END;

			UpdateGenPot(V, NcLoc, Vgn, Vgp); {update generalised potentials}
			Contn(n, Vgn, oldp, mun, Dn, gen, Lang, diss_prob, SRHB); {calc. new elec. density}

			n_contn := n;
			{now use SOR (to speed up) or SUR (to slow down) changes in n and p:}
			FOR i:=0 TO NP+1 DO
				n[i]:=accDens * n[i] + (1-accDens) * oldn[i];
				
			FOR i:=0 TO NP+1 DO
				totn:=totn + n[i];

			del_totn := totn - old_totn;
	    
			IF (ABS(del_totn) > max_del) THEN
			BEGIN
				red_fac := max_del / ABS(del_totn);
				{ writeln('max del (n) = ',max_del, ' red rac = ',red_fac); }
				totn:= 0;
				FOR i:=0 TO NP+1 DO
				BEGIN
					n[i]:= oldn[i] + (n[i]-oldn[i]) * red_fac;
					totn:= totn + n[i];
				END;
				del_totn := del_totn * red_fac
			END;
	    
			Contp(p, Vgp, n, mup, Dp, gen, Lang, diss_prob, SRHB); {calc. new hole density}
	
			p_contp := p;
			
			{now use SOR (to speed up) or SUR (to slow down) changes in n and p:}
			FOR i:=0 TO NP+1 DO
				p[i]:=accDens * p[i] + (1-accDens) * oldp[i];
			{set a maximum if we are overshooting}
			FOR i:=0 TO NP+1 DO
				totp:=totp + p[i];
			del_totp := totp - old_totp;

			IF (ABS(del_totp) > max_del) THEN
			BEGIN
				red_fac := max_del / ABS(del_totp);
				{ writeln('max del (p) = ',max_del, ' red fac = ', red_fac ); }
				test_p := 0;
				totp:=0;
				FOR i:=0 TO NP+1 DO
				BEGIN
					p[i]:= oldp[i] + (p[i]-oldp[i]) * red_fac;
					test_p := test_p + (p[i]-oldp[i]) * red_fac;
					totp:= totp + p[i];
				END;
				{		writeln('test_p = ',test_p);}
				del_totp := del_totp * red_fac
			END;

			max_found_deln := MAX(ABS(max_found_deln), ABS(del_totn));    
			IF (((old_del_totn > 0) AND (del_totn < 0)) OR ((old_del_totn < 0) AND (del_totn > 0))) THEN flip_count_deln := flip_count_deln +1;
			max_found_delp := MAX(ABS(max_found_delp), ABS(del_totp));
			IF (((old_del_totp > 0) AND (del_totp < 0)) OR ((old_del_totp < 0) AND (del_totp > 0))) THEN flip_count_delp := flip_count_delp +1;
	    
			IF ((flip_count_deln > 5) AND (flip_count_delp > 5)) THEN
			BEGIN
				{ writeln('max del halved, new max_del: ',max_del:9); }
				max_del := 0.5 * MAX( max_found_deln, max_found_delp);
				max_found_deln := 0;
				max_found_delp := 0;
				flip_count_deln := 0;
				flip_count_delp := 0;
			END;
	    
			{ writeln('it #',MainIt,' new deln =  ', del_totn:9, ' delp = ',del_totp:9,' ntot=',totn:13,' ptot=',totp:13 ); }

			old_del_totn := del_totn;
			old_del_totp := del_totp;
			old_totn := totn;
			old_totp := totp;
	    
			Solve_Poisson(V, n, p, nid, pid, JVSim[Vcount].UpdateIons, check_Poisson);

			{Update ions if needed:}
			IF JVSim[Vcount].UpdateIons THEN Calc_Ion_Distribution(nion, pion, V); 

			Calc_elec_mob(mun, V, n); {calc. the new mobilities}
			Calc_hole_mob(mup, V, p);
			Calc_Dn(Dn, mun, n); {update the diffusivities}
			Calc_Dp(Dp, mup, p);
			Calc_Langevin_factor(Lang, mun, mup); {Calc. the Langevin recombination strength}
			Calc_Dissociation(diss_prob, gen, Lang, V, mun, mup); {update generation rate}
			Calc_SRH_bulk_prefactor(SRHB, V, n, p); {calc. the Shockley-Read-Hall recombi. prefactor}
			UpdateGenPot(V, NcLoc, Vgn, Vgp); {update generalised potentials}
			Calc_elec_curr(Jn, Vgn, n, Dn, mun); {calc. the current densities}
			Calc_hole_curr(Jp, Vgp, p, Dp, mup);
			Calc_Doping(V, nid, pid); {calc. ionised doping densities}

			{Check for convergence of main loop:}
			delDens:=Norm(Difference(oldn,n_contn, 0, NP+1), 0, NP+1) + Norm(Difference(oldp,p_contp, 0, NP+1), 0, NP+1);
			{delDens = total change in the loop}
			oldJtot:=Jtot; {store the old value of the current}
			Jtot:=(Jn[ROUND(NP/2)]+Jp[ROUND(NP/2)]); {and compute the new value}
			IF (Conv_Var = 1) AND (ABS(Jtot) > 1e-4) {check for convergence of current density}
				THEN Conv_Main:=ABS( (Jtot-oldJtot)/Jtot ) < TolMain
				ELSE Conv_Main:=delDens/(totn + totp) < TolMain;
	    
			MainIt:=MainIt+1
		END; {main loop}
	
		IF JVSim[VCount].Store THEN 
			WITH JVSim[VCount] DO
			BEGIN
				Jint:=Jtot; {the internal current}
				{now correct for Rshunt and Rseries to give us Vext and Jext:}
				IF Rshunt>0 THEN Jext:=Jtot + Va/Rshunt ELSE Jext:=Jtot; {add current through shunt to measured current (=J)}
				Vext:=Va + Jint * Rseries; {add voltage drop across series resistance, use measured current J, not Jtot}
				Use:=check_Poisson AND Conv_Main; {store whether convergerence was achieved}
				quit_Voc:=(Jext>0) AND until_Voc; {use the externally measured JV-curve to see if we're past Voc}
			END;
    
		IF check_Poisson AND Conv_Main THEN
		BEGIN   {Convergence!}
			Recombi(rec,V, n, p, empty, Lang, diss_prob, empty); {calc. rec. without SRH and generation}
			recLan:=Average(rec, h, 0, NP+1); {store average Langevin recombination}
			Recombi(rec,V, n, p, empty, empty, diss_prob, SRHB); {calc. rec. without Langevin and generation}
			recSRH:=Average(rec, h, 0, NP+1); {store average SRH recombination}
			Av_Diss:=Average_Dissociation(diss_prob);
			Jbimo:=0; JSRH_bulk:=0; JSRH_LI:=0; JSRH_RI:=0; Jph:=0;
			Calc_Recombi_current(Jbimo, JSRH_bulk, JSRH_LI, JSRH_RI, Jph, Jn_l, Jn_r, Jp_l, Jp_r); {calc. net recombination}
			WRITELN('At Vint=', Va:4:3,' converged in',MainIt:4,' loop(s), Jint=',Jtot:7:4);
			IF JVSim[VCount].Store THEN {now only store the point in the JV_file if needed}
				WITH JVSim[VCount] DO {note: we take Vint,Vext,Jint,Jext from JVSim[VCount]}
					WRITELN(uitv,Vext:nd,' ', Vint:nd, ' ',Jext:nd,' ',
							Jint:nd,' ',Av_Diss:nd,' ',recLan:nd,' ',recSRH:nd,' ',
							Jbimo:nd,' ',JSRH_bulk:nd,' ',JSRH_LI:nd,' ',JSRH_RI:nd,' ',
							Jph:nd,' ',Jn_l:nd,' ',Jp_l:nd,' ',Jn_r:nd,' ',Jp_r:nd);
			FLUSH(uitv);
		END
		{no convergence:}
		ELSE
		BEGIN
			WRITELN('No convergence for V=',Va:5:3,' Rel. change densities:',
				delDens/(totn + totp):4,' Conv. Poisson: ', Check_Poisson);
			IF (Exit_after_fail = 1) THEN Stop_prog('not converged, quitting');
		END;
        VCount:=VCount + 1;
    END; {loop over voltages}

    IF (Gmax * Gfrac <> 0) AND (V0 <> VL) THEN  {we do have a solar cell}
        Calc_and_Output_Solar_Cell_Parameters(JVExp, JVSim);

    IF UseExpData THEN Compare_Exp_Sim_JV(JVExp, JVSim);
    {note: if Rseries <> 0 then Vext in JVExp and JVSim will be different so we can't do a direct comparison}

    WRITELN('The JV characteristic is written in ', JV_file,'.');
    {Store variables (x, V, etc) for last applied voltage:}
    Write_Vars_To_File(V, n, p, Vgn, Vgp, rec, Lang, diss_prob, SRHB, nid, pid, mun, mup, Gm, Jn, Jp, nion, pion, Va);
    WRITELN('Finished, press enter to exit');

    IF Pause_at_end = 1 THEN READLN; {pause at the end of the program}

    CLOSE(log);
END.
