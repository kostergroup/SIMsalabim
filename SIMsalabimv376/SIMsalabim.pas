PROGRAM SIMsalabim;

{
SIMsalabim: a 1D drift-diffusion simulator 
Copyright (c) 2020 Dr T.S. Sherkar, V.M. Le Corre, M. Koopmans,
F.O.B. Wobben, and Prof. Dr. L.J.A. Koster, University of Groningen
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
email: l.j.a.koster@rug.nl
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
    version = '3.76';   {version, 1.00 = 10-03-2004}
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


VAR i, VCount, MainIt, gridpoint, i1, i2, i11, i22 : INTEGER;

    V, Vgn, Vgp, n, p, nid, pid, empty, n_trap, p_trap, oldn, oldp, Jn, Jp, rec, h, x,
    mun, mup, Dn, Dp, Gm, gen, Lang, SRH, diss_prob, eps, Nt_trap, Pt_trap, nion, pion, Q_trap : vector;

    {for tabulated electron mobility:}
    n_vals, Fn_vals, p_vals, Fp_vals : row; {values of n and F for tabulated electron mobility}
    n_mob_tab, p_mob_tab : Table; {the table of mobilities}
    n_points, p_points, Fn_points, Fp_points : INTEGER; {number of points in F and n direction}

    check_Poisson, Conv_Main, extra_F, extra_ln_mob, Field_dep_G, Use_gen_profile,
    UseLangevin, ImageLowering, Traps, until_Voc, quit_Voc,
    TLsAbsorb, TLsTrap, PreCond, resetNegDens, AutoTidy, UseExpData : BOOLEAN;

    Egap, ni, Jtot, oldJtot, totn, totp, Vt, F, Braun_rec, V0, VL, Av_Diss,
    recLan, recSRH, phi_left, phi_right, phi_left_eff, phi_right_eff, epsi : myReal;

    T, L, eps_r, CB, VB, Nc, n_0, p_0, mun_0, mup_0, beta_n, beta_p,
    gamma_n, gamma_p, W_L, W_R, L_LTL, L_RTL, doping_LTL , doping_RTL,
    mob_LTL , mob_RTL, eps_r_LTL, eps_r_RTL, CB_LTL, VB_LTL, CB_RTL, VB_RTL,
    RoughLeft, RoughRight, SnMaj, SpMaj, SnMin, SpMin, Rshunt, Rseries,
    eta_n, eta_p, Gmax, Gfrac, P0, a, kf, kdirect, Lang_pre, 
    Nt, Pt, Bulk_tr, LTL, LTR, St_L, St_R, Etrap, Ttr, GB_tr, L_GB, Cn, Cp,
    tolPois, maxDelV, accPois, accDens,tolMain, grad, TolRomb, Vpre, 
    Vmin, Vmax, Vstep, Vacc, LowerLimBraun, UpperLimBraun,
    CIM, rms_threshold : myReal; {paramaters from input file}

    Jbimo, JSRH_bulk, JSRH_LI, JSRH_RI, Jph, Jn_l, Jn_r, Jp_l, Jp_r : myReal; {recombination current}

    MaxItPois, MaxItMain, MaxRombIt, mob_n_dep, mob_p_dep, Dn_dep, Dp_dep,
    ThermLengDist, Trtype, Q_charge, Conv_Var, Vdistribution,
    NJV, mob_ion_spec, ion_red_rate, num_GBs, NP, Vscan,
    Pause_at_end : INTEGER; { also form input file}

	rms_mode : Tfitmode; {parameter from input file}

    inv, uitv , log : TEXT; {the input and output files}

    delDens, Va, Vaold : myReal;

    Gen_profile, ExpJV, JV_file, Var_file, n_file, p_file, scPars_file : STRING;

    JVSim, JVExp : JVList; {stores the current-voltage characteristics}
    

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

PROCEDURE Prepare_Log_File(VAR log : TEXT);
VAR PIDStr : STRING;
	dum : myReal;
   	gotit : BOOLEAN;
BEGIN
   	Str(GetProcessID, PIDStr); {put process ID into PIDStr}
   	{if user specifies -PIDlog 1 in command line then we add the PID number to the log file name to make it unique}
   	getRealfromCL('-PIDlog', gotit, dum);
   	IF gotit AND (ROUND(dum) = 1)
		THEN ASSIGN(log, 'log' + '_' + PIDStr + '.txt') {create log file with PID number}
   		ELSE ASSIGN(log, 'log.txt'); {create log file without PID number}

    REWRITE(log);
    WRITELN(log,'Version ', version);
    WRITELN(log,'Size of reals used in simulation: ',SizeOf(myReal),' bytes');
    IF ParamCount > 0 THEN
   		WRITELN(log, 'Values from command line:');
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
		IF (LeftStr(line, 1)='*') OR (LeftStr(line, 1)='') THEN
			outline:=line;
		
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

PROCEDURE Read_Parameters(VAR log : TEXT);
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
    Get_Float(inv, log,'T',T);  {abs. temperature, K}
    Get_Float(inv, log, 'L', L); {Read_Number(inv, L);  {device thickness, m}
    Get_Float(inv, log, 'eps_r', eps_r);  {relative dielectric constant}
    Get_Float(inv, log, 'CB', CB);  {eV, conduction band edge}
    Get_Float(inv, log, 'VB', VB);  {eV, valence band edge}
    Get_Float(inv, log, 'Nc',Nc);  {effective DOS, m^-3}
    Get_Float(inv, log, 'n_0', n_0);  {ionised n-doping density, m^-3}
    Get_Float(inv, log, 'p_0', p_0);  {ionised p-doping density, m^-3}
    Get_Float(inv, log, 'beta_n', beta_n); {(m/V)^0.5, field  activation factor of n_0}
    Get_Float(inv, log, 'beta_p', beta_p); {(m/V)^0.5, field  activation factor of p_0}

{**Mobilities************************************************************************}
    Get_Float(inv, log, 'mun_0', mun_0); {electron zero-field mobility, m^2/Vs}
    Get_Float(inv, log, 'mup_0', mup_0); {hole zere-field mobility, m^2/Vs}
    Get_Integer(inv, log, 'mob_n_dep', mob_n_dep);{dependence of elec mobility, 0 : const. mob, 1 : field-dep., 2 : get mobility from table specified in n_file}
    Get_Integer(inv, log, 'mob_p_dep', mob_p_dep);  {dependence of hole mobility, 0 : const. mob, 1 : field-dep., 2 : get mobility from table specified in p_file}
    Get_Float(inv, log,'gamma_n', gamma_n); {field depedence of mobility, eV/(V/m)^0.5}
    Get_Float(inv, log, 'gamma_p', gamma_p); {Poole-Frenkel form}
    Get_String(inv, log, 'n_file', n_file); {read name of file with tabulated electron mobility}
    Get_String(inv, log, 'p_file', p_file); {read name of file with tabulated hole mobility}
    Get_Integer(inv, log, 'extra_F', dumint); {whether (1) or not (0) to extrapolate F beyond its max. value in n_file or p_file}
    extra_F:=ROUND(dumint) = 1;
    Get_Integer(inv, log, 'extra_ln_mob', dumint); {whether (1) or not (0) to inter-/extrapolate ln(mob) instead of mob in n_file or p_file}
    extra_ln_mob:=ROUND(dumint) = 1;
    Get_Integer(inv, log, 'Dn_dep', dumint); {read Dn_dep}
    Dn_dep:=ROUND(dumint); {0 : Einstein, 1 : value (eta_n) times Einstein}
    Get_Integer(inv, log, 'Dp_dep', dumint); {read Dp_dep}
    Dp_dep:=ROUND(dumint); {0 : Einstein, 1 : value (eta_p) times Einstein}
    Get_Float(inv, log, 'eta_n', eta_n); {multiplication factor Dn; Dn = eta_n Vt mun}
    Get_Float(inv, log, 'eta_p', eta_p); {multiplication factor Dp; Dp = eta_p Vt mup}

{**Contacts**************************************************************************}
	Get_Float(inv, log, 'W_L', W_L); {eV, work function left electrode (= cathode)}
	Get_Float(inv, log, 'W_R', W_R); {eV, work function right electrode (= anode)}
	
	phi_left:=W_L-CB;
	phi_right:=VB-W_R;

    V0:=0.5 * (VB+CB) - W_L;
    VL:=0.5 * (VB+CB) - W_R;
 
	Get_Integer(inv, log, 'ImageLowering', dumint);
	ImageLowering:= (dumint = 1); {whether (1) or not (<>1) to use image-force lowering of barriers.
		    Note: only works for electrons at left and holes at right electrode
		    and only if surface recombination is infinite. }
	Get_Float(inv, log, 'RoughLeft', RoughLeft); {extra barrier lowering due to roughness of cathode}
	Get_Float(inv, log, 'RoughRight', RoughRight); {extra barrier lowering due to roughness of cathode}
	Get_Float(inv, log, 'SnMin', SnMin); {m/s, surface recombination of electrons at anode}
	Get_Float(inv, log, 'SpMin', SpMin); {m/s, surface recombination of holes at cathode}
	Get_Float(inv, log, 'SnMaj', SnMaj); {m/s, surface recombination of electrons at anode, majority carriers}
	Get_Float(inv, log, 'SpMaj', SpMaj); {m/s, surface recombination of holes at cathode, majority carriers}
	Get_Float(inv, log, 'Rshunt', Rshunt); {Ohms m2, shunt resistance. Use negative value for infinite Rshunt}
	Get_Float(inv, log, 'Rseries', Rseries); {Ohms m2, series resistance}

{**Transport Layers******************************************************************}
{**syntax: [name variable]_LTL/RTL for left, resp., right transport layer}
	Get_Float(inv, log, 'L_LTL', L_LTL); {m, thickness left TL}
	Get_Float(inv, log, 'L_RTL', L_RTL); {m, thickness right TL}
	Get_Float(inv, log, 'doping_LTL', doping_LTL);  {m^-3, doping in left TL if >0 p-type doping if <0 n-type doping}
	Get_Float(inv, log, 'doping_RTL', doping_RTL);  {m^-3, doping in right TL if >0 p-type doping if <0 n-type doping}
	Get_Float(inv, log, 'mob_LTL', mob_LTL); {m2/Vs, mobility of left TL}
	Get_Float(inv, log, 'mob_RTL', mob_RTL); {m2/Vs, mobility of right TL}
	Get_Float(inv, log, 'eps_r_LTL', eps_r_LTL); {relative dielectric constant left TL}
	Get_Float(inv, log, 'eps_r_RTL', eps_r_RTL); {relative dielectric constant right TL}
	Get_Float(inv, log, 'CB_LTL', CB_LTL); {eV, conduction band left TL}
	Get_Float(inv, log, 'CB_RTL', CB_RTL); {eV, conduction band right TL}
	Get_Float(inv, log, 'VB_LTL', VB_LTL); {eV, valence left TL}
	Get_Float(inv, log, 'VB_RTL', VB_RTL); {eV, valence right TL}
	Get_Integer(inv, log, 'TLsAbsorb', dumint); {TLsAbsorb, TLs absorb yes(1)/no(0), overrides the profile}
	TLsAbsorb:=ROUND(dumint)=1;
	Get_Integer(inv, log, 'TLsTrap', dumint); {TLsTrap, traps in TLs yes(1)/no(0)}
	TLsTrap:=ROUND(dumint)=1;

{**Ions*******************************************************************}
	Get_Float(inv, log, 'CIM', CIM); {m^-3, concentration of ions}
	Get_Float(inv, log, 'mob_ion_spec', dummy);{mobile ion species: -1: negative, 0: both, 1: positive ions}
	mob_ion_spec := ROUND(dummy);
	Get_Integer(inv, log, 'ion_red_rate', ion_red_rate);{number of voltage steps after which ions redistribute, }
	
{**Generation and recombination******************************************************}
    Get_Float(inv, log, 'Gmax', Gmax);  {maximum generation rate, m^-3/s}
    Get_Float(inv, log, 'Gfrac', Gfrac);
    Gmax:=Gmax*Gfrac;
    Get_String(inv, log, 'Gen_profile', Gen_profile); {name of file generation profile (or 'none')}
    Use_gen_profile:= lowercase(Trim(Gen_profile))<>'none'; {use the profile if Gen_profile isn't 'none'}
    Get_Integer(inv, log, 'Field_dep_G', dumint);  {field-dependent G, true or false}
    Field_dep_G:=(ROUND(dumint) = 1);
    Get_Float(inv, log, 'P0', P0); {0<=P0<1, fraction of quenched excitons that direcltly yield free carriers}
    Get_Float(inv, log, 'a', a); {thermalization length, Braun model used, m}
    Get_Integer(inv, log, 'ThermLengDist', dumint);
    ThermLengDist:=ROUND(dumint);
    Get_Float(inv, log, 'kf', kf); {decay rate of CT state, 1/s}
    Get_Float(inv, log, 'kdirect', kdirect); {m3/s, direct (band-to-band, bimolecular) recombination rate}
    Get_Float(inv, log, 'Lang_pre', Lang_pre); {Langevin prefactor}
    Get_Integer(inv, log, 'UseLangevin', dumint);
    UseLangevin:=(ROUND(dumint) = 1); {Calculate recombination using Langevin equation (1) or direct input (<>1, kdirect is used))}

{**Trapping**************************************************************************}
    {FOR ALL TRAP DENSITIES: IF < 0 ELECTRON TRAPS, IF > 0 HOLE TRAPS}
    {** Bulk traps}
    Get_Float(inv, log,'Bulk_tr', Bulk_tr); {density of electron or hole traps, m^-3}
    IF Bulk_tr >=0
	THEN BEGIN
		Pt:=Bulk_tr;
		Nt:=0;
	END
    ELSE BEGIN
      Pt:=0;
      Nt:=-Bulk_tr;
	END;
    {** Surface traps}
    Get_Float(inv, log, 'LTL', LTL); {m, thickness left trap layer}
    Get_Float(inv, log, 'LTR', LTR); {m, thickness right layer}
    Get_Float(inv, log, 'St_L', St_L); {m^-2, density of left surface traps}
    Get_Float(inv, log, 'St_R', St_R); {m^-2, density of right surface traps}
    {** Grain boundaries}
    Get_Integer(inv, log,'num_GBs',num_GBs); {number of grain boundaries}
    Get_Float(inv, log, 'GB_tr', GB_tr); {m^-2, grain boundary trap density per grain boundary}
    Get_Float(inv, log, 'L_GB', L_GB); {m, thickness of the grain boundary}
    {** traps coefficients}
    Get_Float(inv, log, 'Cn', Cn); {m^3/s, capture coefficient for electrons (put to 0 to exclude SRH)}
    Get_Float(inv, log, 'Cp', Cp); {m^3/s, capture coefficient for holes (put to 0 to exclude SRH) }
    Get_Integer(inv, log, 'Trtype', dumint);
    Trtype:=ROUND(dumint); {trap type electrons 0: single level, 1 : exponential}
    Get_Float(inv, log, 'Etrap', Etrap); {eV, electron trapping level, w.r.t. conduction band}
    Get_Float(inv, log, 'Ttr', Ttr); {K, charac. temperature of exponential distr. traps}
    Get_Integer(inv, log, 'Q_charge', dumint);
    Q_charge:=ROUND(dumint);{sign of filled trap 0: neutral, 1: charged}
	Traps:=(Nt<>0) OR (Pt<>0) OR (LTR>0) OR (LTL>0) OR (num_GBs>0); {are there any traps?}

{**Numerical Parameters**************************************************************}
    Get_Integer(inv, log, 'NP', NP); {number of grid points}
    Get_Float(inv, log, 'tolPois', tolPois); {abs. tolerance of Poisson solver}
	Get_Float(inv, log, 'maxDelV', maxDelV); {maximum change (in Vt) of the potential per loop}
    Get_Float(inv, log, 'accPois', accPois); {SOR/SUR acceleration parameter for Poisson, (0<accPois<2)}
    Get_Float(inv, log, 'accDens', accDens); {SOR/SUR acceleration parameter for densities, (0<accDens<2)}
    Get_Integer(inv, log, 'resetNegDens', dumint);
    resetNegDens:=(dumint = 1); {whether(1) or not(<>1) to reset points with a negative density}
    Get_Float(inv, log, 'tolMain', tolMain); {rel. tolerance of main loop}
    Get_Float(inv, log, 'Conv_VAR', dummy); {1 selects current to be monitored for convergence in main loop}
    Conv_Var:=ROUND(dummy); {<> selects the Slotboom variables (cf. densities)}
    Get_Float(inv, log, 'MaxItPois', dummy); {Max. number of loops Poisson solver}
    MaxItPois:=ROUND(dummy);
    Get_Float(inv, log, 'MaxitMain', dummy); {Max. number of main loops}
    MaxItMain:=ROUND(dummy);
    Get_Float(inv, log, 'grad', grad); {gradient of grid, increase grad for smaller h[1]}
    Get_Float(inv, log, 'TolRomb', TolRomb); {rel. tolerance of Romberg integration}
    Get_Float(inv, log, 'MaxRombIt', dummy);
    MaxRombIt:=ROUND(dummy); {max. # Romberg iterations}
    Get_Float(inv, log, 'LowerLimBraun', LowerLimBraun); {lower limit of integration over distribution Braun}
    Get_Float(inv, log, 'UpperLimBraun', UpperLimBraun); {upper limit}

{**Voltage range of simulation*******************************************************}
    Get_Float(inv, log, 'Vdistribution', dummy);
    Vdistribution:=ROUND(dummy); {type of V distribution, 1=linear, 2=logarithmic}
    Get_Integer(inv, log, 'PreCond', dumint); {Pre-condition in light(1)/dark(0)}
    PreCond:=ROUND(dumint)=1;
    Get_Float(inv, log, 'Vpre', Vpre); {V, pre-conditioned voltage}
    Get_Integer(inv, log, 'Vscan', Vscan); {integer, direction of voltage scan: up > 0, down < 0}
    Get_Float(inv, log, 'Vmin', Vmin); {V, minimum voltage in JV characteristic}
    Get_Float(inv, log, 'Vmax', Vmax); {V, max. voltage in JV}
    Get_Float(inv, log, 'Vstep', Vstep); {V, voltage step}
    Get_Float(inv, log, 'Vacc', Vacc); {accumulation voltage for logarithmic JV, should be outside [Vmin, Vmax]}
    Get_Integer(inv, log, 'NJV', NJV); {Number of JV points, for logarithmic JV}
    IF (Vstep<>0) AND (Vdistribution=1) THEN
        NJV:=TRUNC((Vmax - Vmin)/Vstep + 1e-10) + 1; {needed for setting length of Jdat and Vdat}
                                                         {1e-10 is needed to get right value}
    Get_Integer(inv, log, 'until_Voc', dumint); {if 1 then SIMsalabim stops at Voc}
	until_Voc:=(dumint=1);
	
{**User interface********************************************************************}
    Get_Integer(inv, log, 'Pause_at_end', Pause_at_end); {pause at the end of the simulation yes(1) or no (0)}
	Get_Integer(inv, log, 'AutoTidy', dumint);
	AutoTidy:=dumint = 1;	{if 1, then SIMsalabim will always tidy up the device_parameter file}
    Get_Integer(inv, log, 'UseExpData', dumint);
    UseExpData:=dumint = 1; {if 1 then  SIMsalabim will try to read ExpJV and use it}
    Get_String(inv, log, 'ExpJV', ExpJV); {name of file with experimental JV points}
    Get_String(inv, log, 'rms_mode', dumstr); {lin or log: use J or log(J) in calc. of rms error}
	dumstr:=lowercase(dumstr);
	IF NOT ((dumstr='lin') OR (dumstr='log')) THEN Stop_Prog('rms_mode has to be either lin or log.');
	IF dumstr='lin' THEN rms_mode:=linear ELSE rms_mode:=logarithmic;
    Get_Float(inv, log, 'rms_threshold', rms_threshold); {threshold of fraction converged points in calc. rms error}
    Get_String(inv, log, 'JV_file', JV_file); {name of file with simulated JV points}
    Get_String(inv, log, 'Var_file', Var_file);{name of file with various variables at V=Vmax}
    Get_String(inv, log, 'scPars_file', scPars_file);{solar cell parameter file}

    FLUSH(log);
    CLOSE(inv);
END;

PROCEDURE Check_Parameters;
{performs a number of checks on the parameters. Just to ensure that they are 
valid, consistent, make sense}
BEGIN
{when adding new check, please keep the order in line with the device_parameter file}
	{Check first if Vt and ni have been initialised. This should have happened in the main code}
	IF (Vt=0) OR (ni=0) THEN Stop_Prog('Vt and ni need to be initialised before calling Check_Parameters.');
{checks on general parameters:}
	IF CB >= VB THEN Stop_Prog('CB should be smaller than VB.');
	IF (CB<0) OR (VB<0) THEN Stop_Prog('CB and VB should be positive.');

{checks on mobilities:}
    IF NOT (mob_n_dep IN [0, 1, 2]) THEN Stop_Prog('Invalid mob_dep_n selected.');
    IF NOT (mob_p_dep IN [0, 1, 2]) THEN Stop_Prog('Invalid mob_dep_p selected.');
{checks on contacts:}
	{part of this can only be done after checking the TLs!}
	IF ImageLowering AND ( (SnMin>=0) OR (SpMin>=0) OR (SnMaj>=0) OR (SpMaj>=0) )
		THEN Stop_Prog('Cannot use finite surface recombination and image-force lowering at the same time.');
	{check if left electrode is cathode and right one is anode if using finite surface recombination:}
	IF (SnMin>=0) AND (phi_left>0.5*Egap) THEN Stop_Prog('Left contact is no longer the cathode. Cannot use finite or zero SnMin!');
	IF (SnMaj>=0) AND (phi_left>0.5*Egap) THEN Stop_Prog('Left contact is no longer the cathode. Cannot use finite or zero SnMaj!');
	IF (SpMin>=0) AND (phi_right>0.5*Egap) THEN Stop_Prog('Right contact is no longer the anode. Cannot use finite or zero SpMin!');
	IF (SpMaj>=0) AND (phi_right>0.5*Egap) THEN Stop_Prog('Right contact is no longer the anode. Cannot use finite or zero SpMaj!');
	IF (SnMaj = 0) OR (SpMaj=0) THEN Stop_Prog('Majority carrier surface recombination velocities cannot be zero.');
	IF Rshunt=0 THEN Stop_Prog('Rshunt cannot be zero, use positive (negative) value for finite (infinite) shunt resistance.');
{checks on transport layers:}
	IF L_LTL+L_RTL>=L THEN Stop_Prog('Sum of L_LTL and L_RTL (transport layers) should be less than total thickness L.');
	{more checks on the energy levels of the TLs:}
	{first check left:}
	IF L_LTL>0 
	THEN BEGIN
		IF W_L<CB_LTL THEN Stop_Prog('W_L cannot be smaller than CB_LTL.');
		IF W_L>VB_LTL THEN Stop_Prog('W_L cannot be larger than VB_LTL.');
		IF CB_LTL>= VB_LTL THEN Stop_Prog('CB_LTL must be smaller than VB_LTL');
	END 
	ELSE BEGIN {no TLs:}
		IF W_L<CB THEN Stop_Prog('W_L cannot be smaller than CB.');
		IF W_L>VB THEN Stop_Prog('W_L cannot be larger than VB.');
	END;
	{check the right:}
	IF L_RTL > 0 
	THEN BEGIN
		IF W_R<CB_RTL THEN Stop_Prog('W_R cannot be smaller than CB_RTL.');
		IF W_R>VB_RTL THEN Stop_Prog('W_R cannot be larger than VB_RTL.');
		IF CB_RTL>= VB_RTL THEN Stop_Prog('CB_RTL must be smaller than VB_RTL');
	END
	ELSE BEGIN
		IF W_R<CB THEN Stop_Prog('W_R cannot be smaller than CB.');
		IF W_R>VB THEN Stop_Prog('W_R cannot be larger than VB.');
	END;

{checks on ions:}
	IF NOT(ABS(mob_ion_spec) IN [0,1]) THEN Stop_Prog('Invalid mob_ion_spec selected.');
	{note: pascal set cannot contain negative numbers, hence the ABS()}	
	IF (ion_red_rate<0) THEN Stop_Prog('Scan rate cannot be lower than zero, should be: 0 <= ion_red_rate < NJV');
{checks on generation and recombination parameters}
	IF NOT (ThermLengDist IN [1,2,3,4,5]) THEN Stop_Prog('Invalid ThermLengDist selected.');
	IF (P0>=1) OR (P0<0) THEN Stop_Prog('Invalid value of P0, should be: 0<=P0<1');
    IF (P0<>0) AND (Field_dep_G = FALSE) THEN Stop_Prog('P0 should be zero if not using field dependent generation');
{checks on trapping:}
	IF Traps AND (NOT(Trtype IN [0,1])) OR (NOT(Trtype IN [0,1])) {validity of trap type}
		THEN Stop_Prog('You have selected an invalid e/h trap type.');
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
	IF UseExpData AND until_Voc 
		THEN Stop_Prog('You cannot use until_Voc = 1 and UseExpData = 1 at the same time.');
	IF ((rms_threshold<=0) OR (rms_threshold>1)) AND UseExpData 
		THEN Stop_Prog('rms_threshold has to be larger than 0 but not larger than 1.');
END;

PROCEDURE DefineLayers(VAR eps : vector; VAR i1, i2, i11, i22 : INTEGER);
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

    {Define trap layers:}
    {need i11, i22 to specify trap layer thickness}
    i11:=i1;
    WHILE x[i11+1]<L_LTL+LTL DO INC(i11);
    i22:=i2;
    WHILE x[i22-1]>L-L_RTL-LTR DO DEC(i22);

    {define dielectric constants of the layers:}
    FOR i:=0 TO NP+1 DO eps[i]:=eps_r * eps_0; {first define bulk}

    {now the transport layers}
    FOR i:=0 TO NP+1 DO
    BEGIN
		nid[i]:=0;
		pid[i]:=0;
    END;
    IF i1>0 THEN
		FOR i:=0 TO i1 DO
		BEGIN
			eps[i]:=eps_r_LTL * eps_0;
      IF doping_LTL<=0 THEN
      nid[i]:=-doping_LTL
      else pid[i]:=doping_LTL;
		END;
	IF i2<NP+1 THEN
		FOR i:=i2 TO NP+1 DO
		BEGIN
			eps[i]:=eps_r_RTL * eps_0;
      IF doping_RTL<=0 THEN
        nid[i]:=-doping_RTL
        else pid[i]:=doping_RTL;
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

PROCEDURE Read_Experimental_JV(VAR JVExp : JVList; VAR log : TEXT);
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
	IF NJV < minExpData THEN Stop_Prog('Not enough experimental data points.', FALSE);
	IF NJV = maxExpData THEN Stop_Prog('It looks like the experimental JV file contains more points than I can handle.', FALSE);
	
	{now we can copy the arrays V and J into ExpJV}
	NJV:=i; {number of data points in JV curve}
	SetLength(JVExp, NJV+1); {the regular voltages start at 1}
	FOR i:=1 TO NJV DO {note, we only use part of the JVList record}
	BEGIN
		JVExp[i].V:=V[i];
		JVExp[i].Vext:=V[i];
		JVExp[i].J:=J[i];
		JVExp[i].Jext:=J[i];
		JVExp[i].Use:=TRUE
	END;

	{voltage step, Vmin, Vmax: we will need these later on}
	Vstep:=ABS(JVexp[2].V - JVexp[1].V);
	IF JVexp[2].V > JVexp[1].V THEN Vscan:=1 ELSE Vscan:=-1; {determine scan direction}
	Vmin:=MIN(JVexp[1].V, JVexp[NJV].V); {we use min, max function as Vscan might be -1}
	Vmax:=MAX(JVexp[1].V, JVexp[NJV].V);
	Vdistribution:=1; {force the linear distribution of voltages}

	{check that the data is well-behaved, i.e. a constant step in voltage:}
	FOR i:=1 to NJV-1 DO
		IF ABS(JVexp[i+1].V - JVexp[i].V - Vscan * Vstep) > tolReal {comparing reals, hence the tolerance}
			THEN Stop_Prog('Experimental voltages should have a constant step in voltage.', FALSE);
	
	{now document what happened and write to log file:}
	WRITELN(log);
	WRITELN(log, 'Read experiment JV curve from ', ExpJV,'.');
	WRITELN(log, 'This overrides the voltage distribution that is in ',parameter_file);
	WRITELN(log, 'and any voltage parameters passed via the command line.');
	WRITELN(log, 'Vmin: ',Vmin:6:4,' Vmax: ',Vmax:6:4,' Vstep: ',Vstep:6:4,' Vscan: ',Vscan:2);
	WRITELN(log, 'Vdistribution: 1, so linearly distributed.');
	WRITELN(log);
	
	WRITELN('Read experimental JV curve from ', ExpJV,'.');
END;

PROCEDURE Init_Voltages_and_Tasks(VAR JVData : JVList; VAR log : TEXT);
VAR i, k : INTEGER;
BEGIN
	SetLength(JVData, NJV+1); {if there is a pre-bias, we need an extra data point, so we make the array one longer than NJV}
	IF PreCond THEN {the pre-bias will be stored in point 0}
		WITH JVData[0] DO {pre-bias point}
		BEGIN
			V:=Vpre;
			UpdateIons:=TRUE; {yes, ions are moving. We checked in Read_Parameters that if PreCond then CIM <>0}
			Store:=FALSE; {however, don't store the pre-bias point}
		END;
	
	k:=ORD(NOT PreCond); {so if PreCond=true then k=0, we will use this below}
	{now for the rest of the voltages:}
	FOR i:=1 TO NJV DO {the regular voltages start at 1}
	BEGIN
		WITH JVData[i] DO
		BEGIN
			V:=Applied_voltage(i);
			IF ion_red_rate=0 {means ions are fixed}
				THEN BEGIN
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
		WITH JVData[i] DO
			WRITELN(log,i:3,'  ',V:8:3,'    ',UpdateIons:5,'      ',Store);
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
         n_vals[i]:=Nc * dum;
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
         p_vals[i]:=Nc * dum;
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


PROCEDURE Trap_Distribution(VAR Nt_trap, Pt_trap, Q_trap : vector);
{Trap distribution along the thickness of the device}
VAR i, gb_nr, gb_start_idx, gb_end_idx, cur_idx : INTEGER;
    gb_centre_pos, gb_start_pos, gb_end_pos, L_GB_on_grid, L_ST_L_on_grid, L_ST_R_on_grid : myReal;
    new_gb_found : BOOLEAN;

    {grain_size,dum : myReal;}
BEGIN
	FOR i := 0 TO NP+1 DO {default}
		Q_trap[i] := Q_charge;
	FOR i :=0 TO NP+1 DO {default}
	BEGIN
		Nt_trap[i] := Nt; 
		Pt_trap[i] := Pt;
	END;

	IF LTL > 0 THEN {if there is traps at the left interface}
    BEGIN
        L_ST_L_on_grid := (x[i11+1] + x[i11]) / 2 - (x[i1] + x[i1+1]) / 2;
		FOR i := i1+1 TO i11 DO {at n-layer/absorber interface}
			IF ST_L > 0 
			THEN Pt_trap[i] := ST_L / L_ST_L_on_grid
			ELSE Nt_trap[i] := -St_L / L_ST_L_on_grid
    END;
    
	IF LTR > 0 THEN {if there is traps at the right interface}
    BEGIN	
        L_ST_R_on_grid := (x[i2-1] + x[i2]) / 2 - (x[i22-1] + x[i22]) / 2;  
		FOR i := i22 TO i2-1 DO {at p-layer/absorber interface}
			IF ST_R > 0 
			THEN Pt_trap[i] := ST_R / L_ST_R_on_grid
			ELSE Nt_trap[i] := -St_R / L_ST_R_on_grid
    END;
	
	IF num_GBs > 0 THEN {include grain boundary recombination}
    BEGIN
        cur_idx := i1;
        FOR gb_nr := 1 TO num_GBs DO
        BEGIN 
            new_gb_found := FALSE;
            gb_centre_pos := (L - L_RTL - L_LTL) / (num_GBs + 1) * gb_nr + L_LTL;
            gb_start_pos := gb_centre_pos - L_GB / 2;
            gb_end_pos := gb_centre_pos + L_GB / 2;

            {find the start and end idx of the grain boundary and set trap density}
            WHILE new_gb_found = FALSE DO 
            BEGIN 
                IF (gb_start_pos > x[cur_idx]) AND (gb_start_pos < x[cur_idx+1]) THEN 
                BEGIN 
                    gb_start_idx := cur_idx+1;
                    WHILE gb_end_pos > x[cur_idx] DO 
                    BEGIN
                        gb_end_idx := cur_idx;
                        inc(cur_idx) 
                    END;
                    new_gb_found := TRUE;

                    {put traps in the newly found grain boundary}
                    L_GB_on_grid := (x[gb_end_idx+1] + x[gb_end_idx]) / 2 - (x[gb_start_idx] + x[gb_start_idx-1]) / 2;
                    IF GB_tr > 0 THEN 
                        FOR i := gb_start_idx TO gb_end_idx DO 
                        Pt_trap[i] := GB_tr / L_GB_on_grid
                    ELSE 
                        FOR i := gb_start_idx TO gb_end_idx DO 
                        Nt_trap[i] := -GB_tr / L_GB_on_grid
                END;

                inc(cur_idx);
            END
        END   
    END;

	IF NOT TLsTrap THEN {whether traps in transport layers:}
	BEGIN
		{left transport layer}
		IF i1 > 0 THEN
			FOR i := 0 TO i1 DO
			BEGIN
				Nt_trap[i] := 0;
				Pt_trap[i] := 0;
			END;
		{right transport layer}
		IF i2 < NP+1 THEN
			FOR i:=i2 TO NP+1 DO
			BEGIN
				Nt_trap[i] := 0;
				Pt_trap[i] := 0;
			END;
	END;
END;

PROCEDURE Calc_Ion_Distribution(VAR nion, pion, V : vector);
{calculates the ion distribution in steady state. This is either constant 
(if the ionic species doesn't move), or it follows from requiring that
the associated ionic particle current be zero.} 
VAR nrmn, nrmp, fac, sum_nion, sum_pion, Vtinv : myReal;
	istart, ifinish : INTEGER;
BEGIN
	sum_nion:=0; {sums of ion concentration}
	sum_pion:=0;
	Vtinv:=1/Vt;
	
	{ions are limited to the middle layer (which can take up the whole device i=0...NP+1.}
	IF i1>0 THEN istart:=i1+1 ELSE istart:=0;
	IF i2<NP+1 THEN ifinish:=i2-1 ELSE ifinish:=NP+1;
	
	nion[istart]:=1;
	pion[istart]:=1;

	FOR i:=istart+1 TO ifinish DO
	BEGIN
		fac:=B(Vtinv*(V[i-1]-V[i]))/B(Vtinv*(V[i]-V[i-1]));
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
	
	nrmn:=CIM*(L-L_LTL-L_RTL)/sum_nion;
	nrmp:=CIM*(L-L_LTL-L_RTL)/sum_pion;
	FOR i:=istart TO ifinish DO
	BEGIN
		nion[i]:=nion[i]*nrmn;
		pion[i]:=pion[i]*nrmp
	END;
	{now the total number of ions should be equal to the concentrion CMI times (L-L_LTL-L_RTL)}
END;

PROCEDURE Solve_Poisson(VAR V, n, p, nid, pid : vector; UpdateIons : BOOLEAN; VAR conv : BOOLEAN);

VAR it, i : INTEGER;
    delV, rhs, lower, upper, main : vector;
    fac, ntd, ntd2, ptd, ptd2 : myReal; {n (resp. p) trap dummy - useful for notation}
BEGIN
    FOR i:=1 TO NP DO delV[i]:=1; {init delV}
    delV[0]:=0; {delV=0 at contacts, by definition, since we know V(0,L)}
    delV[NP+1]:=0;
    it:=0;
    fac:=0.5*q*L*L;
    WHILE (Norm(delV, 0, NP+1) > tolPois) AND (it < MaxItPois) DO
    BEGIN
        {first: fill the matrix of Poisson equations, with or without traps}
        IF Traps {there are traps, use different form of matrix}
        THEN FOR i:=1 TO NP DO {filling the matrix en the right-hand side}
                BEGIN
                    {dummy needed for single trap level: }
                    ntd:=(Nc/n[i])*EXP(-Etrap/Vt);
                    ptd:=(Nc/p[i])*EXP(-Etrap/Vt);
                    {dummies for exponential trap distribution: }
                    ntd2:=Nt*Power(n[i]/Nc, T/Ttr);
                    ptd2:=Pt*Power(p[i]/Nc, T/Ttr);
                    {note: factors Trtype, (1-Trtype), etc. select which terms
                    are used to take the traps into account in rhs[i] and main[i]:}
                    rhs[i]:=-eps[i]*h[i-1]*V[i+1] + (h[i]*eps[i-1]+h[i-1]*eps[i])*V[i] - eps[i-1]*h[i]*V[i-1]
							+fac*(h[i]+h[i-1])*h[i]*h[i-1]*(n[i]-p[i]-nid[i]+pid[i] + nion[i]-pion[i]
							+ (1-Trtype)*Nt_trap[i]/(1+ntd)*Q_trap[i] - (1-Trtype)*Pt_trap[i]/(1+ptd)*Q_trap[i]
							+ Trtype*ntd2 - Trtype*ptd2 );
					lower[i]:=eps[i-1]*h[i];
					upper[i]:=eps[i]*h[i-1];
					main[i]:=-(h[i-1]*eps[i]+h[i]*eps[i-1]) - fac*(h[i]+h[i-1])*h[i]*h[i-1]/Vt
							* (n[i] + p[i] + nion[i] + pion[i] + (1-Trtype)*(Q_trap[i]*Nt_trap[i]*ntd)/SQR(1+ntd) + (1-Trtype)*Q_trap[i]*Pt_trap[i]*ptd/SQR(1+ptd) + Trtype*ntd2*T/Ttr + Trtype*ptd2*T/Ttr )
                END
        ELSE FOR i:=1 TO NP DO {filling the matrix and the right-hand side}
                BEGIN  {there are no traps}
                    rhs[i]:=-eps[i]*h[i-1]*V[i+1] - eps[i-1]*h[i]*V[i-1] + (h[i]*eps[i-1]+h[i-1]*eps[i])*V[i]
							+ fac*(n[i]-p[i]-nid[i]+pid[i] + nion[i]-pion[i])*(h[i]+h[i-1])*h[i]*h[i-1];
					lower[i]:=eps[i-1]*h[i];
					upper[i]:=eps[i]*h[i-1];
					main[i]:=-(h[i-1]*eps[i]+h[i]*eps[i-1]) - (n[i]+p[i] + nion[i]+pion[i])*fac*(h[i]+h[i-1])*h[i]*h[i-1]/Vt
                END;
        Tridiag(delV, lower, main, upper, rhs, 1, NP); {solve for delV}
        FOR i:=1 TO NP DO BEGIN
            delV[i]:=accPois*delV[i]; {use SOR/SUR}
			IF ABS(delV[i])>maxDelV*Vt THEN delV[i]:=SIGN(delV[i])*maxDelV*Vt;         
            V[i]:=V[i]+delV[i];  {and add delV to V}
            n[i]:=n[i]*EXP(delV[i]/Vt); {now update the densities: we have to do this}
            p[i]:=p[i]*EXP(-delV[i]/Vt); {in order to conserve Gummel iteration}

        IF UpdateIons THEN
           	CASE mob_ion_spec OF
    			-1 : nion[i]:=nion[i]*EXP(delV[i]/Vt);
    			 0 : BEGIN
    					nion[i]:=nion[i]*EXP(delV[i]/Vt);
    					pion[i]:=pion[i]*EXP(-delV[i]/Vt);
    				END;
    			1 : pion[i]:=pion[i]*EXP(-delV[i]/Vt);
    		END
        END;
        it:=it+1;
    END;
    IF Norm(delV, 0, NP+1) <= tolPois  {finally, check for convergence}
        THEN conv:=TRUE
        ELSE conv:=FALSE
END;

PROCEDURE UpdateGenPot(VAR V, Vgn, Vgp : vector);
VAR i : INTEGER;
BEGIN
	Vgn:=V;
	Vgp:=V;
	{left transport layer, offset generalised potentials:}
	IF i1>0 THEN
		FOR i:=0 TO i1 DO
		BEGIN
			Vgn[i]:=V[i] + CB_LTL - CB;
			Vgp[i]:=V[i] - VB + VB_LTL;
		END;
	{right transport layer, offset generalised potentials:}
	IF i2<NP+1 THEN
		FOR i:=i2 TO NP+1 DO
		BEGIN
			Vgn[i]:=V[i] + CB_RTL - CB;
			Vgp[i]:=V[i] - VB + VB_RTL;
		END;
END;

PROCEDURE Init_Pot_Dens(VAR V, Vgn, Vgp, n, p, empty : vector; Va : myReal);
{init. for V, Vgn,p, n, p at bias voltage Va}
VAR i : INTEGER;
BEGIN
    FOR i:=0 TO NP+1 DO {guess new V in interior of device} {linear interpolation of V}
	    V[i]:=V0 - Va/2 + x[i]*(VL-V0+Va)/L;
	UpdateGenPot(V, Vgn, Vgp); {update generalised potentials}

    FOR i:=0 TO NP+1 DO {guess n, p in interior of device}
    BEGIN
		n[i]:=ni*EXP( (Va * (0.5 -x[i]/L) + Vgn[i])/Vt ); {note: we use Vgn,p, NOT V}
        p[i]:=ni*EXP( (Va * (x[i]/L - 0.5) - Vgp[i])/Vt );
        empty[i]:=0;  {Empty vector, used to separately calculate Lan & SRH Recombination}
     END;
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
            END; {case 5}
    END; {case selector}

    mu[NP+1]:=mu[NP]; {does not have a meaning, should not be used}
    IF i1>0 THEN
		FOR i:=0 TO i1 DO mu[i]:=mob_LTL ;
	IF i2<NP+1 THEN
		FOR i:=i2 TO NP+1 DO mu[i]:=mob_RTL;
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
        2 : BEGIN {hole mobility from table in file n_file}
            FOR i:=0 TO NP DO
                BEGIN
                     conc:=0.5*(n[i] + n[i+1]);
                     {note: n, V are defined on grid, but mu on interleaved mesh!}
                     Field:=ABS(V[i+1]-V[i])/(L*h[i]);
                     IF extra_ln_mob {decide whether to extrapolate ln(mob) or just mob}
                     THEN mu[i]:=EXP(BilinearInterpolation(conc, Field, p_vals, Fp_vals, p_points, Fp_points, p_mob_tab, extra_F))
                     ELSE mu[i]:=BilinearInterpolation(conc, Field, p_vals, Fp_vals, p_points, Fp_points, p_mob_tab, extra_F);
                END; {for loop}
            END; {case 5}
    END; {case selector}
    mu[NP+1]:=mu[NP]; {does not have a meaning, should not be used}

    IF i1>0 THEN
		FOR i:=0 TO i1 DO mu[i]:=mob_LTL ;
	IF i2<NP+1 THEN
		FOR i:=i2 TO NP+1 DO mu[i]:=mob_RTL;
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
    IF UseLangevin {use Langevin formula to calculate bimolecular, direct, band-to-band recombination rate}
    THEN FOR i:=1 TO NP DO
		BEGIN
			rec_mob:=(mob_n[i-1] + mob_n[i]+ mob_p[i-1] + mob_p[i] )/2;
			{we take mob(x=xi)=(mob(x=xi-1/2)+mob(x=xi+1/2))/2}
			Lan[i]:=Lang_pre * q * rec_mob/eps[i];
		END
	ELSE FOR i:=1 TO NP DO {use input value for direct recombination}
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
    kdF:=3*Braun_rec/(4*PI*r*SQR(r))*EXP(-delE/Vt)*Bessel(b);
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

FUNCTION SRH_pre_EXP(E : myReal) : myReal;
{ Computes the Shockley-Read-Hall prefactor for a fixed energy in one grid point}
{ for an exponential energetic distribution of electron traps}
VAR nt0, pt0, Vtrapn, Vtrapp : myReal;
BEGIN
    nt0:=Nc*EXP(-E/Vt);
    pt0:=Nc*EXP( (E-Egap)/Vt );
    Vtrapn:=k*Ttr/q; {thermal voltage with T=Ttr}
    Vtrapp:=k*Ttr/q;
    SRH_pre_EXP:=(EXP(-E/Vtrapn) *Nt/Vtrapn + EXP(-E/Vtrapp)*Pt/Vtrapp )
                /( (n[gridpoint]+nt0)/cn + (p[gridpoint]+pt0)/cp );
END;

PROCEDURE Calc_SRH_Prefactor(VAR srhp : vector; V, n, p : vector);
{Calculates the Shockley-Read-Hall recombination prefactor,
 i.e., R_SRH = srhp * (np - n_int^2), where R_SRH is the SRH recombination rate }
VAR i : INTEGER;
    nt0, pt0 : myReal;
BEGIN
    IF (Cn= 0) OR (Cp = 0) OR (NOT Traps)
    THEN FOR i:=0 TO NP+1 DO srhp[i]:=0 {no SRH recombination}
    ELSE CASE Trtype OF
        0 : BEGIN {single electron trapping level}
                nt0:=Nc*EXP(-Etrap/Vt);
                pt0:=Nc*EXP( (Etrap-Egap)/Vt );
                FOR i:=1 TO NP DO
                BEGIN
                    srhp[i]:=(Nt_trap[i] + Pt_trap[i])/ ( (n[i] + nt0) /Cp
                            + (p[i] + pt0)/Cn)
                END;
            END;
        1 : BEGIN { exponential electron trapping levels}
                FOR i:=1 TO NP DO
                BEGIN 
					gridpoint:=i;
                    srhp[i]:=RombergIntegration(SRH_pre_EXP,0,Egap, TolRomb, MaxRombIt, FALSE)
					{the variables (gridpoint,n,p) are global, so their scope}
					{includes the FUNCTION SRH_pre, so we don't have to pass}
					{them on through the function heading! => we can use RombergIntegration}
                END; {for loop}
            END {case 1: exponential traps}
    END; {CASE selector}
END;

PROCEDURE Recombi(VAR r : vector; V, n, p, g, Lan, dp, srhp : vector);
VAR i : INTEGER;
BEGIN {Langevin and Shockley-Read-Hall recombination with generation term G}
    FOR i:=0 TO NP+1 DO 
        {note that srhp, g, and, lan are empty depending on calculated recombination}
        r[i]:=( (1-dp[i])*Lan[i]+srhp[i] ) * (n[i]*p[i]-SQR(ni)) - g[i] 
END;


PROCEDURE Calc_Recombi_current(VAR Jbimo, JSRH_bulk, JSRH_LI, JSRH_RI, Jph, Jn_l, Jn_r, Jp_l, Jp_r : myReal);
VAR i : INTEGER;
BEGIN {Langevin and Shockley-Read-Hall recombination with generation term G}
    Jbimo:=0;
    JSRH_bulk:=0;
    JSRH_LI:=0;
    JSRH_RI:=0;
    Jph:=0;
    FOR i:=0 TO NP+1 DO Jbimo := Jbimo + q*L*h[i]*( (1-dp[i])*Lang[i]) * (n[i]*p[i]-SQR(ni));
    FOR i:=i11+1 TO i22-1 DO JSRH_bulk := JSRH_bulk + q*L*h[i]*( SRH[i] ) * (n[i]*p[i]-SQR(ni));
    FOR i:=i1+1 TO i11 DO JSRH_LI := JSRH_LI + q*L*h[i]*( SRH[i] ) * (n[i]*p[i]-SQR(ni));
    FOR i:=i22 TO i2-1 DO JSRH_RI := JSRH_RI + q*L*h[i]*( SRH[i] ) * (n[i]*p[i]-SQR(ni));
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
		n[0]:=4*SQR(psi)*Nc*EXP(SQRT(f)-phi_left/Vt);
		n[0]:=n[0]*EXP( RoughLeft * SQRT(q*F0/(4*pi*eps[0]))/(2*Vt)); {now account for influence of roughness}
		n[0]:=MIN(n[0],Nc); {n[0] cannot exceed Nc}
		p[0]:=SQR(ni)/n[0];
	END;
	phi_left_eff:=-Vt*LN(n[0]/Nc); {calculate effective barrier height}

	{now do the same at the anode:}
	F0:=(V[NP+1]-V[NP])/(L*h[NP]); {electric field at anode}
	IF F0>0 THEN {only if field is positive do we have lowering:}
	BEGIN
		f:=q*F0/(4*pi*eps[NP+1]*SQR(Vt));
		psi:=1/f * (1-SQRT(1+2*SQRT(f))) + 1/SQRT(f);
		p[NP+1]:=4*SQR(psi)*Nc*EXP(SQRT(f)-phi_right/Vt);
		p[NP+1]:=p[NP+1]*EXP( RoughRight * SQRT(q*F0/(4*pi*eps[NP+1]))/(2*Vt)); {now account for influence of roughness}
		p[NP+1]:=MIN(p[NP+1],Nc); {p[NP+1] cannot exceed Nc}
		n[NP+1]:=SQR(ni)/p[NP+1];
	END;
	phi_right_eff:=-Vt*LN(p[NP+1]/Nc); {calculate effective barrier height}
END;


PROCEDURE Contn(VAR n : vector; V, p, mu, D, g, Lan, dp, srhp : vector);
VAR i : INTEGER;
    lo, m, u, rhs : vector;
BEGIN

	{Set the boundary conditions for finite surface recombination on the left side of the device}
	IF (SnMaj < 0) THEN
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
		m[0]:=-D[0]*B(mu[0]*(V[0]-V[1])/D[0])/(L*h[0]) - SnMaj;
		u[0]:=D[0]*B(mu[0]*(V[1]-V[0])/D[0])/(L*h[0]);
		rhs[0]:= - SnMaj*NC * EXP(-phi_left/Vt);
	END;

    FOR i:=1 TO NP DO  {continuity eq. in matrix vorm}
    BEGIN
        rhs[i]:=-0.5 * ((1-dp[i])*Lan[i]*SQR(ni)+g[i]+srhp[i]*SQR(ni))
            *SQR(L)*h[i]*h[i-1]*(h[i]+h[i-1]);
        lo[i]:=h[i]*D[i-1]*B(mu[i-1]*(V[i-1]-V[i])/D[i-1]);
        m[i]:=-(h[i-1]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i]) +
            h[i]*D[i-1]*B(mu[i-1]*(V[i]-V[i-1])/D[i-1]))
            -0.5*((1-dp[i])*Lan[i]*p[i]+ p[i]*srhp[i])*SQR(L)*h[i]*h[i-1]*(h[i]+h[i-1]);
         u[i]:=h[i-1]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i]);
    END;

	{Set the boundary conditions on the right side of the device}
	IF (SnMin < 0) THEN
	BEGIN
		{Infinite surface recombination}
		lo[NP+1]:= 0;
		m[NP+1]:= 1;
		u[NP+1]:= 0;
		rhs[NP+1]:= n[NP+1];
	END
	ELSE
	BEGIN
		{Finite surface recombination}
		lo[NP+1]:=-D[NP] * B(mu[NP]*(V[NP]-V[NP+1])/D[NP])/(L*h[NP]);
		m[NP+1]:=D[NP] * B(mu[NP]*(V[NP+1]-V[NP])/D[NP])/(L*h[NP]) + SnMin;
		u[NP+1]:=0;
		rhs[NP+1]:=SnMin*NC * EXP(-(Egap-phi_right)/Vt);
	END;

	Tridiag(n, lo, m, u, rhs, 0, NP+1); {Solve for the new electron densities}
    FOR i:=1 TO NP DO {and check whether n>=0}
        IF n[i] < 0 THEN BEGIN
			IF resetNegDens THEN n[i]:=1 ELSE Stop_Prog('Negative electron concentration encountered!')
	END;
END;


PROCEDURE Contp(VAR p : vector; V, n,  mu, D, g, Lan, dp, srhp : vector);
VAR i : INTEGER;
    lo, m, u, rhs : vector;
BEGIN

	{Set the boundary conditions for finite surface recombination on the left side of the device}
	IF (SpMin < 0) THEN
	BEGIN
		{Infinite surface recombination}
		lo[0]:=0;
		m[0]:=1;
		u[0]:=0;
		rhs[0]:=p[0];
	END
	ELSE
	BEGIN
		{Finite surface recombination}
		lo[0]:= 0;
		m[0]:=D[0] * B(mu[0]*(V[1]-V[0])/D[0])/(L*h[0]) + SpMin;
		u[0]:=-D[0] * B(mu[0]*(V[0]-V[1])/D[0])/(L*h[0]);
		rhs[0]:=-SpMin*NC * EXP(-(Egap-phi_left)/Vt);
	END;

    FOR i:=1 TO NP DO  {continuity eq. in matrix vorm}
    BEGIN
        rhs[i]:=-0.5 * ((1-dp[i])*Lan[i]*SQR(ni)+g[i]+srhp[i]*SQR(ni))
            *SQR(L)*h[i]*h[i-1]*(h[i]+h[i-1]);
        lo[i]:=h[i]*D[i-1]*B(mu[i-1]*(V[i]-V[i-1])/D[i-1]);
        m[i]:=-(h[i-1]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i]) +
            h[i]*D[i-1]*B(mu[i-1]*(V[i-1]-V[i])/D[i-1]))
            -0.5*((1-dp[i])*Lan[i]*n[i]+ n[i]*srhp[i])*SQR(L)*h[i]*h[i-1]*(h[i]+h[i-1]);
         u[i]:=h[i-1]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i]);
    END;

	{Set the boundary conditions on the right side of the device}
	IF (SpMaj < 0) THEN
	BEGIN
		{Infinite surface recombination}
		lo[NP+1]:=0;
		m[NP+1]:=1;
		u[NP+1]:=0;
		rhs[NP+1]:=p[NP+1];
	END
	ELSE
	BEGIN
		{Finite surface recombination}
		lo[NP+1]:=D[NP]*B(mu[NP]*(V[NP+1]-V[NP])/D[NP])/(L*h[NP]);
		m[NP+1]:=-D[NP]*B(mu[NP]*(V[NP]-V[NP+1])/D[NP])/(L*h[NP]) - SpMaj;
		u[NP+1]:=0;
		rhs[NP+1]:=-SpMaj*NC * EXP(-phi_right/Vt);
	END;

	Tridiag(p, lo, m, u, rhs, 0, NP+1); {Solve for the new hole densities}
    FOR i:=1 TO NP DO {and check whether p>=0}
        IF p[i] < 0 THEN BEGIN
			IF resetNegDens THEN p[i]:=1 ELSE Stop_Prog('Negative hole concentration encountered!')
	END;
END;

PROCEDURE Calc_elec_curr(VAR Curr, V, n, D, mu : vector);
{This procedure calculates the electron current density from Selberherr eq. 6.1-39 modified for generalised Einstein relation}
VAR i : INTEGER;
BEGIN
	FOR i:=0 TO NP DO
			Curr[i]:=-q*(n[i+1]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i]) - n[i]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i]))/(L*h[i])
END;

PROCEDURE Calc_hole_curr(VAR Curr, V, p, D, mu : vector);
{This procedure calculates the electron current density from Selberherr eq. 6.1-41 modified for generalised Einstein relation}
VAR i : INTEGER;
BEGIN
	FOR i:=0 TO NP DO
		Curr[i]:=-q*(p[i]*D[i]*B(mu[i]*(V[i+1]-V[i])/D[i]) - p[i+1]*D[i]*B(mu[i]*(V[i]-V[i+1])/D[i]))/(L*h[i])
END;

PROCEDURE Calc_Elec_Trap_Dens(V, n : vector; VAR ntr : vector);
VAR i : INTEGER;
BEGIN
    CASE Trtype OF
        0 : FOR i:=0 TO NP+1 DO {single trap level}
                ntr[i]:=Nt_trap[i]/(1 + (Nc/n[i])*EXP(-Etrap/Vt));
        1 : FOR i:=0 TO NP+1 DO {exponential trap distribution}
                ntr[i]:=Nt_trap[i]*Power(n[i]/Nc, T/Ttr)
    END
END;

PROCEDURE Calc_Hole_Trap_Dens(V, p : vector; VAR ptr : vector);
VAR i : INTEGER;
BEGIN
    CASE Trtype OF
        0 : FOR i:=0 TO NP+1 DO {single trap level}
                ptr[i]:=Pt_trap[i]/(1 + (Nc/p[i])*EXP(-Etrap/Vt));
       1 : FOR i:=0 TO NP+1 DO {exponential trap distribution}
                ptr[i]:=Pt_trap[i]*Power(p[i]/Nc, T/Ttr)
    END
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
		ELSE 
			WITH SCParExp DO 
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
                bracketed:=(JVExp[j].V < JVExp[i].Vext) AND (JVExp[j+1].V > JVExp[i].Vext);
                WHILE (NOT bracketed) AND (j < NJV-1) DO
				BEGIN
					INC(j);
					bracketed:=(JVExp[j].V < JVExp[i].Vext) AND (JVExp[j+1].V > JVExp[i].Vext);
				END;
				
				IF bracketed THEN 
					JVExp[i].Jext:=JVExp[j].J + (JVExp[j+1].J - JVExp[j].J) * (JVExp[i].Vext - JVExp[j].V) / (JVExp[j+1].V - JVExp[j].V)
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
        THEN BEGIN
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
			END;
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
    ELSE BEGIN
		WRITELN('Could not compute a meaningful rms-error.');
		WRITELN('Possible reasons:');
		WRITELN('-not enough simulated points converged, check rms_threshold');
		WRITELN('-the range in experimental currents is not large enough.');
		IF rms_mode = logarithmic 
			THEN WRITELN('-too many voltages were sim. and exp. currents have different signs. Try using rms_mode = lin.');
		WRITELN
	END
END;

PROCEDURE Write_Vars_To_File(VAR V, n, p, Vgn, Vgp, n_trap, p_trap, rec, Lang, diss_prob, SRH, nid, pid, mun, mup, Gm, Jn, Jp, 
								nion, pion : vector; Va : myReal);
{writes all relevant parameters to file. } 
VAR i : INTEGER;
	Evac, Ec, Ev, phin, phip : myReal;
	uitv : TEXT;
BEGIN
	{Create file:}
	ASSIGN(uitv, Var_file);
	REWRITE(uitv);

	Calc_Elec_Trap_Dens(Vgn, n, n_trap); {Calc. the densities of trapped charges}
	Calc_Hole_Trap_Dens(Vgp, p, p_trap);
	Recombi(rec,V, n, p, empty, Lang, diss_prob, SRH); {calc. net recombination}

	{note: we limit the length to 25 (nd) characters, otherwise, if myReal=EXTENDED, Origin cannot read the file!} 
	WRITELN(uitv, 'x V n p Evac Ec Ev phin phip ntrap ptrap nid pid nion pion mun mup rec dp Gm Jn Jp');
	FOR i:=0 TO NP+1 DO
	BEGIN
		Evac:=V[0] - V[i]; {the vacuum level. We take it zero at x=0}
		{use the generalised potentials to take care of the band-offsets}
		Ec:=V[0] - Vgn[i] - CB;
		Ev:=V[0] - Vgp[i] - VB;
		{electron and hole quasi-Fermi levels:}
		phin:=Ec + Vt*LN(n[i]/Nc);
		phip:=Ev - Vt*LN(p[i]/Nc);
		WRITELN(uitv,x[i]:nd,' ',V[i]:nd,' ',n[i]:nd,' ',p[i]:nd,' ',Evac:nd,' ',Ec:nd,' ',Ev:nd,' ',
				phin:nd,' ',phip:nd,' ',n_trap[i]:nd,' ',p_trap[i]:nd,' ',nid[i]:nd,' ',pid[i]:nd,' ',
				nion[i]:nd,' ', pion[i]:nd,' ',mun[i]:nd,' ',mup[i]:nd,' ',rec[i]:nd,' ',
				diss_prob[i]:nd,' ',Gm[i]:nd,' ',Jn[i]:nd,' ',Jp[i]:nd);
	END;
	CLOSE(uitv);

	WRITELN('The values of (V, n, p) for the last voltage are written in ', Var_file,'.');
END;


BEGIN {main program}
    WRITELN('Welcome to SIMsalabim version ',version,'.');
	WRITELN('Copyright (C) 2020 Dr T.S. Sherkar, V.M. Le Corre, M. Koopmans');
	WRITELN('F.O.B. Wobben, and Prof L.J.A. Koster, University of Groningen.');
    WRITELN;
    {if '-h' or '-H' option is given then display some help and exit:}
    IF hasCLoption('-h') THEN DisplayHelpExit;
    IF hasCLoption('-tidy') THEN Tidy_up_parameter_file(TRUE); {tidy up file and exit}
    IF NOT Correct_version_parameter_file THEN Stop_Prog('Version of SIMsalabim and '+parameter_file+' do not match.');

{Initialisation:}
    Prepare_Log_File(log); {open log file}
    Read_Parameters(log); {Read parameters from input file}
	IF AutoTidy THEN Tidy_up_parameter_file(FALSE); {clean up file but don't exit!}
    Egap:=VB-CB; {define band gap}
    Vt:=k*T/q;  {thermal voltage}
    ni:=Nc*EXP(-Egap/(2*Vt)); {equilibrium concentration}
    Check_Parameters; {perform a number of chekcs on the paramters. Note: we need Vt and ni}
    Make_Grid(h, x); {Initialize the grid}
    DefineLayers(eps, i1, i2, i11, i22); {define layers}
    Trap_Distribution(Nt_trap, Pt_trap, Q_trap);
	IF UseExpData THEN Read_Experimental_JV(JVExp, log);
	Init_Voltages_and_Tasks(JVSim, log);

    IF PreCond
		THEN VCount:=0 {Counter of the number of voltages which have been computed}
		ELSE VCount:=1; {index VCount=0 is reserved for the pre-bias voltage (if any)}

    Va:=JVSim[VCount].V; {get applied voltage from the list}
	Vaold:=Va; {Vaold: old value of Va, in 1st loop equal to Va}

    Init_Generation_Profile(Use_gen_profile, Gen_profile, Gm); {init. the Gm array}

    IF mob_n_dep = 2 THEN Init_Elec_Mob_Table(n_vals, Fn_vals, n_points, Fn_points, n_mob_tab); {init. the table with values of mob_n vs. F and n}
    IF mob_p_dep = 2 THEN Init_Hole_Mob_Table(p_vals, Fp_vals, p_points, Fp_points, p_mob_tab); {init. the table with values of mob_p vs. F and p}

	Init_Pot_Dens(V, Vgn, Vgp, n, p, empty, Va); {init. (generalised) potentials and densities}

    ASSIGN(uitv, JV_file);  {create the JV-file}
	REWRITE(uitv);
	WRITELN(uitv,'V   J   P  recLan  recSRH  Jbimo  JSRH_bulk  JSRH_LI  JSRH_RI  Jph  Jn_l  Jp_l  Jn_r  Jp_r  Vext  Jext');
	WRITELN('The calculation has started, please wait.');
	
	quit_Voc:=FALSE;

    WHILE (Vcount<=NJV) AND (NOT quit_Voc) DO  
    BEGIN
        Vaold:=Va; {keep old applied voltage}
        Va:=JVSim[VCount].V; {get new applied voltage from the list}
        FOR i:=0 TO NP+1 DO {guess new V, n, p in interior of device}
			V[i]:=V[i] + (Va-Vaold) * (x[i]-L/2)/L; {new V guessed from old V and linear term proportional to Va-Vaold}
 
        Solve_Poisson(V, n, p, nid, pid, JVSim[VCount].UpdateIons, check_Poisson);
        Calc_Elec_Trap_Dens(V, n, n_trap); {Calc. the densities of trapped charges}
        Calc_elec_mob(mun, V, n); {calc. mobilities with the new V}
        Calc_hole_mob(mup, V, p);
        Calc_Dn(Dn, mun, n); {calc. the diffusivities}
        Calc_Dp(Dp, mup, p);
        Calc_Langevin_factor(Lang, mun, mup);{Calc. the Langevin recombination strength}
        Calc_Dissociation(diss_prob, gen, Lang, V, mun, mup); {calculate the field-dependent generation}
        Calc_SRH_Prefactor(SRH, V, n, p); {calc. the Shockley-Read-Hall recombi. prefactor}
		Calc_Doping(V, nid, pid); {calc. ionised doping densities}

        MainIt:=0;
        delDens:=TolMain + 1;  {to ensure that delDens>TolMain in 1st loop}
        Jtot:=0;
        totn:=0.5;
        totp:=0.5; {this will ensure delDens/(totn + totp)>TolMain in 1st loop}
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

            UpdateGenPot(V, Vgn, Vgp); {update generalised potentials}
            Contn(n, Vgn, oldp, mun, Dn, gen, Lang, diss_prob, SRH); {calc. new elec. density}
            Contp(p, Vgp, oldn, mup, Dp, gen, Lang, diss_prob, SRH); {calc. new hole density}

			{now use SOR (to speed up) or SUR (to slow down) changes in n and p:}
			FOR i:=0 TO NP+1 DO
			BEGIN
				n[i]:=accDens * n[i] + (1-accDens) * oldn[i];
				p[i]:=accDens * p[i] + (1-accDens) * oldp[i];
			END;

			Solve_Poisson(V, n, p, nid, pid, JVSim[Vcount].UpdateIons, check_Poisson);

			{Update ions if needed:}
            IF JVSim[Vcount].UpdateIons THEN Calc_Ion_Distribution(nion, pion, V); 

            Calc_Elec_Trap_Dens(V, n, n_trap); {Calc. the densities of trapped charges}
            Calc_elec_mob(mun, V, n); {calc. the new mobilities}
            Calc_hole_mob(mup, V, p);
            Calc_Dn(Dn, mun, n); {update the diffusivities}
            Calc_Dp(Dp, mup, p);
            Calc_Langevin_factor(Lang, mun, mup); {Calc. the Langevin recombination strength}
            Calc_Dissociation(diss_prob, gen, Lang, V, mun, mup); {update generation rate}
            Calc_SRH_Prefactor(SRH, V, n, p); {calc. the Shockley-Read-Hall recombi. prefactor}
            UpdateGenPot(V, Vgn, Vgp); {update generalised potentials}
			Calc_elec_curr(Jn, Vgn, n, Dn, mun); {calc. the current densities}
			Calc_hole_curr(Jp, Vgp, p, Dp, mup);
			Calc_Doping(V, nid, pid); {calc. ionised doping densities}

            {Check for convergence of main loop:}
            FOR i:=0 TO NP+1 DO
            BEGIN
                totn:=totn + n[i];
                totp:=totp + p[i]
            END;
            delDens:=Norm(Difference(oldn,n, 0, NP+1), 0, NP+1) + Norm(Difference(oldp,p, 0, NP+1), 0, NP+1);
            {delDens = total change in the loop}
            oldJtot:=Jtot; {store the old value of the current}
            Jtot:=(Jn[ROUND(NP/2)]+Jp[ROUND(NP/2)]); {and compute the new value}
            IF Rshunt>0 THEN Jtot:=Jtot + Va/Rshunt; {add leakage current}
            IF (Conv_Var = 1) AND (ABS(Jtot) > 1e-4) {check for convergence of current density}
                THEN Conv_Main:=ABS( (Jtot-oldJtot)/Jtot ) < TolMain
                ELSE Conv_Main:=delDens/(totn + totp) < TolMain;

            MainIt:=MainIt+1
        END;
		
		IF JVSim[VCount].Store THEN 
			WITH JVSim[VCount] DO
			BEGIN
				J:=Jtot;
				{now correct for Rshunt and Rseries to give us Vext and Jext:}
				IF Rshunt>0 THEN Jext:=Jtot + Va/Rshunt ELSE Jext:=Jtot; {add current through shunt to measured current (=J)}
				Vext:=Va + J * Rseries; {add voltage drop across series resistance, use measured current J, not Jtot}
				Use:=check_Poisson AND Conv_Main; {store whether convergerence was achieved}
				quit_Voc:=(Jext>0) AND until_Voc; {use the externally measured JV-curve to see if we're past Voc}
			END;
       
		IF check_Poisson AND Conv_Main THEN
		BEGIN   {Convergence!}
			Recombi(rec,V, n, p, empty, Lang, diss_prob, empty); {calc. rec. without SRH and generation}
			recLan:=Average(rec, h, 0, NP+1); {store average Langevin recombination}
			Recombi(rec,V, n, p, empty, empty, diss_prob, SRH); {calc. rec. without Langevin and generation}
			recSRH:=Average(rec, h, 0, NP+1); {store average SRH recombination}
			Av_Diss:=Average_Dissociation(diss_prob);
			Jbimo:=0; JSRH_bulk:=0; JSRH_LI:=0; JSRH_RI:=0; Jph:=0;
			Calc_Recombi_current(Jbimo, JSRH_bulk, JSRH_LI, JSRH_RI, Jph, Jn_l, Jn_r, Jp_l, Jp_r); {calc. net recombination}
			WRITELN('At V=', Va:4:3,' converged in',MainIt:4,' loop(s), tot. curr. ',Jtot:7:4);
            IF JVSim[VCount].Store THEN {now only store the point in the JV_file if needed}
				WRITELN(uitv, Va:8:6, '  ',Jtot:nd, '  ',Av_Diss:nd,' ',recLan:nd,'  ',recSRH:nd,' ',Jbimo:nd,' ',
						JSRH_bulk:nd,' ',JSRH_LI:nd,' ',JSRH_RI:nd,' ',Jph:nd,' ',Jn_l:nd,' ',Jp_l:nd,' ',Jn_r:nd,
						' ',Jp_r:nd,' ',JVSim[VCount].Vext:nd,' ',JVSim[VCount].Jext:nd);
			FLUSH(uitv);
		END
		{no convergence:}
		ELSE WRITELN('No convergence for V=',Va:5:3,' Rel. change densities:',
                    delDens/(totn + totp):4,' Conv. Poisson: ', Check_Poisson);

		VCount:=VCount + 1;
    END;

    IF (Gmax * Gfrac <> 0) AND (V0 <> VL) THEN  {we do have a solar cell}
		Calc_and_Output_Solar_Cell_Parameters(JVExp, JVSim);  

	IF UseExpData THEN Compare_Exp_Sim_JV(JVExp, JVSim);
	{note: if Rseries <> 0 then Vext in JVExp and JVSim will be different so we can't do a direct comparison}

    WRITELN('The JV characteristic is written in ', JV_file,'.');
    {Store variables (x, V, etc) for last applied voltage:}
	Write_Vars_To_File(V, n, p, Vgn, Vgp, n_trap, p_trap, rec, Lang, diss_prob, SRH, nid, pid, mun, mup, Gm, Jn, Jp, nion, pion, Va);
	WRITELN('Finished, press enter to exit');

	IF Pause_at_end = 1 THEN READLN; {pause at the end of the program}

	CLOSE(log);
END.
