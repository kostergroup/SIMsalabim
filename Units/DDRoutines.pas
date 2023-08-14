unit DDRoutines;
{provides drift-diffusion procedures and functions}

{
SIMsalabim:a 1D drift-diffusion simulator 
Copyright (c) 2021, 2022, 2023 S. Heester, Dr T.S. Sherkar, Dr V.M. Le Corre, 
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
Nijenborgh 4, 9747 AG Groningen, the Netherlands
}

{$MODE OBJFPC} {force OBJFPC mode}

INTERFACE

USES sysutils, 
	 typinfo,
	 DateUtils,
	 Math, 
	 TypesAndConstants,
	 TransferMatrix,
	 InputOutputUtils, 
     NumericalUtils,
     StrUtils,
     DDTypesAndConstants;

CONST DDRoutinesVersion = '4.56'; {version of this unit}

PROCEDURE Print_Welcome_Message(ProgName : TProgram; version : STRING);
{Prints a welcome message, lists the authors and shows the name and verion of the program.}

PROCEDURE Display_Help_Exit(ProgName : TProgram); 
{displays a short help message and exits}

FUNCTION Max_Value_myReal : EXTENDED;
{Determines the largest value that can be stored in myReal. The result, thus,
depends on how myReal was defined. It should be a single, double or extended.}

PROCEDURE Set_Number_Digits(LimOutput : BOOLEAN; BytesInReal : INTEGER);
{sets the number of digits when printing real numbers (in global var nd), 
depending on the size of the floating point type}

PROCEDURE Determine_Name_Parameter_File(VAR parameterFile : STRING);
{Determines which parameter file should be used}

FUNCTION Correct_Version_Parameter_File(ProgName : TProgram; parameterFile, version : STRING) : BOOLEAN;
{checks if the version and name of the program and the parameter file match. Stops program if not!}

PROCEDURE Prepare_Log_File(VAR log : TEXT; MsgStr : ANSISTRING; CONSTREF par : TInputParameters; version : STRING);
{opens a log_file for later use and writes date/time and MsgStr to the log file.}

PROCEDURE Finalize_Log_File(VAR log : TEXT; MsgStr : ANSISTRING);
{writes final comments, date, time, run time and closes log file.}

PROCEDURE Read_Parameters(parameterFile : STRING; VAR msg : ANSISTRING; VAR par : TInputParameters; VAR stv : TStaticVars; ProgName : TProgram);
{Reads-in all the parameters. Some bits are specific to either ZimT or SimSS}

PROCEDURE Check_Parameters(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; ProgName : TProgram);
{performs a number of checks on the parameters. Just to ensure that they are valid, consistent, make sense}
{Some bits are specific to either ZimT or SimSS}

PROCEDURE Make_Grid(VAR k, x : vector; VAR i1, i2 : INTEGER; CONSTREF par : TInputParameters);
{Makes an exponential symmetric grid, for every layer}
{k[i] = (x[i+1] - x[i])/L and initialises the array with x-positions}
{i1 is the last point in the left insulator (or 0 if there isn't any)
i2 is the first point in the right insulator (or NP+1 if there is none)}

PROCEDURE Define_Layers(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Note, stv are not CONSTREF as we need to change them}
{Sets layer dependent properties}

PROCEDURE Init_Generation_Profile(VAR stv : TStaticVars; VAR log : TEXT; CONSTREF par : TInputParameters);
{Inits the generation profile, either constant or from a file. This is the SHAPE of the profile}
{When using a profile from file, a message is written on screen and in the log file.}

PROCEDURE Update_Generation_Profile(org: vector; VAR new : vector; Geff : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Rescale profile such that Geff equals the average of the profile over the length of the absorber. 
The latter equals L if TLsGen or there are no transport layers, and L-L_LTL-L_RTL otherwise. }

PROCEDURE Init_Pot_Dens_Ions_Traps(VAR V, Vgn, Vgp, n, p, nion, pion : vector; VAR f_tb, f_ti : TrapArray; Va : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{init. for V, Vgn,p, n, p, ion densities and trap parameters at bias voltage Va}

PROCEDURE Init_Trap_Distribution(VAR log : TEXT; VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Places all types of traps (bulk and interface) in the device at places deterined by define_layers.}

PROCEDURE Init_nt0_And_pt0(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{inits nt0 and pt0 arrays needed for SRH recombination}
{note: stv are changed here (nt0 and pt0), so they are VAR parameters}

PROCEDURE Extrapolate_Solution(CONSTREF prev, curr : TState; VAR new : TState; AccSols : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Obtains a guess for state new based on an extrapolation of states prev and curr}

PROCEDURE Main_Solver(VAR curr, new : TState; VAR it : INTEGER; VAR conv : BOOLEAN; VAR StatusStr : ANSISTRING; 
					  CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Iteratively solves the Poisson and continuity equations, including traps and ions}
{can be used in steady-state and transient cases}

PROCEDURE Prepare_tJV_File(VAR uitv : TEXT; filename : STRING; transient : BOOLEAN); 
{create a new tJV_file with appropriate heading
after running this, the TEXT file 'uitv' is still open and ready for writing}

PROCEDURE Write_To_tJV_File(VAR uitv : TEXT; CONSTREF CurrState, PrevState : Tstate; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{before running this proc, uitv must be open (by running Prepare_tJV_File). It must be closed in the main code.
This proc writes the (time), voltage, currents, recombination currents to a file that contains the JV-curve}

PROCEDURE Prepare_Var_File(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN); 
{create a new var_file with appropriate heading}

PROCEDURE Write_Variables_To_File(VAR CurrState : TState; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{writes the internal variables (astate) to file par.Var_file. It assumes the file has a header produced by Prepare_Var_File}

PROCEDURE Tidy_Up_Parameter_File(parameterFile : STRING; QuitWhenDone : BOOLEAN);
{This procedure cleans up the parameter file: it makes sure that all the * are aligned, that every parameter has a 
unit or other description/comment and left-aligns any line that starts with **.}

IMPLEMENTATION

VAR RecDum : TRec; {global dummy variable that will be used in Calc_All_Currents, Calc_Recombination_n, and Calc_Recombination_p}
{using a local variable in those functions will cause problems (unpredictable behaviour) as this var is pretty large.}

VAR TimeStart : TDateTime; {stores the date/time at the start of the sim. Used in Prepare_Log_File and Finalize_Log_File}

VAR nd : INTEGER; {used to specify the field_width of floating point type}

PROCEDURE Print_Welcome_Message(ProgName : TProgram; version : STRING);
{Prints a welcome message, lists the authors and shows the name and verion of the program.}
VAR strprogname : STRING;
BEGIN
    Str(ProgName, strprogname); {convert variable ProgName to a string}
    WRITELN('Welcome to ',strprogname,' version ',version,'.');
    WRITELN('Copyright (C) 2020, 2021, 2022, 2023 S. Heester, Dr T.S. Sherkar,'); 
    WRITELN('Dr V.M. Le Corre, Dr M. Koopmans, F. Wobben,');
    WRITELN('and Prof L.J.A. Koster, University of Groningen.');
    WRITELN;
END;

PROCEDURE Display_Help_Exit(ProgName : TProgram); 
{displays a short help message and exits}
VAR strprogname : STRING;
BEGIN
    Str(ProgName, strprogname); {convert variable ProgName to a string}
    strprogname:=LOWERCASE(strprogname); {ensure it's lower case}
    WRITELN('This program can be used with the following options:');
    WRITELN;
    WRITELN('All parameters can be defined in a parameter file');
    WRITELN('or via the command line, which overrides the values in the file.');
    WRITELN('Example: ./',strprogname,' -T 400 -Var_file Var.dat');
    WRITELN;
    WRITELN('To tidy-up the parameter file, use ''-tidy''');
    WRITELN;
    Stop_Prog('''-h''     : displays this help message', EC_Warning);
END;

FUNCTION B(x : myReal) : myReal; {The Bernoulli function, an approximation}
{we use several approximations to the Bernoulli function:
 close to 0: 1st order, then 2nd order, then 4th order. Finally, for large (abs.) x, we use the full expression}
{Taylor series around 0 is: B(x) = 1 - x/2 + x^2/12 -x^4/720 + O(x^6)}
CONST A1 : myReal = 1/12;
	  A2 : myReal = 1/720; 
	  C1 = 0.05; {if |x| < C1, then use 1st order Taylor}
	  C2 = 0.45; {if C1 <= |x| < C2, then use 2nd order Taylor}
	  C3 = 1.5; {if C2 <= |x| < C3, then use 4th order Taylor}
	  {if |x| is even larger, then we use the full expression}
VAR absx : myReal;
BEGIN
    absx:=ABS(x);
    IF absx < C1 {if x close to zero (=most common case!)}
		THEN B:=1 - 0.5*x {then use 1st order Taylor}
		ELSE
		BEGIN {absx >= C1}
			IF absx < C2 THEN B:=1 - x*(0.5 - A1 * x) {use 2nd order}
			ELSE
			BEGIN {absx >= C2}
				IF absx < C3 THEN B:=1 + x*(-0.5 + x*(A1 - x*x*A2)) {use 4th order Taylor}
				ELSE B:=x/(EXP(x)-1); {x > C3, will be a very small value}
			END {absx >= C2}
		END; {absx >= C1}
END;

FUNCTION Average(Vec, h : Vector; istart, ifinish : INTEGER; onGrid : BOOLEAN = FALSE) : myReal;
{calcs the average of a vector (Vec) from index istart -> ifinish. 
As the grid may be non-uniform, it also takes and uses grid spacing h.}
{On-grid: like densities, V, etc. Not on-grid: currents, mobility}
VAR i   : INTEGER;
    sumy, sumh, val : myReal;
BEGIN
    sumy:=0; {sum of the y-values}
    sumh:=0; {sum of grid spacing}
    
	FOR i:=istart TO ifinish-1 DO BEGIN
		IF onGrid THEN
			val:=0.5*(Vec[i+1] + Vec[i]) {we use trapezoidal integration, so interpolate y-values (Vec)}
		ELSE
			val:=Vec[i]; {quantity (Vec) is defined between grid points, so Vec[i] is the value between x[i] and x[i+1]}
		sumy:=sumy + h[i]*val;
		sumh:=sumh + h[i]
	END;
   
    Average := sumy/sumh;
END;

FUNCTION Max_Value_myReal : EXTENDED;
{Determines the largest value that can be stored in myReal. The result, thus,
depends on how myReal was defined. It should be a single, doulbe or extended.}
VAR dummy : EXTENDED;
BEGIN
	dummy:=0;

	{we assume that singles and doubles are available}
	IF TypeInfo(myReal) = TypeInfo(single) THEN dummy:=MaxSingle;
	IF TypeInfo(myReal) = TypeInfo(double) THEN dummy:=MaxDouble;

	{however, extended reals might not be available, so we check:}
{$IFDEF FPC_HAS_TYPE_EXTENDED}
	IF TypeInfo(myReal) = TypeInfo(extended) THEN dummy:=MaxExtended;
{$ENDIF}

	{at this point: if dummy is still zero, then myReal was not the right type!}

{$IFDEF FPC_HAS_TYPE_EXTENDED}
	IF dummy=0 THEN Stop_Prog('Error in function Max_Value_myReal. myReal should be single, double, or extended.', EC_ProgrammingError);
{$ELSE}
	{we're assuming that we do have singles and doubles}
	IF dummy=0 THEN Stop_Prog('Error in function Max_Value_myReal. myReal should be single or double.', EC_ProgrammingError);
{$ENDIF}	
	
	{OK, now we can be sure that myReal is of the right type}
	Max_Value_myReal:=dummy;
END;

PROCEDURE Set_Number_Digits(LimOutput : BOOLEAN; BytesInReal : INTEGER);
{sets the number of digits when printing real numbers (in global var nd), depending on the size of the floating point type}
BEGIN	

	IF LimOutput THEN {limit output, depending on size of floating point type}
		nd:=ROUND(0.333*(4*BytesInReal+20))
	ELSE
		{we use the full floating point type, which depends on the size of the type}
		CASE BytesInReal OF
			4 : nd:=16;
			8 : nd:=24;
			10 : nd:=29;
		OTHERWISE 
			Stop_Prog('Invalid parameter BytesInReal passed to procedure Set_Number_Digits', EC_ProgrammingError)
		END
END;

PROCEDURE Determine_Name_Parameter_File(VAR parameterFile : STRING);
{Determines which parameter file should be used}
VAR key : STRING;
BEGIN
	{copy first parameter (after prog name!) from command line, this works even if there aren't any!}
	key:=TRIM(ParamStr(1));

	{now take default file if the first character is a '-' (=> variable name is coming!) or there parameters at all!}
	IF StartsStr('-', key) OR (ParamCount=0) THEN
		parameterFile:=defaultParameterFile
	ELSE {only in this case do we take the name of the parameter file from the command line}
		parameterFile:=key
END;

FUNCTION Correct_Version_Parameter_File(ProgName : TProgram; parameterFile, version : STRING) : BOOLEAN;
{checks if the version and name of the program and the parameter file match. Stops program if not!}
{check for programe name: we take the variable ProgName and convert this to a lower-case string strprogname. This
will be either 'zimt' or 'simsalabim'. Then we check if either 'zimt' or 'simsalabim' appears (at least once) in the 
parameter file.}
{check for version info:
we do this by checking if there is a line that contains both the
string 'version' and the string that contains the version number of the program (version).}
{Note: neither version- nor progname checking is very strict!}
VAR inp : TEXT;
    found_version, found_progname : BOOLEAN;
    line, strprogname : STRING;
BEGIN
    IF NOT FileExists(parameterFile) {the file with input par. is not found}
        THEN Stop_Prog('Could not find file '+parameterFile+'.', EC_FileNotFound);
    ASSIGN(inp, parameterFile);
    RESET(inp);
    
    Str(ProgName, strprogname); {convert variable ProgName to a string}
    strprogname:=LOWERCASE(strprogname); {force lower case}
    found_version:=FALSE;
    found_progname:=FALSE;
    
    WHILE NOT EOF(inp) DO
    BEGIN
        READLN(inp, line); {read a line from the file}
        line:=LOWERCASE(line); {convert all characters to lower case to make search easier}
        {first: did we find the version number:}
        IF (POS('version', line) > 0) AND (POS(LOWERCASE(version), line) > 0) THEN
            found_version:=TRUE; {we have found the version number and it is correct!}
        {did we find the correct name of the program:}
        IF POS(strprogname, line) > 0 THEN found_progname:=TRUE;
        {POS returns the index of Substr in S, if S contains Substr. In case Substr isn't found, 0 is returned.}
    END;
    CLOSE(inp);
    Correct_Version_Parameter_File:=found_version AND found_progname;
END;

PROCEDURE Prepare_Log_File(VAR log : TEXT; MsgStr : ANSISTRING; CONSTREF par : TInputParameters; version : STRING);
{opens a log_file for later use and writes date/time and MsgStr to the log file.}
BEGIN
    ASSIGN(log, par.log_file); 
    REWRITE(log);
    WRITELN(log,'Version ', version);
    WRITELN(log,'Size of reals used in simulation: ',SizeOf(myReal),' bytes');
    TimeStart:=NOW; {this is a global variable that will also be used in Finalize_Log_File}
    WRITELN(log,'Timestamp at start: ',DateTimeToStr(TimeStart));
    
    {Read_Parameters may have added something to MsgStr (values from the command line):}
    IF Length(MsgStr) > 0 THEN
    BEGIN
		WRITELN(log, 'Values from command line:');
		WRITE(log, MsgStr) {MsgStr contains any messages from Read_Parameter. No need for writeln as it either empty, or end with LineEnding}
    END;
    FLUSH(log);
END;

PROCEDURE Finalize_Log_File(VAR log : TEXT; MsgStr : ANSISTRING);
{writes final comments, date, time, run time and closes log file.}
VAR TimeEnd : TDateTime;
BEGIN
	{first write any messages:}
	IF Length(MsgStr) > 0 THEN WRITE(log, MsgStr); 
	
	{note: TimeStart is a global variable in this unit and gets its value in Prepare_Log_File}
	TimeEnd:=NOW;
	WRITELN(log,'Timestamp at end: ',DateTimeToStr(TimeEnd));
	WRITELN(log,'Total run time: ',SecondSpan(TimeEnd, TimeStart):6:3,' seconds.');
	CLOSE(log)
END;

PROCEDURE Read_Parameters(parameterFile : STRING; VAR msg : ANSISTRING; VAR par : TInputParameters; VAR stv : TStaticVars; ProgName : TProgram);
{Reads-in all the parameters. Some bits are specific to either ZimT or SimSS}
VAR countStart, countFinish : INTEGER;
    dumint : INTEGER; {a dummy integer variable}
    dumstr : STRING; {again, a dummy variable}
    inv : TEXT;
    ZimT, SimSS, UsedSpecialParFile	: BOOLEAN; 
BEGIN
    {use 2 booleans to check if we're using ZimT or SimSS}
    ZimT:= (ProgName=TProgram.ZimT);
    SimSS:= (ProgName=TProgram.SimSS);
    
    IF NOT FileExists(parameterFile) {the file with input par. is not found}
		THEN Stop_Prog('Could not find file '+parameterFile+'.', EC_FileNotFound);
    ASSIGN(inv, parameterFile);
    RESET(inv);
	countStart:=Count_Substring_In_String(LineEnding, msg); {first, count the number of LineEndings in msg-str.}

	WITH par DO BEGIN
{**General**************************************************************************}
		Get_Float(inv, msg,'T',T);  {abs. temperature, K}
		stv.Vt:=k*T/q;  {thermal voltage}
		stv.Vti:=1/stv.Vt; {inverse of Vt, we'll use this a lot!}
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
	
		stv.V0:=0.5 * (VB+CB) - W_L;
		stv.VL:=0.5 * (VB+CB) - W_R;

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
		IF SimSS THEN Get_Integer(inv, msg, 'mob_ion_spec', dumint);{mobile ion species: -1: negative, 0: both, 1: positive ions}
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
		stv.Traps:=(Bulk_tr<>0) OR (St_L<>0) OR (St_R<>0) OR ((GB_tr<>0) AND (num_GBs>0)); {are there any traps?}
		stv.Traps_int:=(St_L<>0) OR (St_R<>0) OR ((GB_tr<>0) AND (num_GBs>0)); {are there any interface traps?}
		stv.Traps_int_poisson:=((St_L<>0) AND (Tr_type_L<>0)) OR ((St_R<>0) AND (Tr_type_R<>0)) OR ((GB_tr<>0) AND (num_GBs>0) AND (Tr_type_B<>0));

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
			Get_Integer(inv, msg, 'NJV', stv.NJV); {Number of JV points, for logarithmic JV}
			{Note: NJV is TStaticVars as we sometimes need to calculate it, so it's not a direc input parameter}
			IF (Vstep<>0) AND (Vdistribution=1) THEN
				stv.NJV:=TRUNC((Vmax - Vmin)/Vstep + 1e-10) + 1; {needed for setting length of Jdat and Vdat}
			{1e-10 is needed to get right value}
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
    {now check if all parameters that were in the command line have been used:}
	{did we use a special parameter file or not?}
	UsedSpecialParFile:=NOT( (ParamCount=0) OR StartsStr('-', TRIM(ParamStr(1))) );

    {we use countStart/Finish to count the number of parameters obtained from the command line}
	{every time we find such a parameter, a LineEnding is added to msg}
    countFinish:=Count_Substring_In_String(LineEnding, msg);
	
	{now the number of arguments (bar a dev par file) needs to be 2 x the difference in the counter:}
    IF ParamCount - ORD(UsedSpecialParFile) <> 2*(countFinish-countStart) THEN
		Stop_Prog('The command line contains invalid arguments, see the manual.', EC_InvalidCLInput);
    msg:=msg + 'Read parameters from ' + parameterFile + LineEnding
END;

PROCEDURE Check_Parameters(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; ProgName : TProgram);
{performs a number of checks on the parameters. Just to ensure that they are valid, consistent, make sense}
{Some bits are specific to either ZimT or SimSS}
VAR ZimT, SimSS	: BOOLEAN;
BEGIN
    {use 2 booleans to check if we're using ZimT or SimSS}
    ZimT:= (ProgName=TProgram.ZimT);
    SimSS:= (ProgName=TProgram.SimSS);
  
    {when adding new check, please keep the order in line with the device_parameter file}
    {Check first if stv.Vt has been initialised. This should have happened in Read_Parameters}
    IF (stv.Vt=0) THEN Stop_Prog('stv.Vt needs to be initialised before calling Check_Parameters.', EC_ProgrammingError);

	WITH par DO BEGIN
{checks on general parameters:}
		IF CB >= VB THEN Stop_Prog('CB should be smaller than VB.', EC_InvalidInput);
		IF (CB<0) OR (VB<0) THEN Stop_Prog('CB and VB should be positive.', EC_InvalidInput);

{checks on mobilities:}
		IF NOT (mob_n_dep IN [0, 1]) THEN Stop_Prog('Invalid mob_dep_n selected.', EC_InvalidInput);
		IF NOT (mob_p_dep IN [0, 1]) THEN Stop_Prog('Invalid mob_dep_p selected.', EC_InvalidInput);
{checks on contacts:}
		{part of this can only be done after checking the TLs!}
		IF Rseries<0 THEN Stop_Prog('Rseries cannot be negative.', EC_InvalidInput);
		IF SimSS AND (Rshunt=0) THEN Stop_Prog('Rshunt cannot be zero, use positive (negative) value for finite (infinite) shunt resistance.', EC_InvalidInput);
{checks on transport layers:}
		IF L_LTL+L_RTL>=L THEN Stop_Prog('Sum of L_LTL and L_RTL (transport layers) should be less than total thickness L.', EC_InvalidInput);
		{more checks on the energy levels of the TLs:}
		{first check left:}
		IF L_LTL>0 THEN
		BEGIN
			IF W_L<CB_LTL THEN Stop_Prog('W_L cannot be smaller than CB_LTL.', EC_InvalidInput);
			IF W_L>VB_LTL THEN Stop_Prog('W_L cannot be larger than VB_LTL.', EC_InvalidInput);
			IF CB_LTL>= VB_LTL THEN Stop_Prog('CB_LTL must be smaller than VB_LTL', EC_InvalidInput);
		END
		ELSE BEGIN {no TLs:}
			IF W_L<CB THEN Stop_Prog('W_L cannot be smaller than CB.', EC_InvalidInput);
			IF W_L>VB THEN Stop_Prog('W_L cannot be larger than VB.', EC_InvalidInput);
		END;
		{check the right:}
		IF L_RTL > 0 THEN
		BEGIN
			IF W_R<CB_RTL THEN Stop_Prog('W_R cannot be smaller than CB_RTL.', EC_InvalidInput);
			IF W_R>VB_RTL THEN Stop_Prog('W_R cannot be larger than VB_RTL.', EC_InvalidInput);
			IF CB_RTL>= VB_RTL THEN Stop_Prog('CB_RTL must be smaller than VB_RTL', EC_InvalidInput);
		END
		ELSE
		BEGIN
			IF W_R<CB THEN Stop_Prog('W_R cannot be smaller than CB.', EC_InvalidInput);
			IF W_R>VB THEN Stop_Prog('W_R cannot be larger than VB.', EC_InvalidInput);
		END;

{checks on traps}
		{check whether there are a possible number of traps (negative not allowed)}
		IF (St_L < 0) OR (St_R < 0) OR (Bulk_tr < 0) OR (GB_tr < 0) THEN Stop_Prog('Negative trap density not allowed.', EC_InvalidInput);
		IF (St_L > 0) AND (L_LTL=0) THEN Stop_Prog('You cannot have interface traps (St_L>0) without a TL (L_LTL=0).', EC_InvalidInput);
		IF (St_R > 0) AND (L_RTL=0) THEN Stop_Prog('You cannot have interface traps (St_R>0) without a TL (L_RTL=0).', EC_InvalidInput);
		IF num_GBs<0 THEN Stop_Prog('The number of grain boundaries (num_GBs) cannot be negative.', EC_InvalidInput);
		IF (num_GBs>0) AND (GB_tr=0) THEN Stop_Prog('Trap density at grain boundaries (GB_tr) must be > 0 if num_GBs>0', EC_InvalidInput);
		
		{Only Cn OR Cp = 0 are allowed, if both are zero, no charge can reach the traps, this makes no sense.}
		IF (Cn = 0) AND (Cp = 0) THEN Stop_Prog('Cn and Cp cannot both be zero, change parameters please.', EC_InvalidInput);	

		{trap types: -1, 0 or 1. Sets in pascal are an ordinal type with a range between 0 and 255, hence the ABS:}
		IF NOT (ABS(Tr_type_L) IN [0, 1]) THEN Stop_Prog('Invalid left interface trap type.', EC_InvalidInput);
		IF NOT (ABS(Tr_type_R) IN [0, 1]) THEN Stop_Prog('Invalid right interface trap type.', EC_InvalidInput);
		IF NOT (ABS(Tr_type_B) IN [0, 1]) THEN Stop_Prog('Invalid bulk trap type.', EC_InvalidInput);        

		{check on the energies of the traps are performed in proc Init_Trap_Distribution}
		
{checks on ions:}
		IF SimSS AND (ion_red_rate<0) THEN Stop_Prog('Ion redistribution rate cannot be lower than zero, should be: 0 <= ion_red_rate < NJV', EC_InvalidInput);
{checks on optics, generation and recombination parameters}
		IF Use_gen_profile = 1 THEN BEGIN
			IF L_TCO < 0 THEN Stop_Prog('L_TCO cannot be negative.', EC_InvalidInput);
			IF L_BE <= 0 THEN Stop_Prog('L_BE must be positive.', EC_InvalidInput);
			IF lambda_min > lambda_max THEN Stop_Prog('lambda_min cannot be larger than lambda_max.', EC_InvalidInput);
			IF lambda_min <= 0 THEN Stop_Prog('lambda_min must be positive.', EC_InvalidInput); 
			{note, now we are sure that lambda_max is also positive!}		
		END;
		IF NOT (ThermLengDist IN [1,2,3,4,5]) THEN Stop_Prog('Invalid ThermLengDist selected.', EC_InvalidInput);
		IF (P0>=1) OR (P0<0) THEN Stop_Prog('Invalid value of P0, should be: 0<=P0<1', EC_InvalidInput);
		IF (P0<>0) AND (Field_dep_G = FALSE) THEN Stop_Prog('P0 should be zero if not using field dependent generation', EC_InvalidInput);
{checks on numerical parameters:}
		IF (NP<=15) OR (NP>Max_NP) THEN Stop_Prog('Invalid number of grid points (NP) selected, must be >=15 and <'+IntToStr(Max_NP)+'.', EC_InvalidInput);
		IF (CurrDiffInt <> 1) AND (CurrDiffInt <> 2) THEN Stop_Prog('CurrDiffInt can only be 1 or 2.', EC_InvalidInput);
		IF maxDelV<=0 THEN Stop_Prog('maxDelV should be positive.', EC_InvalidInput);
		IF maxDelV*stv.Vt <= tolPois THEN Stop_Prog('maxDelV*Vt should be (much) larger than tolPois.', EC_InvalidInput);
		IF tolDens <= 0 THEN Stop_Prog('tolDens must be larger than zero.', EC_InvalidInput);
		IF couplePC < 0 THEN Stop_Prog('couplePC must be non-negative.', EC_InvalidInput);
		{check if values of minAcc and maxAcc makes any sense:}
		IF maxAcc >= 2 THEN Stop_Prog('maxAcc must be smaller than 2.', EC_InvalidInput);  
		IF minAcc <= 0 THEN Stop_Prog('minAcc must be positive.', EC_InvalidInput);  
		IF minAcc > maxAcc THEN Stop_Prog('minAcc cannot be larger than maxAcc.', EC_InvalidInput);  
		IF NOT (FailureMode IN [0,1,2]) THEN Stop_Prog('Invalid FailureMode selected.', EC_InvalidInput);
{checks on voltages, SimSS only:}
		IF SimSS THEN
		BEGIN
			IF NOT (Vdistribution IN [1,2]) THEN Stop_Prog('Invalid voltage distribution selected.', EC_InvalidInput);
			IF ABS(Vscan) <> 1 THEN Stop_Prog('Vscan must be either -1 or 1.', EC_InvalidInput);
			{check if Vmin and Vmax are not too small or large:}
			IF Vmin*stv.Vti < -1.95 * LN(Max_Value_myReal) THEN Stop_Prog('Vmin is too small.', EC_InvalidInput);
			IF Vmax*stv.Vti > 1.95 * LN(Max_Value_myReal) THEN Stop_Prog('Vmax is too large.', EC_InvalidInput);
			IF Vmin > Vmax THEN Stop_Prog('Vmin should not be greater than Vmax.', EC_InvalidInput);
			{now check for redundancy of pre-bias:}
			IF PreCond THEN
			BEGIN
				IF Rseries>0 THEN Warn_User('Pre-bias voltage does not take Rseries into account, so Vpre=Vint.');
				IF ABS(Vpre)*stv.Vti > 1.95 * LN(Max_Value_myReal) THEN Stop_Prog('|Vpre| is too large.', EC_InvalidInput);
				IF ((CNI=0) AND (CPI=0)) THEN Stop_Prog('Do not use a pre-bias without any ions, makes no sense.', EC_InvalidInput);
				IF (Vscan=1) AND (Vpre=Vmin) THEN Stop_Prog('Pre-bias voltage is equal to Vmin, makes no sense.', EC_InvalidInput);
				IF (Vscan=-1) AND (Vpre=Vmax) THEN Stop_Prog('Pre-bias voltage is equal to Vmax, makes no sense.', EC_InvalidInput);
			END;
			IF ((CNI<>0) OR (CPI<>0))AND (NOT (ion_red_rate in [0,1])) AND (Vdistribution =2) THEN 
				Stop_Prog('Do not use Vdistribution=2 with ion_red_rate other than 0 or 1.', EC_InvalidInput);
			IF (Vacc >= Vmin) AND (Vacc <= Vmax) AND (Vdistribution = 2) THEN {Vacc is not valid} 
				Stop_Prog('Invalid Vacc selected, must be outside [Vmin, Vmax].', EC_InvalidInput);
			IF (Vdistribution=1) AND (Vstep <= 0) THEN Stop_Prog('Vstep should be positive.', EC_InvalidInput);	
			IF (ABS(Vmin-Vmax) < 1e-10) AND (Vdistribution=2) {to avoid infinite loop of Va}
			THEN Stop_Prog('Do not use Vdistribution=2 when Vmin = Vmax.', EC_InvalidInput);	
		END;

{checks on user-interface:}
		IF timeout = 0 THEN Stop_Prog('Invalid timeout: either positive number in seconds or negative for unlimited run time.', EC_InvalidInput);
		IF SimSS 
		THEN BEGIN
			IF (Gehp * Gfrac <> 0) AND (stv.V0 <> stv.VL) AND UseExpData AND (rms_mode=logarithmic) {this is a weird combination, warn user}
				THEN Warn_User('You are fitting a solar cell with rms_mode=log.');
			IF UseExpData AND until_Voc THEN Stop_Prog('You cannot use until_Voc = 1 and UseExpData = 1 at the same time.', EC_InvalidInput);
			IF UseExpData AND PreCond THEN Stop_Prog('You cannot use pre-conditioning (PreCond) and UseExpData = 1 at the same time.', EC_InvalidInput);
			IF PreCond AND (ion_red_rate=1) THEN Stop_Prog('You cannot use PreCond and have ion_red_rate=1 at the same time as that defeats the purpose.', EC_InvalidInput);
			IF ((rms_threshold<=0) OR (rms_threshold>1)) AND UseExpData THEN
				Stop_Prog('rms_threshold has to be larger than 0 but not larger than 1.', EC_InvalidInput);
		END;		
		IF SimSS AND (OutputRatio < 0) THEN Stop_Prog('OutputRatio should be 0 (no output) or positive.', EC_InvalidInput); {if zero, then simply no var file output}
		IF ZimT AND (OutputRatio <= 0) THEN Stop_Prog('OutputRatio should be positive.', EC_InvalidInput); {In ZimT it cannot be zero as we NEED to write output as it also limits the output to screen.}
    END
END;

PROCEDURE Make_Sub_Grid(VAR k : vector; istart, ifinish : INTEGER; grad, xstart, xfinish, L : myReal);
{this makes a grid on part of the volume}
VAR i, ending : INTEGER;
    norm : myReal;
BEGIN
    FOR i:=istart TO ifinish DO k[i]:=1;
    {note: we assign a value to k[np+1], but it doesn't have any meaning!}
    ending:=ROUND(0.25*(ifinish-istart)); {defines the exponential part of the grid}
    FOR i:=0 TO ending DO
    BEGIN
        k[i+istart]:=EXP(grad*(1/ending - 1/(i+1)));
        k[ifinish-i]:=k[i+istart]
    END;
    norm:=0;
    FOR i:=istart TO ifinish DO norm:=norm + k[i];
    norm:=(xfinish-xstart)/(L*norm);
    FOR i:=istart TO ifinish DO k[i]:=k[i]*norm;
END;

PROCEDURE Make_Grid(VAR k, x : vector; VAR i1, i2 : INTEGER; CONSTREF par : TInputParameters);
{Makes an exponential symmetric grid, for every layer}
{k[i] = (x[i+1] - x[i])/L and initialises the array with x-positions}
{i1 is the last point in the left insulator (or 0 if there isn't any)
i2 is the first point in the right insulator (or NP+1 if there is none)}
VAR del : myReal;
	i : INTEGER;
BEGIN
   
	WITH par DO 
	BEGIN
		{number of grid point in a layer: at least 5, proportional to thickness of layer:}
		IF L_LTL=0 THEN i1:=0 ELSE i1:=MAX(5, ROUND(NP*L_LTL/L));
		IF L_RTL=0 THEN i2:=NP+1 ELSE i2:=NP+1 - MAX(5, ROUND(NP*L_RTL/L));

		{first make sub-grids in TLs:}
		IF L_LTL>0 THEN Make_Sub_Grid(k, 0, i1, grad, 0, L_LTL, L);
		IF L_RTL>0 THEN Make_Sub_Grid(k, i2, NP+1, grad, L-L_RTL, L, L);
		{now main layer:}
		Make_Sub_Grid(k, i1, i2, grad, L_LTL, L-L_RTL, L);
    
		{and compute the x-positions:}
		x[0]:=0;    
		FOR i:=1 TO NP+1 DO x[i]:=x[i-1] + L*k[i-1];
    
		{Now we revisite the points at the interfaces, if any, and make their spacing uniform}
		{this helps with convergence if there are traps and we're using De Mari's form of the current}
		IF L_LTL>0 THEN 
		BEGIN
			{next, make sure that points around interfaces have the same spacing:}
			del:=(x[i1+2]-x[i1-1])/3;
			k[i1-1]:=del/L; k[i1]:=del/L; k[i1+1]:=del/L;
			x[i1]:=x[i1-1]+del;
			x[i1+1]:=x[i1] + del;
		END;
		
		IF L_RTL>0 THEN 
		BEGIN
			del:=(x[i2+1]-x[i2-2])/3;
			k[i2-2]:=del/L; k[i2-1]:=del/L; k[i2]:=del/L;
			x[i2-1]:=x[i2-2]+del;
			x[i2]:=x[i2-1]+del;
		END;
		
		{note: i1: last point in LTL, i2: first point in RTL}
		{however, x[i1] is a little bit off from L_LTL, likewise for x[i2]}
		{because the interface sits between grid points}

	END; {WITH par statement}
END;

PROCEDURE Define_Layers(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Note, stv are not CONSTREF as we need to change them}
{Sets layer dependent properties}
VAR i : INTEGER;
BEGIN
	{first define layers:}
    WITH stv DO BEGIN
		{define different properties of the layers:}
 
		{first define bulk}
		FOR i:=0 TO par.NP+1 DO 
		BEGIN
			NcLoc[i]:=par.Nc; {bulk value of Nc}
			E_CB[i]:=par.CB;
			E_VB[i]:=par.VB;
			ni[i]:=par.Nc*EXP(-0.5*Vti*(E_VB[i]-E_CB[i])); {equilibrium concentration}
			nid[i]:=par.n_0;
			pid[i]:=par.p_0;
			eps[i]:=par.eps_r * eps_0;
		END;    
      
		{now the transport layers}
 
		{left TL:}
		IF i1>0 THEN
			FOR i:=0 TO i1 DO
			BEGIN
				NcLoc[i]:=par.Nc_LTL;
				E_CB[i]:=par.CB_LTL;
				E_VB[i]:=par.VB_LTL;
				ni[i]:=par.Nc*EXP(-0.5*Vti*(E_VB[i]-E_CB[i])); {equilibrium concentration}				
				nid[i]:=0; {first set nid and pid to zero and then make them non-zero if required}
				pid[i]:=0;
				IF par.doping_LTL<=0 THEN nid[i]:=-par.doping_LTL ELSE pid[i]:=par.doping_LTL;
				eps[i]:=par.eps_r_LTL * eps_0;
			END;
	
		{right TL:}	
		IF i2<par.NP+1 THEN
			FOR i:=i2 TO par.NP+1 DO
			BEGIN
				NcLoc[i]:=par.Nc_RTL;
				E_CB[i]:=par.CB_RTL;
				E_VB[i]:=par.VB_RTL;
				ni[i]:=par.Nc*EXP(-0.5*Vti*(E_VB[i]-E_CB[i])); {equilibrium concentration}				
				nid[i]:=0; {first set nid and pid to zero and then make them non-zero if required}
				pid[i]:=0;
				IF par.doping_RTL<=0 THEN nid[i]:=-par.doping_RTL ELSE pid[i]:=par.doping_RTL;
				eps[i]:=par.eps_r_RTL * eps_0;
			END;	
	  
	END {WITH stv statement}
END;

PROCEDURE Init_Generation_Profile(VAR stv : TStaticVars; VAR log : TEXT; CONSTREF par : TInputParameters);
{Inits the generation profile, either constant, calculated by the transfer matrix unit, or from a file. 
This is the SHAPE of the profile}
{When using a profile from file, a message is written on screen and in the log file.}
VAR a, gr : Row; {a : x-coordinate, gr: generation rate}
    numLines, i, j : INTEGER;
    maxG : myReal;
BEGIN
	CASE par.Use_gen_profile OF
	0:	BEGIN {uniform generation}
			FOR i:=0 TO par.NP+1 DO stv.orgGm[i]:=1; {constant Gehp=generation rate of e-h pairs}
			maxG:=1
		END;
	1:	Calc_TransferMatrix(stv,par); {Calculate the generation profile using the TransferMatrix script.}
	2:	BEGIN {use profile supplied by user}
			WRITELN('Reading generation profile from ',par.Gen_profile);
			WRITELN(log, 'Reading generation profile from ',par.Gen_profile);
			Read_XY_Table(a, gr, par.Gen_profile, 'x Gehp', numLines);
			{a: x-coordinate, gr: corresponding generation rate}

            IF numLines=0 THEN Stop_Prog('The file generation_profile.txt is empty.', EC_InvalidInput);

           { rescale a to ensure that a[numLines-1]=L: }
            FOR i:=0 TO numLines-1 DO a[i]:=a[i] * par.L/a[numLines-1];
			{the input file contains the generation profile, however, the grid points
             in this file may not correspond to our grid here}

            {Now interpolate gr to get genrate on x[i] grid:}
            i:=0; {counter for x[i] grid}
            maxG:=0; {max. of genrate}
            FOR j:=0 TO numLines-2 DO
                WHILE (stv.x[i] < a[j+1]) AND (i<par.NP+1) DO
                    BEGIN
						stv.orgGm[i]:=gr[j] + (stv.x[i]-a[j]) * (gr[j+1]-gr[j])/(a[j+1]-a[j]);
                        maxG:=Max(stv.orgGm[i], maxG); {look for maximum of genrate}
                        {we're using linear interpolation here}
                        i:=i+1
                    END;

            {now rescale genrate array so that the maximum is equal to 1:}
			FOR i:=0 TO par.NP+1 DO
				stv.orgGm[i]:=stv.orgGm[i]/maxG;
        END; {reading generation profile from file}
	END;
	
	{if the insulators (if there are any) don't absorb/generate we need set the generation
	 rate to zero. This overrides the generation profile (if any).}
	IF NOT par.TLsGen THEN
	BEGIN
		{left insulator, offset generalised potentials:}
		IF stv.i1>0 THEN { I would have thought that one shouldn't need this if statement}
			FOR i:=0 TO stv.i1 DO
				stv.orgGm[i]:=0;
		{right insulator, offset generalised potentials:}
		IF stv.i2<par.NP+1 THEN
			FOR i:=stv.i2 TO par.NP+1 DO
				stv.orgGm[i]:=0;
	END;
END;

PROCEDURE Update_Generation_Profile(org: vector; VAR new : vector; Geff : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Rescale profile such that Geff equals the average of the profile over the length of the absorber. 
The latter equals L if TLsGen or there are no transport layers, and L-L_LTL-L_RTL otherwise. }
VAR i : INTEGER;
	Gav : myReal;
BEGIN
	IF Geff<>0 
	THEN BEGIN
		IF par.Use_gen_profile <> 1 THEN {either uniform generation, or from user defined profile}
			BEGIN
			{Rescale such that Geff = Gehp * Gfrac equals the average of the profile over the
			length of the absorber. The latter equals L if TLsGen or there
			are no transport layers, and L-L_LTL-L_RTL otherwise. }
			Gav:=Average(org, stv.h, 0, par.NP+1, TRUE); {so average gen rate of original profile}
			IF NOT par.TLsGen THEN Gav:=Gav * par.L / (par.L-par.L_LTL-par.L_RTL);
			FOR i:=0 TO par.NP+1 DO
				new[i]:=org[i] * (Geff/Gav) {now set new generation rate to correct value:}
			END
		ELSE {calculated using transfer matrix, Geff = Gfrac:}
			FOR i:=0 TO par.NP+1 DO new[i]:=org[i]*Geff; {Generation rate already set in original profile}
		END
	ELSE FOR i:=0 TO par.NP+1 DO new[i]:=0;
END;

PROCEDURE Update_Gen_Pot(V : vector; VAR Vgn, Vgp : vector; CONSTREF stv : TstaticVars; CONSTREF par : TInputParameters);
{updates the generalised potentials (Vgn, Vgp) after potential V has been altered.}
VAR i : INTEGER;
	facDOS : myReal;
BEGIN
	Vgn:=V;
	Vgp:=V;
	{left insulator, offset generalised potentials:}
	IF stv.i1>0 THEN { I would have thought that one shouldn't need this if statement}
		FOR i:=0 TO stv.i1 DO
		BEGIN
			facDOS:=stv.Vt*LN(stv.NcLoc[i]/par.Nc);
			Vgn[i]:=V[i] + par.CB_LTL - par.CB + facDOS;
			Vgp[i]:=V[i] - par.VB + par.VB_LTL - facDOS;
		END;
	{right insulator, offset generalised potentials:}
	IF stv.i2<par.NP+1 THEN
		FOR i:=stv.i2 TO par.NP+1 DO
		BEGIN
			facDOS:=stv.Vt*LN(stv.NcLoc[i]/par.Nc);
			Vgn[i]:=V[i] + par.CB_RTL - par.CB + facDOS;
			Vgp[i]:=V[i] - par.VB + par.VB_RTL - facDOS;
		END;
END;

PROCEDURE Init_Pot_Dens_Ions_Traps(VAR V, Vgn, Vgp, n, p, nion, pion : vector; VAR f_tb, f_ti : TrapArray; Va : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{init. for V, Vgn,p, n, p, ion densities and trap parameters at bias voltage Va}
VAR i, istart, ifinish : INTEGER;
	facDOS : myReal;
BEGIN
    FOR i:=0 TO par.NP+1 DO {guess new V in interior of device} {linear interpolation of V}
	    V[i]:=stv.V0 - Va/2 + stv.x[i]*(stv.VL-stv.V0+Va)/par.L;
	Update_Gen_Pot(V, Vgn, Vgp, stv, par); {update generalised potentials}

    FOR i:=0 TO par.NP+1 DO {guess n, p in interior of device, not critical for n and p}
    BEGIN
		n[i]:=stv.ni[i]*EXP( (Va * (0.5 -stv.x[i]/par.L) + Vgn[i])*stv.Vti ); {note: we use Vgn,p, NOT V}
        p[i]:=stv.ni[i]*EXP( (Va * (stv.x[i]/par.L - 0.5) - Vgp[i])*stv.Vti );
		nion[i]:=0; {just to make sure these have been initialised}
		pion[i]:=0; {later on, we will set them to CNI/CPI}   
    END;
             
    {now we set the boundary conditions on the densities, these are the ones that matter:} 
	{left, i=0:}
	facDOS:=stv.Vt*LN(stv.NcLoc[0]/par.Nc);
	n[0]:=stv.NcLoc[0]*EXP((Vgn[0]-V[0]+par.CB-facDOS-par.W_L)*stv.Vti);
    p[0]:=SQR(stv.ni[0])/n[0]; {use law of mass action for p}

    {right, i=NP+1:}
    facDOS:=stv.Vt*LN(stv.NcLoc[par.NP+1]/par.Nc);
    p[par.NP+1]:=stv.NcLoc[par.NP+1]*EXP((par.W_R+V[par.NP+1]-Vgp[par.NP+1]-par.VB-facDOS)*stv.Vti);
    n[par.NP+1]:=SQR(stv.ni[par.NP+1])/p[par.NP+1];
   
	{the trap occupancies are first set to 
	 zero just to make sure they have a well defined value:}
	FILLCHAR(f_tb, SIZEOF(f_tb), 0);
	FILLCHAR(f_ti, SIZEOF(f_ti), 0);
    
    {Last but not least: the ions}
    {ions might be limited to the middle layer}
	IF (stv.i1>0) AND (NOT par.IonsInTLs) THEN istart:=stv.i1+1 ELSE istart:=0;
	IF (stv.i2<par.NP+1) AND (NOT par.IonsInTLs) THEN ifinish:=stv.i2-1 ELSE ifinish:=par.NP+1;
    FOR i:=istart TO ifinish DO
    BEGIN
		nion[i]:=par.CNI;
		pion[i]:=par.CPI;
    END;
    
END;

PROCEDURE Init_Trap_Distribution_Single_Level(VAR stv : TStaticVars; e : INTEGER; CONSTREF par : TInputParameters);
{Places all types of traps (bulk and interface) in the device at places determined by define_layers.}
VAR i, gb_nr, cur_idx, cwe_dummy : INTEGER;
    gb_centre_pos, inv_int_trap_width : myReal;
    new_gb_found : BOOLEAN;
BEGIN    
	{init the vectors}
	FOR i:=0 TO par.NP+1 DO
	BEGIN
		stv.Nti[i,e] := 0;
		stv.cwe_i[i] := 0;
	END;
	
	{we consider a trap filled when it contains an electron}
	{a donor type trap is charged when empty, and an acceptor trap is charged when filled.}
	CASE par.Tr_type_B OF
		-1 : stv.cwe_b := 0;
		0  : stv.cwe_b := 0; {we will not use this variable in this case, but it should have a value regardless}
		1  : stv.cwe_b := 1;
	END;

	FOR i:=0 TO par.NP+1 DO stv.Ntb[i,e] := par.Bulk_tr * stv.BulkTrapDist[e]; {first put bulk traps everywhere}
	{next: check if we need to remove traps in the TLs:}
	IF NOT par.TLsTrap THEN
	BEGIN {so the TLs (if any!) do not contain traps}
		{left transport layer}
		IF stv.i1 > 0 THEN
			FOR i := 0 TO stv.i1 DO
				stv.Ntb[i,e] := 0;
		{right transport layer}
		IF stv.i2 < par.NP+1 THEN
			FOR i:=stv.i2 TO par.NP+1 DO
				stv.Ntb[i,e] := 0;
	END;
	
	FOR i:=0 TO par.NP+1 DO stv.q_tr_igb[i]:= 1; {this array determines whether trap charge at interfaces and grain boundaries go in the Poisson equation (1: yes, 0:no).}
		
	IF (par.St_L > 0) THEN {if there is traps at the left interface}
	BEGIN
		{input St_L is #traps/area, will convert this to #traps/volume}
		inv_int_trap_width := 2 / (par.L * (0.5 * stv.h[stv.i1-1] + stv.h[stv.i1] + 0.5 * stv.h[stv.i1+1]));
		
		stv.Nti[stv.i1,e] := stv.IntTrapDist[e] * par.St_L * inv_int_trap_width; {#traps/volume}
		stv.Nti[stv.i1+1,e] := stv.IntTrapDist[e] * par.St_L * inv_int_trap_width; {#traps/volume}

		CASE par.Tr_type_L OF
			-1 : cwe_dummy := 0;
			0  : BEGIN cwe_dummy := 0; stv.q_tr_igb[stv.i1]:= 0; stv.q_tr_igb[stv.i1+1]:= 0 END; {we will not use this variable in this case, but it should have a value regardless}
			1  : cwe_dummy := 1;
		END;
		stv.cwe_i[stv.i1] := cwe_dummy;
		stv.cwe_i[stv.i1+1] := cwe_dummy;
	END;
	
	IF (par.St_R > 0) THEN {if there is traps at the right interface}
	BEGIN	
		{input St_R is #traps/area, will convert this to #traps/volume}
		inv_int_trap_width := 2 / (par.L * (0.5*stv.h[stv.i2] + stv.h[stv.i2-1] + 0.5*stv.h[stv.i2-2]));
		
		stv.Nti[stv.i2-1,e] := stv.IntTrapDist[e] * par.St_R * inv_int_trap_width; {#traps/volume}
		stv.Nti[stv.i2,e] := stv.IntTrapDist[e] * par.St_R * inv_int_trap_width; {#traps/volume}

		CASE par.Tr_type_R OF
		  -1 : cwe_dummy := 0;
		  0  : BEGIN cwe_dummy := 0; stv.q_tr_igb[stv.i2]:= 0; stv.q_tr_igb[stv.i2-1]:= 0 END;{we will not use this variable in this case, but it should have a value regardless}
		  1  : cwe_dummy := 1;
		END;
		
		stv.cwe_i[stv.i2-1] := cwe_dummy;
		stv.cwe_i[stv.i2] := cwe_dummy;
	END;

	IF (par.num_GBs > 0) AND (par.GB_tr <> 0) THEN {include grain boundary recombination}
	{If there are grain boundaries (GBs), then things get a bit complicated:
	We want to put the GBs at regular intervals. However, there may (=likely!) not be grid
	points exactly at these positions. Thus, we try to get as close as possible.} 
	BEGIN
		cur_idx := stv.i1;
		CASE par.Tr_type_B OF
			-1 : cwe_dummy := 0;
			0  : BEGIN cwe_dummy := 0; FOR i:=stv.i1+2 TO stv.i2-2 DO stv.q_tr_igb[i]:= 0 END; {we will not use this variable in this case, but it should have a value regardless}
			1  : cwe_dummy := 1;	
		END;
		
		FOR gb_nr := 1 TO par.num_GBs DO
		BEGIN 
			new_gb_found := FALSE;
			gb_centre_pos := (par.L - par.L_RTL - par.L_LTL) / (par.num_GBs + 1) * gb_nr + par.L_LTL;

			{find the start and end idx of the grain boundary and set trap density}
			WHILE new_gb_found = FALSE DO 
			BEGIN 
				IF (gb_centre_pos >= stv.x[cur_idx]) AND (gb_centre_pos < stv.x[cur_idx+1]) THEN 
				BEGIN 
					new_gb_found := TRUE;
					inv_int_trap_width := 2 / (par.L * (0.5 * stv.h[cur_idx-1] + stv.h[cur_idx] + 0.5 * stv.h[cur_idx+1]));
					stv.Nti[cur_idx,e] := stv.IntTrapDist[e] * ABS(par.GB_tr) * inv_int_trap_width;
					stv.Nti[cur_idx+1,e] := stv.IntTrapDist[e] * ABS(par.GB_tr) * inv_int_trap_width;
					
					stv.cwe_i[cur_idx] := cwe_dummy;
					stv.cwe_i[cur_idx+1] := cwe_dummy;
				END;
				INC(cur_idx);
			END
		END   
	END;

	{There cannot be two touching interfaces, as this is conflicting with how the equations are derived.}
	FOR i := 1 TO par.NP DO
		IF (stv.Nti[i-1,e] <> 0) AND (stv.Nti[i,e] <> 0) AND (stv.Nti[i+1,e] <> 0) THEN
			Stop_Prog('There are multiple consecutive interfaces, decrease interface density.', EC_InvalidInput);
	
END;

PROCEDURE Init_Traps_From_File(VAR log : TEXT; VAR Energies, Traps : TrapEnArray; VAR NumLevels : INTEGER; TrapFile : ANSISTRING);
{Inits traps if their energy levels are specified in a file}
{Energies: the energies of the traps. Traps: their relative number (per volume or area). This sums to 1}
{NumLevels: number of trap levels that are specified. Capped to Max_NEtr (in unit TypesAndConstants)}
VAR e : INTEGER;
	sum : myReal;
	x, y : Row;
BEGIN
	WRITELN('Reading trap profile from ',TrapFile);
	WRITELN(log, 'Reading trap profile from ',TrapFile);
	
	{Read data from TrapFile and put them in local vars x and y:}
	Read_XY_Table(x, y, TrapFile, 'E Ntrap', NumLevels);
	
	{check the number of levels:}
	IF NumLevels > Max_NEtr THEN Stop_Prog('There are more trap levels in file '+TrapFile+' than we can handle.'+LineEnding+'See Max_NEtr in TypesAndConstants.', EC_InvalidInput);
	IF NumLevels=0 THEN Stop_Prog('Could not find any trap levels in file '+TrapFile+'.', EC_InvalidInput);

	{now copy data into the proper arrays Traps and Energies:}
	sum:=0;
	FOR e:=1 TO NumLevels DO
	BEGIN
		Energies[e]:=x[e-1]; {note: Energies and Traps start at index 1.}
		Traps[e]:=y[e-1];
		sum:=sum + Traps[e]; 
	END;	

	{now we reschale the total number of Traps (now in sum) to 1}
	FOR e:=1 TO NumLevels DO 
		Traps[e]:=Traps[e]/sum;

END;


PROCEDURE Init_Trap_Distribution(VAR log : TEXT; VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Places all types of traps (bulk and interface) in the device at places determined by define_layers.}
VAR e, NumBulkLevels, NumIntLevels : INTEGER;
	MinETrap, MaxETrap : myReal;
BEGIN
	WITH stv DO {first set all trap-related arrays to zero}
	BEGIN
		FILLCHAR(ETrapBulk, SIZEOF(ETrapBulk), 0); 
		FILLCHAR(ETrapInt, SIZEOF(ETrapInt), 0); 
		FILLCHAR(BulkTrapDist, SIZEOF(BulkTrapDist), 0); 
		FILLCHAR(IntTrapDist, SIZEOF(IntTrapDist), 0); 
	END;

	IF par.BulkTrapFromFile AND (par.Bulk_tr <> 0) THEN {if no traps, then don't bother with the file}
		Init_Traps_From_File(log, stv.ETrapBulk, stv.BulkTrapDist, NumBulkLevels, par.BulkTrapFile)
	ELSE BEGIN {simplest case: take trap energy from device parameters}
		NumBulkLevels:=1;
		stv.ETrapBulk[1]:=par.ETrapSingle;
		stv.BulkTrapDist[1]:=1;
	END;
	
	IF par.IntTrapFromFile AND ((par.St_L <> 0) OR (par.St_R <> 0) OR ((par.GB_tr <> 0) AND (par.num_GBs <> 0))) THEN {if no traps, then don't bother with the file}
		Init_Traps_From_File(log, stv.ETrapInt, stv.IntTrapDist, NumIntLevels, par.IntTrapFile)
	ELSE BEGIN {simplest case: take trap energy from device parameters}
		NumIntLevels:=1;
		stv.ETrapInt[1]:=par.ETrapSingle;
		stv.IntTrapDist[1]:=1;
	END;

	{Now let's check if the trap energies make sense with the CB and VB in the layers in the device}
	{first: the bulk.}
	IF par.Bulk_tr <> 0 THEN 
	BEGIN
		{First: determine the relevant limits, Min/MaxETrap}
		{simplest case for bulk traps, should fall within CB and VB of bulk}
		MinETrap:=par.CB; 
		MaxETrap:=par.VB;
		
		{account for LTL:}
		IF par.TLsTrap AND (par.L_LTL<>0) THEN 
		BEGIN
			MinETrap:=MAX(MinETrap, par.CB_LTL);
			MaxETrap:=MIN(MaxETrap, par.VB_LTL);
		END;
		
		{account for RTL:}
		IF par.TLsTrap AND (par.L_RTL<>0) THEN 
		BEGIN
			MinETrap:=MAX(MinETrap, par.CB_RTL);
			MaxETrap:=MIN(MaxETrap, par.VB_RTL);
		END;
	
		{now check the bulk trap energies and see if they are in (MinETrap, MaxETrap)}
		FOR e:=1 TO NumBulkLevels DO
		BEGIN
			IF stv.ETrapBulk[e] >= MaxETrap THEN Stop_Prog('Found a bulk trap energy that is larger than physically possible.', EC_InvalidInput);
			IF stv.ETrapBulk[e] <= MinETrap THEN Stop_Prog('Found a bulk trap energy that is smaller than physically possible.', EC_InvalidInput);
		END;
	END;
	
	{now check the interface traps and GB traps}
	IF stv.Traps_int THEN 
	BEGIN
		{First: determine the relevant limits, Min/MaxETrap}
		{note: if there are grain boundary traps, then par.CB/VB are the limits and things are easy}
		MinETrap:=par.CB; 
		MaxETrap:=par.VB;
		
		{account for LTL:}
		IF (par.St_L<>0) AND (par.L_LTL<>0) THEN 
		BEGIN
			MinETrap:=MAX(MinETrap, par.CB_LTL);
			MaxETrap:=MIN(MaxETrap, par.VB_LTL);
		END;
		
		{account for RTL:}
		IF (par.St_R<>0) AND (par.L_RTL<>0) THEN 
		BEGIN
			MinETrap:=MAX(MinETrap, par.CB_RTL);
			MaxETrap:=MIN(MaxETrap, par.VB_RTL);
		END;
		
		{now check the interface trap energies and see if they are in (MinETrap, MaxETrap)}
		FOR e:=1 TO NumIntLevels DO
		BEGIN
			IF stv.ETrapInt[e] >= MaxETrap THEN Stop_Prog('Found an interface trap energy that is larger than physically possible.', EC_InvalidInput);
			IF stv.ETrapInt[e] <= MinETrap THEN Stop_Prog('Found an interface trap energy that is smaller than physically possible.', EC_InvalidInput);
		END;
		
	END;

	stv.N_Etr:=MAX(NumBulkLevels, NumIntLevels); {this guarentees that N_Etr is at least 1 but not larger than Max_NEtr}

	FOR e:=1 TO stv.N_Etr DO
		Init_Trap_Distribution_Single_Level(stv, e, par); {Places all types of traps (bulk and interface) in the device at places determined by define_layers.}
END;

PROCEDURE Rescale_Ion_Density(VAR ion : vector; conc : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{normalises an ion distribution to the correct overall concentration conc}
VAR i, istart, ifinish : INTEGER;
	sum, norm : myReal;
BEGIN
	sum:=0;
	{ions might be limited to the middle layer (which can take up the whole device i=0...NP+1.}
	IF (stv.i1>0) AND (NOT par.IonsInTLs) THEN istart:=stv.i1+1 ELSE istart:=0;
	IF (stv.i2<par.NP+1) AND (NOT par.IonsInTLs) THEN ifinish:=stv.i2-1 ELSE ifinish:=par.NP+1;
	
	FOR i:=istart+1 TO ifinish DO {note: we start at istart + 1 as we access ion[i-1]}
		sum:=sum + 0.5*(ion[i]+ion[i-1])*stv.h[i-1]*par.L; {if the grid is non-uniform, we need to take this into account}

	norm:=conc*(stv.x[ifinish]-stv.x[istart])/sum; {Note: length where ions sit depends on TLs and IonsInTLs, use ifinish and istart}
	FOR i:=istart TO ifinish DO {now renormalise the ionic densities such that the total number of ions is correct}
		ion[i]:=ion[i]*norm;
END;

PROCEDURE Calc_Ion_Distribution_Steady_State(VAR nion, pion, V : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the ion distribution in steady state. This is either constant 
(if the ionic species doesn't move), or it follows from requiring that
the associated ionic particle current be zero.} 
VAR fac : myReal;
	i, istart, ifinish : INTEGER;
BEGIN
	
	{ions can be limited to the middle layer (which can take up the whole device i=0...NP+1.}
	IF (stv.i1>0) AND (NOT par.IonsInTLs) THEN istart:=stv.i1+1 ELSE istart:=0;
	IF (stv.i2<par.NP+1) AND (NOT par.IonsInTLs) THEN ifinish:=stv.i2-1 ELSE ifinish:=par.NP+1;

	nion[istart]:=1;
	pion[istart]:=1;

	FOR i:=istart+1 TO ifinish DO
	BEGIN
		fac:=B(stv.Vti*(V[i-1]-V[i]))/B(stv.Vti*(V[i]-V[i-1]));
		{we use the expressions for the ion currents to make sure that they are zero. This yields
		the profile of neg/pos ions.}
		IF par.negIonsMove THEN nion[i]:=nion[i-1]*fac ELSE nion[i]:=1;
		IF par.posIonsMove THEN pion[i]:=pion[i-1]/fac ELSE pion[i]:=1;
	END;

	{nomalize concentrations}
	Rescale_Ion_Density(nion, par.CNI, stv, par);
	Rescale_Ion_Density(pion, par.CPI, stv, par); 
	{now the total number of ions should be equal to the concentrion CNI/CPI times (L-L_LTL-L_RTL)}
END;

PROCEDURE Solve_Neg_Ions(VAR nion : vector; nionPrevTime, V : vector; dti : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{dti: DeltaTInverse, so dti = 1/delt. If dti=0 then we're in steady-state.
If dti > 0 we're using the transient equations (Selberherr 6-4.32)}
VAR i, istart, ifinish : INTEGER;
    lo, m, u, rhs : vector;
BEGIN
	IF (stv.i1>0) AND (NOT par.IonsInTLs) THEN
	WITH stv DO 
	BEGIN {ions cannot move towards the contacts}
		istart:=i1+1;
		i:=istart;
		rhs[i]:=-0.5 * nionPrevTime[i]*dti*SQR(par.L)*h[i]*(h[i]+h[i-1]);
		lo[i]:=0;
		m[i]:=-Vt*par.mobnion*B(Vti*(V[i]-V[i+1]))
				-0.5*dti*SQR(par.L)*h[i]*(h[i]+h[i-1]);
		u[i]:=Vt*par.mobnion*B(Vti*(V[i+1]-V[i]));   
	END
	ELSE BEGIN {ions can move towards the contacts, nion[0]<>0:}
		istart:=0;
		{we need to ensure that the ionic currents into/out of the contacts be zero, so let's do that first:}
		lo[istart]:=0;
		m[istart]:=1;
		u[istart]:=-B(stv.Vti*(V[1]-V[0])) / B(stv.Vti*(V[0]-V[1]));
		rhs[istart]:=0; 
	END;
  
	{at the other end:}
	IF (stv.i2<par.NP+1) AND (NOT par.IonsInTLs) THEN
	WITH stv DO
	BEGIN  {ions cannot move towards the contacts}
		ifinish:=i2-1;
		i:=ifinish;
		rhs[i]:=-0.5 * nionPrevTime[i]*dti*SQR(par.L)*h[i-1]*(h[i]+h[i-1]);
		lo[i]:=Vt*par.mobnion*B(Vti*(V[i-1]-V[i]));
		m[i]:=-Vt*par.mobnion*B(Vti*(V[i]-V[i-1]))
				-0.5*dti*SQR(par.L)*h[i-1]*(h[i]+h[i-1]);
		u[i]:=0;
	END 
	ELSE BEGIN {ions can move towards the contacts, nion[NP+1]<>0:}
		ifinish:=par.NP+1;
		lo[ifinish]:=-B(stv.Vti*(V[par.NP]-V[par.NP+1])) / B(stv.Vti*(V[par.NP+1]-V[par.NP]));
		m[ifinish]:=1;
		u[ifinish]:=0;
		rhs[ifinish]:=0;
	END;
	
	{now set the interior part:}
    FOR i:=istart+1 TO ifinish-1 DO  {continuity eq. in matrix vorm, equivalent to that of n and p, but without generation and recombination}
    WITH stv DO
    BEGIN
        rhs[i]:=-0.5 * nionPrevTime[i]*dti*SQR(par.L)*h[i]*h[i-1]*(h[i]+h[i-1]);
        lo[i]:=h[i]*Vt*par.mobnion*B(Vti*(V[i-1]-V[i]));
        m[i]:=-(h[i-1]*Vt*par.mobnion*B(Vti*(V[i]-V[i+1])) +
            h[i]*Vt*par.mobnion*B(Vti*(V[i]-V[i-1])))
            -0.5*dti*SQR(par.L)*h[i]*h[i-1]*(h[i]+h[i-1]);
        u[i]:=h[i-1]*Vt*par.mobnion*B(Vti*(V[i+1]-V[i]));          
    END;

 	{Solve nion from istart to ifinish:}
	Tridiag(nion, lo, m, u, rhs, istart, ifinish);
	
	{now check if nion is still well behaved, i.e. positive, and has the correct overall value}
	FOR i:=istart TO ifinish DO
		IF (nion[i]<0) THEN 
		BEGIN
			IF par.IgnoreNegDens 
				THEN nion[i]:=-nion[i] ELSE
				Stop_Prog('Negative concentration of negative ions encountered!', EC_NumericalFailure);
		END;
	{make sure the number of ions is preserved, i.e. correct:}
	Rescale_Ion_Density(nion, par.CNI, stv, par);
END;

PROCEDURE Solve_Pos_Ions(VAR pion : vector; pionPrevTime, V : vector; dti : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{dti: DeltaTInverse, so dti = 1/delt. If dti=0 then we're in steady-state.
If dti > 0 we're using the transient equations (Selberherr 6-4.33)}
VAR i, istart, ifinish : INTEGER;
    lo, m, u, rhs : vector;
BEGIN
	IF (stv.i1>0) AND (NOT par.IonsInTLs) THEN
	WITH stv DO
	BEGIN {ions cannot move towards the contacts}
		istart:=i1+1;
		i:=istart;
		rhs[i]:=-0.5 * pionPrevTime[i]*dti*SQR(par.L)*h[i]*(h[i]+h[i-1]);
		lo[i]:=0;
		m[i]:=-Vt*par.mobpion*B(Vti*(V[i+1]-V[i])) 
				-0.5*dti*SQR(par.L)*h[i]*(h[i]+h[i-1]);
		u[i]:=Vt*par.mobpion*B(Vti*(V[i]-V[i+1])); 
	END
	ELSE BEGIN {ions can move towards the contacts, pion[0]<>0:}
		istart:=0;
		{we need to ensure that the ionic currents into/out of the contacts be zero, so let's do that first:}
		lo[istart]:=0;
		m[istart]:=1;
		u[istart]:=-B(stv.Vti*(V[0]-V[1])) / B(stv.Vti*(V[1]-V[0]));
		rhs[istart]:=0;
	END;
  
	{at the other end:}
	IF (stv.i2<par.NP+1) AND (NOT par.IonsInTLs) THEN
	WITH stv DO
	BEGIN  {ions cannot move towards the contacts}
		ifinish:=i2-1;
		i:=ifinish;
		rhs[i]:=-0.5 * pionPrevTime[i]*dti*SQR(par.L)*h[i-1]*(h[i]+h[i-1]);
		lo[i]:=Vt*par.mobpion*B(Vti*(V[i]-V[i-1]));
		m[i]:=-Vt*par.mobpion*B(Vti*(V[i-1]-V[i]))
				-0.5*dti*SQR(par.L)*h[i-1]*(h[i]+h[i-1]);
		u[i]:=0;  
	END 
	ELSE BEGIN {ions can move towards the contacts, pion[NP+1]<>0:}
		ifinish:=par.NP+1;
		lo[ifinish]:=-B(stv.Vti*(V[par.NP+1]-V[par.NP])) / B(stv.Vti*(V[par.NP]-V[par.NP+1])); 
		m[ifinish]:=1;
		u[ifinish]:=0;
		rhs[ifinish]:=0;
	END;	
		
	{now set the interior part:}
    FOR i:=istart+1 TO ifinish-1 DO  {continuity eq. in matrix vorm, equivalent to that of n and p, but without generation and recombination}
    WITH stv DO
    BEGIN
        rhs[i]:=-0.5 * pionPrevTime[i]*dti*SQR(par.L)*h[i]*h[i-1]*(h[i]+h[i-1]);
        lo[i]:=h[i]*Vt*par.mobpion*B(Vti*(V[i]-V[i-1]));
        m[i]:=-(h[i-1]*Vt*par.mobpion*B(Vti*(V[i+1]-V[i])) +
            h[i]*Vt*par.mobpion*B(Vti*(V[i-1]-V[i])))
            -0.5*dti*SQR(par.L)*h[i]*h[i-1]*(h[i]+h[i-1]);
        u[i]:=h[i-1]*Vt*par.mobpion*B(Vti*(V[i]-V[i+1]));         
    END;

 	{Solve pion from istart to ifinish:}
    Tridiag(pion, lo, m, u, rhs, istart, ifinish); {Solve for the new ion densities}

	{now check if pion is still well behaved, i.e. positive, and has the correct overall value}
	FOR i:=istart TO ifinish DO
		IF (pion[i]<0) THEN 
		BEGIN
			IF par.IgnoreNegDens 
				THEN pion[i]:=-pion[i] ELSE
				Stop_Prog('Negative concentration of positive ions encountered!', EC_NumericalFailure);
		END;
	{make sure the number of ions is preserved, i.e. correct:}
	Rescale_Ion_Density(pion, par.CPI, stv, par);
END;


FUNCTION Calc_f_ti_Numer(n, p : vector; e : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters): vector;
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculates the numerator of the fraction of interface traps filled. The numerator and denominator are split because both are re-used in different places in the Poisson and continuity equations.}
VAR i, j : INTEGER;
    f_t	 : myReal; {fraction of filled traps, i.e. traps that contain an electron}
BEGIN
    Calc_f_ti_Numer[0] := 0;
	Calc_f_ti_Numer[par.NP+1] := 0;
    
	FOR i := 1 TO par.NP DO
    BEGIN
		IF stv.Nti[i,e] <> 0 THEN
		BEGIN
			f_t := 0;
			FOR j:=-1 TO 1 DO f_t := f_t + par.Cn * (stv.Nti[i+j,e] * n[i+j]) + par.Cp * (stv.Nti[i+j,e] * stv.pt0i[i+j,e]);
			Calc_f_ti_Numer[i] := f_t;
		END
		ELSE Calc_f_ti_Numer[i] := 0;
    END;   
END;

FUNCTION Calc_Inv_f_ti_Denom(n, p : vector; e : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters): vector;
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculates the inverse of the denominator of the fraction of interface traps filled. 
The numerator and denominator are split because both are re-used in different places in the Poisson and continuity equations.}
VAR i, j : INTEGER;
    denom_trap_occup : myReal;
BEGIN
    Calc_Inv_f_ti_Denom[0] := 0;
    Calc_Inv_f_ti_Denom[par.NP+1] := 0;    
    FOR i := 1 TO par.NP DO
    BEGIN
		IF (stv.Nti[i,e] <> 0) THEN
		BEGIN
			denom_trap_occup:=0;
			FOR j:=-1 TO 1 DO denom_trap_occup := denom_trap_occup + par.Cn * (stv.Nti[i+j,e] * n[i+j] + stv.Nti[i+j,e] * stv.nt0i[i+j,e]) + par.Cp * (stv.Nti[i+j,e] * p[i+j] + stv.Nti[i+j,e] * stv.pt0i[i+j,e]);
			calc_inv_f_ti_denom[i] := 1/denom_trap_occup;
		END
		ELSE Calc_Inv_f_ti_Denom[i] := 0;
    END;
END;

FUNCTION Calc_f_ti(CONSTREF n, p : vector; Old_f_ti : TrapArray; dti : myReal; e : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : vector;
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculates the fraction of interface traps filled with an electron, 
and takes as input the result of Calc_f_ti_Numer and Calc_Inv_f_ti_Denom.}
VAR i : INTEGER;
	a, b : myReal;
	f_ti_numer, f_ti_inv_denom : vector;
BEGIN
	{we have excluded the possibility of having interface traps on the electrodes (= no TLs), see Check_Parameters:}
	Calc_f_ti[0] := 0;
	Calc_f_ti[par.NP+1] := 0;
	
	IF dti=0 THEN {steady-state}
	BEGIN
		f_ti_numer:=Calc_f_ti_Numer(n, p, e, stv, par);
		f_ti_inv_denom:=Calc_Inv_f_ti_Denom(n, p, e, stv, par);	
		FOR i := 1 TO par.NP DO 
			Calc_f_ti[i] := f_ti_numer[i] * f_ti_inv_denom[i]
	END
	ELSE {transient}
		FOR i := 1 TO par.NP DO
			IF stv.Nti[i-1,e] * stv.Nti[i,e] <> 0 THEN
			BEGIN {now we have just crossed an interface}
				a:= par.Cn*n[i] + par.Cn*stv.nt0i[i,e] + par.Cp*stv.pt0i[i,e] + par.Cp*p[i] +
					par.Cn*n[i-1] + par.Cn*stv.nt0i[i-1,e] + par.Cp*stv.pt0i[i-1,e] + par.Cp*p[i-1];
				b:= par.Cn*n[i] + par.Cp*stv.pt0i[i,e] +
					par.Cn*n[i-1] + par.Cp*stv.pt0i[i-1,e];
				
				Calc_f_ti[i]:= b/a + (Old_f_ti[i,e] - b/a)*EXP(-a/dti);
				Calc_f_ti[i-1]:= Calc_f_ti[i];
				{we solve f_tb from d(f_tb)/dt = a f_tb + b, which is a trivial ODE}
			END
			ELSE Calc_f_ti[i] := 0
END;

FUNCTION Calc_f_tb(CONSTREF n, p : vector; Old_f_tb : TrapArray; dti : myReal; e : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : vector;
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
VAR i : INTEGER;
    a, b : myReal;
	{Calculates the fraction of bulk traps filled with an electron.}
BEGIN
	FOR i := 0 TO par.NP+1 DO
    BEGIN
		IF (stv.Ntb[i,e] <> 0) THEN
		BEGIN 
			IF dti=0 THEN {steady-state}
				Calc_f_tb[i]:=(par.Cn*n[i] + par.Cp*stv.pt0b[i,e]) / (par.Cn*(n[i]+stv.nt0b[i,e]) + par.Cp*(p[i]+stv.pt0b[i,e]))
			ELSE BEGIN {transient}
			    b:= par.Cn*n[i] + par.Cp*stv.pt0b[i,e];
			    a:= par.Cn*n[i] + par.Cn*stv.nt0b[i,e] + par.Cp*stv.pt0b[i,e] + par.Cp*p[i];
			    Calc_f_tb[i]:= b/a + (Old_f_tb[i,e] - b/a)*EXP(-a/dti);
				{we solve f_tb from d(f_tb)/dt = a f_tb + b, which is a trivial ODE}
			END
		END
		ELSE {no traps, simply set f_tb to zero}
			Calc_f_tb[i] := 0
    END;
END;

PROCEDURE Calc_Trap_Filling_Charge_Single_Level(VAR f_tb, f_ti : TrapArray; VAR Ntb_charge, Nti_charge : vector; CONSTREF n, p : vector; Old_f_tb, Old_f_ti : TrapArray; dti : myReal; e : INTEGER;
												CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the charge of the (bulk and interface) traps at every gridpoint for a given trap energy.}
VAR i								: INTEGER;
	zeros, f_tb_single, f_ti_single	: vector;
BEGIN
	FOR i:=0 TO par.NP+1 DO zeros[i]:=0;

	{note: f_tb/i is the occupancy of the bulk/interface traps}
	{if there are no traps, or if the traps are neutral (Tr_type_B/I=0), then this doesn't matter}
	
	{bulk traps:}
	IF stv.BulkTrapDist[e] <> 0 THEN
		f_tb_single := Calc_f_tb(n, p, Old_f_tb, dti, e, stv, par)
	ELSE
		f_tb_single := zeros;
	FOR i:=1 TO par.NP DO f_tb[i,e]:= f_tb_single[i];

	IF (stv.BulkTrapDist[e] <> 0) AND (par.Tr_type_B <> 0) THEN
		{bulk trap charge}
		FOR i:=1 TO par.NP DO Ntb_charge[i] := ABS(par.Tr_type_B)*(stv.cwe_b - f_tb[i,e]) * stv.Ntb[i,e] {note that we have f_tb minus charges or 1-f_tb positive charges}
	ELSE
		Ntb_charge := zeros;
	
	{interface traps:}
	IF stv.Traps_int AND (stv.IntTrapDist[e] <> 0) THEN
		f_ti_single := Calc_f_ti(n, p, Old_f_ti, dti, e, stv, par)
	ELSE
		f_ti_single := zeros;
	FOR i:=1 TO par.NP DO f_ti[i,e]:= f_ti_single[i];
		
	IF stv.Traps_int_poisson AND (stv.IntTrapDist[e] <> 0) THEN
		FOR i:=0 TO par.NP+1 DO Nti_charge[i] := stv.q_tr_igb[i]*(stv.cwe_i[i] - f_ti[i,e]) * 0.5 * stv.Nti[i,e] {each side of the interface hosts half the traps}
	ELSE
		Nti_charge := zeros;
END;


PROCEDURE Calc_Trap_Filling_Charge(VAR f_tb, f_ti : TrapArray; VAR Ntb_charge, Nti_charge : vector; CONSTREF n, p : vector; Old_f_tb, Old_f_ti : TrapArray; dti : myReal;
									CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the charge of the (bulk and interface) traps at every gridpoint.}
VAR i, e : INTEGER;
	Ntb_charge_single, Nti_charge_single : vector;
BEGIN
	FOR e:=1 TO stv.N_Etr DO
	BEGIN
		Calc_Trap_Filling_Charge_Single_Level(f_tb, f_ti, Ntb_charge_single, Nti_charge_single, n, p, Old_f_tb, Old_f_ti, dti, e, stv, par);

		IF e=1 THEN
			FOR i:=1 TO par.NP DO
			BEGIN
				Ntb_charge[i]:= Ntb_charge_single[i];
				Nti_charge[i]:= Nti_charge_single[i];
			END
		ELSE
		BEGIN
			FOR i:=1 TO par.NP DO
			BEGIN
				Ntb_charge[i]:= Ntb_charge[i] + Ntb_charge_single[i];
				Nti_charge[i]:= Nti_charge[i] + Nti_charge_single[i];				
			END
		END;		
	END;
END;


PROCEDURE Calc_Linearization_f_t_Single_Level(VAR lin : TLinFt; CONSTREF n, p : vector; dti : myReal; e : INTEGER;
                                  CONSTREF stv		  : TStaticVars; CONSTREF par : TInputParameters); 
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculate the linearization of the filled (with electrons) fraction of interface and bulk 
traps with respect to the potential. This function is written in such a way, 
that the result can be added to the Poisson equation regardless of the presence of traps. 
If traps are not present the elements added to the Poisson equation simply become zero.
A derivation can be found in the docs.}
VAR f_tb_numer, f_tb_inv_denom : myReal;
	zeros, f_ti_numer, f_ti_inv_den	: vector;
	i : INTEGER;
BEGIN	
	FOR i:=0 TO par.NP+1 DO zeros[i]:=0;

	{Bulk traps: if we have them, calculate their linearization in delta V.}
	IF (dti=0) AND (par.bulk_tr <> 0) AND (par.Tr_type_B <> 0) AND (stv.BulkTrapDist[e] <> 0) THEN
		FOR i:=1 TO par.NP DO
		BEGIN
			f_tb_numer := par.Cn*n[i] + par.Cp*stv.pt0b[i,e];
			f_tb_inv_denom := 1 / (par.Cn*(n[i]+stv.nt0b[i,e]) + par.Cp*(p[i]+stv.pt0b[i,e]));
			lin.f_tb_m[i] := par.Cn * n[i+1] * f_tb_inv_denom;
			lin.f_tb_m[i] := lin.f_tb_m[i] - (par.Cn * n[i] - par.Cp * p[i]) * f_tb_numer * SQR(f_tb_inv_denom);
			lin.f_tb_m[i] := par.Tr_type_B*(stv.cwe_b - lin.f_tb_m[i]) * stv.Ntb[i,e];
		END
    ELSE
		lin.f_tb_m:= zeros;
			
	{Linearize f_ti (interface trap occupance fraction) in delV, a derivation can be found in the docs.}
	{Key here is that all linearized elements become zero when they should not affect the Poisson equation}
	IF (dti=0) AND (stv.Traps_int_poisson) AND (stv.IntTrapDist[e] <> 0) THEN
	BEGIN
		FOR i:=1 TO par.NP DO {filling the matrix en the right-hand side}
		BEGIN		
			f_ti_numer := Calc_f_ti_Numer(n, p, e, stv, par);
			f_ti_inv_den := Calc_Inv_f_ti_Denom(n, p, e, stv, par);
		
			IF (stv.Nti[i-1,e] = 0) THEN
				lin.f_ti_lo[i] := 0
			ELSE
			BEGIN
				lin.f_ti_lo[i] := par.Cn * n[i-1] * stv.Nti[i-1,e] * stv.Vti * f_ti_inv_den[i];
				lin.f_ti_lo[i] := lin.f_ti_lo[i] - (par.Cn * n[i-1] - par.Cp * p[i-1]) * stv.Nti[i-1,e] * stv.Vti * f_ti_numer[i] * SQR(f_ti_inv_den[i]);
				lin.f_ti_lo[i] := stv.q_tr_igb[i]*lin.f_ti_lo[i] * 0.5 * stv.Nti[i-1,e];
			END;

			IF (stv.Nti[i,e] = 0) THEN
				lin.f_ti_m[i] := 0
			ELSE
			BEGIN
				{the factor stv.Vti is multiplied in the poisson equation with for this main diagonal}
				lin.f_ti_m[i] := par.Cn * n[i] * stv.Nti[i,e] * f_ti_inv_den[i];
				lin.f_ti_m[i] := lin.f_ti_m[i] - (par.Cn * n[i] - par.Cp * p[i]) * stv.Nti[i,e] * f_ti_numer[i] * SQR(f_ti_inv_den[i]);
				lin.f_ti_m[i] := stv.q_tr_igb[i]*lin.f_ti_m[i] * 0.5 * stv.Nti[i,e];
			END;

			IF (stv.Nti[i+1,e] = 0) THEN
				lin.f_ti_up[i] := 0
			ELSE
			BEGIN
				lin.f_ti_up[i] := par.Cn * n[i+1] * stv.Nti[i+1,e] * stv.Vti * f_ti_inv_den[i];
				lin.f_ti_up[i] := lin.f_ti_up[i] - (par.Cn * n[i+1] - par.Cp * p[i+1]) * stv.Nti[i+1,e] * stv.Vti * f_ti_numer[i] * SQR(f_ti_inv_den[i]);
				lin.f_ti_up[i] := stv.q_tr_igb[i]*lin.f_ti_up[i] * 0.5 * stv.Nti[i+1,e];
			END;
		END;
	END
	ELSE
	BEGIN
		lin.f_ti_up:= zeros;
		lin.f_ti_lo:= zeros;
		lin.f_ti_m:= zeros;
	END;
END;

PROCEDURE Calc_Linearization_f_t_All(VAR lin : TLinFt; CONSTREF n, p : vector; dti : myReal;
                                  CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the charge of the (bulk and interface) traps at every gridpoint.}
VAR i, e : INTEGER;
	lin_single : TLinFt; {this type (a record) stores the linearisation of the trapping terms.}	
BEGIN
	FOR e:=1 TO stv.N_Etr DO
	BEGIN
	    Calc_Linearization_f_t_Single_Level(lin_single, n, p, dti, e, stv, par);
		IF e=1 THEN
			FOR i:=1 TO par.NP DO
			BEGIN
				lin.f_tb_m[i]:= lin_single.f_tb_m[i];
				lin.f_ti_up[i]:= lin_single.f_ti_up[i];
				lin.f_ti_m[i]:= lin_single.f_ti_m[i];
				lin.f_ti_lo[i]:= lin_single.f_ti_lo[i];
			END
		ELSE
		BEGIN
			FOR i:=1 TO par.NP DO
			BEGIN
				lin.f_tb_m[i]:= lin.f_tb_m[i] + lin_single.f_tb_m[i];
				lin.f_ti_up[i]:= lin.f_ti_up[i] + lin_single.f_ti_up[i];
				lin.f_ti_m[i]:= lin.f_ti_m[i] + lin_single.f_ti_m[i];
				lin.f_ti_lo[i]:= lin.f_ti_lo[i] + lin_single.f_ti_lo[i];
			END;
		END;
	END;
END;

PROCEDURE Solve_Poisson(VAR V, n, p, nion, pion	: vector; VAR f_tb, f_ti : TrapArray; VAR Ntb_charge, Nti_charge : vector; CONSTREF Old_f_tb, Old_f_ti : TrapArray;
						VAR conv, coupleIonsPoisson	: BOOLEAN; VAR PoissMsg : STRING; dti : myReal;
						CONSTREF stv				: TStaticVars; CONSTREF par : TInputParameters);
{Solves the Poisson equation, can be used in steady-state and transient simulations}
{Solve_Poisson also modifies the charges (n,p,ions,traps) by estimating the effects of the newly calc'd potential}
VAR it, i : INTEGER;
	sumPre, sumPost, NormDelV : myReal;
    delV, rhs, lower, upper, main : vector;
    fac_m, fac_u, fac_l, fac2, fac3, val, lnr : myReal;
    IonsOK : BOOLEAN;
	lin	  : TLinFt; {this type (a record) stores the linearisation of the trapping terms.}
BEGIN
    FOR i:=1 TO par.NP DO delV[i]:=1; {init delV}
    delV[0]:=0; {delV=0 at contacts, by definition, since we know V(0,L)}
    delV[par.NP+1]:=0;
    it:=0;
    conv:=FALSE;
	PoissMsg:='Poisson solver status:' + LineEnding; {message string}
	lnr:=LN(1 + par.couplePC);

	{if needed, check the total number of ions in volume}
	IF (par.negIonsMove OR par.posIonsMove) AND coupleIonsPoisson THEN
	BEGIN
		sumPre:=0;
		FOR i:=0 TO par.NP DO
			sumPre:=sumPre + nion[i] + pion[i];
	END;
	
    WHILE (NOT conv) AND (it < par.MaxItPois) DO
    BEGIN
		Calc_Trap_Filling_Charge(f_tb, f_ti, Ntb_charge, Nti_charge, n, p, Old_f_tb, Old_f_ti, dti, stv, par);
		Calc_Linearization_f_t_All(lin, n, p, dti, stv, par);

        FOR i:=1 TO par.NP DO {filling the matrix and the right-hand side}
			WITH stv DO
			BEGIN
			    {to properly deal with non-uniform dielectric constants we need an extra term in the Poisson equation, that's where these factors originate.}
				fac_m:= (eps[i]-eps[i-1])/(SQR(h[i-1])*SQR(par.L)) - 2*eps[i]*(1 / (h[i]*(h[i]+h[i-1])*SQR(par.L)) + 1 / (h[i-1]*(h[i]+h[i-1])*SQR(par.L)));
				fac_l:= 2*eps[i] / (h[i-1]*(h[i]+h[i-1])*SQR(par.L)) -(eps[i]-eps[i-1])/(SQR(h[i-1])*SQR(par.L));
				fac_u:= 2*eps[i] / (h[i]*(h[i]+h[i-1])*SQR(par.L));

				rhs[i]:=- fac_l*V[i-1] - fac_u*V[i+1] -fac_m*V[i]
				        + q*(n[i] + nion[i] + pid[i] - p[i] - pion[i] - nid[i] - Nti_charge[i] - Ntb_charge[i]); {add all negative charges and substract all positive charges.}

				lower[i]:= fac_l - q*lin.f_ti_lo[i]; {linearization of f_ti with respect to delV at gridpoint i}
				upper[i]:= fac_u - q*lin.f_ti_up[i]; {linearization of f_ti with respect to delV at gridpoint i}
				{While unintuitive, the main diagonal adds all charges regardless of sign. This is because we don't add the charges but their derivatives to delV, which is also the source of the factor Vti.}
				main[i]:= fac_m - q*Vti*(n[i] + p[i] + nion[i] + pion[i] + lin.f_tb_m[i] + lin.f_ti_m[i]);
			END;

        Tridiag(delV, lower, main, upper, rhs, 1, par.NP); {solve for delV}
	
        FOR i:=1 TO par.NP DO 
        BEGIN
			delV[i]:=SIGN(delV[i])*MIN(par.maxDelV*stv.Vt, ABS(delV[i])); {limit delV to a pre-set max}
            V[i]:=V[i] + delV[i];  {and add delV to V}
			IF dti=0 THEN BEGIN
				{Couple the Poisson to the cont. eqs. and makes convergence a lot easier!}
				val:=SIGN(delV[i]) * MIN(lnr, ABS(delV[i])*stv.Vti);
				fac2:=EXP(val);
				fac3:=1/fac2;
				n[i]:=n[i]*fac2; {now update the densities: we have to do this}
				p[i]:=p[i]*fac3; {in order to conserve Gummel iteration}
				{we also apply this to the ions. If IonsInTLs = 0 then we can still do this even though it's redundant for i<i1 and i>i2}
				IF par.negIonsMove AND coupleIonsPoisson THEN nion[i]:=nion[i]*fac2; {and also apply this to the ions}
				IF par.posIonsMove AND coupleIonsPoisson THEN pion[i]:=pion[i]*fac3
			END
		END; {for loop}

        it:=it+1;
        NormDelV:=Norm_Eucl(delV, 1, par.NP);
        conv:=NormDelV <= par.tolPois  {finally, check for convergence}
    END;
	
	PoissMsg:=PoissMsg +'- delV ='+FloatToStrF(NormDelV, ffGeneral,5,0) + LineEnding;
    
    {OK, now see if we haven't changed the ions too much:}
    IF (par.negIonsMove OR par.posIonsMove) AND coupleIonsPoisson THEN
    BEGIN
		sumPost:=0;
		FOR i:=0 TO par.NP DO
			sumPost:=sumPost + nion[i] + pion[i];
		{note: by now, sumPost cannot be zero as there are ions}
		IonsOK:=(ABS(sumPre-sumPost)/sumPost) < par.tolPois;
		conv:=conv AND IonsOK;
		PoissMsg:=PoissMsg + '- movement of ions acceptable: ' + myBoolStr(IonsOK) + LineEnding;
	END;
	
	PoissMsg:=PoissMsg + '- converged: ' + myBoolStr(conv);
	
END;

PROCEDURE Calc_Elec_Mob(VAR mu : vector; CONSTREF V, n : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
VAR i : INTEGER;
{calculates the elec. mob. on the interleaved mesh, mun[i]=mun at x=i+1/2}
{therefore the field or concentration on x=i+1/2 is needed}
BEGIN
    CASE par.mob_n_dep OF
        0 : FOR i:=0 TO par.NP DO mu[i]:=par.mun_0; { mob. is constant }
        1 : FOR i:=0 TO par.NP DO {field-dep. mob}
                mu[i]:=par.mun_0 * EXP(par.gamma_n*SQRT(ABS(V[i+1]-V[i])/(par.L*stv.h[i])));
		ELSE Stop_Prog('Only very simple mobility models (0 and 1) currently implemented.', EC_ProgrammingError);
    END; {case selector}

    mu[par.NP+1]:=mu[par.NP]; {does not have a meaning, should not be used}

	IF stv.i1>0 THEN
    BEGIN
		FOR i:=0 TO stv.i1-1 DO mu[i]:=par.mob_LTL;
		{remember: i1 is the last point (highest index) in the right TL}
		mu[stv.i1]:=(par.L*stv.h[stv.i1])*par.nu_int_LTL*stv.Vti; {mobility AT the TL/main absorber interface}
	END;
	IF stv.i2<par.NP+1 THEN
	BEGIN
		mu[stv.i2-1]:=(par.L*stv.h[stv.i2-1])*par.nu_int_RTL*stv.Vti; {mobility AT the TL/main absorber interface}
		{remember: i2 is the first point (lowest index) in the right TL}
		FOR i:=stv.i2 TO par.NP+1 DO mu[i]:=par.mob_RTL;
	END
END;

PROCEDURE Calc_Hole_Mob(VAR mu : vector; CONSTREF V, p : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
VAR i : INTEGER;
{calculates the hole mob. on the interleaved mesh, mup[i]=mup at x=i+1/2}
{therefore the field or concentration on x=i+1/2 is needed}
BEGIN
    CASE par.mob_p_dep OF
        0 : FOR i:=0 TO par.NP DO mu[i]:=par.mup_0; { mob. is constant }
        1 : FOR i:=0 TO par.NP DO {field-dep. mob}
                mu[i]:=par.mup_0 * EXP(par.gamma_p*SQRT(ABS(V[i+1]-V[i])/(par.L*stv.h[i])));
		ELSE Stop_Prog('Only very simple mobility models (0 and 1) currently implemented.', EC_ProgrammingError);
    END; {case selector}
    mu[par.NP+1]:=mu[par.NP]; {does not have a meaning, should not be used}

    IF stv.i1>0 THEN
    BEGIN
		FOR i:=0 TO stv.i1-1 DO mu[i]:=par.mob_LTL;
		{remember: i1 is the last point (highest index) in the right TL}
		mu[stv.i1]:=(par.L*stv.h[stv.i1])*par.nu_int_LTL*stv.Vti; {mobility AT the TL/main absorber interface}
	END;
	IF stv.i2<par.NP+1 THEN
	BEGIN
		mu[stv.i2-1]:=(par.L*stv.h[stv.i2-1])*par.nu_int_RTL*stv.Vti; {mobility AT the TL/main absorber interface}
		{remember: i2 is the first point (lowest index) in the right TL}
		FOR i:=stv.i2 TO par.NP+1 DO mu[i]:=par.mob_RTL;
	END

END;

PROCEDURE Calc_Langevin_Factor(VAR Lan : vector; mob_n, mob_p : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculates the Langevin recombination strength. Lan[i] is defined on the
regular grid, i.e., Lan[i]=Lan at x=xi. Note that mobilities are defined on the
interleaved mesh!}
VAR i : INTEGER;
    rec_mob : myReal;
BEGIN
    IF par.UseLangevin THEN {use Langevin formula to calculate bimolecular, direct, band-to-band recombination rate}
    FOR i:=1 TO par.NP DO
    BEGIN
		rec_mob:=(mob_n[i-1] + mob_n[i]+ mob_p[i-1] + mob_p[i] )/2;
		{we take mob(x=xi)=(mob(x=xi-1/2)+mob(x=xi+1/2))/2}
		Lan[i]:=par.Lang_pre * q * rec_mob/stv.eps[i];
    END
    ELSE
        FOR i:=1 TO par.NP DO {use input value for direct recombination}
            Lan[i]:=par.kdirect;

    Lan[0]:=0; {no recombination (or generation) at the contacts}
    Lan[par.NP+1]:=0;
END;


FUNCTION Diss_Prob_Delta(r : myReal; vals : Row) : myReal;
{Calculates the dissociation probability as a function of distance r}
VAR b, kdF, delE, epsi, Vti, Braun_rec, F, kf : myReal;
BEGIN
    Vti:=vals[1]; {local inverse thermal voltage}
    epsi:=vals[2]; {local relative dielectric constant}
    Braun_rec:=vals[3]; {local copy of Braun recombination strength}
    F:=vals[4]; {local copy of electric field}
    kf:=vals[5]; {local copy of kf}
    delE:=q/(4*PI*epsi*r); {binding energy in eV}
    b:=q*ABS(F)*SQR(Vti)/(8*PI*epsi); {scaled ABOLUTE field strength}
    kdF:=3*Braun_rec/(4*PI*r*SQR(r))*EXP(-delE*Vti)*Bessel(b);
    Diss_Prob_Delta:=kdF/(kdF + kf)
END;

FUNCTION Diss_Prob_Gauss(r : myReal; vals : Row) : myReal;
{calculates the diss. prob. as a function of distance r with a Gaussian distribution of a}
VAR a : myReal;
BEGIN
    a:=vals[0]; {charge separation distance in Braun model}
    Diss_Prob_Gauss:=(4/(a*a*a*SQRT(PI)))*SQR(r)*EXP(-SQR(r/a))*Diss_Prob_Delta(r, vals)
END;

FUNCTION Diss_Prob_Exp(r : myReal; vals : Row) : myReal;
{calculates the diss. prob. as a function of distance r with an exponential
distribution of a}
VAR a : myReal;
BEGIN
    a:=vals[0]; {charge separation distance in Braun model}
    Diss_Prob_Exp:=EXP(-r/a)/a * Diss_Prob_Delta(r, vals)
END;

FUNCTION Diss_Prob_SQRrExp(r : myReal; vals : Row) : myReal;
{calculates the diss. prob. as a function of distance r with an r^2 exponential
distribution of a}
VAR a : myReal;
BEGIN
    a:=vals[0]; {charge separation distance in Braun model}
    Diss_Prob_SQRrExp:=SQR(r)*EXP(-r/a)/(2*a*a*a) * Diss_Prob_Delta(r, vals)
END;

FUNCTION Diss_Prob_r4Gauss(r : myReal; vals : Row) : myReal;
{calculates the diss. prob. as a function of distance r with an r^4 * Gaussian
distribution of a}
VAR a : myReal;
BEGIN
    a:=vals[0]; {charge separation distance in Braun model}	
	Diss_Prob_r4Gauss:=8*SQR(r)*SQR(r)*EXP(-SQR(r/a))/(3*Power(a,5)*SQRT(PI)) * Diss_Prob_Delta(r, vals)
END;

PROCEDURE Calc_Dissociation(VAR dp, g : vector; Gm, Lan, V, mob_n, mob_p : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the field-dependent dissociation rate dp and the generation rate,
g[i]=g at x=i, see C. Braun, J. Chem. Phys. 80, p. 4157 (1984)}
VAR i : INTEGER;
    vals : Row; {array with values needed for Romberg integration}
BEGIN
    IF par.Field_dep_G THEN { calc. field dep. G }
    BEGIN
		SetLength(vals, 6); {vals: needed to pass values to MathFuncValues}
		vals[0]:=par.a; {we'll use it (here) to pass the charge separation distance}
		vals[1]:=stv.Vti; {and the inverse thermal voltage} 
        vals[5]:=par.kf; 
        FOR i:=1 TO par.NP DO
        BEGIN
            vals[4]:=(V[i+1]-V[i-1])/(par.L*(stv.h[i]+stv.h[i-1])); {local electric field}
            {F= field on mesh point i, centered difference, global variable}      

            vals[2]:=stv.eps[i]; {copy local dielectric constant into vals as Diss_Prob_Delta needs it}
            vals[3]:=Lan[i]; {Braun recombination strength at x=xi, equal to Langevin (direct)}

            CASE par.ThermLengDist OF  {a = thermalization length}
                1 : dp[i]:=Diss_Prob_Delta(par.a, vals); {delta-function distribution}
                2 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_Gauss, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
                3 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_Exp, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
                4 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_SQRrExp, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
                5 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_r4Gauss, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
            END; {case selector}
            {total free-carrier yield is sum of direct generation (P0) and the field dependent part (1-P0)*(dp)}
            g[i]:=(par.P0 + (1-par.P0)*dp[i]) * Gm[i]
        END {for loop}
    END { calc. field dep. G }
    ELSE FOR i:=1 TO par.NP DO BEGIN dp[i]:=0; g[i]:=Gm[i] END; {G is constant}
    dp[0]:=0;
    dp[par.NP+1]:=0;
    g[0]:=0;   {no generation at the contacts}
    g[par.NP+1]:=0
END;


PROCEDURE Init_nt0_And_pt0(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{inits nt0 and pt0 arrays (detrapping rates) needed for SRH recombination}
{note: stv are changed here (nt0 and pt0), so they are VAR parameters}
VAR i, e : INTEGER;
BEGIN
    WITH stv DO BEGIN
		FOR e:=1 TO stv.N_Etr DO BEGIN
			nt0b[0,e] := 0;
			nt0b[par.NP+1,e] := 0;
			pt0b[0,e] := 0;
			pt0b[par.NP+1,e] := 0;
			nt0i[0,e] := 0;
			nt0i[par.NP+1,e] := 0;
			pt0i[0,e] := 0;
			pt0i[par.NP+1,e] := 0;
			
			FOR i := 1 TO par.NP DO
			BEGIN
				IF stv.BulkTrapDist[e] > 0 THEN BEGIN
				nt0b[i,e]:=NcLoc[i]*EXP((E_CB[i]-stv.ETrapBulk[e])*Vti);
				pt0b[i,e]:=NcLoc[i]*EXP((stv.ETrapBulk[e]-E_VB[i])*Vti);
				END
				ELSE BEGIN
					nt0b[i,e]:= 0;
					pt0b[i,e]:= 0;
				END;
				IF stv.IntTrapDist[e] > 0 THEN BEGIN
				nt0i[i,e]:=NcLoc[i]*EXP((E_CB[i]-stv.ETrapInt[e])*Vti);
				pt0i[i,e]:=NcLoc[i]*EXP((stv.ETrapInt[e]-E_VB[i])*Vti);
				END
				ELSE BEGIN
					nt0i[i,e]:= 0;
					pt0i[i,e]:= 0;
				END;
			END;
		END
	END
END;

PROCEDURE Calc_Recombination_n_Single_Level(VAR Rn : TRec; dti : myReal; CONSTREF e : INTEGER; CONSTREF n, p, dp, Lan : vector; CONSTREF f_tb, f_ti : TrapArray; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculate all recombination processes and their contribution to the continuity equation for electrons. 
For a derivation see the doc files.}
VAR i, j : INTEGER;
	f_ti_inv_denom, zeros	: vector;
	sum_aj, sum_bj, sum_cj, sum_dj,
	ci, ai_min, ai, ai_plus,
	numer, denom, a, b, c1, g, g0, g1, h, h0, h1 : myReal;
BEGIN	
	{initialize a vector with zeros for all processes that do not contribute}
	FOR i:=0 TO par.NP+1 DO zeros[i]:=0;

	{make sure points on the electrodes are zero}
	{Note: extended reals otherwise get initialised as 'NAN'}
	Rn.direct[0]:=0; Rn.direct[par.NP+1]:=0;
	Rn.bulk[0]:=0; Rn.bulk[par.NP+1]:=0;
	Rn.int[0]:=0; Rn.int[par.NP+1]:=0;
	
	{direct recombination}
	IF e = 1 THEN
		FOR i:=1 TO par.NP DO BEGIN
			Rn.dir_cont_rhs[i]:= (1-dp[i])*Lan[i]*SQR(stv.ni[i]);
			Rn.dir_cont_m[i]:= (1-dp[i])*Lan[i]*p[i];
			Rn.direct[i]:= Rn.dir_cont_m[i]*n[i] - Rn.dir_cont_rhs[i];
		END
	ELSE BEGIN
		Rn.dir_cont_rhs:= zeros;
		Rn.dir_cont_m:= zeros;
		Rn.direct:= zeros;
	END;
	
		
  	{calculate bulk recombination and its linearization in n}
	IF (par.bulk_tr > 0) AND (stv.BulkTrapDist[e] <> 0) THEN
		FOR i:=1 TO par.NP DO
		BEGIN
			IF dti=0 THEN {steady-state}
			BEGIN
				numer:= par.Cn*par.Cp*(n[i]*p[i] - stv.nt0b[i,e]*stv.pt0b[i,e]);
				denom:= (par.Cn*(n[i]+stv.nt0b[i,e]) + par.Cp*(p[i]+stv.pt0b[i,e]));
				Rn.bulk[i]:= stv.Ntb[i,e] * (numer / denom);
				Rn.bulk_cont_m[i]:= stv.Ntb[i,e] * (denom*par.Cn*par.Cp*p[i] - numer*par.Cn) / SQR(denom);
				Rn.bulk_cont_rhs[i]:= Rn.bulk[i] - n[i]*Rn.bulk_cont_m[i];
			END
			ELSE BEGIN {transient}
				{We solve the integral of the ODE we solved for calculating the trap occupance.}
				b:= par.Cn*n[i] + par.Cp*stv.pt0b[i,e];
			    a:= par.Cn*n[i] + par.Cn*stv.nt0b[i,e] + par.Cp*stv.pt0b[i,e] + par.Cp*p[i];
				g:= (par.Cn*n[i]*stv.Ntb[i,e] + par.Cn*stv.nt0b[i,e]*stv.Ntb[i,e]);
				h:= par.Cn*n[i]*stv.Ntb[i,e];
				c1:= (f_tb[i,e] - b/a);
				
				Rn.bulk_cont_m[i]:=0;
				Rn.bulk_cont_rhs[i]:=((h - g*b/a) / dti + g*c1* EXP(-a/dti)/a - g*c1/a)*dti;
				Rn.bulk[i]:= Rn.bulk_cont_rhs[i];			
			END
		END
	ELSE
    BEGIN
		Rn.bulk:= zeros;
		Rn.bulk_cont_m:= zeros;
		Rn.bulk_cont_rhs:= zeros;
	END;	

	{Interface traps, first init to zero:}
	Rn.int:= zeros;
	Rn.int_cont_lo:= zeros;
	Rn.int_cont_m:= zeros;
	Rn.int_cont_up:= zeros;
	Rn.int_cont_rhs:= zeros;
	
	{calculate interface recombination and its linearization in n}
	{for the derivation check docs}
	IF (stv.Traps_int) AND (stv.IntTrapDist[e] <> 0) THEN
		IF dti=0 THEN
		BEGIN {steady-state!}
			f_ti_inv_denom := Calc_Inv_f_ti_Denom(n, p, e, stv, par);
	   
			FOR i:=1 TO par.NP DO
			BEGIN
				IF (stv.Nti[i,e] > 0) THEN
				BEGIN {there are interface traps at this grid point}
					sum_aj := 0;
					sum_bj := 0;
					sum_cj := 0;
					sum_dj := 0;
					FOR j:=-1 TO 1 DO
					BEGIN
						sum_aj := sum_aj + par.Cn * stv.Nti[i+j,e] * n[i+j];
						sum_bj := sum_bj + par.Cp * stv.Nti[i+j,e] * p[i+j];
						sum_cj := sum_cj + par.Cn * stv.Nti[i+j,e] * stv.nt0i[i+j,e];
						sum_dj := sum_dj + par.Cp * stv.Nti[i+j,e] * stv.pt0i[i+j,e];
					END;
					ci := par.Cn * stv.Nti[i,e] * stv.nt0i[i,e];
					ai_min := par.Cn * stv.Nti[i-1,e] ;
					ai := par.Cn * stv.Nti[i,e] ;
					ai_plus := par.Cn * stv.Nti[i+1,e];
				
					{Interface recombination as calculated with the current n and p}
					Rn.int[i] := (ai*n[i] * (sum_bj + sum_cj) - ci * (sum_aj + sum_dj)) * f_ti_inv_denom[i];

					{Calculate the partial derivative of recombination to n[i-1], n[i], n[i+1]}
					Rn.int_cont_lo[i] := -ci * ai_min / f_ti_inv_denom[i];
					Rn.int_cont_lo[i] := Rn.int_cont_lo[i] - ai_min * (ai*n[i] * (sum_bj + sum_cj) - ci * (sum_aj + sum_dj));
					Rn.int_cont_lo[i] := Rn.int_cont_lo[i] * SQR(f_ti_inv_denom[i]);
		
					Rn.int_cont_up[i] := -ci * ai_plus / f_ti_inv_denom[i];
					Rn.int_cont_up[i] := Rn.int_cont_up[i] - ai_plus * (ai*n[i] * (sum_bj + sum_cj) - ci * (sum_aj + sum_dj));
					Rn.int_cont_up[i] := Rn.int_cont_up[i] * SQR(f_ti_inv_denom[i]);
				
					Rn.int_cont_m[i] := (-ci + sum_bj + sum_cj) * ai  / (f_ti_inv_denom[i]);
					Rn.int_cont_m[i] := Rn.int_cont_m[i] - ai * (ai*n[i] * (sum_bj + sum_cj) - ci * (sum_aj + sum_dj));
					Rn.int_cont_m[i] := Rn.int_cont_m[i] * SQR(f_ti_inv_denom[i]);

					{The right hand side of the continuity equation contains the recombination term, but because we linearize in n we add the linearization terms as well.}
					Rn.int_cont_rhs[i] := Rn.int[i] - n[i-1] * Rn.int_cont_lo[i] - n[i+1] * Rn.int_cont_up[i] - n[i] * Rn.int_cont_m[i];
				END
			END; {for loop over grid points}
		END {steady-state!}
		ELSE {so transient}
		BEGIN
			FOR i :=1 TO par.NP DO
			BEGIN
				IF stv.Nti[i-1,e] * stv.Nti[i,e] <> 0 THEN
				BEGIN {now we have just crossed an interface}
					{We solve the integral of the ODE we solved for calculating the trap occupance.}
					a:= par.Cn*n[i] + par.Cn*stv.nt0i[i,e] + par.Cp*stv.pt0i[i,e] + par.Cp*p[i] +
					par.Cn*n[i-1] + par.Cn*stv.nt0i[i-1,e] + par.Cp*stv.pt0i[i-1,e] + par.Cp*p[i-1];
					b:= par.Cn*n[i] + par.Cp*stv.pt0i[i,e] +
					par.Cn*n[i-1] + par.Cp*stv.pt0i[i-1,e];

					g0:= (par.Cn*n[i-1]*stv.Nti[i-1,e] + par.Cn*stv.nt0i[i-1,e]*stv.Nti[i-1,e]);				  
					g1:= (par.Cn*n[i]*stv.Nti[i,e] + par.Cn*stv.nt0i[i,e]*stv.Nti[i,e]);
					h0:= par.Cn*n[i-1]*stv.Nti[i-1,e];
					h1:= par.Cn*n[i]*stv.Nti[i,e];					
					c1:= (f_ti[i,e] - b/a);
					
					Rn.int_cont_rhs[i-1]:=((h0 - g0*b/a) / dti + g0*c1* EXP(-a/dti)/a - g0*c1/a)*dti;
					Rn.int_cont_rhs[i]:=((h1 - g1*b/a) / dti + g1*c1* EXP(-a/dti)/a - g1*c1/a)*dti;

					Rn.int[i-1]:= Rn.int_cont_rhs[i-1];
					Rn.int[i]:= Rn.int_cont_rhs[i];
				END
			END
		END {transient}
END;

PROCEDURE Calc_Recombination_p_Single_Level(VAR Rp : TRec; dti : myReal; CONSTREF e : INTEGER; CONSTREF n, p, dp, Lan : vector; CONSTREF f_tb, f_ti : TrapArray; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculate all recombination processes and their contribution to the continuity equation for holes. 
For a derivation see the doc files.}
VAR i, j : INTEGER;
	f_ti_inv_denom, zeros : vector;
	sum_aj, sum_bj, sum_cj, sum_dj,
	di, bi_min, bi, bi_plus,
	numer, denom, a, b, c1, g, g0, g1, h, h0, h1 : myReal;
BEGIN
	{initialize a vector with zeros for all processes that do not contribute}
	FOR i:=0 TO par.NP+1 DO zeros[i]:=0;

	{make sure points on the electrodes are zero}
	{Note: extended reals otherwise get initialised as 'NAN'}
	Rp.direct[0]:=0; Rp.direct[par.NP+1]:=0;
	Rp.bulk[0]:=0; Rp.bulk[par.NP+1]:=0;
	Rp.int[0]:=0; Rp.int[par.NP+1]:=0;

	{direct recombination}
	IF e = 1 THEN
		FOR i:=1 TO par.NP DO BEGIN
			Rp.dir_cont_rhs[i]:= (1-dp[i])*Lan[i]*SQR(stv.ni[i]);
			Rp.dir_cont_m[i]:= (1-dp[i])*Lan[i]*n[i];
			Rp.direct[i]:= Rp.dir_cont_m[i]*p[i] - Rp.dir_cont_rhs[i];
		END
	ELSE BEGIN
		Rp.dir_cont_rhs:= zeros;
		Rp.dir_cont_m:= zeros;
		Rp.direct:= zeros;
	END;


  	{calculate bulk recombination and its linearization in p}
	IF (par.bulk_tr > 0) AND (stv.BulkTrapDist[e] <> 0) THEN
		FOR i:=1 TO par.NP DO
		BEGIN
			IF dti=0 THEN {steady-state}
			BEGIN
				numer:= par.Cn*par.Cp*(n[i]*p[i] - stv.nt0b[i,e]*stv.pt0b[i,e]);
				denom:= (par.Cn*(n[i]+stv.nt0b[i,e]) + par.Cp*(p[i]+stv.pt0b[i,e]));
				Rp.bulk[i]:= stv.Ntb[i,e] * (numer / denom);
				Rp.bulk_cont_m[i]:= stv.Ntb[i,e] * (denom*par.Cp*par.Cn*n[i] - numer*par.Cp) / SQR(denom);
				Rp.bulk_cont_rhs[i]:= Rp.bulk[i] - p[i]*Rp.bulk_cont_m[i];
			END
			ELSE BEGIN {transient}
				{We solve the integral of the ODE we solved for calculating the trap occupance.}
				b:= par.Cp*p[i] + par.Cn*stv.nt0b[i,e]; 
			    a:= par.Cn*n[i] + par.Cn*stv.nt0b[i,e] + par.Cp*stv.pt0b[i,e] + par.Cp*p[i];
				g:= (par.Cp*p[i]*stv.Ntb[i,e] + par.Cp*stv.pt0b[i,e]*stv.Ntb[i,e]);
				h:= par.Cp*p[i]*stv.Ntb[i,e];
				c1:= (1-f_tb[i,e] - b/a);
				
				Rp.bulk_cont_m[i]:=0;
				Rp.bulk_cont_rhs[i]:= ((h - g*b/a) / dti + g*c1* EXP(-a/dti)/a - g*c1/a)*dti;
				Rp.bulk[i]:= Rp.bulk_cont_rhs[i];		
			END
		END
	ELSE
    BEGIN
		Rp.bulk:= zeros;
		Rp.bulk_cont_m:= zeros;
		Rp.bulk_cont_rhs:= zeros;
	END;	

	{Interface traps, first init to zero:}
	Rp.int:= zeros;
	Rp.int_cont_lo:= zeros;
	Rp.int_cont_m:= zeros;
	Rp.int_cont_up:= zeros;
	Rp.int_cont_rhs:= zeros;
		
	{calculate interface recombination and its linearization in p}
	{for the derivation check Interface_trap_derivation.pdf in the doc files}
	IF (stv.Traps_int) AND (stv.IntTrapDist[e] <> 0) THEN
		IF dti=0 THEN
		BEGIN {steady-state!}
			f_ti_inv_denom := Calc_Inv_f_ti_Denom(n, p, e, stv, par);
	
			FOR i:=1 TO par.NP DO
			BEGIN
				IF (stv.Nti[i,e] > 0) THEN 
				BEGIN {there are interface traps at this grid point}
					sum_aj := 0;
					sum_bj := 0;
					sum_cj := 0;
					sum_dj := 0;
					FOR j:=-1 TO 1 DO
					BEGIN
						sum_aj := sum_aj + par.Cn * stv.Nti[i+j,e] * n[i+j];
						sum_bj := sum_bj + par.Cp * stv.Nti[i+j,e] * p[i+j];
						sum_cj := sum_cj + par.Cn * stv.Nti[i+j,e] * stv.nt0i[i+j,e];
						sum_dj := sum_dj + par.Cp * stv.Nti[i+j,e] * stv.pt0i[i+j,e];
					END;
					di := par.Cp * stv.Nti[i,e] * stv.pt0i[i,e];
					bi_min := par.Cp * stv.Nti[i-1,e];
					bi := par.Cp * stv.Nti[i,e];
					bi_plus := par.Cp * stv.Nti[i+1,e] ;

					{Interface recombination as calculated with the current n and p}
					Rp.int[i] := (bi * p[i] * (sum_aj + sum_dj) - di * (sum_bj + sum_cj)) / (sum_aj+sum_bj+sum_cj+sum_dj);
	
					{Calculate the partial derivative of recombination to p[i-1], p[i], p[i+1]}
					Rp.int_cont_lo[i] := -di * bi_min / f_ti_inv_denom[i];
					Rp.int_cont_lo[i] := Rp.int_cont_lo[i] - bi_min * (bi * p[i] * (sum_aj + sum_dj) - di * (sum_bj + sum_cj));
					Rp.int_cont_lo[i] := Rp.int_cont_lo[i] * SQR(f_ti_inv_denom[i]);

					Rp.int_cont_up[i] := -di * bi_plus / f_ti_inv_denom[i];
					Rp.int_cont_up[i] := Rp.int_cont_up[i] - bi_plus * (bi * p[i] * (sum_aj + sum_dj) - di * (sum_bj + sum_cj));
					Rp.int_cont_up[i] := Rp.int_cont_up[i] * SQR(f_ti_inv_denom[i]);

					Rp.int_cont_m[i] := (-di + sum_aj + sum_dj) * bi  / (f_ti_inv_denom[i]);
					Rp.int_cont_m[i] := Rp.int_cont_m[i] - bi * (bi * p[i] * (sum_aj + sum_dj) - di * (sum_bj + sum_cj));
					Rp.int_cont_m[i] := Rp.int_cont_m[i] * SQR(f_ti_inv_denom[i]);

					{The right hand side of the continuity equation contains the recombination term, but because we linearize in n we add the linearization terms as well.}			
					Rp.int_cont_rhs[i] := Rp.int[i] - p[i-1] * Rp.int_cont_lo[i] - p[i+1] * Rp.int_cont_up[i] - p[i] * Rp.int_cont_m[i];	    
				END
			END; {for loop over grid points}
		END {steady-state!}
		ELSE {so transient}
		BEGIN
			FOR i := 1 TO par.NP DO
			BEGIN
				IF stv.Nti[i-1,e] * stv.Nti[i,e] <> 0 THEN
				BEGIN {now we have just crossed an interface}
				{We solve the integral of the ODE we solved for calculating the trap occupance.}
					a:= par.Cn*n[i] + par.Cn*stv.nt0i[i,e] + par.Cp*stv.pt0i[i,e] + par.Cp*p[i] +
					par.Cn*n[i-1] + par.Cn*stv.nt0i[i-1,e] + par.Cp*stv.pt0i[i-1,e] + par.Cp*p[i-1];
					b:= par.Cp*p[i] + par.Cn*stv.nt0i[i,e] +
					par.Cp*p[i-1] + par.Cn*stv.nt0i[i-1,e];
					
					g0:= (par.Cp*p[i-1]*stv.Nti[i-1,e] + par.Cp*stv.pt0i[i-1,e]*stv.Nti[i-1,e]);
					g1:= (par.Cp*p[i]*stv.Nti[i,e] + par.Cp*stv.pt0i[i,e]*stv.Nti[i,e]);
					h0:= par.Cp*p[i-1]*stv.Nti[i-1,e];
					h1:= par.Cp*p[i]*stv.Nti[i,e];
					c1:= (1-f_ti[i,e] - b/a);
					
					Rp.int_cont_rhs[i-1]:=((h0 - g0*b/a) / dti + g0*c1* EXP(-a/dti)/a - g0*c1/a)*dti;
					Rp.int_cont_rhs[i]:=((h1 - g1*b/a) / dti + g1*c1* EXP(-a/dti)/a - g1*c1/a)*dti;
					Rp.int[i-1]:= Rp.int_cont_rhs[i-1];
					Rp.int[i]:= Rp.int_cont_rhs[i];
				END
			END
		END {transient}
END;

PROCEDURE Calc_Recombination_n(VAR Rn : TRec; dti : myReal; CONSTREF n, p, dp, Lan : vector; f_tb, f_ti : TrapArray; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculate all recombination processes and their contribution to the continuity equation for electrons. 
For a derivation see the doc files.}
VAR i, e : INTEGER;
BEGIN
	FOR e:=1 TO stv.N_Etr DO
	BEGIN
		{Note: we use global variable RecDum : TRec as a dummy variable. It is too larger to be a local variable.}
		
		Calc_Recombination_n_Single_Level(RecDum, dti, e, n, p, dp, Lan, f_tb, f_ti, stv, par);
		IF e=1 THEN
			FOR i:=1 TO par.NP DO
			BEGIN
				Rn.direct[i]:= RecDum.direct[i];
				Rn.bulk[i]:= RecDum.bulk[i];
				Rn.int[i] := RecDum.int[i];
				Rn.dir_cont_rhs[i]:= RecDum.dir_cont_rhs[i];
				Rn.dir_cont_m[i]:= RecDum.dir_cont_m[i];
				Rn.bulk_cont_rhs[i]:= RecDum.bulk_cont_rhs[i];
				Rn.bulk_cont_m[i]:= RecDum.bulk_cont_m[i];
				Rn.int_cont_lo[i]:= RecDum.int_cont_lo[i];
				Rn.int_cont_up[i]:= RecDum.int_cont_up[i];
				Rn.int_cont_m[i]:= RecDum.int_cont_m[i];
				Rn.int_cont_rhs[i]:= RecDum.int_cont_rhs[i];
			END
		ELSE
		BEGIN
			FOR i:=1 TO par.NP DO
			BEGIN			
				Rn.direct[i]:= Rn.direct[i] + RecDum.direct[i];
				Rn.bulk[i]:= Rn.bulk[i] + RecDum.bulk[i];
				Rn.int[i] := Rn.int[i] + RecDum.int[i];
				Rn.dir_cont_rhs[i]:= Rn.dir_cont_rhs[i] + RecDum.dir_cont_rhs[i];
				Rn.dir_cont_m[i]:= Rn.dir_cont_m[i] + RecDum.dir_cont_m[i];
				Rn.bulk_cont_rhs[i]:= Rn.bulk_cont_rhs[i] + RecDum.bulk_cont_rhs[i];
				Rn.bulk_cont_m[i]:= Rn.bulk_cont_m[i] + RecDum.bulk_cont_m[i];
				Rn.int_cont_lo[i]:= Rn.int_cont_lo[i] + RecDum.int_cont_lo[i];
				Rn.int_cont_up[i]:= Rn.int_cont_up[i] + RecDum.int_cont_up[i];
				Rn.int_cont_m[i]:= Rn.int_cont_m[i] + RecDum.int_cont_m[i];
				Rn.int_cont_rhs[i]:= Rn.int_cont_rhs[i] + RecDum.int_cont_rhs[i];
			END
		END;
	END;
END;

PROCEDURE Calc_Recombination_p(VAR Rp : TRec; dti : myReal; CONSTREF n, p, dp, Lan : vector; f_tb, f_ti : TrapArray; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculate all recombination processes and their contribution to the continuity equation for electrons. 
For a derivation see the doc files.}
VAR i, e : INTEGER;
BEGIN
	FOR e:=1 TO stv.N_Etr DO
	BEGIN
		{Note: we use global variable RecDum : TRec as a dummy variable. It is too larger to be a local variable.}
		Calc_Recombination_p_Single_Level(RecDum, dti, e, n, p, dp, Lan, f_tb, f_ti, stv, par);
		IF e=1 THEN
			FOR i:=1 TO par.NP DO
			BEGIN
				Rp.direct[i]:= RecDum.direct[i];
				Rp.bulk[i]:= RecDum.bulk[i];
				Rp.int[i] := RecDum.int[i];
				Rp.dir_cont_rhs[i]:= RecDum.dir_cont_rhs[i];
				Rp.dir_cont_m[i]:= RecDum.dir_cont_m[i];
				Rp.bulk_cont_rhs[i]:= RecDum.bulk_cont_rhs[i];
				Rp.bulk_cont_m[i]:= RecDum.bulk_cont_m[i];
				Rp.int_cont_lo[i]:= RecDum.int_cont_lo[i];
				Rp.int_cont_up[i]:= RecDum.int_cont_up[i];
				Rp.int_cont_m[i]:= RecDum.int_cont_m[i];
				Rp.int_cont_rhs[i]:= RecDum.int_cont_rhs[i];
			END
		ELSE
		BEGIN
			FOR i:=1 TO par.NP DO
			BEGIN			
				Rp.direct[i]:= Rp.direct[i] + RecDum.direct[i];
				Rp.bulk[i]:= Rp.bulk[i] + RecDum.bulk[i];
				Rp.int[i] := Rp.int[i] + RecDum.int[i];
				Rp.dir_cont_rhs[i]:= Rp.dir_cont_rhs[i] + RecDum.dir_cont_rhs[i];
				Rp.dir_cont_m[i]:= Rp.dir_cont_m[i] + RecDum.dir_cont_m[i];
				Rp.bulk_cont_rhs[i]:= Rp.bulk_cont_rhs[i] + RecDum.bulk_cont_rhs[i];
				Rp.bulk_cont_m[i]:= Rp.bulk_cont_m[i] + RecDum.bulk_cont_m[i];
				Rp.int_cont_lo[i]:= Rp.int_cont_lo[i] + RecDum.int_cont_lo[i];
				Rp.int_cont_up[i]:= Rp.int_cont_up[i] + RecDum.int_cont_up[i];
				Rp.int_cont_m[i]:= Rp.int_cont_m[i] + RecDum.int_cont_m[i];
				Rp.int_cont_rhs[i]:= Rp.int_cont_rhs[i] + RecDum.int_cont_rhs[i];
			END
		END;
	END;
END;

PROCEDURE Cont_Eq_Elec(VAR n : vector; nPrevTime, V, Jn, p, mu, g, Lan, dp : vector; VAR f_tb, f_ti : TrapArray; VAR Rn : TRec;
				CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters;  dti : myReal = 0);
{dti: DeltaTInverse, so dti = 1/delt. If dti=0 then we're using the steady-state equations
(Selberherr 6-1.73) as they correspond to an infinite timestep.
If dti > 0 we're using the transient equations (Selberherr 6-4.32)}
VAR i, NP		  : INTEGER;
	nmin, fac	  : myReal;
	ResetDens	  : BOOLEAN;
    lo, m, u, rhs : vector;
BEGIN
	NP:=par.NP; {copy number of grid points into local variable NP}
	
	{Set the boundary conditions for finite surface recombination on the left side of the device}
	IF (par.Sn_L < 0) OR (Jn[0]>0) THEN
	BEGIN
		{Infinite surface recombination}
		lo[0]:= 0;
		m[0]:=1;
		u[0]:=0;
		rhs[0]:=stv.NcLoc[0] * EXP(-(par.W_L - stv.E_CB[0])*stv.Vti);
	END
	ELSE
	BEGIN
		{Finite surface recombination}
		lo[0]:= 0;
		m[0]:=-mu[0]*stv.Vt*B((V[0]-V[1])*stv.Vti)/(par.L*stv.h[0]) - par.Sn_L;
		u[0]:=mu[0]*stv.Vt*B((V[1]-V[0])*stv.Vti)/(par.L*stv.h[0]);
		rhs[0]:= - par.Sn_L*stv.NcLoc[0] * EXP((stv.E_CB[0] - par.W_L)*stv.Vti);
	END;

	{Calculate recombination and its contribution to the continuity equation} 
	Calc_Recombination_n(Rn, dti, n, p, dp, Lan, f_tb, f_ti, stv, par);

    {now set the interior part:}
    FOR i:=1 TO NP DO  {continuity eq. in matrix vorm}
		WITH stv DO {ni and h are static variables!}
		BEGIN
			fac := 0.5*SQR(par.L)*h[i]*h[i-1]*(h[i]+h[i-1]); {repeats often in the equations}
			
			rhs[i]:=- fac * (g[i] +nPrevTime[i]*dti)
        			- fac*Rn.dir_cont_rhs[i] {direct / Langevin recombination}
		        	+ fac*Rn.bulk_cont_rhs[i] {the part of Rn_bulk that does not depend on n_new}
			        + fac*Rn.int_cont_rhs[i]; {the part of Rn_int that does not depend on n_new}
			
			lo[i]:=  h[i]*mu[i-1]*stv.Vt*B((V[i-1]-V[i])*stv.Vti) +
         		   - fac*Rn.int_cont_lo[i]; {the part of Rn_int that depends on n[i-1]}
			
			m[i]:=- (h[i-1]*mu[i]*stv.Vt*B((V[i]-V[i+1])*stv.Vti) +
				    h[i]*mu[i-1]*stv.Vt*B((V[i]-V[i-1])*stv.Vti))
				  - fac*dti
			      - fac*Rn.dir_cont_m[i] {direct / Langevin recombination}
			      - fac*Rn.bulk_cont_m[i] {the part of Rn_bulk that depends on n[i]}
                  - fac*Rn.int_cont_m[i]; {the part of Rn_int that depends on n[i]}
			
			u[i]:=  h[i-1]*mu[i]*stv.Vt*B((V[i+1]-V[i])*stv.Vti)
			      - fac*Rn.int_cont_up[i]; {the part of Rn_int that depends on n[i+1]}
		END;

  	{Set the boundary conditions on the right side of the device}
	IF (par.Sn_R < 0) OR (Jn[NP]>0) THEN
	BEGIN
		{Infinite surface recombination}
		lo[NP+1]:= 0;
		m[NP+1]:= 1;
		u[NP+1]:= 0;
		rhs[NP+1]:=stv.NcLoc[NP+1] * EXP(-(par.W_R - stv.E_CB[NP+1])*stv.Vti);
	END
	ELSE
	BEGIN
		{Finite surface recombination}
		lo[NP+1]:=-mu[NP]*stv.Vt * B((V[NP]-V[NP+1])*stv.Vti)/(par.L*stv.h[NP]);
		m[NP+1]:=mu[NP]*stv.Vt * B((V[NP+1]-V[NP])*stv.Vti)/(par.L*stv.h[NP]) - par.Sn_R;
		u[NP+1]:=0;
		rhs[NP+1]:=-par.Sn_R*stv.NcLoc[NP+1] * EXP((stv.E_CB[NP+1] - par.W_R)*stv.Vti);
	END;

    Tridiag(n, lo, m, u, rhs, 0, NP+1); {Solve for the new electron densities}

	{now check if n is still well behaved, i.e. positive.
	They're sometimes negative, for example close to the cathode (and if recombination is strong)
	the very large electron density will make the hole density very small, virtually zero.
	Due to finite numerical accuracy, the hole density can then get negative (but still |p| is small).
	So, we reset them to a small (pmin) positive value.}

	{find minimal values of n at the contacts:}
	nmin:=MIN(n[0], n[NP+1]);
	ResetDens:=FALSE;
	FOR i:=1 TO NP DO
		IF (n[i]<=0)
		THEN BEGIN
			ResetDens:=TRUE;
			n[i]:=nmin
		END;

	IF ResetDens AND (NOT par.IgnoreNegDens) THEN Stop_Prog('Negative electron concentration encountered!' , EC_NumericalFailure)
END;


PROCEDURE Cont_Eq_Holes(VAR p : vector; pPrevTime, V, Jp, n, mu, g, Lan, dp : vector; VAR f_tb, f_ti : TrapArray; VAR Rp : TRec;
				CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; dti : myReal = 0);
{dti: DeltaTInverse, so dti = 1/delt. If dti=0 then we're using the steady-state equations
(Selberherr 6-1.74) as they correspond to an infinite timestep.
If dti > 0 we're using the transient equations (Selberherr 6-4.33)}
VAR i, NP : INTEGER;
	pmin, fac : myReal;
	ResetDens : BOOLEAN;
    lo, m, u, rhs : vector;
BEGIN
    NP:=par.NP; {copy number of grid points into local variable NP}

    {Set the boundary conditions for finite surface recombination on the left side of the device}
	IF (par.Sp_L < 0) OR (Jp[0]>0) THEN
	BEGIN
		{Infinite surface recombination}
		lo[0]:= 0;
		m[0]:=  1;
		u[0]:= 0; 
		rhs[0]:=stv.NcLoc[0] *EXP((par.W_L - stv.E_VB[0])*stv.Vti);
	END
	ELSE
	BEGIN
		{Finite surface recombination}
		lo[0]:= 0;
		m[0]:=mu[0]*stv.Vt * B((V[1]-V[0])*stv.Vti)/(par.L*stv.h[0]) - par.Sp_L;
		u[0]:=-mu[0]*stv.Vt * B((V[0]-V[1])*stv.Vti)/(par.L*stv.h[0]);
		rhs[0]:=-par.Sp_L*stv.NcLoc[0] * EXP((par.W_L - stv.E_VB[0])*stv.Vti);
	END;

	{Calculate recombination and its contribution to the continuity equation} 	
	Calc_Recombination_p(Rp, dti, n, p, dp, Lan, f_tb, f_ti, stv, par);

    {now do the interior points:}
    FOR i:=1 TO NP DO  {continuity eq. in matrix vorm}
		WITH stv DO {ni and h are static variables!}
		BEGIN
			fac := 0.5*SQR(par.L)*h[i]*h[i-1]*(h[i]+h[i-1]); {repeats often in the equations}
			
			rhs[i]:=- fac * (g[i] +pPrevTime[i]*dti)
			        - fac*Rp.dir_cont_rhs[i] {direct / Langevin recombination}
			        + fac*Rp.bulk_cont_rhs[i] {the part of Rp_bulk that is independent of p_new}
         			+ fac*Rp.int_cont_rhs[i]; {the part of Rp_int that is independent of p_new}

			lo[i]:=  h[i]*mu[i-1]*stv.Vt*B((V[i]-V[i-1])*stv.Vti)
			       - fac*Rp.int_cont_lo[i]; {the part of Rp_int that depends on p[i-1]}
			
			m[i]:=- (h[i-1]*mu[i]*stv.Vt*B((V[i+1]-V[i])*stv.Vti) +
				    h[i]*mu[i-1]*stv.Vt*B((V[i-1]-V[i])*stv.Vti))
				  - fac*dti
			      - fac*Rp.dir_cont_m[i] {direct / Langevin recombination}
			      - fac*Rp.bulk_cont_m[i] {part of Rp_bulk that depends on p[i]}
			      - fac*Rp.int_cont_m[i]; {the part of Rp_int that depends on p[i]}

			
			u[i]:=  h[i-1]*mu[i]*stv.Vt*B((V[i]-V[i+1])*stv.Vti)
			      - fac*Rp.int_cont_up[i]; {the part of Rp_int that depends on p[i+1]}
		END;

    {Set the boundary conditions on the right side of the device}
	IF (par.Sp_R < 0) OR (Jp[NP]>0) THEN
	BEGIN
		{Infinite surface recombination}
		lo[NP+1]:= 0;
		m[NP+1]:= 1;
		u[NP+1]:= 0;
		rhs[NP+1]:=stv.NcLoc[NP+1] * EXP(-(stv.E_VB[NP+1]-par.W_R)*stv.Vti);
	END
	ELSE
	BEGIN
		{Finite surface recombination}
		lo[NP+1]:=mu[NP]*stv.Vt*B((V[NP+1]-V[NP])*stv.Vti)/(par.L*stv.h[NP]);
		m[NP+1]:=-mu[NP]*stv.Vt*B((V[NP]-V[NP+1])*stv.Vti)/(par.L*stv.h[NP]) - par.Sp_R;
		u[NP+1]:=0;
		rhs[NP+1]:=-par.Sp_R*stv.NcLoc[NP+1] * EXP((par.W_R - stv.E_VB[NP+1])*stv.Vti);
	END;

    {Solve for all grid points, so including i=0 and i=NP+1:}
    Tridiag(p, lo, m, u, rhs, 0, NP+1); {Solve for the new hole densities}

    {now check if p is still well behaved, i.e. positive.
	They're sometimes negative, for example close to the cathode (and if recombination is strong)
	the very large electron density will make the hole density very small, virtually zero.
	Due to finite numerical accuracy, the hole density can then get negative (but still |p| is small).
	So, we reset them to a small (pmin) positive value.}

	{find minimal values of n at the contacts:}
	pmin:=MIN(p[0], p[NP+1]);
	ResetDens:=FALSE;
	FOR i:=1 TO NP DO
		IF (p[i]<=0)
		THEN BEGIN
			ResetDens:=TRUE;
			p[i]:=pmin
		END;

	IF ResetDens AND (NOT par.IgnoreNegDens) THEN Stop_Prog('Negative hole concentration encountered!', EC_NumericalFailure)

END;

PROCEDURE Calc_Displacement_Curr(VAR JD : vector; V, VPrevTime : vector; dti : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure calculates the displacement current}
VAR i : INTEGER;
BEGIN
	FOR i:=0 TO par.NP DO
		JD[i]:=stv.eps[i] * (V[i+1]-V[i]-VPrevTime[i+1]+VPrevTime[i]) * dti / (par.L*stv.h[i]);
	JD[par.NP+1]:=JD[par.NP]; {doesn't have a physical meaning though}
END;

PROCEDURE Calc_Curr_Diff(sn : ShortInt; istart, ifinish : INTEGER; VAR J : vector; V, dens, mu, Rint : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure calculates the current density in differential form, see Selberherr eq. 6.1-39 or 6.1-41}
{sn denotes the sign of the carrier, so -1 for electrons, +1 for holes}
VAR i, e	  : INTEGER;
	int_traps : BOOLEAN;
BEGIN
	IF ABS(sn)<>1 THEN Stop_Prog('Incorrect sn passed to Calc_Curr_Diff', EC_ProgrammingError);
	IF (istart<0) OR (ifinish>=par.NP+1) THEN Stop_Prog('Incorrect istart and/or ifinish passed to Calc_Curr_Diff.', EC_ProgrammingError);

	{the current is only non-zero between istart and ifinish, so first init the zero's:}
	FOR i:=0 TO istart-1 DO J[i]:=0;
	FOR i:=ifinish+1 TO par.NP DO J[i]:=0;
	
	{now actually calc the current:}
	FOR i:=istart TO ifinish DO
	BEGIN
		J[i]:=sn*q*mu[i]*stv.Vt*(dens[i+1]*B(sn*(V[i]-V[i+1])*stv.Vti) - dens[i]*B(sn*(V[i+1]-V[i])*stv.Vti))/(par.L*stv.h[i]);
		int_traps := FALSE;
		FOR e:=1 TO stv.N_Etr DO
			IF (stv.Nti[i,e] <> 0) AND (stv.Nti[i+1,e] <> 0) THEN int_traps:= TRUE;
		IF int_traps THEN
			J[i]:=J[i] - sn*0.5*q*par.L*stv.h[i]*(Rint[i] - Rint[i+1]);
	END;

	{last point as this is in the output (Var_file)}
	J[par.NP+1]:=J[par.NP]; {doesn't have a physical meaning: J[NP+1] is current between NP+1 and NP+2}
END;

PROCEDURE Calc_Curr_Int(sn : ShortInt; istart, ifinish : INTEGER; dti : myReal; VAR J : vector; V, dens, olddens, mu, g : vector; 
						CONSTREF Rec : TRec; CONSTREF stv : TStaticVars; 
						CONSTREF par : TInputParameters);
{Calculates the current density in integral form, see De Mari, solid-state elec. vol 11 p.33 (68) eq. 15}
{istart and ifinish exclude the electrodes, so istart>=1, ifinish<=NP}
{sn denotes the sign of the carrier, so -1 for electrons, +1 for holes}
VAR i, e : INTEGER;
    K, single_int, double_int, int_U, int_U_old, dx : myReal;
    U : vector;
	int_traps : BOOLEAN;
BEGIN
    {first check a few things:}
    IF ABS(sn)<>1 THEN Stop_Prog('Incorrect sn passed to Calc_Curr_Int', EC_ProgrammingError);
    IF (istart<0) OR (ifinish>=par.NP+1) THEN Stop_Prog('Incorrect istart and/or ifinish passed to Calc_Curr_Int.', EC_ProgrammingError);

    single_int:=0;
    double_int:=0;
    int_U:=0;

    {the current is only non-zero between istart and ifinish, so first init the zero's:}
    FOR i:=0 TO istart-1 DO J[i]:=0;
    FOR i:=ifinish+1 TO par.NP DO J[i]:=0;

    {now we compute the various integrals. Note: some variables are on-grid (V, dens, Rnet.direct/bulk), but 
    others are not: Rnet.int is defined between 2 grid points. We approximate the integral in either case 
    by value(grid point i) * grid spacing between i and i+1. This works as we're using a grid with a uniform
    spacing near the interfaces (see Make_Grid)} 
	
    i:=istart;
    U[i]:=Rec.direct[i] + Rec.bulk[i] + Rec.int[i] - g[i] {recombination - generation in grid point i}
			    + dti * (dens[i]-olddens[i]); {change in density also contributes}  
   
    {now we do trapezoidal integration for single and double int:}
    FOR i:=istart+1 TO ifinish DO
    BEGIN
        U[i]:=Rec.direct[i] + Rec.bulk[i] + Rec.int[i] - g[i] {recombination - generation in grid point i}
				   + dti * (dens[i]-olddens[i]); {change in density also contributes}  
        dx:=stv.x[i]-stv.x[i-1]; {grid spacing between x[i] and x[i-1]}
        int_U_old:=int_U; {keep the old one}
        int_U:=int_U + 0.5*dx*(U[i-1] + U[i]);
        single_int:=single_int + 0.5*dx*(EXP(sn*V[i-1]*stv.Vti)/mu[i-1] + EXP(sn*V[i]*stv.Vti)/mu[i]);
        double_int:=double_int + 0.5*dx*(EXP(sn*V[i-1]*stv.Vti)*int_U_old/mu[i-1] + EXP(sn*V[i]*stv.Vti)*int_U/mu[i])    
    END;
    K:=(stv.Vt*(sn*dens[ifinish+1]*EXP(sn*V[ifinish+1]*stv.Vti) - sn*dens[istart]*EXP(sn*V[istart]*stv.Vti)) 
		- sn*double_int)/single_int;

    J[istart]:=q*K;
    FOR i:=istart+1 TO ifinish DO
        J[i]:=J[i-1] + sn*q*par.L*stv.h[i]*U[i];
	
	{now correct for interface recombination. This represents the current THROUGH the interface traps.}
	FOR i:=1 TO par.NP+1 DO 
	BEGIN
		int_traps := FALSE;
		FOR e:=1 TO stv.N_Etr DO
			IF (stv.Nti[i,e] <> 0) AND (stv.Nti[i+1,e] <> 0) THEN int_traps:= TRUE;
		IF int_traps THEN
			J[i]:=J[i-1] + sn*0.5*q*par.L*stv.h[i]*(Rec.int[i] + Rec.int[i+1])
	END;
	
    J[par.NP+1]:=J[par.NP]; {doesn't have a physical meaning: J[NP+1] is current between NP+1 and NP+2}
END;

PROCEDURE Calc_All_Currents(VAR new : TState; CONSTREF curr : TState; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters); 
{calculates all the currents for state new}
VAR mob, Jsum, GenDum : vector;
	i, istart, ifinish : INTEGER;
BEGIN
	WITH new DO 
	BEGIN
		{Note: we use global variable RecDum : TRec as a dummy variable. It is too larger to be a local variable.}
		
		{first do the electrons}
		CASE par.CurrDiffInt OF
			1 : Calc_Curr_Diff(-1, 0, par.NP, Jn, Vgn, n, mun, Rn.int, stv, par); {only needs part of Rn with interface recombination}
			2 : Calc_Curr_Int(-1, 0, par.NP, dti, Jn, Vgn, n, curr.n, mun, gen, Rn, stv, par); {needs full Rn and curr.n}
		END;	

		{now do the holes}
		CASE par.CurrDiffInt OF
			1 : Calc_Curr_Diff(1, 0, par.NP, Jp, Vgp, p, mup, Rp.int, stv, par);{only needs part of Rp with interface recombination}
			2 : Calc_Curr_Int(1, 0, par.NP, dti, Jp, Vgp, p, curr.p, mup, gen, Rp, stv, par); {needs full Rp and curr.p}
		END;	

		{ions can be limited to the middle layer (which can take up the whole device i=0...NP+1.}
		{we calc the current between points istart and ifinish:}
		IF (stv.i1>0) AND (NOT par.IonsInTLs) THEN istart:=stv.i1+1 ELSE istart:=0;
		IF (stv.i2<par.NP+1) AND (NOT par.IonsInTLs) THEN ifinish:=stv.i2-2 ELSE ifinish:=par.NP;

		FILLCHAR(RecDum, SIZEOF(RecDum), 0); {set all fields of RecDum to zero as ions don't have generation/recombination}
		FILLCHAR(GenDum, SIZEOF(GenDum), 0); {likewise for the generation rate of ions}

		IF par.negIonsMove AND (dti<>0) THEN {dti=0, then we're in steady-state => ionic currents are zero!}
		BEGIN 
			FOR i:=0 TO par.NP+1 DO mob[i]:=par.mobnion;
			CASE par.CurrDiffInt OF
				1 : Calc_Curr_Diff(-1, istart, ifinish, Jnion, V, nion, mob, RecDum.int, stv, par);		
				2 : Calc_Curr_Int(-1, istart, ifinish, dti, Jnion, V, nion, curr.nion, mob, GenDum, RecDum, stv, par);
			END	
		END
		ELSE FILLCHAR(Jnion, SIZEOF(Jnion), 0); {set ionic current to zero}
	
		IF par.posIonsMove AND (dti<>0) THEN {dti=0, then we're in steady-state => ionic currents are zero!}
		BEGIN 
			FOR i:=0 TO par.NP+1 DO mob[i]:=par.mobpion;
			CASE par.CurrDiffInt OF
				1 : Calc_Curr_Diff(1, istart, ifinish, Jpion, V, pion, mob, RecDum.int, stv, par);
				2 : Calc_Curr_Int(1, istart, ifinish, dti, Jpion, V, pion, curr.pion, mob, GenDum, RecDum, stv, par);
			END
		END
		ELSE FILLCHAR(Jpion, SIZEOF(Jpion), 0); {set ionic current to zero}
		
		{lastly, the displacement current:}
		IF (dti<>0)
			THEN Calc_Displacement_Curr(JD, V, curr.V, dti, stv, par) {calc. displacement current}
			ELSE FILLCHAR(JD, SIZEOF(JD), 0); {if dti=0 => steady-state, so zero.}
		
		{now compute the total summed current:}
		FOR i:=0 TO par.NP+1 DO
			Jsum[i]:=Jn[i] + Jp[i] + Jnion[i] + Jpion[i] + JD[i];
			
		{calc the total interal current density}
		Jint:=Average(Jsum, stv.h, 0, par.NP+1);	
		{and its rms error:}
		errJ:=Calc_RMS_Error(Jsum, Jint, 0, par.NP+1)

	END;
END; 

PROCEDURE Extrapolate_Solution(CONSTREF prev, curr : TState; VAR new : TState; AccSols : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Obtains a guess for state new based on an extrapolation of states prev and curr}
VAR i : INTEGER;
	rV, rt : myReal;
BEGIN
	{what we will doing depends on the number of accepted solutions (AccSols)}
	
	CASE AccSols OF
		0 : ; {here we might add the initial guess for the first solution}
		1 : IF new.Vint <> curr.Vint THEN {all we can do is extrapolate the electric field}
				FOR i:=1 TO par.NP DO {only loop over interior points, not the contacts}
					new.V[i]:=curr.V[i] + 1*(stv.x[i]/par.L - 0.5)*(new.Vint-curr.Vint);
	OTHERWISE BEGIN
		{so at least 2 accepted solutions are available}
		{first calculate the factors (rV,G,t) by which the voltage, generation and time have increased:}
		IF curr.Vint <> prev.Vint THEN
			rV:=(new.Vint - curr.Vint)/(curr.Vint - prev.Vint)
		ELSE rV:=0;
				IF curr.tijd <> prev.tijd THEN
			rt:=(new.tijd - curr.tijd)/(curr.tijd - prev.tijd)
		ELSE rt:=0; {In SimSS this will always be zero, in ZimT mostly not}
		IF (rV<>0) AND (curr.Gehp = prev.Gehp) THEN 
			FOR i:=1 TO par.NP DO {only loop over interior points, not the contacts}
			BEGIN
				{use linear extrapolation. Simple and it works!}
				new.V[i]:=curr.V[i] + rV*(curr.V[i]-prev.V[i]); 

				{for the densities, we have to ensure that they are: 0< n <= NcLoc:}
				{If n<=0, then we set n to MinFloat, the smallest positive real:}
				new.n[i]:=MAX(MinFloat, MIN(stv.NcLoc[i], curr.n[i] + rV*(curr.n[i]-prev.n[i])));
				new.p[i]:=MAX(MinFloat, MIN(stv.NcLoc[i], curr.p[i] + rV*(curr.p[i]-prev.p[i])));

				{now IF there are moving ions:}
				IF new.UpdateIons THEN BEGIN
					IF par.negIonsMove THEN	new.nion[i]:=MAX(0, curr.nion[i] + rV*(curr.nion[i]-prev.nion[i])); {these need to be non-negative}
					IF par.posIonsMove THEN new.pion[i]:=MAX(0, curr.pion[i] + rV*(curr.pion[i]-prev.pion[i]))
				END; {ions}
			
			END {for loop interior grid points}
		
		ELSE {rV=0 and/or (curr.Gehp<>prev.Gehp). In either case, simply extrapolate using the time difference}
			FOR i:=1 TO par.NP DO {only loop over interior points, not the contacts}
			BEGIN
				{use linear extrapolation. Simple and it works!}
				new.V[i]:=curr.V[i] + rt*(curr.V[i]-prev.V[i]);
			
				{for elec/holes and ions, we essentially extrapolate the quasi-fermi level:}
				{for the densities, we have to ensure that they are: 0< n <= NcLoc:}
				new.n[i]:=MIN(stv.NcLoc[i], curr.n[i]*power(curr.n[i]/prev.n[i], rt));
				new.p[i]:=MIN(stv.NcLoc[i], curr.p[i]*power(curr.p[i]/prev.p[i], rt));
			
				{now IF there are moving ions:}
				IF new.UpdateIons THEN BEGIN
					IF par.negIonsMove AND (prev.nion[i]<>0) THEN new.nion[i]:=curr.nion[i]*power(abs(curr.nion[i]/prev.nion[i]), rt); {these need to be non-negative}
					IF par.posIonsMove AND (prev.nion[i]<>0) THEN new.pion[i]:=curr.pion[i]*power(abs(curr.pion[i]/prev.pion[i]), rt)
				END; {ions}
		
			END {for loop interior grid points}
	
	END {otherwise}
	END {case statement}	
	
END;

FUNCTION Deterimine_Convergence_Densities(CONSTREF deln, delp, delnion, delpion, n, p, nion, pion : vector; 
							UpdateIons : BOOLEAN; accDens : myReal; VAR ConvMsg : STRING; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : BOOLEAN;
{Determines whether the densities have converged}
VAR totRelChange : myReal;
BEGIN
	totRelChange:=0; {total relative change of all charges}
	
	{now compute relative changes of densities:}
	totRelChange:=totRelChange + Norm_Eucl(deln, 0, par.NP+1) / Norm_Eucl(n, 0, par.NP+1);
	totRelChange:=totRelChange + Norm_Eucl(delp, 0, par.NP+1) / Norm_Eucl(p, 0, par.NP+1);
	
	{for ions, only if they move (note: if they move, their densities is non-zero!}
	IF par.negIonsMove THEN
		totRelChange:=totRelChange + Norm_Eucl(delnion, 0, par.NP+1) / Norm_Eucl(nion, 0, par.NP+1);
	IF par.posIonsMove THEN
		totRelChange:=totRelChange + Norm_Eucl(delpion, 0, par.NP+1) / Norm_Eucl(pion, 0, par.NP+1);

	ConvMsg:='- relative change of densities: ' + FloatToStrF(totRelChange, ffExponent,5,0) + LineEnding; {our message string}

	{so, converged if change is small enough, we take the acceleration into accout as well:}
	Deterimine_Convergence_Densities:=totRelChange <= accDens*par.tolDens
END;

PROCEDURE Main_Solver(VAR curr, new : TState; VAR it : INTEGER; VAR conv : BOOLEAN; VAR StatusStr : ANSISTRING; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Iteratively solves the Poisson and continuity equations, including traps and ions}
{can be used in steady-state and transient cases}
VAR i, MaxIt : INTEGER;
	check_Poisson, coupleIonsPoisson, convDensities : BOOLEAN;
	oldn, oldp, oldnion, oldpion, deln, delp, delnion, delpion : vector; {we need this to monitor the iteration loops}
	PoissMsg, ConvMsg : STRING;
	accDens : myReal;
BEGIN
	{apply new bias to electrodes:}
	new.V[0]:=stv.V0 - 0.5 *new.Vint;
	new.V[par.NP+1]:=stv.VL + 0.5 *new.Vint;
	Update_Gen_Pot(new.V, new.Vgn, new.Vgp, stv, par); {init. generalised potentials}
	it:=0;

	IF new.dti=0 {max number of iterations depends on whether we're in steady-state(dti=0) or not}
		THEN MaxIt:=par.MaxItSS 
		ELSE MaxIt:=par.MaxItTrans;

	coupleIonsPoisson:=TRUE; {signfies whether ion density can be modified by Poisson solver}
	
	REPEAT WITH new DO
	BEGIN
		INC(it);
		{perform 1 iteration, calc mobilities, densities, etc:}
		Calc_Elec_Mob(mun, V, n, stv, par); {calc. the new mobilities}
		Calc_Hole_Mob(mup, V, p, stv, par);
		Calc_Langevin_Factor(Lang, mun, mup, stv, par);{Calc. the Langevin recombination strength}
		Calc_Dissociation(diss_prob, gen, Gm, Lang, V, mun, mup, stv, par); {update generation rate}

		oldn:=n; {keep the old densities to monitor the progress}
		oldp:=p; {and so we can calculate their time derivatives. we need those to cacl. the current}
		oldnion:=nion;
		oldpion:=pion;

		{note: Solve_Poisson also modifies the charges (n,p,ions,traps) by estimating the effects of the newly calc'd potential}
		Solve_Poisson(V, n, p, nion, pion, f_tb, f_ti, Ntb_charge, Nti_charge, curr.f_tb, curr.f_ti, check_Poisson, coupleIonsPoisson, PoissMsg, dti, stv, par); 
		Update_Gen_Pot(V, Vgn, Vgp, stv, par); {update generalised potentials}
		{note: pass (new) p to Cont_Eq_Elec nor (new) n to Cont_Eq_Holes. This is needed so we can also do t=0!}
		Cont_Eq_Elec(n, curr.n, Vgn, Jn, p, mun, gen, Lang, diss_prob, curr.f_tb, curr.f_ti, Rn, stv, par, dti); {calc. new elec. density}
		Cont_Eq_Holes(p, curr.p, Vgp, Jp, n, mup, gen, Lang, diss_prob, curr.f_tb, curr.f_ti, Rp, stv, par, dti); {calc. new hole density}
		{note: transient ion solvers cannot (as yet) do steady-state, as that yields all densities=0!}
		IF UpdateIons THEN 
		BEGIN
			IF dti=0 {dti=0 => steady-state}
				THEN Calc_Ion_Distribution_Steady_State(nion, pion, V, stv, par) {use steady-state proc for ions}
			ELSE BEGIN {use transient versions:}
				IF par.negIonsMove THEN Solve_Neg_Ions(nion, curr.nion, V, dti, stv, par); {update neg ions}
				IF par.posIonsMove THEN Solve_Pos_Ions(pion, curr.pion, V, dti, stv, par) {update neg ions}
			END
		END;

		{now use SUR to get updated densities:}
		accDens:=par.maxAcc - it/MaxIt * (par.maxAcc - par.minAcc);
		FOR i:=0 TO par.NP+1 DO
		BEGIN
			deln[i]:=accDens*(n[i]-oldn[i]); {use SUR to get change in n}
			delp[i]:=accDens*(p[i]-oldp[i]); {use SUR to get change in p}
			n[i]:=oldn[i] + deln[i]; {now compute new densities, with damping (SUR)}
			p[i]:=oldp[i] + delp[i];

			{now do the same for the ions, but only if they need updating:}
			IF UpdateIons THEN BEGIN
				delnion[i]:=accDens*(nion[i]-oldnion[i]);
				delpion[i]:=accDens*(pion[i]-oldpion[i]);
				nion[i]:=oldnion[i] + delnion[i];
				pion[i]:=oldpion[i] + delpion[i]
			END
		END;

		{now check if the charge densities (n,p,nion,pion) have converged:}
		convDensities:=Deterimine_Convergence_Densities(deln, delp, delnion, delpion, n, p, nion, pion, UpdateIons, accDens, ConvMsg, stv, par);
			
		{if there are ions: Until the first time that convDensities, we have coupled the Poisson solver
		 and the ion solver by allowing the Poisson solver to modify the ions densities on the fly.
		 Now we see if the convergence is also OK if we forbid the Poisson solver to change the ion densities.}
		IF convDensities AND (par.negIonsMove OR par.posIonsMove) AND coupleIonsPoisson THEN
		BEGIN {force main solver to keep iterating, but now don't touch ions in Poisson solver}
			coupleIonsPoisson:=FALSE; {no longer update ions inside Poisson solver}
			convDensities:=FALSE {reset convDensities to make sure we'll keep iterating}
		END;	
		
		conv:=convDensities AND check_Poisson; {so convergence=true if index is positive}
		
		{now check if time is up!}
		IF (par.timeOut > 0) AND (SecondSpan(NOW, TimeStart) > par.timeOut) THEN 
			Stop_Prog('SIMsalabim terminates due to time-out.', EC_TimeOut)
		
	END; {WITH new statement}

	UNTIL conv OR (it = MaxIt); 
	
	{now compute the currents:}
	Calc_All_Currents(new, curr, stv, par); {calcs vectors Jn, Jp, Jnion, Jpion, JD and overall current Jint plus its rms error} 

	{finally, compute the effects of series and shunt resistance:}
	WITH new DO BEGIN
		IF par.Rshunt>0 THEN Jext:=Jint + Vint/par.Rshunt ELSE Jext:=Jint; {note: infinite Rshunt (no shunt) means Rshunt<0}
		Vext:=Vint + Jext*par.Rseries
	END;
	
	{now construct string to report on our progress:}
	StatusStr:='Overall convergence: '+myBoolStr(conv) + LineEnding;
	StatusStr:=StatusStr + 'Iterations perfomed: '+IntToStr(it) + LineEnding;
	StatusStr:=StatusStr + PoissMsg + LineEnding;
	StatusStr:=StatusStr + ConvMsg;
	
END;

PROCEDURE Prepare_tJV_File(VAR uitv : TEXT; filename : STRING; transient : BOOLEAN); 
{create a new tJV_file with appropriate heading
after running this, the TEXT file 'uitv' is still open and ready for writing}
BEGIN
	ASSIGN(uitv, filename);
	REWRITE(uitv); {rewrite old file (if any) or create new one}
    {write header, in the simulation we'll simply output the variables, but not change this header:}
	IF transient THEN WRITE(uitv,' t');
	WRITE(uitv,' Vext Jext errJ P Jphoto Jdir ');
	IF transient THEN WRITE(uitv, 'JBulkSRHn JBulkSHRp JIntLeftn JIntLeftp JIntRightn JintRightp ') 
		ELSE WRITE(uitv, 'JBulkSRH JIntLeft JIntRight ');
	WRITE(uitv,'JminLeft JminRight JShunt'); 
	IF transient THEN WRITELN(uitv,' Jdndt Jdpdt Jnion Jpion JD') ELSE WRITELN(uitv);
END;

PROCEDURE Write_To_tJV_File(VAR uitv : TEXT; CONSTREF CurrState, PrevState : Tstate; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{before running this proc, uitv must be open (by running Prepare_tJV_File). It must be closed in the main code.
This proc writes the (time), voltage, currents, recombination currents to a file that contains the JV-curve}
VAR JminLeft, JminRight, J1, J2 : myReal;
	NP : INTEGER;

	{first: 2 short functions that will aid in compact notation:}
	FUNCTION Ave(Vec : vector; onGrid : BOOLEAN = FALSE) : myReal;
	{local version of Average to calc the average over the full interval 0..NP+1 without having to specify the bounds}
	BEGIN
		Ave:=Average(Vec, stv.h, 0, par.NP+1, onGrid);
	END;
	
	FUNCTION EquiCurr(Vec : vector) : myReal;
	{local function to calc the Equivalent Current (A/m2) of a volume generation or recombination process}
	BEGIN
		EquiCurr:=q*Ave(Vec, FALSE)*par.L
	END;
	
BEGIN
    WITH CurrState DO 
    BEGIN
        IF transient THEN WRITE(uitv,tijd:nd,' ');  

		{now determine the minority currents left and right}
		NP:=par.NP; {local copy of number of grid points}
		IF n[0]<p[0] THEN JminLeft:=Jn[0] ELSE JminLeft:=Jp[0];
		IF n[NP+1]<p[NP+1] THEN JminRight:=Jn[NP+1] ELSE JminRight:=Jp[NP+1];   
        
        WRITE(uitv,Vext:nd,' ',Jext:nd,' ',errJ:nd,' ',Ave(diss_prob):nd,' ',EquiCurr(gen):nd,' ',
			EquiCurr(Rn.direct):nd,' ',EquiCurr(Rn.bulk):nd);
		IF transient THEN WRITE(uitv,' ',EquiCurr(Rp.bulk):nd); {transient => SRH Bulk n and p might be different!}

		{LEFT interface recombination currents:}
		IF stv.i1>0 THEN {interface rec at left interface for electrons}
			J1:=0.5*q*par.L*stv.h[stv.i1]*(Rn.int[stv.i1] + Rn.int[stv.i1+1])
		ELSE J1:=0; 
		WRITE(uitv,' ',J1:nd); {in steady-state J1 is equal to J2}
		IF transient THEN BEGIN
			IF stv.i1>0 THEN 
				J2:=0.5*q*par.L*stv.h[stv.i1]*(Rp.int[stv.i1] + Rp.int[stv.i1+1])
			ELSE J2:=0; {for transient, we also need the interface rec for holes}
			WRITE(uitv,' ',J2:nd);
		END;
		
		{RIGHT interface recombination currents:}
		IF stv.i2<par.NP+1 THEN {interface rec at right interface for electrons}
			J1:=0.5*q*par.L*stv.h[stv.i2]*(Rn.int[stv.i2-1] + Rn.int[stv.i2])
		ELSE J1:=0; 
		WRITE(uitv,' ',J1:nd); {in steady-state J1 is equal to J2}
		IF transient THEN BEGIN
			IF stv.i2<par.NP+1 THEN 
				J2:=0.5*q*par.L*stv.h[stv.i2]*(Rp.int[stv.i2-1] + Rp.int[stv.i2]) 
			ELSE J2:=0; {for transient, we also need the interface rec for holes}
			WRITE(uitv,' ',J2:nd);
		END;
	
		WRITE(uitv,' ',JminLeft:nd,' ',JminRight:nd,' ',Jext-Jint:nd);

        IF transient 
			THEN WRITELN(uitv,' ',q*par.L*dti*(Ave(n, TRUE)-Ave(PrevState.n, TRUE)):nd,' ',q*par.L*dti*(Ave(p, TRUE)-Ave(PrevState.p, TRUE)):nd,' ',Ave(Jnion):nd,' ',Ave(Jpion):nd,' ',Ave(JD):nd) 
			ELSE WRITELN(uitv);
    END; {with astate}
    FLUSH(uitv);
END;

PROCEDURE Prepare_Var_File(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN); {create a new var_file with appropriate heading}
VAR uitv : TEXT;
	e : INTEGER;
BEGIN
	ASSIGN(uitv, par.Var_file);
	REWRITE(uitv); {rewrite old file (if any) or create new one}
    {write header, in the simulation we'll simply output the variables, but not change this header:}
    WRITE(uitv, ' x V Evac Ec Ev phin phip n p nion pion ');
    
    {header depends on how many trap levels we have:}
    FOR e:=1 TO stv.N_Etr DO WRITE(uitv, 'ftb',e,' '); {filling of bulk trap level}
    FOR e:=1 TO stv.N_Etr DO WRITE(uitv, 'fti',e,' '); {filling interface trap level}

    {now write the rest:}
    WRITE(uitv, 'mun mup Gehp Gfree Rdir BulkSRHn BulkSRHp IntSRHn IntSRHp Jn Jp Jtot');
    IF transient 
		THEN WRITELN(uitv,' Jnion Jpion JD time') {add time! the ion & displacement currents are zero if not transient!}
		ELSE WRITELN(uitv,' Vext'); {add Vext so we can identify the different voltages}
    CLOSE(uitv);
END;

PROCEDURE Write_Variables_To_File(VAR CurrState : TState; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{writes the internal variables (astate) to file par.Var_file. It assumes the file has a header produced by Prepare_Var_File}
VAR i, e : INTEGER;
	uitv : TEXT;
	Evac, Ec, Ev, Jtot, phin, phip : myReal;
BEGIN
    ASSIGN(uitv, par.Var_file);
    APPEND(uitv);
	
    FOR i:=0 TO par.NP+1 DO
    WITH CurrState DO BEGIN
		Evac:=V[0] - V[i]; {the vacuum level. We take it zero at x=0}
		{use the generalised potentials to take care of the band-offsets}
		{but we need to correct for the effect of the DOS on Vgn/p:}
		Ec:=V[0] - Vgn[i] - par.CB + stv.Vt*LN(stv.NcLoc[i]/par.Nc);
		Ev:=V[0] - Vgp[i] - par.VB - stv.Vt*LN(stv.NcLoc[i]/par.Nc);
		Jtot:=Jn[i]+Jp[i]+JD[i]+Jnion[i]+Jpion[i]; {total current on grid point}
		{electron and hole quasi-Fermi levels:}
		phin:=Ec + stv.Vt*LN(n[i]/stv.NcLoc[i]);
		phip:=Ev - stv.Vt*LN(p[i]/stv.NcLoc[i]); 

        WRITE(uitv, stv.x[i]:nd,' ',V[i]:nd,' ',
				Evac:nd,' ',Ec:nd,' ',Ev:nd,' ',phin:nd,' ',phip:nd,' ', {band diagram}
				{all charged species:}
				n[i]:nd,' ',p[i]:nd,' ',nion[i]:nd,' ', pion[i]:nd,' ');
		{traps are seperate, per level:}
		FOR e:=1 TO stv.N_Etr DO {first bulk}
			WRITE(uitv, f_tb[i,e]:nd,' ');
		
		FOR e:=1 TO stv.N_Etr DO {then interface traps}
			WRITE(uitv, f_ti[i,e]:nd,' ');
		
		{continue with the rest:}
		WRITE(uitv,
				{transport:}
				mun[i]:nd,' ',mup[i]:nd,' ',
				{generation:}
				Gm[i]:nd,' ',gen[i]:nd,' ',
				{recombination:}
				Rn.direct[i]:nd,' ',Rn.bulk[i]:nd,' ',Rp.bulk[i]:nd,' ',Rn.int[i]:nd,' ',Rp.int[i]:nd,' ',
				{current densities:}
				Jn[i]:nd,' ',Jp[i]:nd,' ',Jtot:nd);	        
        IF transient 
			THEN WRITELN(uitv,' ',Jnion[i]:nd,' ',Jpion[i]:nd,' ',JD[i]:nd,' ',tijd:nd)
			ELSE WRITELN(uitv,' ',Vext:nd)
    END;
    CLOSE(uitv);
END;

PROCEDURE Tidy_Up_Parameter_File(parameterFile : STRING; QuitWhenDone : BOOLEAN);
{This procedure cleans up the parameter file: it makes sure that all the * are aligned, that every parameter has a 
unit or other description/comment and left-aligns any line that starts with **. It does this by reading the original device
parameter file line by line and writing corrected lines to a temp file. Once the correct position of the descriptions 
(starting with *) has been found, it uses the temp file to create the new, tidy parameter file. Lastly, the temp file is
removed and the program exits.}
VAR inp, outp : TEXT;
    line, dumstr, outline : STRING;
    max_pos, linecount, pos_asterix, i : INTEGER;
BEGIN
    {open the original parameterFile}
    ASSIGN(inp, parameterFile);
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
				Stop_Prog('See line number '+IntToStr(linecount)+'.', EC_DevParCorrupt);
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
    {and start writing into the real parameterFile:}
    ASSIGN(outp, parameterFile);
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
    
    {now delete old parameterFile and replace with temp_file}
    ERASE(inp); {the file should be assigned with Assign, but not opened with Reset or Rewrite}
    RENAME(outp, parameterFile); {file must be assigned, but not opened}
    IF QuitWhenDone 
		THEN Stop_Prog('Replaced old '+parameterFile+' with a tidy version.', EC_Warning)
		ELSE WRITELN('Tidied-up ',parameterFile,'.');
END;


BEGIN 

END.
