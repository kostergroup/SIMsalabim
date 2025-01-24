unit DDRoutines;
{provides drift-diffusion procedures and functions}

{
SIMsalabim:a 1D drift-diffusion simulator 
Copyright (c) 2021, 2022, 2023, 2024 S. Heester, Dr T.S. Sherkar, Dr V.M. Le Corre, 
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

CONST DDRoutinesVersion = '5.14'; {version of this unit}

{now check to see if the versions of the units match that of this code:}
{$IF (TransferMatrixVersion <> DDRoutinesVersion) OR (DDTypesAndConstantsVersion <> DDRoutinesVersion)} 
	{$STOP Wrong version of one or more units!}
{$ENDIF}


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

FUNCTION Correct_Version_Parameter_File(parameterFile, version : STRING; CheckProgName : BOOLEAN = FALSE; ProgName : TProgram = ZimT) : BOOLEAN;
{checks if the version and the parameter file match. Stops program if not! If CheckProgName, then we also check this.}

PROCEDURE Prepare_Log_File(VAR log : TEXT; MsgStr : ANSISTRING; CONSTREF par : TInputParameters; version : STRING);
{opens a logFile for later use and writes date/time and MsgStr to the log file.}

PROCEDURE Finalize_Log_File(VAR log : TEXT; MsgStr : ANSISTRING);
{writes final comments, date, time, run time and closes log file.}

PROCEDURE Read_Parameters(parameterFile : STRING; VAR msg : ANSISTRING; VAR par : TInputParameters; VAR stv : TStaticVars; ProgName : TProgram);
{Reads-in all the parameters. Some bits are specific to either ZimT or SimSS}

PROCEDURE Check_Parameters(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; ProgName : TProgram);
{performs a number of checks on the parameters. Just to ensure that they are valid, consistent, make sense}
{Some bits are specific to either ZimT or SimSS}

PROCEDURE Make_Grid(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Makes an exponential symmetric grid, for every layer}
{h[i] = (x[i+1] - x[i])/Ltot and initialises the array with x-positions}
{we also define an array stv.lid: the Layer ID. This allows us to determine in which layer a grid point falls.}
{and 2 arrays (i0 and i1) that contain the first and last grid point of each layer}

PROCEDURE Define_Layers(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Note, stv are not CONSTREF as we need to change them}
{Sets layer dependent properties}

PROCEDURE Init_Generation_Profile(VAR stv : TStaticVars; VAR log : TEXT; CONSTREF par : TInputParameters);
{Inits the generation profile, either constant, calculated by the transfer matrix unit, or from a file. 
This is the SHAPE of the profile. Also inits stv.Lgen: the sum of the lengths/thicknesses of all layers that can generate elec-hole pairs.}
{When using a profile from file, a message is written on screen and in the log file.}

PROCEDURE Update_Generation_Profile(org: vector; VAR new : vector; G_frac : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Rescales generation profile (org) by factor G_frac to obtain new profile}

PROCEDURE Init_Pot_Dens_Ions_Traps(VAR V, Vgn, Vgp, n, p, nion, pion : vector; VAR f_tb, f_ti, f_ti_numer, f_ti_inv_denom : TTrapArray; Va : myReal; VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{init. for V, Vgn,p, n, p, ion densities at bias voltage Va. Also sets lengths of f_tb/i arrays}

PROCEDURE Init_Trapping(VAR log : TEXT; VAR stv : TStaticVars; CONSTREF par : TInputParameters); 
{Inits all variables needed for trapping and SRH recombination}

PROCEDURE Extrapolate_Solution(CONSTREF prev, curr : TState; VAR new : TState; AccSols : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Obtains a guess for state new based on an extrapolation of states prev and curr}

PROCEDURE Main_Solver(VAR curr, new : TState; VAR it : INTEGER; VAR conv : BOOLEAN; VAR StatusStr : ANSISTRING; 
					  CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Iteratively solves the Poisson and continuity equations, including traps and ions}
{can be used in steady-state and transient cases}

PROCEDURE Prepare_tJV_File(VAR uitv : TEXT; filename : STRING; transient : BOOLEAN; CONSTREF stv : TStaticVars); 
{create a new tJV_file with appropriate heading
after running this, the TEXT file 'uitv' is still open and ready for writing}

PROCEDURE Write_To_tJV_File(VAR uitv : TEXT; CONSTREF CurrState, PrevState : Tstate; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{before running this proc, uitv must be open (by running Prepare_tJV_File). It must be closed in the main code.
This proc writes the (time), voltage, currents, recombination currents to a file that contains the JV-curve}

PROCEDURE Prepare_Var_File(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN); 
{create a new var_file with appropriate heading}

PROCEDURE Write_Variables_To_File(VAR CurrState : TState; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{writes the internal variables (astate) to file par.varFile. It assumes the file has a header produced by Prepare_Var_File}

PROCEDURE Tidy_Up_Parameter_Files(parameterFile : STRING; QuitWhenDone : BOOLEAN; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure cleans up the parameter files: the main file and the files with the parameters of each layer}

FUNCTION Copy_State(CONSTREF a : TState; CONSTREF par : TInputParameters) : TState;
{Returns a copy of state a. We need this because of the dynamic arrays in TState.}

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
    WRITELN('Copyright (C) 2020, 2021, 2022, 2023, 2024, S. Heester, Dr T.S. Sherkar,'); 
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
    WRITELN('Example: ./',strprogname,' -T 400 -varFile Var.dat');
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

FUNCTION Correct_Version_Parameter_File(parameterFile, version : STRING; CheckProgName : BOOLEAN = FALSE; ProgName : TProgram = ZimT) : BOOLEAN;
{checks if the version and the parameter file match. Stops program if not! If CheckProgName, then we also check this.}
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
    Correct_Version_Parameter_File:=found_version AND (NOT CheckProgName OR found_progname);
END;

PROCEDURE Prepare_Log_File(VAR log : TEXT; MsgStr : ANSISTRING; CONSTREF par : TInputParameters; version : STRING);
{opens a logFile for later use and writes date/time and MsgStr to the log file.}
BEGIN
    ASSIGN(log, par.logFile); 
    REWRITE(log);
    WRITELN(log,'Version ', version);
    WRITELN(log,'Size of reals used in simulation: ',SizeOf(myReal),' bytes');
    TimeStart:=NOW; {this is a global variable that will also be used in Finalize_Log_File}
    WRITELN(log,'Timestamp at start: ',DateTimeToStr(TimeStart));
    
    {Read_Parameters may have added something to MsgStr (values from the command line):}
    IF Length(MsgStr) > 0 THEN
    BEGIN
		WRITELN(log, 'Reading in of parameters:');
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


PROCEDURE Read_Layer_Parameters(VAR layerPar : TLayerParameters; layerNumber : INTEGER; VAR msg : ANSISTRING; VAR stv : TStaticVars);
{Reads-in parameters of layer layerNumber}

VAR dumint : INTEGER; {a dummy integer variable}
    dumstr, CLpre : STRING; 
    inv : TEXT;
BEGIN
    
    IF NOT FileExists(layerPar.layerFile) {the file with input par. is not found}
		THEN Stop_Prog('Could not find file '+layerPar.layerFile+'.', EC_FileNotFound);
    
    {now check if the version is correct. For the main simulation setup file, we do this in the main simss/zimt code.}
    IF NOT Correct_Version_Parameter_File(layerPar.layerFile, DDRoutinesVersion) THEN Stop_Prog('Version of SIMsalabim and '+layerPar.layerFile+' do not match.', EC_DevParCorrupt);
    
    ASSIGN(inv, layerPar.layerFile);
    RESET(inv);
	
	WRITELN('Reading parameters of layer ',layerNumber,' from file ',layerPar.layerFile);
	
	CLpre:='l' + IntToStr(layerNumber) + '.'; {this string (l1, l2, etc) is used as a prefix in the command line}
	
	WITH layerPar DO BEGIN
{**General**************************************************************************}
		Get_Float(inv, msg, 'L', L, CLpre); {device thickness, m}
		Get_Float(inv, msg, 'eps_r', eps_r, CLpre);  {relative dielectric constant}
		Get_Float(inv, msg, 'E_c', E_c, CLpre);  {eV, conduction band edge}
		Get_Float(inv, msg, 'E_v', E_v, CLpre);  {eV, valence band edge}
		Get_Float(inv, msg, 'N_c',N_c, CLpre);  {effective DOS, m^-3}
		Get_Float(inv, msg, 'N_D', N_D, CLpre);  {ionised n-doping density, m^-3}
		Get_Float(inv, msg, 'N_A', N_A, CLpre);  {ionised p-doping density, m^-3}

{**Mobilities************************************************************************}
		Get_Float(inv, msg, 'mu_n', mu_n, CLpre); {electron zero-field mobility, m^2/Vs}
		Get_Float(inv, msg, 'mu_p', mu_p, CLpre); {hole zere-field mobility, m^2/Vs}
		Get_Integer(inv, msg, 'mobnDep', mobnDep, CLpre);{dependence of elec mobility, 0 : const. mob, 1 : field-dep}
		Get_Integer(inv, msg, 'mobpDep', mobpDep, CLpre);  {dependence of hole mobility, 0 : const. mob, 1 : field-dep}
		Get_Float(inv, msg,'gamma_n', gamma_n, CLpre); {field depedence of mobility, eV/(V/m)^0.5}
		Get_Float(inv, msg, 'gamma_p', gamma_p, CLpre); {field depedence of mobility, eV/(V/m)^0.5}

{**Interface-layer-to-right***********************************************************}
		Get_Float(inv, msg, 'nu_int_n', nu_int_n, CLpre); {m/s, interface transfer velocity for electrons}
		Get_Float(inv, msg, 'nu_int_p', nu_int_p, CLpre); {m/s, interface transfer velocity for holes}
		Get_Float(inv, msg, 'N_t_int', N_t_int, CLpre); {m^-2, interface trap density}
		Get_Float(inv, msg, 'E_t_int', E_t_int, CLpre); {eV, energy level of traps at interface}
		Get_String(inv, msg, 'intTrapFile', intTrapFile, CLpre); {name of file with interface trap energy profile (or 'none'). If specified, overrides E_t_int}
		intTrapFromFile:= lowercase(Trim(intTrapFile)) <> 'none'; {use the profile if intTrapFile isn't 'none'}	
		Get_Integer(inv, msg, 'intTrapType', intTrapType, CLpre); {Trap type for the interface to the right: -1: acceptor, 0: neutral, 1: donor}	
		stv.Traps_int_poisson:=stv.Traps_int_poisson OR (N_t_int<>0) AND (intTrapType<>0) AND (layerNumber<stv.NLayers); {note: last layer does not specify interface traps}
		Get_Float(inv, msg, 'C_n_int', C_n_int, CLpre); {m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)}
		Get_Float(inv, msg, 'C_p_int', C_p_int, CLpre); {m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)}

{**Ions*******************************************************************}
		Get_Float(inv, msg, 'N_anion', N_anion, CLpre); {m^-3, concentration of negative ions}
		Get_Float(inv, msg, 'N_cation', N_cation, CLpre); {m^-3, concentration of positive ions}
		Get_Float(inv, msg, 'mu_anion', mu_anion, CLpre);{mobility of negative ions}
		Get_Float(inv, msg, 'mu_cation', mu_cation, CLpre);{mobility of negative ions}
		Get_Integer(inv, msg, 'ionsMayEnter', dumint, CLpre); {may ions enter from other layers? yes(1) or no(<>1)}
		ionsMayEnter:=dumint=1;

{**Generation and recombination******************************************************}
		Get_Float(inv, msg, 'G_ehp', G_ehp, CLpre); {m^-3 s^-1, generation rate of electron-hole pairs in this layer}
		Get_Integer(inv, msg, 'layerGen', dumint, CLpre); {does this layer generate electron/hole pairs? yes(1) or no (0)')}
		layerGen:=dumint=1;
		Get_String(inv, msg, 'nkLayer', nkLayer, CLpre); {name of file with n and k of this layer}
		Get_Integer(inv, msg, 'fieldDepG', dumint, CLpre);  {field-dependent G, true or false}
		fieldDepG:=(ROUND(dumint) = 1);
		Get_Float(inv, msg, 'P0', P0, CLpre); {0<=P0<1, fraction of quenched excitons that directly yield free carriers}
		Get_Float(inv, msg, 'a', a, CLpre); {thermalization length, Braun model used, m}
		Get_Integer(inv, msg, 'thermLengDist', thermLengDist, CLpre);
		Get_Float(inv, msg, 'k_f', k_f, CLpre); {decay rate of CT state, 1/s}
		Get_Float(inv, msg, 'k_direct', k_direct, CLpre); {m3/s, direct (band-to-band, bimolecular) recombination rate}
		Get_Float(inv, msg, 'preLangevin', preLangevin, CLpre); {Langevin prefactor}
		Get_Integer(inv, msg, 'useLangevin', dumint, CLpre);
		useLangevin:=(ROUND(dumint) = 1); {Calculate recombination using Langevin equation (1) or direct input (<>1, k_direct is used))}

{**Bulk trapping*********************************************************************}
		Get_Float(inv, msg,'N_t_bulk', N_t_bulk, CLpre); {m^-3, trap density (in bulk)}
		Get_Float(inv, msg, 'C_n_bulk', C_n_bulk, CLpre); {m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)}
		Get_Float(inv, msg, 'C_p_bulk', C_p_bulk, CLpre); {m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)}
		Get_Float(inv, msg, 'E_t_bulk', E_t_bulk, CLpre); {eV, energy level of all traps}
		Get_String(inv, msg, 'bulkTrapFile', bulkTrapFile, CLpre); {name of file with bulk trap energy profile (or 'none'). If specified, overrides E_t_bulk}
		bulkTrapFromFile:= lowercase(Trim(bulkTrapFile))<>'none'; {use the profile if bulkTrapFile isn't 'none'}	
		Get_Integer(inv, msg, 'bulkTrapType', bulkTrapType, CLpre); {Trap type of bulk and grain boundary traps: -1: acceptor, 0: neutral, 1: donor}	

    END; {WITH par statement}
    
    CLOSE(inv);

    {now write to screen and log (via msg) that we've read from the file:}
    dumstr:='Read parameters for layer '+IntToStr(layerNumber)+' from '+layerPar.layerFile+LineEnding;
	msg:=msg + dumstr

END;

FUNCTION ObtainWorkFunction(VAR WF : myReal; dumstr : STRING; CONSTREF lyr : TLayerParameters; CONSTREF stv : TStaticVars) : BOOLEAN;
{computes the work function of electrode from dumstr: 1) direct input (convert to numerical value), or 2) semi-flat band 
(dumstr = 'sfb'), then WF follows from net doping at electrode. Returns TRUE if successful.}
VAR netDoping : myReal;
	code : WORD;
BEGIN
	ObtainWorkfunction:=FALSE;
	
	IF LowerCase(dumstr) = 'sfb' THEN {we will have to compute WF from doping}
	BEGIN 
		netDoping:=lyr.N_A - lyr.N_D; {be careful, this can be zero!}
		IF netDoping > 0 THEN {efficively, p-doped}
			WF:=lyr.E_v - stv.Vt*LN(lyr.N_c/netDoping);
		IF netDoping < 0 THEN {so efficively n-doped:}
			WF:=lyr.E_c + stv.Vt*LN(-lyr.N_c/netDoping); 
		IF netDoping = 0 THEN {intrinsic, or fuly compensated}
			WF:=0.5 * (lyr.E_c + lyr.E_v);
		ObtainWorkfunction:=TRUE
	END
	ELSE BEGIN {WF's value is simply given}
		VAL(dumstr, WF, code); {if conversion successful, then code=0}
		ObtainWorkfunction:=code=0
	END

END;

PROCEDURE Read_Parameters(parameterFile : STRING; VAR msg : ANSISTRING; VAR par : TInputParameters; VAR stv : TStaticVars; ProgName : TProgram);
{Reads-in all the parameters. Some bits are specific to either ZimT or SimSS}
VAR countStart, countFinish, i : INTEGER;
    dumint : INTEGER; {a dummy integer variable}
    dumstr, dumstr2 : ANSISTRING; {again, dummy variables}
    inv : TEXT;
    ZimT, SimSS, UsedSpecialParFile, FoundLayer	: BOOLEAN; 

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

{**Layers****************************************************************************}
		{now we try to read the parameter files for each layer.}
		stv.NLayers:=1;
		SetLength(lyr, stv.NLayers + 1); {l is the array that contains all layer parameters. First: length is 2 as Layer 0 is reserved for the TCO and layer 1 is the first real layer.}
		Get_String(inv, msg, 'l1', lyr[1].layerFile); {parameter file for layer 1, mandatory!}

		{this is the only part of the parameter files that is flexible: we'll try to see if there are more layers!}
		REPEAT
			readln(inv, dumstr);
			dumstr:=DelSpace(Trim(dumstr)); {remove all whitespace}
			{dumstr is either 1) a layer so starts with 'lXX', 2) blank/whitespace, or 3) starts with *}
			dumstr2:='l' + IntToStr(stv.NLayers+1) + '=';
			FoundLayer:=LeftStr(LowerCase(dumstr), LENGTH(dumstr2)) = dumstr2;
			IF FoundLayer THEN TRY
				dumstr:=ExtractWord(2, dumstr, ['=', '*']); {the file name sits between = and a *}
				INC(stv.NLayers);
				SetLength(lyr, stv.NLayers+1); {this adds the layer in the array lyr}
				lyr[stv.Nlayers].layerFile:=Trim(dumstr);
			EXCEPT
				Stop_Prog('Error reading name of layer '+IntToStr(stv.NLayers)+' in file '+parameterFile+'.', EC_DevParCorrupt);
			END
		UNTIL NOT FoundLayer;

		{now see if the command line specifies any other names of the layers:}
		FOR i:=2 TO stv.NLayers DO {exclude layer 1 as we treated that one as per usual}
		BEGIN
			dumstr:='-l' + IntToStr(i);
			getStringfromCL(dumstr, FoundLayer, dumstr2); {try to get it from command line}
			IF FoundLayer THEN 
			BEGIN
				lyr[i].layerFile:=Trim(dumstr2);
				msg:=msg + 'l' + IntToStr(i) + ' = '+ lyr[i].layerFile + LineEnding
			END
		END;

		{read the parameters of all layers:}
		stv.Traps_int_poisson:=FALSE; {init this field, this can become TRUE if a layer contains interface traps}
		FOR i:=1 to stv.NLayers DO
			Read_Layer_Parameters(lyr[i], i, msg, stv); 

{**Contacts**************************************************************************}
		Get_String(inv, msg, 'W_L', dumstr); {eV, work function left electrode (= cathode), or 'sfb'}
		IF NOT ObtainWorkFunction(W_L, dumstr, lyr[1], stv) THEN Stop_Prog('Could not convert value of W_L to a numerical value.', EC_InvalidInput);	
		Get_String(inv, msg, 'W_R', dumstr); {eV, work function right electrode (= anode), or 'sfb'}
		IF NOT ObtainWorkFunction(W_R, dumstr, lyr[stv.NLayers], stv) THEN Stop_Prog('Could not convert value of W_R to a numerical value.', EC_InvalidInput);

		stv.V0:=0.5*(-W_L + W_R);
		stv.VL:=0.5*(W_L - W_R);
		
		Get_Float(inv, msg, 'S_n_L', S_n_L); {m/s, surface recombination of electrons at the left electrode.}
		Get_Float(inv, msg, 'S_p_L', S_p_L); {m/s, surface recombination of holes at the left electrode.}
		Get_Float(inv, msg, 'S_n_R', S_n_R); {m/s, surface recombination of electrons at the right electrode.}
		Get_Float(inv, msg, 'S_p_R', S_p_R); {m/s, surface recombination of holes at the right electrode.}
		Get_Float(inv, msg, 'R_shunt', R_shunt); {Ohms m2, shunt resistance. Use negative value for infinite R_shunt}
		Get_Float(inv, msg, 'R_series', R_series); {Ohms m2, series resistance}

{**Optics****************************************************************************}
		IF SimSS THEN Get_Float(inv, msg, 'G_frac', G_frac);
		Get_String(inv, msg, 'genProfile', genProfile); {name of file generation profile (or lStr+'none')}
		CASE lowercase(Trim(genProfile)) OF
			'none' : Use_gen_profile := 0; {Uniform generation}
			'calc' : Use_gen_profile := 1; {Calculate generation profile using the Transfermatrix method}
		ELSE
			Use_gen_profile := 2; {Use an user-defined generation profile}
		END;
		Get_Float(inv, msg, 'L_TCO', L_TCO); {m, thickness of the TCO. Set to 0 if layer is not used}
		Get_Float(inv, msg, 'L_BE', L_BE); {m, thickness of back electrode, must be >0}
		Get_String(inv, msg, 'nkSubstrate', nkSubstrate); {name of file with n,k values of substrate}
		Get_String(inv, msg, 'nkTCO', nkTCO); {name of file with n,k values of TCO}
		Get_String(inv, msg, 'nkBE', nkBE); {name of file with n,k values of back electrode}
		Get_String(inv, msg, 'spectrum', spectrum); {name of file that contains the spectrum}
		Get_Float(inv, msg, 'lambda_min', lambda_min); {m, lower bound wavelength}
		Get_Float(inv, msg, 'lambda_max', lambda_max); {m, upper bound wavelength}

{**Numerical Parameters**************************************************************}
		Get_Integer(inv, msg, 'NP', NP); {number of grid points}
		Get_Float(inv, msg, 'tolPois', tolPois); {abs. tolerance of Poisson solver}
		Get_Float(inv, msg, 'maxDelV', maxDelV); {maximum change (in Vt) of the potential per loop}
		Get_Integer(inv, msg, 'maxItPois', maxItPois); {Max. number of loops Poisson solver}
		Get_Integer(inv, msg, 'maxItSS', maxItSS); {max. number it. steady-state loops}
		IF ZimT THEN Get_Integer(inv, msg, 'maxItTrans', maxItTrans); {max. number it. transient solver}
		Get_Integer(inv, msg, 'currDiffInt', currDiffInt); {Calc. current from differential (1) or integral (2) expression}
		Get_Float(inv, msg, 'tolDens', tolDens); {relative tolerance of density solver}
		Get_Float(inv, msg, 'couplePC', couplePC); {>= 0, coupling between Poisson equation and continuity equations}
		Get_Float(inv, msg, 'minAcc', minAcc); {>0, min. acceleration parameter}
		Get_Float(inv, msg, 'maxAcc', maxAcc); {<2, max. acceleration parameter}
		Get_Integer(inv, msg, 'ignoreNegDens', dumint);
		ignoreNegDens:= dumint=1; {whether(1) or not(<>1) to ignore negative densities}
		Get_Integer(inv, msg, 'failureMode', failureMode); {how treat failed (t,V,G) points: 0: stop, 1: ignore, 2: skip}
		Get_Float(inv, msg, 'grad', grad); {gradient of grid, increase grad for smaller h[1]}
		IF ZimT THEN Get_Float(inv, msg, 'tolVint', tolVint); {V, tolerance in internal voltage (Vint)}

{**Voltage range of simulation*******************************************************}
		IF SimSS {this entire block is only relevant to SimSS}
		THEN BEGIN
			Get_Integer(inv, msg, 'Vdist', Vdist); {type of V distribution, 1=linear, 2=logarithmic}
			Get_Integer(inv, msg, 'preCond', dumint); {Pre-condition in light(1)/dark(0)}
			preCond:=ROUND(dumint)=1;
			Get_Float(inv, msg, 'Vpre', Vpre); {V, pre-conditioned voltage}
			Get_Integer(inv, msg, 'fixIons', dumint); {fix ions at first applied voltage? yes(1) or no (0).}
			fixIons:= (dumint=1);			
			Get_Integer(inv, msg, 'Vscan', Vscan); {integer, direction of voltage scan: up = 1, down = -1}
			Get_Float(inv, msg, 'Vmin', Vmin); {V, minimum voltage in JV characteristic}
			Get_Float(inv, msg, 'Vmax', Vmax); {V, max. voltage in JV}
			Get_Float(inv, msg, 'Vstep', Vstep); {V, voltage step}
			Get_Float(inv, msg, 'Vacc', Vacc); {accumulation voltage for logarithmic JV, should be outside [Vmin, Vmax]}
			Get_Integer(inv, msg, 'NJV', stv.NJV); {Number of JV points, for logarithmic JV}
			{Note: NJV is TStaticVars as we sometimes need to calculate it, so it's not a direc input parameter}
			IF (Vstep<>0) AND (Vdist=1) THEN
				stv.NJV:=TRUNC((Vmax - Vmin)/Vstep + 1e-10) + 1; {needed for setting length of Jdat and Vdat}
			{1e-10 is needed to get right value}
			Get_Integer(inv, msg, 'untilVoc', dumint); {if 1 then SimSS stops at Voc}
			untilVoc:= dumint;
		END;

{**User interface********************************************************************}
		Get_Float(inv, msg, 'timeout', timeout); {s, max run time. use negative value for unlimited run time.}
		Get_Integer(inv, msg, 'pauseAtEnd', dumint);
		pauseAtEnd:=dumint = 1;  {pause at the end of the simulation yes(1) or no (0)}
		Get_Integer(inv, msg, 'autoTidy', dumint);
		autoTidy:=dumint = 1;	{if 1, then we will always tidy up the device_parameter file}
		IF SimSS 
		THEN BEGIN
			Get_Integer(inv, msg, 'useExpData', dumint);
			useExpData:=dumint = 1; {if 1 then  SimSS will try to read expJV and use it}
			Get_String(inv, msg, 'expJV', expJV); {name of file with experimental JV points}
			Get_String(inv, msg, 'fitMode', dumstr); {lin or log: use J or log(J) in calc. of fit error}
			dumstr:=lowercase(dumstr);
			IF NOT ((dumstr='lin') OR (dumstr='log')) THEN Stop_Prog('fitMode has to be either lin or log.', EC_InvalidInput);
			IF dumstr='lin' THEN fitMode:=linear ELSE fitMode:=logarithmic;
			Get_Float(inv, msg, 'fitThreshold', fitThreshold); {threshold of fraction converged points in calc. fit error}
		END;		
		IF ZimT THEN
		BEGIN
			Get_Integer(inv, msg, 'autoStop', dumint);
			autoStop:= dumint=1; {stop ZimT if change of system stops chaning, yes(1) or no (<>1).	}
			Get_String(inv, msg, 'tVGFile', tVGFile); {name of file that specifies time t, voltage V and gen. rate G}
			Get_String(inv, msg, 'tJFile', tJFile); {name of file with (t, J, V, G)}
		END;
		IF SimSS THEN Get_String(inv, msg, 'JVFile', JVFile); {name of file with simulated JV points}
		Get_String(inv, msg, 'varFile', varFile); {name of file with internal variables}
		Get_Integer(inv, msg, 'limitDigits', dumint); {if 1, then number of digits in output is limited}
		limitDigits:=dumint = 1;
		Get_Integer(inv, msg, 'outputRatio', outputRatio); {output (ZimT: J to screen and) variables to var_file every outputRatio timesteps/voltages}	
		IF SimSS THEN StoreVarFile:=outputRatio>0;
		IF ZimT THEN StoreVarFile:=lowercase(Trim(varFile))<>'none'; {only store var_file if varFile isn't 'none'}    
		IF SimSS THEN Get_String(inv, msg, 'scParsFile', scParsFile); {name of file with solar cell parameters}
		Get_String(inv, msg, 'logFile', logFile); { name of log file}
    END; {WITH par statement}

    CLOSE(inv);
    WRITELN('Read simulation setup from ',parameterFile);
    {now check if all parameters that were in the command line have been used:}
	{did we use a special parameter file or not?}
	UsedSpecialParFile:=NOT( (ParamCount=0) OR StartsStr('-', TRIM(ParamStr(1))) );

    {we use countStart/Finish to count the number of parameters obtained from the command line}
	{every time we find such a parameter, a LineEnding is added to msg}
    countFinish:=Count_Substring_In_String(LineEnding, msg);

	{now the number of arguments (bar a dev par file) needs to be 2 x the difference in the counter:}
    IF ParamCount - ORD(UsedSpecialParFile) <> 2*(countFinish-countStart-stv.NLayers) THEN
		Stop_Prog('The command line contains invalid arguments, see the manual.', EC_InvalidCLInput);
    msg:=msg + 'Read simulation setup from ' + parameterFile + LineEnding
END;

PROCEDURE Check_Parameters_Layer(CONSTREF layerPar : TLayerParameters; layerNumber : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{performs a number of checks on the parameters of an individual layer}
VAR dumStr : STRING;
BEGIN
	dumStr:='Error in layer ' + IntToStr(layerNumber) +': ';
	
	WITH layerPar DO BEGIN
{checks on general parameters:}
		IF L <=0 THEN Stop_Prog(dumStr+'Thickness L should be positive.', EC_InvalidInput);
		IF (N_D < 0) OR (N_A < 0) THEN Stop_Prog(dumStr+'Doping densities cannot be negative.', EC_InvalidInput);
		IF E_c >= E_v THEN Stop_Prog(dumStr+'E_c should be smaller than E_v.', EC_InvalidInput);
		IF (E_c<0) OR (E_v<0) THEN Stop_Prog(dumStr+'E_c and E_v should be positive.', EC_InvalidInput);
		
{checks on mobilities:}
		IF (mu_n<=0) OR (mu_p<=0) THEN Stop_Prog(dumStr+'Mobilities mu_n and mu_p must be positive.', EC_InvalidInput); 
		IF NOT (mobnDep IN [0, 1]) THEN Stop_Prog(dumStr+'Invalid mob_dep_n selected.', EC_InvalidInput);
		IF NOT (mobpDep IN [0, 1]) THEN Stop_Prog(dumStr+'Invalid mob_dep_p selected.', EC_InvalidInput);

{checks on interface parameters:}
		IF (nu_int_n <= 0) OR (nu_int_p <= 0) THEN Stop_Prog(dumStr+'Interface transfer velocity nu_int_n/p must be positive.', EC_InvalidInput);
		
{checks on ions:}
		IF (N_anion <0) OR (N_cation<0) THEN Stop_Prog(dumStr+'Ionic concentrations cannot be negative.', EC_InvalidInput);
		IF (mu_anion<0) OR (mu_cation<0) THEN Stop_Prog(dumStr+'Ion mobilities cannot be negative.', EC_InvalidInput);
		IF ((N_anion>0) OR (N_cation>0)) AND NOT ionsMayEnter THEN Stop_Prog('If a layer contains ions, then ionsMayEnter (of that layer) must be 1.', EC_InvalidInput);
		IF ionsMayEnter AND (mu_anion*mu_cation = 0) THEN Stop_Prog('If ionsMayEnter = 1, then the ion mobilities cannot be 0.', EC_InvalidInput);
		
{checks on generation and recombination parameters}
		IF G_ehp < 0 THEN Stop_Prog('G_ehp cannot be negative.', EC_InvalidInput);
		IF (P0>=1) OR (P0<0) THEN Stop_Prog(dumStr+'Invalid value of P0, should be: 0<=P0<1', EC_InvalidInput);
		IF (P0<>0) AND (fieldDepG = FALSE) THEN Stop_Prog(dumStr+'P0 should be zero if not using field dependent generation', EC_InvalidInput);
		IF NOT (thermLengDist IN [1,2,3,4,5]) THEN Stop_Prog(dumStr+'Invalid thermLengDist selected.', EC_InvalidInput);
		IF (k_f<=0) AND fieldDepG THEN Stop_Prog(dumStr+'k_f must be positive.', EC_InvalidInput);
		
{checks on bulk traps}
		{check whether there are a possible number of traps (negative not allowed)}
		IF N_t_bulk < 0 THEN Stop_Prog(dumStr+'Negative bulk trap density not allowed.', EC_InvalidInput);
		IF (N_t_bulk > 0) AND (bulkTrapFile = 'none') AND ((E_t_bulk > E_v) OR (E_t_bulk < E_c)) THEN Stop_Prog(dumStr+'E_t_bulk must fall within E_c and E_v.', EC_InvalidInput);
		{Only Cn OR Cp = 0 are allowed, if both are zero, no charge can reach the traps, this makes no sense.}
		IF (N_t_bulk>0) AND (C_n_bulk = 0) AND (C_p_bulk = 0) THEN Stop_Prog(dumStr+'C_n_bulk and C_p_bulk cannot BOTH be zero, change parameters please.', EC_InvalidInput);	
		{trap types: -1, 0 or 1. Sets in pascal are an ordinal type with a range between 0 and 255, hence the ABS:}
		IF NOT (ABS(bulkTrapType) IN [0, 1]) THEN Stop_Prog(dumStr+'Invalid bulk trap type.', EC_InvalidInput);        
		{more checks on the energies of the traps are performed in proc Init_Trap_Distribution}

	END {with statement}
END;

PROCEDURE Check_Parameters(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; ProgName : TProgram);
{performs a number of checks on the parameters. Just to ensure that they are valid, consistent, make sense}
{Some bits are specific to either ZimT or SimSS}
VAR ZimT, SimSS	: BOOLEAN;
	i : INTEGER;
BEGIN
  
    {use 2 booleans to check if we're using ZimT or SimSS}
    ZimT:= (ProgName=TProgram.ZimT);
    SimSS:= (ProgName=TProgram.SimSS);
  
    {when adding new check, please keep the order in line with the device_parameter file}
    {Check first if stv.Vt has been initialised. This should have happened in Read_Parameters}
    IF (stv.Vt=0) THEN Stop_Prog('stv.Vt needs to be initialised before calling Check_Parameters.', EC_ProgrammingError);

{now we check all parameters in 3 blocks: 1) general, 2) layer specific, 3) at interfaces}

	WITH par DO BEGIN
{checks on contacts:}
		{part of this can only be done after checking the TLs!}
		IF R_series<0 THEN Stop_Prog('R_series cannot be negative.', EC_InvalidInput);
		IF R_shunt=0 THEN Stop_Prog('R_shunt cannot be zero, use positive (negative) value for finite (infinite) shunt resistance.', EC_InvalidInput);

{checks on energy levels:}
		IF W_L < lyr[1].E_c THEN Stop_Prog('W_L cannot be smaller than E_c of leftmost layer.',EC_InvalidInput);
		IF W_L > lyr[1].E_v THEN Stop_Prog('W_L cannot be larger than E_v of leftmost layer.', EC_InvalidInput);
		IF W_R < lyr[stv.NLayers].E_c THEN Stop_Prog('W_R cannot be smaller than E_c of rightmost layer.', EC_InvalidInput);
		IF W_R > lyr[stv.NLayers].E_v THEN Stop_Prog('W_R cannot be larger than E_v of rightmost layer.', EC_InvalidInput);
		
{checks on optics, generation and recombination parameters}
		IF Use_gen_profile = 1 THEN BEGIN
			IF L_TCO < 0 THEN Stop_Prog('L_TCO cannot be negative.', EC_InvalidInput);
			IF L_BE <= 0 THEN Stop_Prog('L_BE must be positive.', EC_InvalidInput);
			IF lambda_min > lambda_max THEN Stop_Prog('lambda_min cannot be larger than lambda_max.', EC_InvalidInput);
			IF lambda_min <= 0 THEN Stop_Prog('lambda_min must be positive.', EC_InvalidInput); 
			{note, now we are sure that lambda_max is also positive!}		
		END;

{checks on numerical parameters:}
		IF (NP<=15) OR (NP>Max_NP) THEN Stop_Prog('Invalid number of grid points (NP) selected, must be >=15 and <='+IntToStr(Max_NP)+'.', EC_InvalidInput);
		IF MinGridPointsPerLayer * stv.NLayers > NP THEN Stop_Prog('Not enough grid points (NP) to put the minimum number of grid points per layer.', EC_InvalidInput);
		IF (currDiffInt <> 1) AND (currDiffInt <> 2) THEN Stop_Prog('currDiffInt can only be 1 or 2.', EC_InvalidInput);
		IF maxDelV<=0 THEN Stop_Prog('maxDelV should be positive.', EC_InvalidInput);
		IF maxDelV*stv.Vt <= tolPois THEN Stop_Prog('maxDelV*Vt should be (much) larger than tolPois.', EC_InvalidInput);
		IF tolDens <= 0 THEN Stop_Prog('tolDens must be larger than zero.', EC_InvalidInput);
		IF couplePC < 0 THEN Stop_Prog('couplePC must be non-negative.', EC_InvalidInput);
		{check if values of minAcc and maxAcc makes any sense:}
		IF maxAcc >= 2 THEN Stop_Prog('maxAcc must be smaller than 2.', EC_InvalidInput);  
		IF minAcc <= 0 THEN Stop_Prog('minAcc must be positive.', EC_InvalidInput);  
		IF minAcc > maxAcc THEN Stop_Prog('minAcc cannot be larger than maxAcc.', EC_InvalidInput);  
		IF NOT (failureMode IN [0,1,2]) THEN Stop_Prog('Invalid failureMode selected.', EC_InvalidInput);

{checks on voltages, SimSS only:}
		IF SimSS THEN
		BEGIN
			IF NOT (Vdist IN [1,2]) THEN Stop_Prog('Invalid voltage distribution selected.', EC_InvalidInput);
			IF ABS(Vscan) <> 1 THEN Stop_Prog('Vscan must be either -1 or 1.', EC_InvalidInput);
			{check if Vmin and Vmax are not too small or large:}
			IF Vmin*stv.Vti < -1.95 * LN(Max_Value_myReal) THEN Stop_Prog('Vmin is too small.', EC_InvalidInput);
			IF Vmax*stv.Vti > 1.95 * LN(Max_Value_myReal) THEN Stop_Prog('Vmax is too large.', EC_InvalidInput);
			IF Vmin > Vmax THEN Stop_Prog('Vmin should not be greater than Vmax.', EC_InvalidInput);
			{now check for redundancy of pre-bias:}
			IF preCond THEN
			BEGIN
				IF R_series>0 THEN Warn_User('Pre-bias voltage does not take R_series into account, so Vpre=Vint.');
				IF ABS(Vpre)*stv.Vti > 1.95 * LN(Max_Value_myReal) THEN Stop_Prog('|Vpre| is too large.', EC_InvalidInput);
				IF (Vscan=1) AND (Vpre=Vmin) THEN Stop_Prog('Pre-bias voltage is equal to Vmin, makes no sense.', EC_InvalidInput);
				IF (Vscan=-1) AND (Vpre=Vmax) THEN Stop_Prog('Pre-bias voltage is equal to Vmax, makes no sense.', EC_InvalidInput);
			END;
			IF (Vacc >= Vmin) AND (Vacc <= Vmax) AND (Vdist = 2) THEN {Vacc is not valid} 
				Stop_Prog('Invalid Vacc selected, must be outside [Vmin, Vmax].', EC_InvalidInput);
			IF (Vdist=1) AND (Vstep <= 0) THEN Stop_Prog('Vstep should be positive.', EC_InvalidInput);	
			IF (ABS(Vmin-Vmax) < 1e-10) AND (Vdist=2) {to avoid infinite loop of Va}
			THEN Stop_Prog('Do not use Vdist=2 when Vmin = Vmax.', EC_InvalidInput);	
		END;

{checks on user-interface:}
		IF timeout = 0 THEN Stop_Prog('Invalid timeout: either positive number in seconds or negative for unlimited run time.', EC_InvalidInput);
		IF SimSS 
		THEN BEGIN
			IF (G_frac <> 0) AND (stv.V0 <> stv.VL) AND useExpData AND (fitMode=logarithmic) {this is a weird combination, warn user}
				THEN Warn_User('You are fitting a solar cell with fitMode=log.');
   			IF useExpData AND (untilVoc <> 0) THEN Stop_Prog('You cannot use untilVoc = 1 and useExpData = 1 at the same time.', EC_InvalidInput);
			IF useExpData AND preCond THEN Stop_Prog('You cannot use pre-conditioning (preCond) and useExpData = 1 at the same time.', EC_InvalidInput);
			IF ((fitThreshold<=0) OR (fitThreshold>1)) AND useExpData THEN
				Stop_Prog('fitThreshold has to be larger than 0 but not larger than 1.', EC_InvalidInput);
		END;		
		IF SimSS AND (outputRatio < 0) THEN Stop_Prog('outputRatio should be 0 (no output) or positive.', EC_InvalidInput); {if zero, then simply no var file output}
		IF ZimT AND (outputRatio <= 0) THEN Stop_Prog('outputRatio should be positive.', EC_InvalidInput); {In ZimT it cannot be zero as we NEED to write output as it also limits the output to screen.}
    END;
    
{now check the individual layers:}
	FOR i:=1 TO stv.NLayers DO
		Check_Parameters_Layer(par.lyr[i], i, stv, par);    
    
{now we check what happens at the interfaces:}    
    FOR i:=1 TO stv.NLayers-1 DO 
    WITH par DO BEGIN
		IF lyr[i].E_c > lyr[i+1].E_v THEN Stop_Prog('Invalid band alignment between layers '+IntToStr(i)+' and '+IntToStr(i+1)+'.', EC_InvalidInput);
		IF lyr[i].E_v < lyr[i+1].E_c THEN  Stop_Prog('Invalid band alignment between layers '+IntToStr(i)+' and '+IntToStr(i+1)+'.', EC_InvalidInput);
		IF (lyr[i].N_t_int > 0) AND NOT (ABS(lyr[i].intTrapType) IN [0, 1]) THEN Stop_Prog('Invalid right interface trap type.', EC_InvalidInput);
		{the trap energies (if any) are checked in proc Init_Trap_Distribution}
    END;
 
END;

PROCEDURE Make_Sub_Grid(VAR k : vector; istart, ifinish : INTEGER; grad, xstart, xfinish, L : myReal);
{this makes a grid on part of the volume}
VAR i, ending : INTEGER;
    norm : myReal;
BEGIN
    FOR i:=istart TO ifinish DO k[i]:=1;
    {note: we assign a value to k[np+1], but it doesn't have any meaning!}
    ending:=ROUND(0.5*(ifinish-istart)); {defines the exponential part of the grid}
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

PROCEDURE Make_Grid(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Makes an exponential symmetric grid, for every layer}
{h[i] = (x[i+1] - x[i])/Ltot and initialises the array with x-positions}
{we also define an array stv.lid: the Layer ID. This allows us to determine in which layer a grid point falls.}
{and 2 arrays (i0 and i1) that contain the first and last grid point of each layer}
VAR del, xstart, LThinLayers : myReal;
	i, j, NumPoints, NumThinLayers, NumGridPointsProportional, istart, ifinish : INTEGER;
BEGIN
    SETLENGTH(stv.i0, stv.NLayers+1); {i0: array with index (i) of start in each layer}
    SETLENGTH(stv.i1, stv.NLayers+1); {i1: array with index (i) of end in each layer}
	
	stv.Ltot:=0; {Ltot: total device length, so sum over layers.}	
	FOR i:=1 TO stv.NLayers DO
		stv.Ltot:=stv.Ltot + par.lyr[i].L;
	
	{first and last grid point in a layer:}
	istart:=0;
	ifinish:=0; 

	xstart:=0; {x-coordinate of layer}
	
	WITH stv DO BEGIN
		
		{First: count how many layers are 'thin'. 
		Thin = their thickness is such that they would get fewer than MinGridPointsPerLayer IF we were to assign the number of grid points
		purely based on the thickness of each layer.}
		NumThinLayers:=0;
		LThinLayers:=0; {the total thickness of all 'thin' layers}
		FOR i:=1 TO NLayers DO 
			IF par.NP*par.lyr[i].L/Ltot < MinGridPointsPerLayer THEN BEGIN
				LThinLayers:=LThinLayers + par.lyr[i].L;
				INC(NumThinLayers);
			END;
	
		{next: we will give each layer a number of grid points based on their thickness, but no fewer than MinGridPointsPerLayer}
		{and we will call Make_Sub_Grid to assign the grid spacing (h)}
		FOR i:=1 TO NLayers DO BEGIN
			NumGridPointsProportional:=FLOOR((par.NP-NumThinLayers*MinGridPointsPerLayer) * (par.lyr[i].L/(Ltot-LThinLayers)));
			NumPoints:=MAX(MinGridPointsPerLayer, NumGridPointsProportional);
			
			i0[i]:=istart;
			ifinish:=istart + NumPoints;
			i1[i]:=ifinish;
			
			{now make sure that last point is on right electrode:}
			IF i=NLayers THEN ifinish:=par.NP+1; {the last layer might thus get an extra point, fine!}
			
			{now make grid in this layer:}
			Make_Sub_Grid(stv.h, istart, ifinish, par.grad, xstart, xstart + par.lyr[i].L, Ltot);
	
			{lid: layer ID, so it tells you the number of the layer for a specific grid point:}
			FOR j:=istart TO ifinish DO
				lid[j]:=i;	
				
			{move to the next layer:}
			istart:=ifinish + 1; 
			xstart:=xstart + par.lyr[i].L
		END;
	
		{now the last layer will have 1 extra point at the end, so set i1 to par.NP+1:}
		i1[NLayers]:=par.NP+1; 
	
		{At the interfaces, it is beneficial (numerically!) to make sure the spacing is the same:}
		FOR i:=1 TO NLayers-1 DO BEGIN
			j:=i0[i+1]; {this is the first grid point in the layer to the right of the interface}
			del:=(h[j] + h[j-1] + h[j-2])/3; 
			h[j]:=del; h[j-1]:=del; h[j-2]:=del {we make the spacing astride the interface the same}
		END;
		
		{Now we know the grid spacing (h), we can compute the x-positions:}
		x[0]:=0;    
		FOR i:=1 TO par.NP+1 DO x[i]:=x[i-1] + Ltot*h[i-1];
	
	END {with stv statement}

END;


PROCEDURE Define_Layers(VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{Note, stv are not CONSTREF as we need to change them}
{Sets layer dependent properties}
VAR i, j : INTEGER;
BEGIN
	FOR j:=1 TO stv.NLayers DO {loop over all layers}
		WITH par.lyr[j] DO
			FOR i:=stv.i0[j] TO stv.i1[j] DO BEGIN
				stv.NcLoc[i]:=N_c;
				stv.E_CB[i]:=E_c;
				stv.E_VB[i]:=E_v;
				stv.ni[i]:=N_c*EXP(-0.5*stv.Vti*(E_v-E_c)); {equilibrium concentration}
				stv.nid[i]:=N_D;
				stv.pid[i]:=N_A;
				stv.eps[i]:=eps_r * eps_0			
			END			
END;

PROCEDURE Init_Ionic_Region(VAR IonRegion : TIonicRegions; sn : ShortInt; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Ions can, depending on parameters, move in/out multiple adjacent layers. 
This proc determines which layers from a region of adjacent layers where ions (pos or neg = sn) can move}
VAR start, finish, j, NumRegions : INTEGER;
	FoundStart, FoundEnd : BOOLEAN;
	TotIons : myReal;
BEGIN
	IF ABS(sn)<>1 THEN Stop_Prog('Init_Ionic_Region called with invalid sn (= sign of ions, +-1).',EC_ProgrammingError);
	start:=0;
	finish:=0;
	FoundStart:=FALSE; {true if we find the start of a region with ions}
	FoundEnd:=FALSE; {true if we find the end of a region with ions}
	TotIons:=0; {total density of ions in a region}
	NumRegions:=0; {the number of regions with ions}

	FOR j:=1 TO stv.NLayers DO {loop over the layers}
	BEGIN
		
		{first: can there be ions in this layer?}
		IF par.lyr[j].ionsMayEnter THEN 
		BEGIN {there may be ions in this layer}
			finish:=stv.i1[j];
			IF sn = 1 THEN
				TotIons:=TotIons + par.lyr[j].L * par.lyr[j].N_cation
			ELSE
				TotIons:=TotIons + par.lyr[j].L * par.lyr[j].N_anion;
	
			IF NOT FoundStart THEN		
			BEGIN {we have found a NEW region with ions!}
				FoundStart:=TRUE;
				start:=stv.i0[j]; 
			END
		END
		ELSE 
		BEGIN {so there can be no ions in this layer}
			FoundStart:=FALSE;
			FoundEnd:=FALSE;
			TotIons:=0
		END;
		
		{now check if we are at the end of a region:}
		IF j=stv.NLayers THEN
			FoundEnd:=TRUE
		ELSE
			FoundEnd:=NOT par.lyr[j+1].ionsMayEnter;
			
		{we have found a full new region if we have the start, end, and some ions:}
		IF FoundStart AND FoundEnd AND (TotIons<>0) THEN
		BEGIN
			INC(NumRegions);
			SetLength(IonRegion, NumRegions); {this adds the region to the array}
			WITH IonRegion[NumRegions-1] DO {now copy the values to this record}
			BEGIN
				istart:=start;
				ifinish:=finish;
				AvC:=TotIons/(stv.x[ifinish] - stv.x[istart])
			END
		END;

	END; {loop over layers}
	
END;

PROCEDURE Init_Generation_Profile(VAR stv : TStaticVars; VAR log : TEXT; CONSTREF par : TInputParameters);
{Inits the generation profile, either constant, calculated by the transfer matrix unit, or from a file. 
Also inits stv.Lgen: the sum of the lengths/thicknesses of all layers that can generate elec-hole pairs.}
{When using a profile from file, a message is written on screen and in the log file.}
VAR a, gr : Row; {a : x-coordinate, gr: generation rate}
    numLines, i, j : INTEGER;
BEGIN
	{first init stv.Lgen, by looping over all layers:}
	stv.Lgen:=0;
	FOR j:=1 TO stv.NLayers DO
		IF par.lyr[j].layerGen THEN {if a layer then generate elec-hole pairs from absorbed light, then add thickness to Lgen}
			stv.Lgen:=stv.Lgen + par.lyr[j].L;
	
	
	CASE par.Use_gen_profile OF
		0:	FOR i:=0 TO par.NP+1 DO {uniform generation, so orgGm simply follows from G_ehp of layer}
				stv.orgGm[i]:=par.lyr[stv.lid[i]].G_ehp;	
		1:	Calc_TransferMatrix(stv,log,par); {Calculate the generation profile using the TransferMatrix script.}
		2:	BEGIN {use profile supplied by user}
				WRITELN('Reading generation profile from ',par.genProfile);
				WRITELN(log, 'Reading generation profile from ',par.genProfile);
				Read_XY_Table(a, gr, par.genProfile, 'x G_ehp', numLines);
				{a: x-coordinate, gr: corresponding generation rate}

				IF numLines=0 THEN Stop_Prog('The file '+par.genProfile+' appears to be empty.', EC_InvalidInput);

				{ rescale a to ensure that a[numLines-1]=L: }
				FOR i:=0 TO numLines-1 DO a[i]:=a[i] * stv.Ltot/a[numLines-1];
				{the input file contains the generation profile, however, the grid points
				in this file may not correspond to our grid here}

				{Now interpolate gr to get genrate on x[i] grid:}
				i:=0; {counter for x[i] grid}
				FOR j:=0 TO numLines-2 DO
					WHILE (stv.x[i] < a[j+1]) AND (i<par.NP+1) DO
						BEGIN
							IF par.lyr[stv.lid[i]].layerGen THEN
								stv.orgGm[i]:=gr[j] + (stv.x[i]-a[j]) * (gr[j+1]-gr[j])/(a[j+1]-a[j])
							ELSE
								stv.orgGm[i]:=0; {this layer might absorb, but it does not generate electron-hole pairs}
							{we're using linear interpolation here}
							i:=i+1
						END;  
			END {case 2, reading generation profile from file}
	END	
END;

PROCEDURE Update_Generation_Profile(org: vector; VAR new : vector; G_frac : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Rescales generation profile (org) by factor G_frac to obtain new profile}
VAR i : INTEGER;
BEGIN
	FOR i:=0 TO par.NP+1 DO 
		new[i]:=org[i]*G_frac {Generation rate already set in original profile}
END;

PROCEDURE Update_Gen_Pot(V : vector; VAR Vgn, Vgp : vector; CONSTREF stv : TstaticVars; CONSTREF par : TInputParameters);
{updates the generalised potentials (Vgn, Vgp) after potential V has been altered.}
VAR i, j : INTEGER;
	facDOS : myReal;
BEGIN
	{we loop over the layers:}
	FOR j:=1 TO stv.NLayers DO
	WITH par.lyr[j] DO
		FOR i:=stv.i0[j] TO stv.i1[j] DO 
		BEGIN
			facDOS:=stv.Vt*LN(stv.NcLoc[i]/stv.NcLoc[0]);
			Vgn[i]:=V[i] + E_c - stv.E_CB[0] + facDOS;
			Vgp[i]:=V[i] + E_v - stv.E_VB[0] - facDOS
		END		
END;

PROCEDURE Init_Trap_Filling_Arrays(VAR f_tb, f_ti, f_ti_numer, f_ti_inv_denom : TTrapArray; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{sets the lengths of the trap filling arrays and makes all elements zero.}
VAR e, i, j, intLevels : INTEGER;
BEGIN
	
	{first the bulk:}
	FOR j:=1 TO stv.NLayers DO
		FOR e:=1 TO stv.Ntb[j].NLevels DO {loop over trap levels}
			FOR i:=MAX(1,stv.i0[j]) TO MIN(par.NP,stv.i1[j]) DO {loop over grid points in this layer, but as always, exclude the electrodes}
			BEGIN
				SETLENGTH(f_tb[i], stv.Ntb[j].NLevels + 1); {set length of this part of the array, +1 as we start counting at 1}
				FILLCHAR(f_tb[i,0], LENGTH(f_tb[i]) * SIZEOF(f_tb[i,0]), 0)
			END;
	
	{the number of interface traps in a layer can be different at its left and right interface.}
	{So, we take the number of levels equal to the max of all interfaces:}
	intLevels:=0;
	FOR j:=1 TO stv.NLayers-1 DO
		intLevels:=MAX(intLevels, stv.Nti[j].NLevels);
	
	{now init length:}
	FOR i:=0 TO par.NP+1 DO
	BEGIN
		SETLENGTH(f_ti[i], intLevels+1); 
		SETLENGTH(f_ti_numer[i], intLevels+1); 
		SETLENGTH(f_ti_inv_denom[i], intLevels+1); 
		FILLCHAR(f_ti[i,0], LENGTH(f_ti[i]) * SIZEOF(f_ti[i,0]), 0);
		FILLCHAR(f_ti_numer[i,0], LENGTH(f_ti_numer[i]) * SIZEOF(f_ti_numer[i,0]), 0);
		FILLCHAR(f_ti_inv_denom[i,0], LENGTH(f_ti_inv_denom[i]) * SIZEOF(f_ti_inv_denom[i,0]), 0);
	END
END;

PROCEDURE Init_Pot_Dens_Ions_Traps(VAR V, Vgn, Vgp, n, p, nion, pion : vector; VAR f_tb, f_ti, f_ti_numer, f_ti_inv_denom : TTrapArray; Va : myReal; VAR stv : TStaticVars; CONSTREF par : TInputParameters);
{init. for V, Vgn,p, n, p, ion densities at bias voltage Va. Also sets lengths of f_tb/i arrays}
VAR i, j : INTEGER;
	locEF : myReal;
BEGIN
    FOR i:=0 TO par.NP+1 DO {guess new V in interior of device} {linear interpolation of V}
	    V[i]:=stv.V0 - Va/2 + stv.x[i]*(stv.VL-stv.V0+Va)/stv.Ltot;
	Update_Gen_Pot(V, Vgn, Vgp, stv, par); {update generalised potentials}

    FOR i:=0 TO par.NP+1 DO {guess n, p in interior of device}
    BEGIN
		locEF:=par.W_L - stv.x[i]*(par.W_L - par.W_R)/stv.Ltot; {estimate of local EF, simple interpolation of W_L to W_R}
		{this also sets the boundary conditions on n and p, they matter, the rest is not critical}
		n[i]:=stv.NcLoc[i]*EXP(-stv.Vti*(locEF - stv.E_CB[i]));
		p[i]:=SQR(stv.ni[i])/n[i];
		nion[i]:=0; {just to make sure these have been initialised}
		pion[i]:=0; {later on, we will set them to the correct values}   
    END;
   
    {Last but not least: the ions}
    {identify regions (can be multiple adjacent layers) that contain either pos or neg ions:}
	Init_Ionic_Region(stv.NegIonRegion, -1, stv, par);
	Init_Ionic_Region(stv.PosIonRegion, 1, stv, par);
	
	{now we know the ionic regions, we must indicate (per LAYER) if a layer contains pos/neg ions. Those ions can come from adjacent layers!}
	{first: negative ions:}
	FOR j:=0 TO LENGTH(stv.NegIonRegion)-1 DO 
		FOR i:=stv.lid[stv.NegIonRegion[j].istart] TO stv.lid[stv.NegIonRegion[j].ifinish] DO	
			{neg/posIons only true if layer is 1) part of an Ionic Region, and 2) the average concentration > 0}
			par.lyr[i].negIons:=stv.NegIonRegion[j].AvC>0;
	{and postive ones:}
	FOR j:=0 TO LENGTH(stv.PosIonRegion)-1 DO 
		FOR i:=stv.lid[stv.PosIonRegion[j].istart] TO stv.lid[stv.PosIonRegion[j].ifinish] DO	
			{neg/posIons only true if layer is 1) part of an Ionic Region, and 2) the average concentration > 0}
			par.lyr[i].posIons:=stv.PosIonRegion[j].AvC>0;
	
	{now give the arrays noin/pion and mu_n/p_ion the right (starting) values:}
	FOR i:=0 TO par.NP+1 DO
	BEGIN
		nion[i]:=par.lyr[stv.lid[i]].N_anion;
		pion[i]:=par.lyr[stv.lid[i]].N_cation;
		stv.mu_n_ion[i]:=par.lyr[stv.lid[i]].mu_anion;
		stv.mu_p_ion[i]:=par.lyr[stv.lid[i]].mu_cation
	END;
	
	Init_Trap_Filling_Arrays(f_tb, f_ti, f_ti_numer, f_ti_inv_denom, stv, par); {inits their lengths and sets all elements to zero}

END;


PROCEDURE Init_Traps_From_File(VAR log : TEXT; VAR dist : TTrapDistLayer; TrapFile : ANSISTRING; trapType : INTEGER);
{Inits traps if their energy levels are specified in a file}
VAR e, NumLevels : INTEGER;
	x, y : Row;
BEGIN
	WRITELN('Reading trap profile from ',TrapFile);
	WRITELN(log, 'Reading trap profile from ',TrapFile);
	
	{Read data from TrapFile and put them in local vars x and y:}
	Read_XY_Table(x, y, TrapFile, 'E Ntrap', NumLevels);
	
	{check the number of levels:}
	IF NumLevels=0 THEN Stop_Prog('Could not find any trap levels in file '+TrapFile+'.', EC_InvalidInput);

	{now copy data into dist}
	WITH dist DO BEGIN	
		NLevels:=NumLevels;
		IF trapType = 1 THEN cwe:=1; {otherwise it remains zero, as initialised}
		SETLENGTH(en, NLevels + 1); {note: +1 as we start counting at index 1} 
		SETLENGTH(Nt, NLevels + 1); {note: +1 as we start counting at index 1} 
		SETLENGTH(nt0, NLevels + 1); {note: +1 as we start counting at index 1} 
		SETLENGTH(pt0, NLevels + 1); {note: +1 as we start counting at index 1} 
		FOR e:=1 TO NumLevels DO
		BEGIN
			en[e]:=x[e-1]; {note: en and Nt start at index 1, but x and y start at 0}
			Nt[e]:=y[e-1];
			IF Nt[e]=0 THEN Stop_Prog('Found a trap level with zero density in file '+TrapFile+' level number '+IntToStr(e), EC_InvalidInput)
		END
	END {with dist}
END;


PROCEDURE Init_Trapping(VAR log : TEXT; VAR stv : TStaticVars; CONSTREF par : TInputParameters); 
{Inits all variables needed for trapping and SRH recombination}
{note: stv are changed here (Ntb and Nti), so they are VAR parameters}
VAR e, j : INTEGER;
BEGIN
	WITH stv DO 
	BEGIN
		{set length of arrays that contain energies and density of traps, we start numbering at 1:}
		SETLENGTH(Ntb, NLayers+1);
		SETLENGTH(Nti, NLayers+1); 

		{and set all to zero}
		FILLCHAR(Ntb[0], SIZEOF(Ntb[0]) * LENGTH(Ntb), 0); 
		FILLCHAR(Nti[0], SIZEOF(Nti[0]) * LENGTH(Nti), 0)	
	END; {with stv}

	{First: init, check bulk trapping:}
	FOR j:=1 TO stv.NLayers DO {loop over all layers}
	BEGIN
		{first: init bulk traps:}
		WITH par.lyr[j] DO BEGIN
			IF bulkTrapFromFile THEN
				Init_Traps_From_File(log, stv.Ntb[j], bulkTrapFile, bulkTrapType)			
			ELSE IF N_t_bulk <> 0 THEN 
			WITH stv.Ntb[j] DO BEGIN {if not from file, we take 1 level}
				NLevels:=1;
				{now set cwe: charge when empty}
				IF bulkTrapType = 1 THEN cwe:=1; {otherwise it remains zero, as initialised}
				SETLENGTH(en, 2); {note: 2 as we start counting at index 1} 
				SETLENGTH(Nt, 2); {note: 2 as we start counting at index 1} 
				SETLENGTH(nt0, 2); {note: 2 as we start counting at index 1} 
				SETLENGTH(pt0, 2); {note: 2 as we start counting at index 1} 
				en[1]:=E_t_bulk;
				Nt[1]:=N_t_bulk
			END;			
		END; {with statement}

		{check validity of levels and init nt0 and pt0}
		WITH stv.Ntb[j] DO
			FOR e:=1 to NLevels DO 
			BEGIN {loop over all energy levels, if any}
				{ALL bulk trap levels of a layers should fall in E_c and E_v}				
				IF (Nt[e]>0) AND ((en[e] >= par.lyr[j].E_v) OR (en[e] <= par.lyr[j].E_c)) {only check if there are traps!}
					THEN Stop_Prog('Found a bulk trap energy in layer '+IntToStr(j)+' that does not fall in the band gap.', EC_InvalidInput);
		
				{init nt0 and pt0, needed for SRH recombination}
				nt0[e]:=par.lyr[j].N_c * EXP((par.lyr[j].E_c-en[e])*stv.Vti);
				pt0[e]:=par.lyr[j].N_c * EXP((en[e]-par.lyr[j].E_v)*stv.Vti)
			END
		
	END; {for loop over all layers for bulk traps}

	{now the interface traps, so loop until NLayers-1}
	FOR j:=1 TO stv.NLayers-1 DO 
	BEGIN
		
		{first: init the traps:}
		WITH par.lyr[j] DO BEGIN
			IF intTrapFromFile THEN
				Init_Traps_From_File(log, stv.Nti[j], intTrapFile, intTrapType)			
			ELSE IF N_t_int <> 0 THEN
			WITH stv.Nti[j] DO BEGIN{if not from file, we take 1 level}
				NLevels:=1;
				{now set cwe: charge when empty}
				IF intTrapType = 1 THEN cwe:=1; {otherwise it remains zero, as initialised}
				SETLENGTH(en, 2); {note: 2 as we start counting at index 1} 
				SETLENGTH(Nt, 2); {note: 2 as we start counting at index 1} 
				en[1]:=E_t_int;
				Nt[1]:=N_t_int
			END
				
		END; {with statement}

		{check validity of levels }		
		WITH stv.Nti[j] DO
			FOR e:=1 to NLevels DO {loop over all energy levels, if any}
				{ALL interface trap levels of a layers should fall in E_c and E_v of BOTH layers adjacent this interface}
				IF (Nt[e]>0) AND ((en[e] >= MIN(par.lyr[j].E_v, par.lyr[j+1].E_v)) OR (en[e] <= MAX(par.lyr[j].E_c, par.lyr[j+1].E_c))) {only check if there are traps!}
					THEN Stop_Prog('Found interface trap energy between layers '+IntToStr(j)+' and '+IntToStr(j+1)+' that does not fall between the E_c and E_v of those layers.', EC_InvalidInput);

		
	END; {for loop over all layers -1 for interface traps}

	{now we have all the parameters of the interface traps, we can init nt0_L/R and pt0_L/R}
	FOR j:=1 TO stv.NLayers-1 DO 
	BEGIN
		{each layer has 2 interfaces, left and right and we need to set the correct length:}
		SETLENGTH(stv.Nti[j].nt0_R, stv.Nti[j].NLevels+1);
		SETLENGTH(stv.Nti[j+1].nt0_L, stv.Nti[j].NLevels+1); {we copy the length of the array in next layer}
		SETLENGTH(stv.Nti[j].pt0_R, stv.Nti[j].NLevels+1);
		SETLENGTH(stv.Nti[j+1].pt0_L, stv.Nti[j].NLevels+1);

		FOR e:=1 TO stv.Nti[j].NLevels DO
			WITH stv.Nti[j] DO
			BEGIN
				{init nt0 and pt0, at left and right side of this layer:}
				nt0_R[e]:=par.lyr[j].N_c * EXP((par.lyr[j].E_c-en[e])*stv.Vti);
				pt0_R[e]:=par.lyr[j].N_c * EXP((en[e]-par.lyr[j].E_v)*stv.Vti);
				{at the left side, we work with layer j+1:}
				stv.Nti[j+1].nt0_L[e]:=par.lyr[j+1].N_c * EXP((par.lyr[j+1].E_c-en[e])*stv.Vti);
				stv.Nti[j+1].pt0_L[e]:=par.lyr[j+1].N_c * EXP((en[e]-par.lyr[j+1].E_v)*stv.Vti)
			END
	END
	
END;


PROCEDURE Rescale_Ion_Density(VAR ion : vector; istart, ifinish : INTEGER; conc : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{normalises an ion distribution to the correct overall concentration conc}
VAR i : INTEGER;
	sum, norm : myReal;
BEGIN
	sum:=0;
	
	FOR i:=istart+1 TO ifinish DO {note: we start at istart + 1 as we access ion[i-1]}
		sum:=sum + 0.5*(ion[i]+ion[i-1])*stv.h[i-1]*stv.Ltot; {if the grid is non-uniform, we need to take this into account}

	norm:=conc*(stv.x[ifinish]-stv.x[istart])/sum; 
	FOR i:=istart TO ifinish DO {now renormalise the ionic densities such that the total number of ions is correct}
		ion[i]:=ion[i]*norm;
END;

PROCEDURE Calc_Ion_Distribution_Steady_State(VAR ion, V : vector; sn : ShortInt; IonRegion : TIonicRegions; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the ion distribution in steady state. This is either constant 
(if the ionic species doesn't move), or it follows from requiring that
the associated ionic particle current be zero.} 
VAR fac : myReal;
	i, j: INTEGER;
BEGIN
	IF ABS(sn)<>1 THEN Stop_Prog('Calc_Ion_Distribution_Steady_State called with invalid sn (= sign of ions, +-1).',EC_ProgrammingError);
	
	{we loop over all regions where that contain ions, so this can be multiple layers!}
	FOR j:=0 TO LENGTH(IonRegion)-1 DO 
		WITH IonRegion[j] DO
		BEGIN
			ion[istart]:=AvC;
			FOR i:=istart+1 TO ifinish DO
			BEGIN
				IF sn=1 THEN
					fac:=B(stv.Vti*(V[i]-V[i-1])) / B(stv.Vti*(V[i-1]-V[i])) 
				ELSE
					fac:=B(stv.Vti*(V[i-1]-V[i])) / B(stv.Vti*(V[i]-V[i-1]));
				{we use the expressions for the ion currents to make sure that they are zero. This yields
				the profile of neg/pos ions.}
				IF AvC>0 THEN ion[i]:=ion[i-1]*fac ELSE ion[i]:=0
			END;
			
			{nomalize concentrations}
			IF AvC>0 THEN Rescale_Ion_Density(ion, istart, ifinish, IonRegion[j].AvC, stv, par);
	
		END {with stv.IonRegion[j]}
	
END;

PROCEDURE Solve_Neg_Ions(VAR nion : vector; nionPrevTime, V : vector; dti : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{dti: DeltaTInverse, so dti = 1/delt. If dti=0 then we're in steady-state.
If dti > 0 we're using the transient equations (Selberherr 6-4.32)}
VAR i, j, istart, ifinish : INTEGER;
    lo, m, u, rhs : vector;
BEGIN
	{we loop over the regions that can contain ions. A region can be multiple adjacent layers.
	Per such region, we solve the drift-diffusion equation and compute the ionic distribution.}
	
	FOR j:=0 TO LENGTH(stv.NegIonRegion)-1 DO BEGIN {loop over all regions that can contain ions:}
		istart:=stv.NegIonRegion[j].istart;
		ifinish:=stv.NegIonRegion[j].ifinish;
		
		{first: we do the left interface: either at the contact or at some adjacent layer}
		IF istart<>0 THEN
		WITH stv DO 
		BEGIN {there's a layer left of this region that is impenetrable to ions}
			i:=istart; {i is just shorthand notation for istart}
			rhs[i]:=-0.5 * nionPrevTime[i]*dti*SQR(stv.Ltot)*h[i]*(h[i]+h[i-1]);
			lo[i]:=0;
			m[i]:=-Vt*mu_n_ion[i]*B(Vti*(V[i]-V[i+1])) -0.5*dti*SQR(Ltot)*h[i]*(h[i]+h[i-1]);
			u[i]:=Vt*mu_n_ion[i]*B(Vti*(V[i+1]-V[i]));   
		END
		ELSE BEGIN {ions can move towards the contacts, nion[0]<>0:}
			{we need to ensure that the ionic currents into/out of the contacts be zero, so let's do that first:}
			lo[0]:=0;
			m[0]:=1;
			u[0]:=-B(stv.Vti*(V[1]-V[0])) / B(stv.Vti*(V[0]-V[1]));
			rhs[0]:=0 
		END;
	
		{now we do the right interface. Again: either at the contact or an adjacent layer}
		IF ifinish <> par.NP+1 THEN
		WITH stv DO
		BEGIN  {ions cannot move towards the contacts}
			i:=ifinish;  {i is just shorthand notation for ifinish}
			rhs[i]:=-0.5 * nionPrevTime[i]*dti*SQR(Ltot)*h[i-1]*(h[i]+h[i-1]);
			lo[i]:=Vt*mu_n_ion[i-1]*B(Vti*(V[i-1]-V[i]));
			m[i]:=-Vt*mu_n_ion[i-1]*B(Vti*(V[i]-V[i-1]))-0.5*dti*SQR(Ltot)*h[i-1]*(h[i]+h[i-1]);
			u[i]:=0
		END 
		ELSE BEGIN {ions can move towards the contacts, nion[NP+1]<>0:}
			lo[ifinish]:=-B(stv.Vti*(V[par.NP]-V[par.NP+1])) / B(stv.Vti*(V[par.NP+1]-V[par.NP]));
			m[ifinish]:=1;
			u[ifinish]:=0;
			rhs[ifinish]:=0
		END;

		{now set the interior part:}
		FOR i:=istart+1 TO ifinish-1 DO  {continuity eq. in matrix vorm, equivalent to that of n and p, but without generation and recombination}
		WITH stv DO
		BEGIN
			rhs[i]:=-0.5 * nionPrevTime[i]*dti*SQR(Ltot)*h[i]*h[i-1]*(h[i]+h[i-1]);
			lo[i]:=h[i]*Vt*mu_n_ion[i-1]*B(Vti*(V[i-1]-V[i]));
			m[i]:=-(h[i-1]*Vt*mu_n_ion[i]*B(Vti*(V[i]-V[i+1])) +
				h[i]*Vt*mu_n_ion[i-1]*B(Vti*(V[i]-V[i-1])))
				-0.5*dti*SQR(Ltot)*h[i]*h[i-1]*(h[i]+h[i-1]);
			u[i]:=h[i-1]*Vt*mu_n_ion[i]*B(Vti*(V[i+1]-V[i]));          
		END;

		{Solve nion from istart to ifinish:}
		Tridiag(nion, lo, m, u, rhs, istart, ifinish);
	
		{now check if nion is still well behaved, i.e. positive, and has the correct overall value}
		FOR i:=istart TO ifinish DO
			IF (nion[i]<0) THEN 
			BEGIN
				IF par.ignoreNegDens 
					THEN nion[i]:=-nion[i] 
					ELSE Stop_Prog('Negative concentration of negative ions encountered!', EC_NumericalFailure)
			END;
		
		{make sure the number of ions is preserved, i.e. correct:}
		Rescale_Ion_Density(nion, istart, ifinish, stv.NegIonRegion[j].AvC, stv, par)
	
	END; {loop over all regions that can contain ions}
	
END;

PROCEDURE Solve_Pos_Ions(VAR pion : vector; pionPrevTime, V : vector; dti : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{dti: DeltaTInverse, so dti = 1/delt. If dti=0 then we're in steady-state.
If dti > 0 we're using the transient equations (Selberherr 6-4.33)}
VAR i, j, istart, ifinish : INTEGER;
    lo, m, u, rhs : vector;
BEGIN
	{we loop over the regions that can contain ions. A region can be multiple adjacent layers.
	Per such region, we solve the drift-diffusion equation and compute the ionic distribution.}
	
	FOR j:=0 TO LENGTH(stv.PosIonRegion)-1 DO BEGIN {loop over all regions that can contain ions:}
		istart:=stv.PosIonRegion[j].istart;
		ifinish:=stv.PosIonRegion[j].ifinish;
		
		{first: we do the left interface: either at the contact or at some adjacent layer}
		IF istart<>0 THEN
		WITH stv DO 
		BEGIN {there's a layer left of this region that is impenetrable to ions}
			i:=istart; {i is just shorthand notation for istart}
			rhs[i]:=-0.5 * pionPrevTime[i]*dti*SQR(stv.Ltot)*h[i]*(h[i]+h[i-1]);
			lo[i]:=0;
			m[i]:=-Vt*mu_p_ion[i]*B(Vti*(V[i+1]-V[i])) -0.5*dti*SQR(Ltot)*h[i]*(h[i]+h[i-1]);
			u[i]:=Vt*mu_p_ion[i]*B(Vti*(V[i]-V[i+1]));   
		END
		ELSE BEGIN {ions can move towards the contacts, pion[0]<>0:}
			{we need to ensure that the ionic currents into/out of the contacts be zero, so let's do that first:}
			lo[0]:=0;
			m[0]:=1;
			u[0]:=-B(stv.Vti*(V[0]-V[1])) / B(stv.Vti*(V[1]-V[0]));
			rhs[0]:=0 
		END;
	
		{now we do the right interface. Again: either at the contact or an adjacent layer}
		IF ifinish <> par.NP+1 THEN
		WITH stv DO
		BEGIN  {ions cannot move towards the contacts}
			i:=ifinish;  {i is just shorthand notation for ifinish}
			rhs[i]:=-0.5 * pionPrevTime[i]*dti*SQR(Ltot)*h[i-1]*(h[i]+h[i-1]);
			lo[i]:=Vt*mu_p_ion[i-1]*B(Vti*(V[i]-V[i-1]));
			m[i]:=-Vt*mu_p_ion[i-1]*B(Vti*(V[i-1]-V[i]))-0.5*dti*SQR(Ltot)*h[i-1]*(h[i]+h[i-1]);
			u[i]:=0
		END 
		ELSE BEGIN {ions can move towards the contacts, pion[NP+1]<>0:}
			lo[ifinish]:=-B(stv.Vti*(V[par.NP+1]-V[par.NP])) / B(stv.Vti*(V[par.NP]-V[par.NP+1]));
			m[ifinish]:=1;
			u[ifinish]:=0;
			rhs[ifinish]:=0
		END;

		{now set the interior part:}
		FOR i:=istart+1 TO ifinish-1 DO  {continuity eq. in matrix vorm, equivalent to that of n and p, but without generation and recombination}
		WITH stv DO
		BEGIN
			rhs[i]:=-0.5 * pionPrevTime[i]*dti*SQR(Ltot)*h[i]*h[i-1]*(h[i]+h[i-1]);
			lo[i]:=h[i]*Vt*mu_p_ion[i-1]*B(Vti*(V[i]-V[i-1]));
			m[i]:=-(h[i-1]*Vt*mu_p_ion[i]*B(Vti*(V[i+1]-V[i])) +
				h[i]*Vt*mu_p_ion[i-1]*B(Vti*(V[i-1]-V[i])))
				-0.5*dti*SQR(Ltot)*h[i]*h[i-1]*(h[i]+h[i-1]);
			u[i]:=h[i-1]*Vt*mu_p_ion[i]*B(Vti*(V[i]-V[i+1]));          
		END;

		{Solve pion from istart to ifinish:}
		Tridiag(pion, lo, m, u, rhs, istart, ifinish);
	
		{now check if nion is still well behaved, i.e. positive, and has the correct overall value}
		FOR i:=istart TO ifinish DO
			IF (pion[i]<0) THEN 
			BEGIN
				IF par.ignoreNegDens 
					THEN pion[i]:=-pion[i] 
					ELSE Stop_Prog('Negative concentration of positive ions encountered!', EC_NumericalFailure)
			END;
		
		{make sure the number of ions is preserved, i.e. correct:}
		Rescale_Ion_Density(pion, istart, ifinish, stv.PosIonRegion[j].AvC, stv, par)
	
	END; {loop over all regions that can contain ions}
	
END;


FUNCTION Calc_f_ti_Numer(CONSTREF n, p : vector; ii, j, e : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : myReal;
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculates the numerator of the fraction of interface traps filled. 
The numerator and denominator are split because both are re-used in different places in the Poisson and continuity equations.
ii: grid point left of interface, i.e. interface sits between ii and ii+1}
BEGIN
	Calc_f_ti_Numer:=stv.Nti[j].Nt[e] * (par.lyr[j].C_n_int*(n[ii]+n[ii+1]) + par.lyr[j].C_p_int*(stv.Nti[j].pt0_R[e]+stv.Nti[j+1].pt0_L[e]))	
END;

FUNCTION Calc_Inv_f_ti_Denom(CONSTREF n, p : vector; ii, j, e : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : myReal;
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculates the inverse of the denominator of the fraction of interface traps filled. 
The numerator and denominator are split because both are re-used in different places in the Poisson and continuity equations.
ii: grid point left of interface, i.e. interface sits between ii and ii+1}
VAR denom_trap_occup : myReal;
BEGIN
	denom_trap_occup:=par.lyr[j].C_n_int*stv.Nti[j].Nt[e]*(n[ii] + n[ii+1] + stv.Nti[j].nt0_R[e] + stv.Nti[j+1].nt0_L[e]) 
					+ par.lyr[j].C_p_int*stv.Nti[j].Nt[e]*(p[ii] + p[ii+1] + stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]);

	{finally, we need the inverse:}
	Calc_Inv_f_ti_Denom:=1/denom_trap_occup
END;


FUNCTION Calc_f_tb(n, p, Cn, Cp, nt0, pt0, Old_f_tb, dti : myReal) : myReal;
{Computes f_tb, the fill level of a bulk trap level, in a single grid point}
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
VAR a, b : myReal;
BEGIN
	IF dti=0 THEN {steady-state}
		Calc_f_tb:=(Cn*n + Cp*pt0) / (Cn*(n+nt0) + Cp*(p+pt0))
	ELSE BEGIN {transient}
	    b:= Cn*n + Cp*pt0;
	    a:= Cn*n + Cn*nt0 + Cp*pt0 + Cp*p;
	    Calc_f_tb:= b/a + (Old_f_tb - b/a)*EXP(-a/dti)
		{we solve f_tb from d(f_tb)/dt = a f_tb + b, which is a trivial ODE}
	END
END;

PROCEDURE Calc_Trap_Filling_Charge(VAR f_tb, f_ti, f_ti_numer, f_ti_inv_denom : TTrapArray; VAR Ntb_charge, Nti_charge : vector; CONSTREF n, p : vector; Old_f_tb, Old_f_ti : TTrapArray; dti : myReal;
									CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the filling (f_tb/i) and charge (Ntb/i_charge) of the (bulk and interface) traps at every gridpoint.}
{See Koopmans and Koster, Sol. RRL 2022, 2200560.}
VAR i, ii, j, e : INTEGER;
	f_ti_level_numer, f_ti_level_inv_denom : myReal;
	a, b, int_trap_width  : myReal;
BEGIN
	
	{First, the bulk traps, so loop over all layers:}
	FILLCHAR(Ntb_charge, SIZEOF(Ntb_charge), 0); {reset array to zero}

	FOR j:=1 TO stv.NLayers DO
		FOR e:=1 TO stv.Ntb[j].NLevels DO {loop over trap levels}
			FOR i:=MAX(1,stv.i0[j]) TO MIN(par.NP,stv.i1[j]) DO {loop over grid points in this layer, but as always, exclude the electrodes}
			BEGIN
				{first: f_tb: the fraction of traps that contains an electron}
				f_tb[i,e]:=Calc_f_tb(n[i], p[i], par.lyr[j].C_n_bulk, par.lyr[j].C_p_bulk, stv.Ntb[j].nt0[e], stv.Ntb[j].pt0[e], Old_f_tb[i,e], dti);
				{next the charge in the level:}
				Ntb_charge[i]:=Ntb_charge[i] + ABS(par.lyr[j].bulkTrapType) * (stv.Ntb[j].cwe	- f_tb[i,e]) * stv.Ntb[j].Nt[e]
				{note that we have f_tb minus charges or 1-f_tb positive charges}			
			END;
		
	
	{next: interfaces}
	
	{first, init all arrays:}
	FILLCHAR(Nti_charge, SIZEOF(Nti_charge), 0); {reset array to zero}

	{now f_ti arrays set to zero:}
	FOR i:=0 TO par.NP+1 DO
	BEGIN
		FILLCHAR(f_ti[i,0], SIZEOF(f_ti[i,0]) * LENGTH(f_ti[i]), 0); {reset array to zero}
		FILLCHAR(f_ti_numer[i,0], SIZEOF(f_ti_numer[i,0]) * LENGTH(f_ti_numer[i]), 0); {reset array to zero}
		FILLCHAR(f_ti_inv_denom[i,0], SIZEOF(f_ti_inv_denom[i,0]) * LENGTH(f_ti_inv_denom[i]), 0); {reset array to zero}
	END;	
	
	FOR j:=1 TO stv.NLayers-1 DO {loop over interfaces. Each interface involves 2 points: i1[j] and i1[j]+1 (in the adjacent layer)}
		FOR e:=1 TO stv.Nti[j].NLevels DO {loop over trap levels}
		BEGIN
			IF dti=0 THEN
			BEGIN {steady-state}
				{note: the interface sits between i1[j] and i1[j]+1}
				ii:=stv.i1[j]; {ii: i interface}

				{calc the numerator and inverse of denominator of single interface level:}
				f_ti_level_numer:=Calc_f_ti_Numer(n, p, ii, j, e, stv, par);	
				f_ti_level_inv_denom:=Calc_Inv_f_ti_Denom(n, p, ii, j, e, stv, par);
				{now fill the arrays with the correct values, just left and right of each interface:}
				f_ti[ii,e]:=f_ti_level_numer * f_ti_level_inv_denom;
				f_ti[ii+1,e]:=f_ti[ii,e]; {copy from value just left of interface}
				f_ti_numer[ii,e]:=f_ti_level_numer;
				f_ti_numer[ii+1,e]:=f_ti_numer[ii,e]; {copy from value just left of interface}
				f_ti_inv_denom[ii,e]:=f_ti_level_inv_denom;
				f_ti_inv_denom[ii+1,e]:=f_ti_inv_denom[ii,e]; {copy from value just left of interface}
			END
			ELSE 
			BEGIN {transient}
				ii:=stv.i1[j]; {ii: i interface, so ii is last point before interface}
			
				a:=par.lyr[j].C_n_int*(n[ii] + n[ii+1] + stv.Nti[j].nt0_R[e] + stv.Nti[j+1].nt0_L[e])
					+ par.lyr[j].C_p_int*(p[ii] + p[ii+1] + stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]);
				b:=par.lyr[j].C_n_int*(n[ii] + n[ii+1]) + par.lyr[j].C_p_int*(stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]);
						
				f_ti[ii,e]:= b/a + (Old_f_ti[ii,e] - b/a)*EXP(-a/dti);
				f_ti[ii+1,e]:=f_ti[ii,e];
				{we solve f_ti from d(f_ti)/dt = a f_ti + b, which is a trivial ODE}
			END;

			{next the charge in the level:}

			{input N_t_int is #traps/area, will convert this to #traps/volume by computing the width of the interface:}
			int_trap_width:=0.5 * (stv.Ltot*(0.5*stv.h[ii+1] + stv.h[ii] + 0.5*stv.h[ii-1]));
					
			{each side of the interface hosts half the traps, so we get an extra factor of 0.5:}
			Nti_charge[ii]:=Nti_charge[ii] + ABS(par.lyr[j].intTrapType) * (stv.Nti[j].cwe	- f_ti[ii,e]) * 0.5 * stv.Nti[j].Nt[e] / int_trap_width;
			Nti_charge[ii+1]:=Nti_charge[ii]
			{note that we have f_ti minus charges or 1-f_ti positive charges}			
		END
END;


PROCEDURE Calc_Linearization_f_tb(VAR lin : TLinFt; CONSTREF n, p : vector; dti : myReal; 
                                  CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters); 
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculate the linearization of the filled (with electrons) fraction of bulk 
traps with respect to the potential. This function is written in such a way, 
that the result can be added to the Poisson equation regardless of the presence of traps. 
If traps are not present the elements added to the Poisson equation simply become zero.
See Koopmans and Koster, Sol. RRL 2022, 2200560.}
VAR f_tb_numer, f_tb_inv_denom, f_tb_m_level : myReal;
	i, j, e : INTEGER;
BEGIN	
	{Note: before calling this routine, we set all fields in lin to zero}

	{Bulk traps: if we have them, calculate their linearization in delta V.}
	FOR j:=1 TO stv.NLayers DO
		IF (dti=0) AND (par.lyr[j].bulkTrapType <> 0) THEN
			FOR e:=1 TO stv.Ntb[j].NLevels DO {loop over trap levels}
				FOR i:=MAX(1,stv.i0[j]) TO MIN(par.NP,stv.i1[j]) DO {loop over grid points in this layer, but as always, exclude the electrodes}
				BEGIN	
					f_tb_numer:=par.lyr[j].C_n_bulk*n[i] + par.lyr[j].C_p_bulk * stv.Ntb[j].pt0[e];
					f_tb_inv_denom:=1 / (par.lyr[j].C_n_bulk*(n[i]+stv.Ntb[j].nt0[e]) + par.lyr[j].C_p_bulk*(p[i] + stv.Ntb[j].pt0[e]));
					f_tb_m_level:=par.lyr[j].C_n_bulk * n[i+1] * f_tb_inv_denom;
					f_tb_m_level:=f_tb_m_level - (par.lyr[j].C_n_bulk * n[i] - par.lyr[j].C_p_bulk * p[i]) * f_tb_numer * SQR(f_tb_inv_denom);
					f_tb_m_level:=par.lyr[j].bulkTrapType*(stv.Ntb[j].cwe - f_tb_m_level) * stv.Ntb[j].Nt[e];					
					lin.f_tb_m[i]:=lin.f_tb_m[i] + f_tb_m_level {add the f_tb_m of this level to the overal sum}
				END
END;


PROCEDURE Calc_Linearization_f_ti(VAR lin : TLinFt; CONSTREF n, p : vector; dti : myReal;
                                  CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters); 
{This routine is not visible outside this unit, so it is not listed in INTERFACE.}
{Calculate the linearization of the filled (with electrons) fraction of interface
traps with respect to the potential. This function is written in such a way, 
that the result can be added to the Poisson equation regardless of the presence of traps. 
If traps are not present the elements added to the Poisson equation simply become zero.
See Koopmans and Koster, Sol. RRL 2022, 2200560.}
VAR j, e, ii, q_tr : INTEGER;
	dum, f_ti_level_numer, f_ti_level_inv_denom : myReal;
BEGIN	
			
	{Linearize f_ti (interface trap occupance fraction) in delV, a derivation can be found in the docs.}
	{Key here is that all linearized elements become zero when they should not affect the Poisson equation}

	{Note: before calling this routine, we set all fields in lin to zero}

	FOR j:=1 TO stv.NLayers-1 DO
		IF (dti=0) AND (stv.Traps_int_poisson) THEN
			FOR e:=1 TO stv.Nti[j].NLevels DO {loop over trap levels}
			BEGIN	
				ii:=stv.i1[j]; {ii: i interface}
				q_tr:=ABS(par.lyr[j].intTrapType);
				{calc the numerator and inverse of denominator of single interface level:}
				f_ti_level_numer:=Calc_f_ti_Numer(n, p, ii, j, e, stv, par);	
				f_ti_level_inv_denom:=Calc_Inv_f_ti_Denom(n, p, ii, j, e, stv, par);

				{lower diagonal element:}
				dum:=par.lyr[j].C_n_int * n[ii-1] * stv.Nti[j].Nt[e] * stv.Vti * f_ti_level_inv_denom;
				dum:=dum - (par.lyr[j].C_n_int * n[ii-1] - par.lyr[j].C_p_int * p[ii-1]) * stv.Nti[j].Nt[e] * stv.Vti * f_ti_level_numer * SQR(f_ti_level_inv_denom);
				lin.f_ti_lo[ii]:=lin.f_ti_lo[ii] + q_tr*dum * 0.5 * stv.Nti[j].Nt[e];
				
				{the factor stv.Vti is multiplied in the poisson equation with for this main diagonal}
				{main diagonal:}
				dum:= par.lyr[j].C_n_int * n[ii] * stv.Nti[j].Nt[e] * f_ti_level_inv_denom;
				dum:=dum - (par.lyr[j].C_n_int * n[ii] - par.lyr[j].C_p_int * p[ii]) * stv.Nti[j].Nt[e] * f_ti_level_numer * SQR(f_ti_level_inv_denom);
				lin.f_ti_m[ii]:=lin.f_ti_m[ii] + q_tr*dum * 0.5 * stv.Nti[j].Nt[e];
			
				{upper diagonal element:}
				dum:=par.lyr[j].C_n_int * n[ii+1] * stv.Nti[j].Nt[e] * stv.Vti * f_ti_level_inv_denom;
				dum:=dum - (par.lyr[j].C_n_int * n[ii+1] - par.lyr[j].C_p_int * p[ii+1]) * stv.Nti[j].Nt[e] * stv.Vti * f_ti_level_numer * SQR(f_ti_level_inv_denom);
				lin.f_ti_up[ii]:=lin.f_ti_up[ii] + q_tr*dum * 0.5 * stv.Nti[j].Nt[e]
			
			END

END;

PROCEDURE Calc_Linearization_f_t_All(VAR lin : TLinFt; CONSTREF n, p : vector; dti : myReal;
                                  CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{calculates the charge of the (bulk and interface) traps at every gridpoint.}
BEGIN
	FILLCHAR(lin, SIZEOF(lin), 0); {reset linearisation}
	
	{first, the bulk}
	Calc_Linearization_f_tb(lin, n, p, dti, stv, par);
	
	{next, the interfaces:}
	Calc_Linearization_f_ti(lin, n, p, dti, stv, par)

END;

PROCEDURE Solve_Poisson(VAR V, n, p, nion, pion	: vector; VAR f_tb, f_ti, f_ti_numer, f_ti_inv_denom : TTrapArray; VAR Ntb_charge, Nti_charge : vector; CONSTREF Old_f_tb, Old_f_ti : TTrapArray;
						VAR conv, coupleIonsPoisson	: BOOLEAN; VAR PoissMsg : STRING; dti : myReal;
						CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Solves the Poisson equation, can be used in steady-state and transient simulations}
{Solve_Poisson also modifies the charges (n,p,ions,traps) by estimating the effects of the newly calc'd potential}
VAR it, i, j : INTEGER;
	sumPre, sumPost, NormDelV : myReal;
    delV, rhs, lower, upper, main : vector;
    fac_m, fac_u, fac_l, fac2, fac3, val, lnr : myReal;
    IonsOK : BOOLEAN;
	lin : TLinFt; {this type (a record) stores the linearisation of the trapping terms.}
BEGIN
    FOR i:=1 TO par.NP DO delV[i]:=1; {init delV}
    delV[0]:=0; {delV=0 at contacts, by definition, since we know V(0,L)}
    delV[par.NP+1]:=0;
    it:=0;
    conv:=FALSE;
	PoissMsg:='Poisson solver status:' + LineEnding; {message string}
	lnr:=LN(1 + par.couplePC);

	{if needed, check the total number of ions in volume}
	IF coupleIonsPoisson THEN 
    BEGIN
		sumPre:=0;
		FOR j:=1 TO stv.NLayers DO {loop over layers}
			IF (par.lyr[j].negIons OR par.lyr[j].posIons) AND coupleIonsPoisson THEN {if ions are moving in this layer, then}
				FOR i:=stv.i0[j] TO stv.i1[j] DO {sum over ion concentrations over all grid points in this layer}
					sumPre:=sumPre + nion[i] + pion[i]
	END;

    WHILE (NOT conv) AND (it < par.maxItPois) DO
    BEGIN

		Calc_Trap_Filling_Charge(f_tb, f_ti, f_ti_numer, f_ti_inv_denom, Ntb_charge, Nti_charge, n, p, Old_f_tb, Old_f_ti, dti, stv, par);
		Calc_Linearization_f_t_All(lin, n, p, dti, stv, par);

        FOR i:=1 TO par.NP DO {filling the matrix and the right-hand side}
			WITH stv DO
			BEGIN
			    {to properly deal with non-uniform dielectric constants we need an extra term in the Poisson equation, that's where these factors originate.}
				fac_m:= (eps[i]-eps[i-1])/(SQR(h[i-1])*SQR(Ltot)) - 2*eps[i]*(1 / (h[i]*(h[i]+h[i-1])*SQR(Ltot)) + 1 / (h[i-1]*(h[i]+h[i-1])*SQR(Ltot)));
				fac_l:= 2*eps[i] / (h[i-1]*(h[i]+h[i-1])*SQR(Ltot)) -(eps[i]-eps[i-1])/(SQR(h[i-1])*SQR(Ltot));
				fac_u:= 2*eps[i] / (h[i]*(h[i]+h[i-1])*SQR(Ltot));

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
				IF par.lyr[stv.lid[i]].negIons AND coupleIonsPoisson THEN nion[i]:=nion[i]*fac2; {and also apply this to the ions}
				IF par.lyr[stv.lid[i]].posIons AND coupleIonsPoisson THEN pion[i]:=pion[i]*fac3
			END
		END; {for loop}

        it:=it+1;
        NormDelV:=Norm_Eucl(delV, 1, par.NP);
        conv:=NormDelV <= par.tolPois  {finally, check for convergence}

    END;
	
	PoissMsg:=PoissMsg +'- delV ='+FloatToStrF(NormDelV, ffGeneral,5,0) + LineEnding;
  
    {OK, now see if we haven't changed the ions too much, so again we sum their concentrations}
    IF coupleIonsPoisson THEN 
    BEGIN 
		sumPost:=0;
		FOR j:=1 TO stv.NLayers DO {loop over layers}
			IF (par.lyr[j].negIons OR par.lyr[j].posIons) THEN {if ions are moving in this layer, then}
				FOR i:=stv.i0[j] TO stv.i1[j] DO {sum over ion concentrations over all grid points in this layer}
					sumPost:=sumPost + nion[i] + pion[i];

		{note: by now, sumPost cannot be zero as there are ions}
		IonsOK:=(ABS(sumPre-sumPost)/sumPost) < par.tolPois;
		conv:=conv AND IonsOK;
		PoissMsg:=PoissMsg + '- movement of ions acceptable: ' + myBoolStr(IonsOK) + LineEnding;
	END;

	PoissMsg:=PoissMsg + '- Poisson solver converged: ' + myBoolStr(conv)

END;

PROCEDURE Calc_Elec_Mob(VAR mu : vector; CONSTREF V, n : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
VAR i, j : INTEGER;
{calculates the elec. mob. on the interleaved mesh, mun[i]=mun at x=i+1/2}
{therefore the field or concentration on x=i+1/2 is needed}
BEGIN
	{we loop over the layers:}
	FOR j:=1 TO stv.NLayers DO
		WITH par.lyr[j] DO
			CASE mobnDep OF
				0 : FOR i:=stv.i0[j] TO stv.i1[j] DO mu[i]:=mu_n; { mob. is constant }
				1 : FOR i:=stv.i0[j] TO stv.i1[j] DO {field-dep. mob}
						mu[i]:=mu_n * EXP(gamma_n*SQRT(ABS(V[i+1]-V[i])/(stv.Ltot*stv.h[i])));
				ELSE Stop_Prog('Only very simple mobility models (0 and 1) currently implemented.', EC_ProgrammingError);
			END; {case selector}
		
	{now we need to re-assess the interfaces and use nu_int_n to compute the mobility:}
	FOR j:=1 TO stv.NLayers-1 DO {note: -1!}
		mu[stv.i1[j]]:=(stv.Ltot*stv.h[stv.i1[j]])*par.lyr[j].nu_int_n*stv.Vti; 		
		
END;

PROCEDURE Calc_Hole_Mob(VAR mu : vector; CONSTREF V, p : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
VAR i, j : INTEGER;
{calculates the hole mob. on the interleaved mesh, mun[i]=mun at x=i+1/2}
{therefore the field or concentration on x=i+1/2 is needed}
BEGIN
	{we loop over the layers:}
	FOR j:=1 TO stv.NLayers DO
		WITH par.lyr[j] DO
			CASE mobpDep OF
				0 : FOR i:=stv.i0[j] TO stv.i1[j] DO mu[i]:=mu_p; { mob. is constant }
				1 : FOR i:=stv.i0[j] TO stv.i1[j] DO {field-dep. mob}
						mu[i]:=mu_p * EXP(gamma_p*SQRT(ABS(V[i+1]-V[i])/(stv.Ltot*stv.h[i])));
				ELSE Stop_Prog('Only very simple mobility models (0 and 1) currently implemented.', EC_ProgrammingError);
			END; {case selector}
		
	{now we need to re-assess the interfaces and use nu_int_p to compute the mobility:}
	FOR j:=1 TO stv.NLayers-1 DO {note: -1!}
		mu[stv.i1[j]]:=(stv.Ltot*stv.h[stv.i1[j]])*par.lyr[j].nu_int_p*stv.Vti; 		
		
END;

PROCEDURE Calc_Langevin_Factor(VAR Lan : vector; mob_n, mob_p : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculates the Langevin recombination strength. Lan[i] is defined on the
regular grid, i.e., Lan[i]=Lan at x=xi. Note that mobilities are defined on the
interleaved mesh!}
VAR i, j : INTEGER;
    rec_mob : myReal;
BEGIN
	{we loop over the layers:}
	FOR j:=1 TO stv.NLayers DO
		WITH par.lyr[j] DO 
		BEGIN
			IF useLangevin THEN {use Langevin formula to calculate bimolecular, direct, band-to-band recombination rate}
			FOR i:=MAX(1, stv.i0[j]) TO stv.i1[j] DO {why the MAX? we cannot have i=1, that's why}
			BEGIN
				rec_mob:=0.5 * (mob_n[i-1] + mob_n[i]+ mob_p[i-1] + mob_p[i] );
				{we take mob(x=xi)=(mob(x=xi-1/2)+mob(x=xi+1/2))/2}
				Lan[i]:=preLangevin * q * rec_mob/stv.eps[i]
			END
			ELSE
				FOR i:=stv.i0[j] TO stv.i1[j] DO
					Lan[i]:=k_direct			
		END;

    Lan[0]:=0; {no recombination (or generation) at the contacts}
    Lan[par.NP+1]:=0;
END;


FUNCTION Diss_Prob_Delta(r : myReal; vals : Row) : myReal;
{Calculates the dissociation probability as a function of distance r}
VAR b, kdF, delE, epsi, Vti, Braun_rec, F, k_f : myReal;
BEGIN
    Vti:=vals[1]; {local inverse thermal voltage}
    epsi:=vals[2]; {local relative dielectric constant}
    Braun_rec:=vals[3]; {local copy of Braun recombination strength}
    F:=vals[4]; {local copy of electric field}
    k_f:=vals[5]; {local copy of k_f}
    delE:=q/(4*PI*epsi*r); {binding energy in eV}
    b:=q*ABS(F)*SQR(Vti)/(8*PI*epsi); {scaled ABOLUTE field strength}
    kdF:=3*Braun_rec/(4*PI*r*SQR(r))*EXP(-delE*Vti)*Bessel(b);
    Diss_Prob_Delta:=kdF/(kdF + k_f)
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
VAR i, j : INTEGER;
    vals : Row; {array with values needed for Romberg integration}
BEGIN
	{we loop over the layers:}
	FOR j:=1 TO stv.NLayers DO
		WITH par.lyr[j] DO
			IF fieldDepG THEN { calc. field dep. G }
			BEGIN
				SetLength(vals, 6); {vals: needed to pass values to MathFuncValues}
				vals[0]:=a; {we'll use it (here) to pass the charge separation distance}
				vals[1]:=stv.Vti; {and the inverse thermal voltage} 
				vals[5]:=k_f; 
				{we exclude the contacts (0, par.NP+1) in the loop:}
				FOR i:=MAX(1,stv.i0[j]) TO MIN(par.NP,stv.i1[j]) DO
				BEGIN
					vals[4]:=(V[i+1]-V[i-1])/(stv.x[i+1]-stv.x[i-1]); {local electric field}
					{F= field on mesh point i, centered difference, global variable}      

					vals[2]:=stv.eps[i]; {copy local dielectric constant into vals as Diss_Prob_Delta needs it}
					vals[3]:=Lan[i]; {Braun recombination strength at x=xi, equal to Langevin (direct)}

					CASE thermLengDist OF  {a = thermalization length}
						1 : dp[i]:=Diss_Prob_Delta(a, vals); {delta-function distribution}
						2 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_Gauss, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
						3 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_Exp, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
						4 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_SQRrExp, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
						5 : WITH par DO dp[i]:=RombergIntegrationValues(@Diss_Prob_r4Gauss, vals, LowerLimBraun*a, UpperLimBraun*a, TolRomb, MaxRombIt, FALSE);
					END; {case selector}
					{total free-carrier yield is sum of direct generation (P0) and the field dependent part (1-P0)*(dp)}
					g[i]:=(P0 + (1-P0)*dp[i]) * Gm[i]
				END {for loop}
			END { calc. field dep. G }
			ELSE FOR i:=stv.i0[j] TO stv.i1[j] DO BEGIN dp[i]:=0; g[i]:=Gm[i] END; {G is constant}
      
    dp[0]:=0;
    dp[par.NP+1]:=0;
    g[0]:=0;   {no generation at the contacts}
    g[par.NP+1]:=0
END;


PROCEDURE Calc_Recombination_n(VAR Rn : TRec; dti : myReal; CONSTREF n, p, dp, Lan : vector; f_tb, f_ti : TTrapArray; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculate all recombination processes and their contribution to the continuity equation for electrons.}
{See Koopmans and Koster, Sol. RRL 2022, 2200560.} 
VAR i, ii, j, e : INTEGER;
	iiw,
	sum_aj, sum_bj, sum_cj, sum_dj,
	ciL, ciR, ai,
	f_ti_level_inv_denom, numer, denom, a, b, c1, g, h, g0, g1, h0, h1, 
	dum1, dum2 : myReal; {used to store partial results}
BEGIN 
	FILLCHAR(Rn, SIZEOF(Rn), 0); {first: reset all fields to 0}

	{direct recombination:}
	FOR i:=1 TO par.NP DO BEGIN {note: no net rec/gen on the electrodes}
		Rn.dir_cont_rhs[i]:= (1-dp[i])*Lan[i]*SQR(stv.ni[i]);
		Rn.dir_cont_m[i]:= (1-dp[i])*Lan[i]*p[i];
		Rn.direct[i]:= Rn.dir_cont_m[i]*n[i] - Rn.dir_cont_rhs[i]
	END;

	{recombination via bulk traps}
	FOR j:=1 TO stv.NLayers DO
		FOR e:=1 TO stv.Ntb[j].NLevels DO {loop over trap levels}
			FOR i:=MAX(1,stv.i0[j]) TO MIN(par.NP,stv.i1[j]) DO {loop over grid points in this layer, but as always, exclude the electrodes}
				IF dti=0 THEN {steady-state}
				BEGIN
					numer:=par.lyr[j].C_n_bulk*par.lyr[j].C_p_bulk*(n[i]*p[i] - stv.Ntb[j].nt0[e]*stv.Ntb[j].pt0[e]);
					denom:=(par.lyr[j].C_n_bulk*(n[i]+stv.Ntb[j].nt0[e]) + par.lyr[j].C_p_bulk*(p[i]+stv.Ntb[j].pt0[e]));
					dum1:=stv.Ntb[j].Nt[e] * (numer / denom);
					Rn.bulk[i]:=Rn.bulk[i] + dum1;
					dum2:=stv.Ntb[j].Nt[e] * (denom*par.lyr[j].C_n_bulk*par.lyr[j].C_p_bulk*p[i] - numer*par.lyr[j].C_n_bulk) / SQR(denom);
					Rn.bulk_cont_m[i]:=Rn.bulk_cont_m[i] + dum2;
					Rn.bulk_cont_rhs[i]:=Rn.bulk_cont_rhs[j] + dum1 - n[i]*dum2
				END
				ELSE BEGIN {transient}
					{We solve the integral of the ODE we solved for calculating the trap occupance.}
					b:=par.lyr[j].C_n_bulk*n[i] + par.lyr[j].C_p_bulk*stv.Ntb[j].pt0[e];
					a:=par.lyr[j].C_n_bulk*n[i] + par.lyr[j].C_n_bulk*stv.Ntb[j].nt0[e] + par.lyr[j].C_p_bulk*stv.Ntb[j].pt0[e] + par.lyr[j].C_p_bulk*p[i];
					g:=par.lyr[j].C_n_bulk*n[i]*stv.Ntb[j].Nt[e] + par.lyr[j].C_n_bulk*stv.Ntb[j].nt0[e]*stv.Ntb[j].Nt[e];
					h:=par.lyr[j].C_n_bulk*n[i]*stv.Ntb[j].Nt[e];
					c1:=f_tb[i,e] - b/a;
		
					Rn.bulk_cont_m[i]:=0;
					Rn.bulk_cont_rhs[i]:=Rn.bulk_cont_rhs[i] + ((h - g*b/a) / dti + g*c1* EXP(-a/dti)/a - g*c1/a)*dti;
					Rn.bulk[i]:=Rn.bulk_cont_rhs[i]			
				END;

	{recombination via interface traps:}
	
	FOR j:=1 TO stv.NLayers-1 DO {loop over interfaces. Each interface involves 2 points: i1[j] and i1[j]+1 (in the adjacent layer)}
	BEGIN
		{note: the interface sits between i1[j] and i1[j]+1}
		ii:=stv.i1[j]; {ii: i interface}	
		iiw:=2 / (stv.Ltot*(0.5*stv.h[ii+1] + stv.h[ii] + 0.5*stv.h[ii-1])); {inverse of width of the interface}

		IF dti=0 THEN {steady-state}
		BEGIN
			FOR e:=1 TO stv.Nti[j].NLevels DO {loop over trap levels}
			BEGIN 
				f_ti_level_inv_denom:=Calc_Inv_f_ti_Denom(n, p, ii, j, e, stv, par);
				sum_aj:=par.lyr[j].C_n_int * stv.Nti[j].Nt[e] * (n[ii] + n[ii+1]);
				sum_bj:=par.lyr[j].C_p_int * stv.Nti[j].Nt[e] * (p[ii] + p[ii+1]);
				sum_cj:=par.lyr[j].C_n_int * stv.Nti[j].Nt[e] * (stv.Nti[j].nt0_R[e] + stv.Nti[j+1].nt0_L[e]); 
				sum_dj:=par.lyr[j].C_p_int * stv.Nti[j].Nt[e] * (stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]); 

				ciR:=par.lyr[j].C_n_int * stv.Nti[j].Nt[e] * stv.Nti[j].nt0_R[e]; 
				ciL:=par.lyr[j].C_n_int * stv.Nti[j].Nt[e] * stv.Nti[j+1].nt0_L[e]; 
		
				ai:=par.lyr[j].C_n_int * stv.Nti[j].Nt[e]; {again, different from v4.57 as there's only 1 ai}
				
				{Interface recombination as calculated with the current n and p}
				Rn.int[ii]:=Rn.int[ii] + iiw*(ai*n[ii] * (sum_bj + sum_cj) - ciR * (sum_aj + sum_dj)) * f_ti_level_inv_denom;
				Rn.int[ii+1]:=Rn.int[ii+1] + iiw*(ai*n[ii+1] * (sum_bj + sum_cj) - ciL* (sum_aj + sum_dj)) * f_ti_level_inv_denom;

				{Calculate the partial derivative of recombination to n[ii-1], n[ii], n[ii+1]}
				dum1:=ciR * ai / f_ti_level_inv_denom;
				dum1:=dum1 - ai * (ai*n[ii] * (sum_bj + sum_cj) - ciR * (sum_aj + sum_dj));
				dum2:=ciL * ai / f_ti_level_inv_denom;
				dum2:=dum2 - ai * (ai*n[ii+1] * (sum_bj + sum_cj) - ciL * (sum_aj + sum_dj));
				Rn.int_cont_lo[ii]:=Rn.int_cont_lo[ii] + iiw*dum1 * SQR(f_ti_level_inv_denom);
				Rn.int_cont_lo[ii+1]:=Rn.int_cont_lo[ii+1] + iiw*dum2 * SQR(f_ti_level_inv_denom);
		
				dum1:=-ciR* ai / f_ti_level_inv_denom;
				dum1:=dum1 - ai * (ai*n[ii] * (sum_bj + sum_cj) - ciR * (sum_aj + sum_dj));
				dum2:=-ciL * ai / f_ti_level_inv_denom;
				dum2:=dum2 - ai * (ai*n[ii+1] * (sum_bj + sum_cj) - ciL * (sum_aj + sum_dj));
				Rn.int_cont_up[ii]:= Rn.int_cont_up[ii] + iiw*dum1 * SQR(f_ti_level_inv_denom);
				Rn.int_cont_up[ii+1]:=Rn.int_cont_up[ii+1] + iiw*dum2 * SQR(f_ti_level_inv_denom);
				
				dum1:=(-ciR + sum_bj + sum_cj) * ai  / f_ti_level_inv_denom;
				dum1:=dum1 - ai * (ai*n[ii] * (sum_bj + sum_cj) - ciR * (sum_aj + sum_dj));
				dum2:=(-ciL + sum_bj + sum_cj) * ai  / f_ti_level_inv_denom;
				dum2:=dum2 - ai * (ai*n[ii+1] * (sum_bj + sum_cj) - ciL * (sum_aj + sum_dj));
				Rn.int_cont_m[ii]:=Rn.int_cont_m[ii] + iiw*dum1 * SQR(f_ti_level_inv_denom);
				Rn.int_cont_m[ii+1]:=Rn.int_cont_m[ii+1] + iiw*dum2 * SQR(f_ti_level_inv_denom);
			
			END; {steady-state, loop over trap levels}
			
			{The right hand side of the continuity equation contains the recombination term, but because we linearize in n we add the linearization terms as well.}
			Rn.int_cont_rhs[ii]:=Rn.int[ii] - n[ii-1] * Rn.int_cont_lo[ii] - n[ii+1] * Rn.int_cont_up[ii] - n[ii] * Rn.int_cont_m[ii];
			Rn.int_cont_rhs[ii+1]:=Rn.int[ii+1] - n[ii] * Rn.int_cont_lo[ii+1] - n[ii+2] * Rn.int_cont_up[ii+1] - n[ii+1] * Rn.int_cont_m[ii+1]	
		END {steady-state}
		ELSE 
		BEGIN {transient}
			FOR e:=1 TO stv.Nti[j].NLevels DO {loop over trap levels}
			BEGIN 
				{We solve the integral of the ODE we solved for calculating the trap occupance.}
				a:=par.lyr[j].C_n_int * (n[ii] + n[ii+1] + stv.Nti[j].nt0_R[e] + stv.Nti[j+1].nt0_L[e]) +
				   par.lyr[j].C_p_int * (p[ii] + p[ii+1] + stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]);
				b:=par.lyr[j].C_n_int*(n[ii] + n[ii+1]) + par.lyr[j].C_p_int*(stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]);
				g0:=par.lyr[j].C_n_int*stv.Nti[j].Nt[e]*(n[ii] + stv.Nti[j].nt0_R[e]);			  
				g1:=par.lyr[j].C_n_int*stv.Nti[j].Nt[e]*(n[ii+1] + stv.Nti[j+1].nt0_L[e]);			  

				h0:=par.lyr[j].C_n_int*n[ii]*stv.Nti[j].Nt[e];
				h1:=par.lyr[j].C_n_int*n[ii+1]*stv.Nti[j].Nt[e];
				
				c1:= (f_ti[ii,e] - b/a);
			
				Rn.int_cont_rhs[ii]:=Rn.int_cont_rhs[ii] + iiw*((h0 - g0*b/a) / dti + g0*c1* EXP(-a/dti)/a - g0*c1/a)*dti;
				Rn.int_cont_rhs[ii+1]:=Rn.int_cont_rhs[ii+1] + iiw*((h1 - g1*b/a) / dti + g1*c1* EXP(-a/dti)/a - g1*c1/a)*dti;

				Rn.int[ii]:=Rn.int_cont_rhs[ii];
				Rn.int[ii+1]:=Rn.int_cont_rhs[ii+1]	
			END;

		END
		
	END; {loop over interfaces}

END;

PROCEDURE Calc_Recombination_p(VAR Rp : TRec; dti : myReal; CONSTREF n, p, dp, Lan : vector; f_tb, f_ti : TTrapArray; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Calculate all recombination processes and their contribution to the continuity equation for holes.}
{See Koopmans and Koster, Sol. RRL 2022, 2200560.}
VAR i, ii, j, e : INTEGER;
	iiw,
	sum_aj, sum_bj, sum_cj, sum_dj,
	diL, diR, bi, 
	f_ti_level_inv_denom, numer, denom, a, b, c1, g, h, g0, g1, h0, h1, 
	dum1, dum2 : myReal; {used to store partial results}
BEGIN 
	FILLCHAR(Rp, SIZEOF(Rp), 0); {first: reset all fields to 0}
	
	{direct recombination:}
	FOR i:=1 TO par.NP DO BEGIN {note: no net rec/gen on the electrodes}
		Rp.dir_cont_rhs[i]:= (1-dp[i])*Lan[i]*SQR(stv.ni[i]);
		Rp.dir_cont_m[i]:= (1-dp[i])*Lan[i]*n[i];
		Rp.direct[i]:= Rp.dir_cont_m[i]*p[i] - Rp.dir_cont_rhs[i]
	END;
	
	{recombination via bulk traps}
	FOR j:=1 TO stv.NLayers DO
		FOR e:=1 TO stv.Ntb[j].NLevels DO {loop over trap levels}
			FOR i:=MAX(1,stv.i0[j]) TO MIN(par.NP,stv.i1[j]) DO {loop over grid points in this layer, but as always, exclude the electrodes}
				IF dti=0 THEN {steady-state}
				BEGIN
					numer:=par.lyr[j].C_n_bulk*par.lyr[j].C_p_bulk*(n[i]*p[i] - stv.Ntb[j].nt0[e]*stv.Ntb[j].pt0[e]);
					denom:=(par.lyr[j].C_n_bulk*(n[i]+stv.Ntb[j].nt0[e]) + par.lyr[j].C_p_bulk*(p[i]+stv.Ntb[j].pt0[e]));
					dum1:=stv.Ntb[j].Nt[e] * (numer / denom);
					Rp.bulk[i]:=Rp.bulk[i] + dum1;
					dum2:=stv.Ntb[j].Nt[e] * (denom*par.lyr[j].C_n_bulk*par.lyr[j].C_p_bulk*n[i] - numer*par.lyr[j].C_p_bulk) / SQR(denom);
					Rp.bulk_cont_m[i]:=Rp.bulk_cont_m[i] + dum2;
					Rp.bulk_cont_rhs[i]:=Rp.bulk_cont_rhs[j] + dum1 - p[i]*dum2
				END
				ELSE BEGIN {transient}
					{We solve the integral of the ODE we solved for calculating the trap occupance.}
					b:=par.lyr[j].C_p_bulk*p[i] + par.lyr[j].C_n_bulk*stv.Ntb[j].nt0[e];
					a:=par.lyr[j].C_n_bulk*n[i] + par.lyr[j].C_n_bulk*stv.Ntb[j].nt0[e] + par.lyr[j].C_p_bulk*stv.Ntb[j].pt0[e] + par.lyr[j].C_p_bulk*p[i];
					g:=par.lyr[j].C_p_bulk*p[i]*stv.Ntb[j].Nt[e] + par.lyr[j].C_p_bulk*stv.Ntb[j].pt0[e]*stv.Ntb[j].Nt[e];
					h:=par.lyr[j].C_p_bulk*p[i]*stv.Ntb[j].Nt[e];
					c1:=1-f_tb[i,e] - b/a;
				
					Rp.bulk_cont_m[i]:=0;
					Rp.bulk_cont_rhs[i]:=Rp.bulk_cont_rhs[i] + ((h - g*b/a) / dti + g*c1* EXP(-a/dti)/a - g*c1/a)*dti;
					Rp.bulk[i]:=Rp.bulk_cont_rhs[i]			
				END;

	{recombination via interface traps:}
	
	FOR j:=1 TO stv.NLayers-1 DO {loop over interfaces. Each interface involves 2 points: i1[j] and i1[j]+1 (in the adjacent layer)}
	BEGIN
		{note: the interface sits between i1[j] and i1[j]+1}
		ii:=stv.i1[j]; {ii: i interface}	
		iiw:=2 / (stv.Ltot*(0.5*stv.h[ii+1] + stv.h[ii] + 0.5*stv.h[ii-1])); {inverse of width of the interface}
		
		IF dti=0 THEN {steady-state}
		BEGIN
			FOR e:=1 TO stv.Nti[j].NLevels DO {loop over trap levels}
			BEGIN 
				f_ti_level_inv_denom:=Calc_Inv_f_ti_Denom(n, p, ii, j, e, stv, par);
				sum_aj:=par.lyr[j].C_n_int * stv.Nti[j].Nt[e] * (n[ii] + n[ii+1]);
				sum_bj:=par.lyr[j].C_p_int * stv.Nti[j].Nt[e] * (p[ii] + p[ii+1]);
				sum_cj:=par.lyr[j].C_n_int * stv.Nti[j].Nt[e] * (stv.Nti[j].nt0_R[e] + stv.Nti[j+1].nt0_L[e]); 
				sum_dj:=par.lyr[j].C_p_int * stv.Nti[j].Nt[e] * (stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]); 
			
				diR:=par.lyr[j].C_p_int * stv.Nti[j].Nt[e] * stv.Nti[j].pt0_R[e]; 
				diL:=par.lyr[j].C_p_int * stv.Nti[j].Nt[e] * stv.Nti[j+1].pt0_L[e]; 
				
				bi:=par.lyr[j].C_p_int * stv.Nti[j].Nt[e]; {again, different from v4.57 as there's only 1 bi}
				
				{Interface recombination as calculated with the current n and p}
				Rp.int[ii]:=Rp.int[ii] + iiw*(bi*p[ii] * (sum_aj + sum_dj) - diR * (sum_bj + sum_cj)) / (sum_aj+sum_bj+sum_cj+sum_dj);
				Rp.int[ii+1]:=Rp.int[ii+1] + iiw*(bi*p[ii+1] * (sum_aj + sum_dj) - diL * (sum_bj + sum_cj)) / (sum_aj+sum_bj+sum_cj+sum_dj);
	
				{Calculate the partial derivative of recombination to p[ii-1], p[ii], p[ii+1]}
				dum1:=-diR * bi / f_ti_level_inv_denom;
				dum1:=dum1 - bi * (bi*p[ii] * (sum_aj + sum_dj) - diR * (sum_bj + sum_cj));
				dum2:=-diL * bi / f_ti_level_inv_denom;
				dum2:=dum2 - bi * (bi*p[ii+1] * (sum_aj + sum_dj) - diL * (sum_bj + sum_dj));
				Rp.int_cont_lo[ii]:=Rp.int_cont_lo[ii] + iiw*dum1 * SQR(f_ti_level_inv_denom);
				Rp.int_cont_lo[ii+1]:=Rp.int_cont_lo[ii+1] + iiw*dum2 * SQR(f_ti_level_inv_denom);

				dum1:=-diR * bi / f_ti_level_inv_denom;
				dum1:=dum1 - bi * (bi*n[ii] * (sum_aj + sum_dj) - diR * (sum_bj + sum_cj));
				dum2:=-diL * bi / f_ti_level_inv_denom;
				dum2:=dum2 - bi * (bi*p[ii+1] * (sum_aj + sum_dj) - diL * (sum_bj + sum_cj));
				Rp.int_cont_up[ii]:= Rp.int_cont_up[ii] + iiw*dum1 * SQR(f_ti_level_inv_denom);
				Rp.int_cont_up[ii+1]:=Rp.int_cont_up[ii+1] + iiw*dum2 * SQR(f_ti_level_inv_denom);
				
				dum1:=(-diR + sum_aj + sum_dj) * bi  / f_ti_level_inv_denom;
				dum1:=dum1 - bi * (bi*p[ii] * (sum_aj + sum_dj) - diR * (sum_bj + sum_cj));
				dum2:=(-diL + sum_aj + sum_dj) * bi  / f_ti_level_inv_denom;
				dum2:=dum2 - bi * (bi*p[ii+1] * (sum_aj + sum_dj) - diL * (sum_bj + sum_cj));
				Rp.int_cont_m[ii]:=Rp.int_cont_m[ii] + iiw*dum1 * SQR(f_ti_level_inv_denom);
				Rp.int_cont_m[ii+1]:=Rp.int_cont_m[ii+1] + iiw*dum2 * SQR(f_ti_level_inv_denom)		
			END; {steady-state, loop over trap levels}
			
			{The right hand side of the continuity equation contains the recombination term, but because we linearize in n we add the linearization terms as well.}
			Rp.int_cont_rhs[ii]:=Rp.int[ii] - p[ii-1] * Rp.int_cont_lo[ii] - p[ii+1] * Rp.int_cont_up[ii] - p[ii] * Rp.int_cont_m[ii];
			Rp.int_cont_rhs[ii+1]:=Rp.int[ii+1] - p[ii] * Rp.int_cont_lo[ii+1] - p[ii+2] * Rp.int_cont_up[ii+1] - p[ii+1] * Rp.int_cont_m[ii+1]	
		END {steady-state}
		ELSE 
		BEGIN {transient}
			FOR e:=1 TO stv.Nti[j].NLevels DO {loop over trap levels}
			BEGIN 
				{We solve the integral of the ODE we solved for calculating the trap occupance.}
				a:=par.lyr[j].C_n_int * (n[ii] + n[ii+1] + stv.Nti[j].nt0_R[e] + stv.Nti[j+1].nt0_L[e]) +
				   par.lyr[j].C_p_int * (p[ii] + p[ii+1] + stv.Nti[j].pt0_R[e] + stv.Nti[j+1].pt0_L[e]);
				b:=par.lyr[j].C_p_int*(p[ii] + p[ii+1]) + par.lyr[j].C_n_int*(stv.Nti[j].nt0_R[e] + stv.Nti[j+1].nt0_L[e]);

				g0:=par.lyr[j].C_p_int*stv.Nti[j].Nt[e]*(p[ii] + stv.Nti[j].pt0_R[e]);			  
				g1:=par.lyr[j].C_p_int*stv.Nti[j].Nt[e]*(p[ii+1] + stv.Nti[j+1].pt0_L[e]);			  

				h0:=par.lyr[j].C_p_int*p[ii]*stv.Nti[j].Nt[e];
				h1:=par.lyr[j].C_p_int*p[ii+1]*stv.Nti[j].Nt[e];
				
				c1:= (1-f_ti[ii,e] - b/a);
					
				Rp.int_cont_rhs[ii]:=Rp.int_cont_rhs[ii] + iiw*((h0 - g0*b/a) / dti + g0*c1* EXP(-a/dti)/a - g0*c1/a)*dti;
				Rp.int_cont_rhs[ii+1]:=Rp.int_cont_rhs[ii+1] + iiw*((h1 - g1*b/a) / dti + g1*c1* EXP(-a/dti)/a - g1*c1/a)*dti;

				Rp.int[ii]:=Rp.int_cont_rhs[ii];
				Rp.int[ii+1]:=Rp.int_cont_rhs[ii+1]	
			END;

		END
		
	END; {loop over interfaces}

END;

PROCEDURE Cont_Eq_Elec(VAR n : vector; nPrevTime, V, Jn, p, mu, g, Lan, dp : vector; VAR f_tb, f_ti : TTrapArray; VAR Rn : TRec;
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
	IF (par.S_n_L < 0) OR (Jn[0]>0) THEN
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
		m[0]:=-mu[0]*stv.Vt*B((V[0]-V[1])*stv.Vti)/(stv.Ltot*stv.h[0]) - par.S_n_L;
		u[0]:=mu[0]*stv.Vt*B((V[1]-V[0])*stv.Vti)/(stv.Ltot*stv.h[0]);
		rhs[0]:= - par.S_n_L*stv.NcLoc[0] * EXP((stv.E_CB[0] - par.W_L)*stv.Vti);
	END;

	{Calculate recombination and its contribution to the continuity equation} 
	Calc_Recombination_n(Rn, dti, n, p, dp, Lan, f_tb, f_ti, stv, par);

    {now set the interior part:}
    FOR i:=1 TO NP DO  {continuity eq. in matrix vorm}
		WITH stv DO {ni and h are static variables!}
		BEGIN
			fac := 0.5*SQR(Ltot)*h[i]*h[i-1]*(h[i]+h[i-1]); {repeats often in the equations}
			
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
	IF (par.S_n_R < 0) OR (Jn[NP]>0) THEN
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
		lo[NP+1]:=-mu[NP]*stv.Vt * B((V[NP]-V[NP+1])*stv.Vti)/(stv.Ltot*stv.h[NP]);
		m[NP+1]:=mu[NP]*stv.Vt * B((V[NP+1]-V[NP])*stv.Vti)/(stv.Ltot*stv.h[NP]) - par.S_n_R;
		u[NP+1]:=0;
		rhs[NP+1]:=-par.S_n_R*stv.NcLoc[NP+1] * EXP((stv.E_CB[NP+1] - par.W_R)*stv.Vti);
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
 
	IF ResetDens AND (NOT par.ignoreNegDens) THEN Stop_Prog('Negative electron concentration encountered!' , EC_NumericalFailure)
END;


PROCEDURE Cont_Eq_Holes(VAR p : vector; pPrevTime, V, Jp, n, mu, g, Lan, dp : vector; VAR f_tb, f_ti : TTrapArray; VAR Rp : TRec;
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
	IF (par.S_p_L < 0) OR (Jp[0]>0) THEN
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
		m[0]:=mu[0]*stv.Vt * B((V[1]-V[0])*stv.Vti)/(stv.Ltot*stv.h[0]) - par.S_p_L;
		u[0]:=-mu[0]*stv.Vt * B((V[0]-V[1])*stv.Vti)/(stv.Ltot*stv.h[0]);
		rhs[0]:=-par.S_p_L*stv.NcLoc[0] * EXP((par.W_L - stv.E_VB[0])*stv.Vti);
	END;

	{Calculate recombination and its contribution to the continuity equation} 	
	Calc_Recombination_p(Rp, dti, n, p, dp, Lan, f_tb, f_ti, stv, par);

    {now do the interior points:}
    FOR i:=1 TO NP DO  {continuity eq. in matrix vorm}
		WITH stv DO {ni and h are static variables!}
		BEGIN
			fac := 0.5*SQR(Ltot)*h[i]*h[i-1]*(h[i]+h[i-1]); {repeats often in the equations}
			
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
	IF (par.S_p_R < 0) OR (Jp[NP]>0) THEN
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
		lo[NP+1]:=mu[NP]*stv.Vt*B((V[NP+1]-V[NP])*stv.Vti)/(stv.Ltot*stv.h[NP]);
		m[NP+1]:=-mu[NP]*stv.Vt*B((V[NP]-V[NP+1])*stv.Vti)/(stv.Ltot*stv.h[NP]) - par.S_p_R;
		u[NP+1]:=0;
		rhs[NP+1]:=-par.S_p_R*stv.NcLoc[NP+1] * EXP((par.W_R - stv.E_VB[NP+1])*stv.Vti);
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

	IF ResetDens AND (NOT par.ignoreNegDens) THEN Stop_Prog('Negative hole concentration encountered!', EC_NumericalFailure)

END;

PROCEDURE Calc_Displacement_Curr(VAR JD : vector; V, VPrevTime : vector; dti : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure calculates the displacement current}
VAR i : INTEGER;
BEGIN
	FOR i:=0 TO par.NP DO
		JD[i]:=stv.eps[i] * (V[i+1]-V[i]-VPrevTime[i+1]+VPrevTime[i]) * dti / (stv.Ltot*stv.h[i]);
	JD[par.NP+1]:=JD[par.NP]; {doesn't have a physical meaning though}
END;

PROCEDURE Calc_Curr_Diff(sn : ShortInt; istart, ifinish : INTEGER; VAR J : vector; V, dens, mu, Rint : vector; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure calculates the current density in differential form, see Selberherr eq. 6.1-39 or 6.1-41}
{sn denotes the sign of the carrier, so -1 for electrons, +1 for holes}
VAR i : INTEGER;
BEGIN
	IF ABS(sn)<>1 THEN Stop_Prog('Incorrect sn passed to Calc_Curr_Diff', EC_ProgrammingError);
	IF (istart<0) OR (ifinish>=par.NP+1) THEN Stop_Prog('Incorrect istart and/or ifinish passed to Calc_Curr_Diff.', EC_ProgrammingError);

	{the current is only non-zero between istart and ifinish, so first init the zero's:}
	FOR i:=0 TO istart-1 DO J[i]:=0;
	FOR i:=ifinish+1 TO par.NP DO J[i]:=0;
	
	{now actually calc the current:}
	FOR i:=istart TO ifinish DO
		J[i]:=sn*q*mu[i]*stv.Vt*(dens[i+1]*B(sn*(V[i]-V[i+1])*stv.Vti) - dens[i]*B(sn*(V[i+1]-V[i])*stv.Vti))/(stv.Ltot*stv.h[i]);

	{now correct for interface recombination. This represents the current THROUGH the interface traps.}
	{do NOT be tempted to combine this in the previous for loop!!}
	FOR i:=istart+1 TO ifinish DO
		IF (i = stv.i1[stv.lid[i]]) THEN {we're crossing an interface}
			J[i]:=J[i-1] + sn*0.5*q*stv.Ltot*stv.h[i]*(Rint[i] + Rint[i+1]);
	
	{last point as this is in the output (varFile)}
	J[par.NP+1]:=J[par.NP]; {doesn't have a physical meaning: J[NP+1] is current between NP+1 and NP+2}
END;

PROCEDURE Calc_Curr_Int(sn : ShortInt; istart, ifinish : INTEGER; dti : myReal; VAR J : vector; V, dens, olddens, mu, g : vector; 
						CONSTREF Rec : TRec; CONSTREF stv : TStaticVars; 
						CONSTREF par : TInputParameters);
{Calculates the current density in integral form, see De Mari, solid-state elec. vol 11 p.33 (68) eq. 15}
{istart and ifinish exclude the electrodes, so istart>=1, ifinish<=NP}
{sn denotes the sign of the carrier, so -1 for electrons, +1 for holes}
VAR i : INTEGER;
    K, single_int, double_int, int_U, int_U_old, dx : myReal;
    U : vector;
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
        J[i]:=J[i-1] + sn*q*stv.Ltot*stv.h[i]*U[i];

	{now correct for interface recombination. This represents the current THROUGH the interface traps.}
	{do NOT be tempted to combine this in the previous for loop!!}
	FOR i:=istart+1 TO ifinish DO
		IF (i = stv.i1[stv.lid[i]]) THEN {we're crossing an interface}
			J[i]:=J[i-1] + sn*0.5*q*stv.Ltot*stv.h[i]*(Rec.int[i] + Rec.int[i+1]);

	{last point as this is in the output (varFile)}	
    J[par.NP+1]:=J[par.NP]; {doesn't have a physical meaning: J[NP+1] is current between NP+1 and NP+2}
END;

PROCEDURE Calc_All_Currents(VAR new : TState; CONSTREF curr : TState; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters); 
{calculates all the currents for state new}
VAR Jsum, GenDum : vector;
	i : INTEGER;
BEGIN
	WITH new DO 
	BEGIN
		{Note: we use global variable RecDum : TRec as a dummy variable. It is too larger to be a local variable.}
	
		{first do the electrons}
		CASE par.currDiffInt OF
			1 : Calc_Curr_Diff(-1, 0, par.NP, Jn, Vgn, n, mun, Rn.int, stv, par); {only needs part of Rn with interface recombination}
			2 : Calc_Curr_Int(-1, 0, par.NP, dti, Jn, Vgn, n, curr.n, mun, gen, Rn, stv, par); {needs full Rn and curr.n}
		END;	

		{now do the holes}
		CASE par.currDiffInt OF
			1 : Calc_Curr_Diff(1, 0, par.NP, Jp, Vgp, p, mup, Rp.int, stv, par);{only needs part of Rp with interface recombination}
			2 : Calc_Curr_Int(1, 0, par.NP, dti, Jp, Vgp, p, curr.p, mup, gen, Rp, stv, par); {needs full Rp and curr.p}
		END;	

		{first, set all ionic currents to zero:}
		FILLCHAR(Jnion, SIZEOF(Jnion), 0); 
		FILLCHAR(Jpion, SIZEOF(Jpion), 0);
	
		IF dti<>0 THEN {only in transient simulations can the ionic currents be non-zero}
		BEGIN
			FILLCHAR(RecDum, SIZEOF(RecDum), 0); {set all fields of RecDum to zero as ions don't have generation/recombination}
			FILLCHAR(GenDum, SIZEOF(GenDum), 0); {likewise for the generation rate of ions}
			
			{Next, we loop over the ionic regions to see what's there}
			{first, the negative ions}
			FOR i:=0 TO LENGTH(stv.NegIonRegion)-1 DO
				WITH stv.NegIonRegion[i] DO
					IF AvC > 0 THEN
						CASE par.currDiffInt OF
							1 : Calc_Curr_Diff(-1, istart, ifinish-1, Jnion, V, nion, stv.mu_n_ion, RecDum.int, stv, par);		
							2 : Calc_Curr_Int(-1, istart, ifinish-1, dti, Jnion, V, nion, curr.nion, stv.mu_n_ion, GenDum, RecDum, stv, par);
						END; {case}	
				
			{next: positive ions}
			FOR i:=0 TO LENGTH(stv.PosIonRegion)-1 DO
				WITH stv.PosIonRegion[i] DO
					IF AvC > 0 THEN
						CASE par.currDiffInt OF
							1 : Calc_Curr_Diff(1, istart, ifinish-1, Jpion, V, pion, stv.mu_p_ion, RecDum.int, stv, par);		
							2 : Calc_Curr_Int(1, istart, ifinish-1, dti, Jpion, V, pion, curr.pion, stv.mu_p_ion, GenDum, RecDum, stv, par);
						END

		END; {transients}
		
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
					new.V[i]:=curr.V[i] + 1*(stv.x[i]/stv.Ltot - 0.5)*(new.Vint-curr.Vint);
	OTHERWISE BEGIN
		{so at least 2 accepted solutions are available}
		{first calculate the factors (rV,G,t) by which the voltage, generation and time have increased:}
		IF curr.Vint <> prev.Vint THEN
			rV:=(new.Vint - curr.Vint)/(curr.Vint - prev.Vint)
		ELSE rV:=0;
				IF curr.tijd <> prev.tijd THEN
			rt:=(new.tijd - curr.tijd)/(curr.tijd - prev.tijd)
		ELSE rt:=0; {In SimSS this will always be zero, in ZimT mostly not}
		IF (rV<>0) AND (curr.G_frac = prev.G_frac) THEN 
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
					IF par.lyr[stv.lid[i]].negIons THEN	new.nion[i]:=MAX(0, curr.nion[i] + rV*(curr.nion[i]-prev.nion[i])); {these need to be non-negative}
					IF par.lyr[stv.lid[i]].posIons THEN new.pion[i]:=MAX(0, curr.pion[i] + rV*(curr.pion[i]-prev.pion[i]))
				END; {ions}
			
			END {for loop interior grid points}
		
		ELSE {rV=0 and/or (curr.G_ehp<>prev.G_ehp). In either case, simply extrapolate using the time difference}
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
					IF par.lyr[stv.lid[i]].negIons AND (prev.nion[i]<>0) THEN new.nion[i]:=curr.nion[i]*power(abs(curr.nion[i]/prev.nion[i]), rt); {these need to be non-negative}
					IF par.lyr[stv.lid[i]].posIons AND (prev.nion[i]<>0) THEN new.pion[i]:=curr.pion[i]*power(abs(curr.pion[i]/prev.pion[i]), rt)
				END; {ions}
		
			END {for loop interior grid points}
	
	END {otherwise}
	END {case statement}	

END;

FUNCTION Deterimine_Convergence_Densities(CONSTREF deln, delp, delnion, delpion, n, p, nion, pion : vector; 
							UpdateIons : BOOLEAN; accDens : myReal; VAR ConvMsg : STRING; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : BOOLEAN;
{Determines whether the densities have converged}
VAR totRelChange : myReal;
	j : INTEGER;
BEGIN
	totRelChange:=0; {total relative change of all charges}

	{now compute relative changes of densities:}
	totRelChange:=totRelChange + Norm_Eucl(deln, 0, par.NP+1) / Norm_Eucl(n, 0, par.NP+1);
	totRelChange:=totRelChange + Norm_Eucl(delp, 0, par.NP+1) / Norm_Eucl(p, 0, par.NP+1);

	{now we add the relative changes of the ions}
	IF UpdateIons THEN
		FOR j:=1 TO stv.NLayers DO 
		BEGIN
			{we simply loop over all layers, check if there are ions (NormIons<>0). If so, then we add the relative change to totRelChange}
			IF par.lyr[j].negIons THEN
				totRelChange:=totRelChange + Norm_Eucl(delnion, stv.i0[j], stv.i1[j]) / Norm_Eucl(nion, stv.i0[j], stv.i1[j]);
			IF par.lyr[j].posIons THEN
				totRelChange:=totRelChange + Norm_Eucl(delpion, stv.i0[j], stv.i1[j]) / Norm_Eucl(pion, stv.i0[j], stv.i1[j])
		END;

	ConvMsg:='- relative change of densities: ' + FloatToStrF(totRelChange, ffExponent,5,0) + LineEnding; {our message string}

	{so, converged if change is small enough, we take the acceleration into accout as well:}
	Deterimine_Convergence_Densities:=totRelChange <= accDens*par.tolDens
END;

PROCEDURE Main_Solver(VAR curr, new : TState; VAR it : INTEGER; VAR conv : BOOLEAN; VAR StatusStr : ANSISTRING; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{Iteratively solves the Poisson and continuity equations, including traps and ions}
{can be used in steady-state and transient cases}
VAR i, MaxIt : INTEGER;
	check_Poisson, coupleIonsPoisson, convDensities, AnyMovingIons : BOOLEAN;
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
		THEN MaxIt:=par.maxItSS 
		ELSE MaxIt:=par.maxItTrans;

	{are there ANY moving ions in any of the layers?}
	AnyMovingIons:=FALSE;
	FOR i:=1 TO stv.NLayers DO
		AnyMovingIons:=AnyMovingIons OR par.lyr[i].negIons OR par.lyr[i].posIons;

	coupleIonsPoisson:=AnyMovingIons; {signfies whether ion density can be modified by Poisson solver}

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
		Solve_Poisson(V, n, p, nion, pion, f_tb, f_ti, f_ti_numer, f_ti_inv_denom, Ntb_charge, Nti_charge, curr.f_tb, curr.f_ti, check_Poisson, coupleIonsPoisson, PoissMsg, dti, stv, par); 
		Update_Gen_Pot(V, Vgn, Vgp, stv, par); {update generalised potentials}
		{note: pass (new) p to Cont_Eq_Elec nor (new) n to Cont_Eq_Holes. This is needed so we can also do t=0!}
		Cont_Eq_Elec(n, curr.n, Vgn, Jn, p, mun, gen, Lang, diss_prob, curr.f_tb, curr.f_ti, Rn, stv, par, dti); {calc. new elec. density}
		Cont_Eq_Holes(p, curr.p, Vgp, Jp, n, mup, gen, Lang, diss_prob, curr.f_tb, curr.f_ti, Rp, stv, par, dti); {calc. new hole density}

		{note: transient ion solvers cannot (as yet) do steady-state, as that yields all densities=0!}
		IF UpdateIons THEN 
		BEGIN
			IF dti=0 THEN BEGIN {dti=0 => steady-state} 
				Calc_Ion_Distribution_Steady_State(nion, V, -1, stv.NegIonRegion, stv, par); {use steady-state proc for negative ions}
				Calc_Ion_Distribution_Steady_State(pion, V, 1, stv.PosIonRegion, stv, par) {use steady-state proc for positive ions}
			END
			ELSE BEGIN {use transient versions:}
				Solve_Neg_Ions(nion, curr.nion, V, dti, stv, par); {update neg ions}
				Solve_Pos_Ions(pion, curr.pion, V, dti, stv, par) {update pos ions}
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
		IF convDensities AND AnyMovingIons AND coupleIonsPoisson THEN
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
		IF par.R_shunt>0 THEN Jext:=Jint + Vint/par.R_shunt ELSE Jext:=Jint; {note: infinite R_shunt (no shunt) means R_shunt<0}
		Vext:=Vint + Jext*par.R_series
	END;

	{now construct string to report on our progress:}
	StatusStr:='Overall convergence: '+myBoolStr(conv) + LineEnding;
	StatusStr:=StatusStr + 'Iterations perfomed: '+IntToStr(it) + LineEnding;
	StatusStr:=StatusStr + PoissMsg + LineEnding;
	StatusStr:=StatusStr + ConvMsg;

END;

PROCEDURE Prepare_tJV_File(VAR uitv : TEXT; filename : STRING; transient : BOOLEAN; CONSTREF stv : TStaticVars); 
{create a new tJV_file with appropriate heading
after running this, the TEXT file 'uitv' is still open and ready for writing}
VAR j : INTEGER;
BEGIN
	ASSIGN(uitv, filename);
	REWRITE(uitv); {rewrite old file (if any) or create new one}
    {write header, in the simulation we'll simply output the variables, but not change this header:}
	IF transient THEN WRITE(uitv,' t');
	WRITE(uitv,' Vext Jext errJ Jint ');
	
	{next we show a break down of the photo- and recombination currents for each layer/interface:}
	FOR j:=1 TO stv.NLayers DO WRITE(uitv, 'JphotoL',IntToStr(j),' ');
	FOR j:=1 TO stv.NLayers DO WRITE(uitv, 'JdirL',IntToStr(j),' ');
	IF transient THEN BEGIN {trapping currents for electrons and holes may be different, so show BOTH}
		FOR j:=1 TO stv.NLayers DO WRITE(uitv, 'JbulkElecL',IntToStr(j),' ');
		FOR j:=1 TO stv.NLayers DO WRITE(uitv, 'JbulkHolesL',IntToStr(j),' ');
		FOR j:=1 TO stv.NLayers-1 DO WRITE(uitv, 'JintElecL',IntToStr(j),'L',IntToStr(j+1),' ');		
		FOR j:=1 TO stv.NLayers-1 DO WRITE(uitv, 'JintHolesL',IntToStr(j),'L',IntToStr(j+1),' ')		
	
	END
	ELSE BEGIN {so steady-state}
		FOR j:=1 TO stv.NLayers DO WRITE(uitv, 'JbulkL',IntToStr(j),' ');
		FOR j:=1 TO stv.NLayers-1 DO WRITE(uitv, 'JintL',IntToStr(j),'L',IntToStr(j+1),' ')
	END;
	
	WRITE(uitv,'JminLeft JminRight JShunt'); 
	IF transient THEN WRITELN(uitv,' Jnion Jpion JD') ELSE WRITELN(uitv)
END;

PROCEDURE Write_To_tJV_File(VAR uitv : TEXT; CONSTREF CurrState, PrevState : Tstate; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{before running this proc, uitv must be open (by running Prepare_tJV_File). It must be closed in the main code.
This proc writes the (time), voltage, currents, recombination currents to a file that contains the JV-curve}
VAR JminLeft, JminRight : myReal;
	j, NP : INTEGER;

	{first: 2 short functions that will aid in compact notation:}
	FUNCTION Ave(Vec : vector; onGrid : BOOLEAN = FALSE) : myReal;
	{local version of Average to calc the average over the full interval 0..NP+1 without having to specify the bounds}
	BEGIN
		Ave:=Average(Vec, stv.h, 0, par.NP+1, onGrid);
	END;
	
	FUNCTION EquiCurr(Vec : vector; i0, i1 : INTEGER) : myReal;
	{local function to calc the Equivalent Current (A/m2) of a volume generation or recombination process}
	BEGIN
		EquiCurr:=q*Average(Vec, stv.h, i0, i1, FALSE)*(stv.x[i1]-stv.x[i0])
	END;
	
BEGIN
    WITH CurrState DO 
    BEGIN
        IF transient THEN WRITE(uitv,tijd:nd,' ');  

		{now determine the minority currents left and right}
		NP:=par.NP; {local copy of number of grid points}
		IF n[0]<p[0] THEN JminLeft:=Jn[0] ELSE JminLeft:=Jp[0];
		IF n[NP+1]<p[NP+1] THEN JminRight:=Jn[NP+1] ELSE JminRight:=Jp[NP+1];   
        
        WRITE(uitv,Vext:nd,' ',Jext:nd,' ',errJ:nd,' ',Jint:nd,' ');
        
        {now we look at the photo- and recombination currents in each layer/interface:}
        {first: photo current and direct recombination:}
        FOR j:=1 TO stv.NLayers DO WRITE(uitv, EquiCurr(gen,stv.i0[j],stv.i1[j]):nd,' ');
        FOR j:=1 TO stv.NLayers DO WRITE(uitv, EquiCurr(Rn.direct,stv.i0[j],stv.i1[j]):nd,' ');
        {the other recombination currents can be different in steady-state v transient:}
        IF transient THEN BEGIN
         	FOR j:=1 TO stv.NLayers DO WRITE(uitv, EquiCurr(Rn.bulk,stv.i0[j],stv.i1[j]):nd,' ');       
         	FOR j:=1 TO stv.NLayers DO WRITE(uitv, EquiCurr(Rp.bulk,stv.i0[j],stv.i1[j]):nd,' ');  
         	FOR j:=1 TO stv.NLayers-1 DO WRITE(uitv, 0.5*q*stv.Ltot*stv.h[stv.i1[j]]*(Rn.int[stv.i1[j]] + Rn.int[stv.i1[j]+1]):nd,' ') ;      	     
         	FOR j:=1 TO stv.NLayers-1 DO WRITE(uitv, 0.5*q*stv.Ltot*stv.h[stv.i1[j]]*(Rp.int[stv.i1[j]] + Rp.int[stv.i1[j]+1]):nd,' ')       	     
        END
        ELSE BEGIN
        	FOR j:=1 TO stv.NLayers DO WRITE(uitv, EquiCurr(Rn.bulk,stv.i0[j],stv.i1[j]):nd,' ');
        	FOR j:=1 TO stv.NLayers-1 DO WRITE(uitv, 0.5*q*stv.Ltot*stv.h[stv.i1[j]]*(Rp.int[stv.i1[j]] + Rp.int[stv.i1[j]+1]):nd,' ')
        END;
	
		WRITE(uitv,JminLeft:nd,' ',JminRight:nd,' ',Jext-Jint:nd);

        IF transient 
			THEN WRITELN(uitv,' ',Ave(Jnion):nd,' ',Ave(Jpion):nd,' ',Ave(JD):nd) 
			ELSE WRITELN(uitv);
    END; {with astate}
    FLUSH(uitv)
END;

PROCEDURE Prepare_Var_File(CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN); {create a new var_file with appropriate heading}
VAR uitv : TEXT;
BEGIN
	ASSIGN(uitv, par.varFile);
	REWRITE(uitv); {rewrite old file (if any) or create new one}
    {write header, in the simulation we'll simply output the variables, but not change this header:}
    WRITE(uitv, ' x V Evac Ec Ev phin phip n p ND NA anion cation ntb nti mun mup G_ehp Gfree Rdir BulkSRHn BulkSRHp IntSRHn IntSRHp Jn Jp Jint');

    IF transient 
		THEN WRITELN(uitv,' Jnion Jpion JD lid time') {add time! the ion & displacement currents are zero if not transient!}
		ELSE WRITELN(uitv,' lid Vext'); {add Vext so we can identify the different voltages}
    CLOSE(uitv);
END;

PROCEDURE Write_Variables_To_File(VAR CurrState : TState; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters; transient : BOOLEAN);
{writes the internal variables (astate) to file par.varFile. It assumes the file has a header produced by Prepare_Var_File}
VAR i, j, e : INTEGER;
	uitv : TEXT;
	Evac, Ec, Ev, Jtot, phin, phip, ntb, nti : myReal;
BEGIN
    ASSIGN(uitv, par.varFile);
    APPEND(uitv);

    FOR i:=0 TO par.NP+1 DO
    WITH CurrState DO BEGIN
		Evac:=V[0] - V[i]; {the vacuum level. We take it zero at x=0}
		Ec:=Evac - stv.E_CB[i];
		Ev:=Evac - stv.E_VB[i];
		Jtot:=Jn[i]+Jp[i]+JD[i]+Jnion[i]+Jpion[i]; {total current on grid point}
		{electron and hole quasi-Fermi levels:}
		phin:=Ec + stv.Vt*LN(n[i]/stv.NcLoc[i]);
		phip:=Ev - stv.Vt*LN(p[i]/stv.NcLoc[i]); 

		{calculate number of trapped electrons in bulk/interface traps:}
		j:=stv.lid[i]; {indicates the current layer}
		ntb:=0;

		IF (i>0) AND (i<=par.NP) THEN {as always, exclude the electrodes}
			FOR e:=1 TO stv.Ntb[j].NLevels DO
				ntb:=ntb + stv.Ntb[j].Nt[e] * f_tb[i,e];

		nti:=0;
		{for interfaces: check if we're about to cross an interface:}
		IF i=stv.i1[j] THEN
			FOR e:=1 TO stv.Nti[j].NLevels DO
				nti:=nti + stv.Nti[j].Nt[e] * f_ti[i,e];			

        WRITE(uitv, stv.x[i]:nd,' ',V[i]:nd,' ',
				Evac:nd,' ',Ec:nd,' ',Ev:nd,' ',phin:nd,' ',phip:nd,' ', {band diagram}
				{all charged species:}
				n[i]:nd,' ',p[i]:nd,' ',stv.nid[i]:nd,' ',stv.pid[i]:nd,' ',nion[i]:nd,' ', pion[i]:nd,' ',
				{trapped electrons in bulk and interface traps:}
				ntb:nd,' ',nti:nd,' ',
				{transport:}
				mun[i]:nd,' ',mup[i]:nd,' ',
				{generation:}
				Gm[i]:nd,' ',gen[i]:nd,' ',
				{recombination:}
				Rn.direct[i]:nd,' ',Rn.bulk[i]:nd,' ',Rp.bulk[i]:nd,' ',Rn.int[i]:nd,' ',Rp.int[i]:nd,' ',
				{current densities:}
				Jn[i]:nd,' ',Jp[i]:nd,' ',Jtot:nd);	        
        IF transient 
			THEN WRITELN(uitv,' ',Jnion[i]:nd,' ',Jpion[i]:nd,' ',JD[i]:nd,' ',stv.lid[i],' ',tijd:nd)
			ELSE WRITELN(uitv,' ',stv.lid[i],' ',Vext:nd)
    END;
    CLOSE(uitv);
END;

PROCEDURE Tidy_Up_File(FileName : STRING);
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

PROCEDURE Tidy_Up_Parameter_Files(parameterFile : STRING; QuitWhenDone : BOOLEAN; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure cleans up the parameter files: the main file and the files with the parameters of each layer}
{It makes sure that all the * are aligned, that every parameter has a 
unit or other description/comment and left-aligns any line that starts with **. It does this by reading the original device
parameter file line by line and writing corrected lines to a temp file. Once the correct position of the descriptions 
(starting with *) has been found, it uses the temp file to create the new, tidy parameter file. Lastly, the temp file is
removed and the program exits.}
VAR i : INTEGER;	
BEGIN
	{first, do the main paramter file:}
	Tidy_Up_File(parameterFile);
	
	{next, loop over files that contain the parameters of the layers:}
	FOR i:=1 TO stv.NLayers DO
		Tidy_Up_File(par.lyr[i].layerFile);
		
	IF QuitWhenDone THEN Stop_Prog('Done cleaning parameter files.', EC_Warning)
END;

FUNCTION Copy_State(CONSTREF a : TState; CONSTREF par : TInputParameters) : TState;
{Returns a copy of state a. We need this because of the dynamic arrays in TState.}
VAR b : TState;
	i : INTEGER;
BEGIN
	b:=a; {first copy everything}
	
	{the dynamic arrays in TState need to be copied explicitely to make a real copy. Otherwise, we only copy the pointer!!}
	WITH b DO
		FOR i:=0 TO par.NP+1 DO 
		BEGIN
			{note: SetLength generates a new bit of memory to store a unique state b}
			SETLENGTH(f_tb[i], LENGTH(f_tb[i]));
			SETLENGTH(f_ti[i], LENGTH(f_ti[i]));
			SETLENGTH(f_ti_numer[i], LENGTH(f_ti_numer[i]));
			SETLENGTH(f_ti_inv_denom[i], LENGTH(f_ti_inv_denom[i]))
		END;
	Copy_State:=b
END;

BEGIN 

END.
