PROGRAM ZimT;
{Zimulates Transients and I like cinnamon}
{based on SIMIS 1.09 (which is based on SIMsalabim 3.17)}

{
ZimT:a transient 1D drift-diffusion simulator 
Copyright (c) 2021, 2022, 2023, 2024, 2025, 2026, S. Heester, Dr T.S. Sherkar, Dr V.M. Le Corre, 
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



{If we're using windows we need to make it a console application:}
{$IFDEF WINDOWS}
	{$APPTYPE CONSOLE}
{$ENDIF}

{$MODE DELPHI} {force DELPHI mode}
{$modeSwitch exceptions+} {we need this in Read_tVG to use RAISE}

{$UNITPATH ../Units/} {first tell compiler where our own units are located}

USES {our own, generic ones:}
	 TypesAndConstants,
     InputOutputUtils, 
     NumericalUtils,
     {our drift-diffusion stuff:}
	 DDTypesAndConstants,
	 TransferMatrix,
	 DDRoutines,
	 {and normal units:}
	 Math,  {for min and max functions}
	 StrUtils, {for DelSpace}
     SysUtils, {for TRIM function in Init_Elec_Mob_Table}
     Crt; {for textcolor}


CONST
    ProgName = TProgram.ZimT;  
    version = '5.36';  


{first: check if the compiler is new enough, otherwise we can't check the version of the code}
{$IF FPC_FULLVERSION < 30200} {30200 is 3.2.0}
	{$STOP FPC VERSION SHOULD BE AT LEAST 3.2.0}
{$ENDIF}
{now check to see if the versions of the units match that of this code:}
{$IF (DDRoutinesVersion <> version) OR (DDTypesAndConstantsVersion <> version) OR (TransferMatrixVersion <> version)} 
	{$STOP Wrong version of one or more units!}
{$ENDIF}


VAR parameterFile : ShortString;

	MainIt, CountAcceptedSolutions, CounttVGPoints, CountStatic : INTEGER;

	prev, curr, new, ref : TState; {store the previous point in time, the current one and the new}
	{prev: solved, stored, done, 2 time steps ago
	curr: solved, stored, done, 1 time step ago
	new: to be solved, latest time that was read
	ref: only used for tracking, not for solving, but we need it to be global as we need to update it in the main loop and use it in the solver}

	stv : TStaticVars; {all variables that are calculated at the start of the simulation and then remain constant}

	par : TInputParameters; {all input parameters}

	ResVmn, ResVmx, Vmn, Vmx : myReal;

    inv, uitv, log : TEXT; {the input and output files}
   
    dumstr : STRING;

    conv, foundtVG, keepGoing, staticSystem, acceptNewSolution, foundFirstRef : BOOLEAN;
  
    MsgStr : ANSISTRING = ''; {Ansistrings have no length limit, init string to ''}

	StatusStr : ANSISTRING; 


PROCEDURE Read_tVG(VAR astate : TState; old_tijd : myReal; VAR inv, log : TEXT; VAR foundtVG, foundFirstRef : BOOLEAN);
{try to read a line of tijd, Vext, G_frac}
VAR varline, orgline, parstr, msg : STRING;
	SimOC : BOOLEAN;
BEGIN
	foundtVG:=FALSE;
	IF NOT EOF(inv)
	THEN {there is something in the file}
	TRY {try to read it, but be careful}
		READLN(inv, orgline); {read a line of variables, we will store the orignal line (orgline)}
		varline:=DelWhite1(orgline); {returns a copy of str with all white spaces (ASCII code 9,..13, and 32) reduced to 1 space}
		varline:=TrimLeft(varline); {remove any spaces on the left}
		IF varline<>'' THEN
		BEGIN {OK, we might have something!}
			{Copy2SpaceDel, strutils: Deletes and returns all characters in a string till the first space character (not included).}
			parstr:=Copy2SpaceDel(varline); {contains the first parameter}
			{ConvertStrToFloat, InputOutputUtils: Convert a string to a floating-point value.}
			IF NOT ConvertStrToFloat(parstr, astate.tijd) THEN {first part contains the time}
				RAISE Exception.Create(''); {Raise empty exception, message is set in the EXCEPT block}
		
			{check the new time (astate.tijd) and compute the 1/timestep (=dti)}
			IF (astate.tijd > old_tijd) THEN {simply progress in time}
				astate.dti:=1/(astate.tijd-old_tijd) {dti: inverse of time step} 
			ELSE BEGIN {not going forward in time}
				IF astate.tijd=0 THEN 
					astate.dti:=0 {we're either staying in steady-state or going back to it. timestep = infinity -> dti=0}
				ELSE Stop_Prog_Finalize_Log(log, 'Time steps must either be positive, or you should go back to t=0 (steady-state).', EC_InvalidInput) {going back in time (but not to steady-state) or simply keeping the same time is not allowed!}
			END;
			
			parstr:=Copy2SpaceDel(varline); {contains the second parameter}
			IF LowerCase(parstr)='oc' {now check if we should simulate at open-circuit, or just at some specific value}
			THEN SimOC:=TRUE
			ELSE WITH astate DO 
			BEGIN
				SimOC:=FALSE;
				IF NOT ConvertStrToFloat(parstr, Vext) THEN
					RAISE Exception.Create(''); {Raise empty exception, message is set in the EXCEPT block}
				{Note: ZimT takes the voltage in the tVGFile (in parstr) to be the external voltage, Vext}
				{if simtype=1, then Vint=Vext.}
				{if simytype=2,3, then we're simulating Voc, so Vint&Vext will need to be solved, thus the exact value doesn't matter}
				{if simtype=4, then we only know Vext and need to solve for Vint.}
				{therefore, for all simtypes, we simply set Vext=Vint for new:} 		
				Vint:=Vext; 
				{check if V is not too large or small:}
				IF Vext*stv.Vti < -1.95 * LN(Max_Value_myReal) THEN Stop_Prog_Finalize_Log(log, 'V is too small.', EC_InvalidInput);
				IF Vext*stv.Vti > 1.95 * LN(Max_Value_myReal) THEN Stop_Prog_Finalize_Log(log, 'V is too large.', EC_InvalidInput);
			END;

			parstr:=Copy2SpaceDel(varline); {contains the third parameter}
			IF NOT ConvertStrToFloat(parstr, astate.G_frac) THEN
				RAISE Exception.Create(''); {Raise empty exception, message is set in the EXCEPT block}

			parstr:=Copy2SpaceDel(varline); {contains the fourth parameter track}
			IF NOT ((parstr = '0') OR (parstr = '1') OR (parstr = '2') OR (parstr = '3')) THEN {we only accept 0-3}
				RAISE Exception.Create(''); {Raise empty exception, message is set in the EXCEPT block}
			astate.Track:=StrToInt(parstr);

			IF NOT foundFirstRef AND (astate.Track=1) THEN
				foundFirstRef:=TRUE; {set to true if astate.Track=1 for the first time}
			{we need to check if the first non-zero Track is equal to 1:}
			IF NOT foundFirstRef AND (astate.Track IN [2,3]) THEN
				Stop_Prog_Finalize_Log(log, 'Error in tVGFile: before setting Track=2 or 3, you need to have at least one Track=1', EC_InvalidInput);
			
			{once we're here, foundtVG is true!}
			foundtVG:=TRUE;
		
			{now determine which kind of simulation we have to do:}
			astate.SimType:=1; {default: we know the internal voltage and it's equal to the external one (R_series=0)}
			IF SimOC AND (astate.tijd=0) THEN astate.SimType:=2; {steady-state, open-circuit}
			IF SimOC AND (astate.tijd>0) THEN astate.SimType:=3; {transient, open-circuit}
			IF (NOT SimOC) AND (par.R_series>0) THEN astate.SimType:=4; {steady-state or transient, not open-circuit, R_series important}
		END;
	EXCEPT {reading didn't work, raise exception}
		msg:='Error while reading from file ' + par.tVGFile + LineEnding;
		msg:=msg + 'Offending line: ' + LineEnding;
		msg:=msg + orgline + LineEnding;
		WRITELN(orgline);
		Stop_Prog_Finalize_Log(log, msg + 'See Reference Manual for details.', EC_InvalidInput);
	END 
END;

PROCEDURE Open_and_Read_tVG_file(VAR inv, log : TEXT; VAR new : TState; CONSTREF par : TInputParameters); 
{open tVG file, read header and first time/voltage/G_ehp|G_frac}
VAR foundHeader : BOOLEAN;
BEGIN
	{open input file with times, voltage and generation rate}
	IF NOT FileExists(par.tVGFile) {the file with input par. is not found}
        THEN Stop_Prog_Finalize_Log(log, 'Could not find file '+par.tVGFile, EC_FileNotFound);
	ASSIGN(inv, par.tVGFile);
	RESET(inv);

	WRITELN('Reading t, Vext, G_frac Track from file ',par.tVGFile);
	{now try to read from input file until we find 't Vext G_frac'}
	REPEAT
		READLN(inv, dumstr);
		foundHeader:=LeftStr(DelWhite(LowerCase(dumstr)),16) ='tvextg_fractrack'; {cut relevant part of string}
		{DelWhite removes all white spaces (incl. tabs, etc.) from a string, see myUtils}
	UNTIL foundHeader OR EOF(inv);
		
	IF NOT foundHeader THEN Stop_Prog_Finalize_Log(log, 'Could not find correct header ''t Vext G_frac Track'' in ' + par.tVGFile + '.', EC_InvalidInput);

	Read_tVG(new, 0, inv, log, foundtVG, foundFirstRef); {try to read a line of t, Va, G_ehp}
	IF NOT(foundtVG) OR (new.tijd<>0) THEN Stop_Prog_Finalize_Log(log, 'tVGFile did not specify steady-state (t=0)', EC_InvalidInput)
END;

FUNCTION Residual_Current_Voltage(Vint : myReal; VAR curr, new : TState; VAR conv : BOOLEAN; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : myReal;
{Calculates the difference between the simulated current or voltage and the targeted ones, i.e. the residual.
It takes a guess for Vint and then solves the equations. This yields the simulated current or voltage at that Vint. 
If we're trying to solve for Voc (SimType 2 or 3), then we want the current to be zero. So the target current is 0
and the residual is the simulated Jext. 
If SimType = 4, then we know the target (requested) Vext (which we will store locally), but we do not know the
Vint that will yield this Vext: Series resistance will make Vext <> Vint. So, we take the input (guess) for Vint and solve
the equations to get Vext. The residual, in this case, is then equal to the difference between the new Vext and the 
one we're aiming for (target).}
VAR it : INTEGER;
	VextTarget : myReal;
BEGIN
	VextTarget:=new.Vext; {store Vext. This is the Vext we are TRYING to get (only if simtype=4)}
	new.Vint:=Vint; {apply the guessed internal voltage (Vint)}
	Main_Solver(curr, new, it, conv, StatusStr, stv, par); {run the solver}
	CASE new.SimType OF
		2,3 :  {steady-state, resp. transient open-circuit}
				Residual_Current_Voltage:=new.Jext; {residual current is simply equal to device current}
		4 	:  {SS or transient, not open-circuit, R_series significant}
				BEGIN
					Residual_Current_Voltage:=new.Vext - VextTarget;	{what we should be calculating is the difference between the target Vext and the realised one}
					new.Vext:=VextTarget; {restore Vext in state new}
				END
		ELSE Stop_Prog('Invalid case in function ResCurr', EC_ProgrammingError);
	END; {case}
END;


PROCEDURE Bracket_device_voltage(VAR Vmin, Vmax, Vguess, ResJmin, ResJmax : myReal; VAR curr, new : TState;
								 CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure finds bracketing voltages around Vint, based on a guess (Vguess). Vmin and Vmax bracket Vint, unless unsuccessful}
{See Numerical recipes section 9.1 on Bracketing and Bisection}
CONST ntry = 50; {max. number of iterations}
	  factor = 1.6; {multiplication factor for enlarging interval}
VAR
	it : INTEGER;
	conv, success : BOOLEAN;
BEGIN
	{start with a guess for Vmin,max around Vguess}
	Vmin:=Vguess - 20*par.tolVint;
	Vmax:=Vguess + 20*par.tolVint;
	{calc. the corresponding currents:}
	ResJMin:=Residual_Current_Voltage(Vmin, curr, new, conv, stv, par);
	ResJMax:=Residual_Current_Voltage(Vmax, curr, new, conv, stv, par);

	{we try to find Vmin, Vmax such that ResJmin * ResJmax < 0}
	success:=ResJmin*ResJmax < 0; {note: we don't check if ResCurr converged...}
	it:=0;
	WHILE (NOT success) AND (it<ntry) DO
	BEGIN
		INC(it); 
		IF ABS(ResJmin)<ABS(ResJmax) THEN BEGIN
			Vmin:=Vmin + factor*(Vmin-Vmax);
			ResJMin:=Residual_Current_Voltage(Vmin, curr, new, conv, stv, par);
		END
		ELSE BEGIN
			Vmax:=Vmax + factor*(Vmax-Vmin);
			ResJMax:=Residual_Current_Voltage(Vmax, curr, new, conv, stv, par);
		END;
		success:=ResJmin*ResJmax < 0
	END;

	IF NOT success THEN Stop_Prog('Could not find bracketing values for Vint at time ' + FloatToStrF(new.tijd, ffGeneral,10,0), EC_NumericalFailure);
END;

PROCEDURE Find_Vint_Brent(a, b, fa, fb: myReal; VAR curr, new: TState; VAR conv: BOOLEAN; 
                         CONSTREF stv: TStaticVars; CONSTREF par: TInputParameters);
{Uses Brent's algorithm to find the value of Vint that makes Residual_Current_Voltage = 0}
VAR
    c, d, e, p, q, r, s, tol, m : myReal;
    fc : myReal;
    ItVint : INTEGER;
    converged : BOOLEAN;
    EPS : myReal; {machine floating-point precision, now a variable based on myReal type}
BEGIN
	{Set the maximum number of iterations based on what would be needed for bisection}
	IF (b - a) <= 0 THEN Stop_Prog('Invalid range for Brent''s method', EC_NumericalFailure);

    {Set EPS based on the size of myReal type}
    CASE SizeOf(myReal) OF
        4: EPS := 1.19E-7;  {SINGLE precision ~ 2^-23}
        8: EPS := 2.22E-16; {DOUBLE precision ~ 2^-52}
        10: EPS := 1.08E-19; {EXTENDED precision ~ 2^-63}
        ELSE EPS := 3.0E-8; {fallback value in case of unexpected type size}
    END;
        
    {Check if the root is bracketed}
    IF fa * fb > 0 THEN Stop_Prog('Root must be bracketed for Brent''s method', EC_NumericalFailure);
    c := a; fc := fa;
    
    {If f(b) is closer to zero than f(a), swap a and b}
    IF ABS(fa) < ABS(fb) THEN
    BEGIN
        a := b; b := c; c := a;
        fa := fb; fb := fc; fc := fa;
    END;
    
    d := 0; e := 0; {Initialize these variables}
    ItVint := 0;
    converged := FALSE;
    
    REPEAT
        INC(ItVint);
        
        {Make sure c remains the endpoint with function value having opposite sign from f(b)}
        IF fb * fc > 0 THEN
        BEGIN
            c := a; fc := fa;
            d := b - a; e := d;
        END;
        
        {Ensure |f(c)| >= |f(b)|}
        IF ABS(fc) < ABS(fb) THEN
        BEGIN
            a := b; b := c; c := a;
            fa := fb; fb := fc; fc := fa;
        END;
        
        {Convergence check}
        tol := 2.0 * EPS * ABS(b) + 0.5 * par.tolVint;
        m := 0.5 * (c - b);
        
        {Check if we've converged or found exact root}
        converged := (ABS(m) <= tol) OR (ABS(fb) < EPS);
        IF NOT converged THEN
        BEGIN
            {Decide whether to use bisection or interpolation}
            IF (ABS(e) >= tol) AND (ABS(fa) > ABS(fb)) THEN
            BEGIN
                {Try inverse quadratic interpolation}
                s := fb / fa;
                
                IF a = c THEN
                BEGIN
                    {Linear interpolation (secant method)}
                    p := 2.0 * m * s;
                    q := 1.0 - s;
                END
                ELSE
                BEGIN
                    {Inverse quadratic interpolation}
                    q := fa / fc;
                    r := fb / fc;
                    p := s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                    q := (q - 1.0) * (r - 1.0) * (s - 1.0);
                END;
                
                {Check bounds}
                IF p > 0 THEN q := -q;
                p := ABS(p);
                
                {Accept interpolation if it's reasonable}
                IF (2.0 * p < MIN(3.0 * m * q - ABS(tol * q), ABS(e * q))) THEN
                BEGIN
                    e := d; d := p / q;
                END
                ELSE
                BEGIN
                    {Fall back to bisection}
                    d := m; e := d;
                END;
            END
            ELSE
            BEGIN
                {Use bisection}
                d := m; e := d;
            END;
            
            {Save previous values}
            a := b; fa := fb;
            
            {Update b to be our new guess}
            IF ABS(d) > tol THEN
                b := b + d
            ELSE
            BEGIN
                {Ensure we move at least tol in the direction of the root}
                IF m > 0 THEN b := b + tol
                ELSE b := b - tol;
            END;
            
            {Calculate function value at new position}
            fb := Residual_Current_Voltage(b, curr, new, conv, stv, par);
        END;
    UNTIL converged OR (ItVint >= MaxItVintBrent);
    
    {Set the result}
    new.Vint := b;
    
    {If we exited due to maximum iterations, show warning}
    IF ItVint >= MaxItVintBrent THEN
        WRITELN('Warning: Brent''s method reached maximum iterations without converging.');
END;

FUNCTION States_Are_Equal(CONSTREF A, B : TState; CONSTREF par : TInputParameters) : BOOLEAN;
{Checks if two states are equal, based on their Vint, G_frac and Jint. This is used for checking if the system is static (not changing anymore). We use the same criteria as for accepting a solution, i.e. the difference in Vint and Jint should be smaller than 1, and G_frac should be exactly the same.}
VAR
	VSame : BOOLEAN;
	
	FUNCTION Rel_Diff_Vec(CONSTREF v1, v2 : vector) : myReal;
	VAR norm : myReal;
	{calc the difference between 2 vectors, normalized to their norm. If the norm is zero, then function result is zero}
	BEGIN
		norm:=0.5*(Norm_Eucl(v1, 0, par.NP+1) + Norm_Eucl(v2, 0, par.NP+1));
		IF norm=0 THEN
			Rel_Diff_Vec:=0
		ELSE
			Rel_Diff_Vec:=Norm_Eucl(Difference(v1, v2, 0, par.NP+1), 0, par.NP+1)/norm;
	END;
	
BEGIN

	{first thing to check: G_frac. No need to check Vext as we will check the potentials any way}
	IF A.G_frac <> B.G_frac THEN {if the experimental conditions are different, then the states are not the same}
		States_Are_Equal:=FALSE
	ELSE {light intensity is identical, but Vext and the internal variables might be different.}
	BEGIN
		{First: compare the internal potentials, so if Vext is different, this will show up here}
		VSame:=Norm_Eucl(Difference(A.V, B.V, 0, par.NP+1), 0, par.NP+1) < par.tolPois; {note: we always compare potentials in Volt, so it's absolute!}

		IF VSame THEN {only if potentials are the same, is there a point in doing any of the rest!}
		BEGIN

			{now we need to compare the rest of the internal variables, but this depends on par.convVar}	
			CASE par.convVar OF
			{note: fpc by default does shortcut evaluation of boolean expressions. So if the first condition (cond_1) in a list of
			cond_1 AND cond_2 AND ... is false, then the rest will not be evaluated. This saves time!}	
				1 : States_Are_Equal:=(Rel_Diff_Vec(A.n, B.n)<par.tolDens) AND (Rel_Diff_Vec(A.p, B.p)<par.tolDens)
										AND (Rel_Diff_Vec(A.nion, B.nion)<par.tolDens) AND (Rel_Diff_Vec(A.pion, B.pion)<par.tolDens);
				2 : States_Are_Equal:=ABS(A.Jint-B.Jint) <= MAX(par.tolCurr, MIN(A.errJ,B.errJ));
				3 : States_Are_Equal:=(ABS(A.Jint-B.Jint) <= MAX(par.tolCurr, MIN(A.errJ,B.errJ)))
										OR ((Rel_Diff_Vec(A.n, B.n)<par.tolDens) AND (Rel_Diff_Vec(A.p, B.p)<par.tolDens)
										AND (Rel_Diff_Vec(A.nion, B.nion)<par.tolDens) AND (Rel_Diff_Vec(A.pion, B.pion)<par.tolDens));
				4 : States_Are_Equal:=(ABS(A.Jint-B.Jint) <= MAX(par.tolCurr, MIN(A.errJ,B.errJ)))
										AND (Rel_Diff_Vec(A.n, B.n)<par.tolDens) AND (Rel_Diff_Vec(A.p, B.p)<par.tolDens)
										AND (Rel_Diff_Vec(A.nion, B.nion)<par.tolDens) AND (Rel_Diff_Vec(A.pion, B.pion)<par.tolDens)						
			END {of case}
		END
	END
END;

BEGIN {main program}
	Print_Welcome_Message(ProgName, version); {Welcomes the user, and shows the name and version of the program and the authors}

    {if '-h' or '-H' option is given then display some help and exit:}
    IF hasCLoption('-h') THEN Display_Help_Exit(ProgName);
    Determine_Name_Parameter_File(parameterFile); {either default or user-specified file with all the parameters}
    IF NOT Correct_Version_Parameter_File(parameterFile, version, TRUE, ProgName) THEN Stop_Prog('Version of ZimT and '+parameterFile+' do not match.', EC_DevParCorrupt);

{Initialisation:}
    Read_Parameters(parameterFile, MsgStr, par, stv, ProgName); {Read parameters from input file}
    Prepare_Log_File(log, MsgStr, par, version); {open log file}
    Check_Parameters(log, stv, par, ProgName); {perform a number of chekcs on the paramters. Note: we need Vt}
    Set_Number_Digits(par.limitDigits, SizeOf(myReal)); {limits number of digits in floating point}
    IF par.autoTidy THEN Tidy_Up_Parameter_Files(parameterFile, FALSE, stv, par); {clean up file but don't exit!}
    
    Make_Grid(stv, par); {Initialize the grid}
    Define_Layers(stv, par); {define layers: Note, stv are not CONSTREF as we need to change them}
	Init_Trapping(log, stv, par); {Inits all variables needed for trapping and SRH recombination}

	Open_and_Read_tVG_file(inv, log, new, par); {open tVG file, read header and first time/voltage/G_ehp}

	WITH new DO Init_Pot_Dens_Ions_Traps(V, Vgn, Vgp, n, p, nion, pion, f_tb, f_ti, f_ti_numer, f_ti_inv_denom, Vint, stv, par); {init. (generalised) potentials and densities}

	new.UpdateIons:=TRUE; {in ZimT this is always true as we don't artificially fix the ions like we can in SimSS}
	{just to make sure these are initialised!}
	curr:=new;
	prev:=curr;
	ref:=curr; {ref is used as a reference for tracking, but we need to initialise it as well}

	Prepare_tJV_File(uitv, par.tJFile, TRUE, stv, par);   {create the tJV-file}
	IF par.StoreVarFile THEN Prepare_Var_File(stv, par, TRUE); {create a new var_file with appropriate heading}

	WRITELN('The calculation has started, please wait.');

    {Init all parameters and start solving for steady-state (tijd=0)}
	Init_Generation_Profile(stv, log, par); {init. the stv.orgGm array. This is the SHAPE of the profile}
	Update_Generation_Profile(stv.orgGm, new.Gm, new.G_frac, stv, par); {update current Gm array to correct value}
	CountAcceptedSolutions:=0; {counts the number of solutions that were accepted}
	CounttVGPoints:=0; {count number of transient (t,V,G) points}
	countStatic:=0;
	conv:=FALSE;
	keepGoing:=foundtVG; {=true at this point!}
	foundFirstRef:=FALSE; {indicates whether we've found the first ref state in the tVGFile, so a 'Track=1'}

	WHILE keepGoing DO
	BEGIN {main loop over times and t,V,G from input file:}
		Update_Generation_Profile(stv.orgGm, new.Gm, new.G_frac, stv, par); {update current Gm array to correct value}
		Extrapolate_Solution(prev, curr, new, CountAcceptedSolutions, stv, par);
	
		IF new.SimType = 1 
		THEN Main_Solver(curr, new, MainIt, conv, StatusStr, stv, par) {easy: Vint = Vext, Vext is in input tVGFile}
		ELSE BEGIN {SimType > 1: Voc or R_series <> 0: in both cases, we don't know Vint}
			{first find Vmn and Vmx, voltages that bracket the device voltage new.Vint}
			Bracket_device_voltage(Vmn, Vmx, curr.Vint, ResVmn, ResVmx, curr, new, stv, par);
			
			{Then use Brent's algorithm to find Vint precisely}
			Find_Vint_Brent(Vmn, Vmx, ResVmn, ResVmx, curr, new, conv, stv, par);
		END; {SimType > 1}

		IF conv THEN BEGIN
			INC(CountAcceptedSolutions);
			acceptNewSolution:=TRUE
		END
		ELSE BEGIN {o dear, now what?}
			{put error messages in log file:}
			WRITELN(log);
			WRITELN(log, 'Messages from main solver at time = ',FloatToStrF(new.tijd, ffGeneral,10,0),':');
			WRITELN(log, StatusStr);
			FLUSH(log);
			{now assess whether we accept the new solution, or skip it, or quit:}
			CASE par.failureMode OF
				0 : Stop_Prog_Finalize_Log (log, 'Convergence failed at time = ' + FloatToStrF(new.tijd, ffGeneral,10,0)+ '. Maybe try smaller time steps.', EC_ConvergenceFailedHalt);
				1 : acceptNewSolution:=TRUE;
				2 : IF curr.Track <> 1 THEN {'normal case'}
						acceptNewSolution:=(new.dti=0) {is true if steady-state, false otherwise}
					ELSE {current state is the reference state (Track=1), so we're extra careful:}
						Stop_Prog_Finalize_Log (log, 'Convergence failed and Track=1 at time = ' + FloatToStrF(new.tijd, ffGeneral,10,0)+ '.', EC_ConvergenceFailedHalt);						
			END;
			{if we get here, then conv=false, but we did not halt the program, so set the ExitCode:}
			ExitCode:=EC_ConvergenceFailedNotHalt
		END;

		IF acceptNewSolution {OK, new solution is good (even if conv might be FALSE)}
		THEN BEGIN
			prev:=Copy_State(curr, par); {we move the current solution to the previous one}
			curr:=Copy_State(new, par); {and we keep the newest solution}
			{output:}
			Write_To_tJV_File(uitv, curr, prev, stv, par, TRUE)
		END
		ELSE BEGIN 
			WRITELN(log, 'Skipping (t,Vext,G_ehp) point at time ',FloatToStrF(new.tijd, ffGeneral,10,0)); 
			FLUSH(log) 
		END;

		IF CounttVGPoints MOD par.outputRatio = 0 THEN {output to screen and var_file every outputRatio timesteps}
        BEGIN
			IF NOT Conv THEN TextColor(LightRed) ELSE TextColor(LightGray); {Reset to default font colour: it's not white, it's light grey!}
			WITH new DO WRITELN('Time:',tijd:12,' Vext: ',Vext:6:4,' G_frac:',G_frac:11,' Jext:',Jext  :7:2,' +- ',errJ:4:2);
			IF par.StoreVarFile AND acceptNewSolution THEN Write_Variables_To_File(curr, stv, par, TRUE)
        END;
	
		INC(CounttVGPoints);
		
		{Now we have to assess if we have to do automatic stopping, based on curr.Track:}
		
		CASE curr.Track OF {if 0, we don't do anything}
			1 : IF Conv THEN 
					ref:=Copy_State(curr, par); {if we're tracking, we need to update the ref state as well}	
			2 : IF Conv OR (par.failureMode = 1) THEN 
					IF States_Are_Equal(ref, curr, par) 
						THEN INC(CountStatic)
						ELSE CountStatic:=0; {reset counter}
			3 :	IF Conv OR (par.failureMode = 1) THEN 
				BEGIN
					IF States_Are_Equal(ref, curr, par) 
						THEN INC(CountStatic)
						ELSE CountStatic:=0; {reset counter}				
					
					ref:=Copy_State(curr, par) {the current state becomes the new reference}
				END
		END;

		staticSystem:=(CountStatic >= MinCountStatic);

		{now see if there's more in the tVGFile:}
		Read_tVG(new, curr.tijd, inv, log, foundtVG, foundFirstRef); {try to read a line of t, Va, G_ehp}

		keepGoing:=foundtVG AND NOT staticSystem
    END; {main loop}

	MsgStr:=''; {reset this string and fill it with stuff we want to show on screen and in the log file}

	IF staticSystem AND foundtVG THEN
	BEGIN
		ExitCode:=EC_AutoStop; {this sets the exitcode to the correct one (EC_AutoStop) without halting the code (like Stop_Prog would do)}
		MsgStr:=MsgStr + 'Stopped ZimT as the simulation has stabilised.' + LineEnding
	END;

	TextColor(LightGray); {whatever happened, switch back to normal colour}
    CLOSE(uitv);

	MsgStr:=MsgStr + IntToStr(CounttVGPoints) + ' (t,V,G)-points' + LineEnding;
	IF CountAcceptedSolutions < CounttVGPoints THEN
		MsgStr:=MsgStr + IntToStr(CounttVGPoints-CountAcceptedSolutions) + ' did not converge.' + LineEnding;

	MsgStr:=MsgStr + 'Output written in file(s) ' + par.tJFile;
	IF par.StoreVarFile THEN MsgStr:=MsgStr + ' and ' + par.varFile;
	MsgStr:=MsgStr + '.' + LineEnding;

	WRITELN(MsgStr); {write messages on screen}
	Finalize_Log_File(log, MsgStr); {writes final comments, date, time, run time and closes log file.}

    WRITELN('Finished, press enter to exit');
    IF par.pauseAtEnd THEN READLN {pause at the end of the program}

END.
