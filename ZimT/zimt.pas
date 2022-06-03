PROGRAM ZimT;
{Zimulates Transients and I like cinnamon}
{based on SIMIS 1.09 (which is based on SIMsalabim 3.17)}

{
ZimT:a transient 1D drift-diffusion simulator 
Copyright (c) 2020, 2021, 2022 Dr T.S. Sherkar, Dr V.M. Le Corre, M. Koopmans,
F. Wobben, and Prof. Dr L.J.A. Koster, University of Groningen
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

{$MODE DELPHI} {force DELPHI mode}


{$UNITPATH ../Units/} {first tell compiler where our own units are located}

USES {our own, generic ones:}
	 TypesAndConstants,
     InputOutputUtils, 
     NumericalUtils,
     {our drift-diffusion stuff:}
	 DDTypesAndConstants,
	 DDRoutines,
	 {and normal units:}
	 Math,  {for min and max functions}
	 StrUtils, {for DelSpace}
     SysUtils, {for TRIM function in Init_Elec_Mob_Table}
     Crt; {for textcolor}


CONST
    ProgName = TProgram.ZimT;  
    version = '4.35';  


{first: check if the compiler is new enough, otherwise we can't check the version of the code}
{$IF FPC_FULLVERSION < 30200} {30200 is 3.2.0}
	{$STOP FPC VERSION SHOULD BE AT LEAST 3.2.0}
{$ENDIF}
{now check to see if the versions of the units match that of this code:}
{$IF (DDRoutinesVersion <> version) OR (DDTypesAndConstantsVersion <> version)} 
	{$STOP Wrong version of the DD units!}
{$ENDIF}


VAR MainIt, ItVint, MaxItVint, CountNotConv, CounttVGPoints, CountStatic : INTEGER;

	prev, curr, new : TState; {store the previous point in time, the current one and the new}
	{prev: solved, stored, done, 2 time steps ago
	curr: solved, stored, done, 1 time step ago
	new: to be solved, latest time that was read}

	stv : TStaticVars; {all variables that are calculated at the start of the simulation and then remain constant}

	par : TInputParameters; {all input parameters}

	ResVint, ResVmn, Vmn, Vmx : myReal;

    inv, uitv, log : TEXT; {the input and output files}
   
    dumstr : STRING;

    conv, foundtVG, keepGoing, staticSystem, acceptNewSolution : BOOLEAN;
  
    MsgStr : ANSISTRING = ''; {Ansistrings have no length limit, init string to ''}

	StatusStr : ANSISTRING; 

PROCEDURE Read_tVG(VAR astate : TState; old_tijd : myReal; VAR inv : TEXT; VAR foundtVG : BOOLEAN);
{try to read a line of tijd, Va, Gehp}
VAR varline, orgline, parstr : STRING;
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
			{StrToFloat, sysutils: Convert a string to a floating-point value.}
			astate.tijd:=StrToFloat(parstr); {first part contains the time}
			IF old_tijd<>astate.tijd
				THEN astate.dti:=1/(astate.tijd-old_tijd) {dti: inverse of time step}
				ELSE astate.dti:=0; {old_tijd=tijd, so we're looking at steady-state, so infinite time step}
			IF astate.dti<0 THEN Stop_Prog('Negative time steps are not allowed');
			parstr:=Copy2SpaceDel(varline); {contains the second parameter}
			IF LowerCase(parstr)='oc' {now check if we should simulate at open-circuit, or just at some specific value}
			THEN SimOC:=TRUE
			ELSE WITH astate DO 
			BEGIN
				SimOC:=FALSE;
				Vext:=StrToFloat(parstr);
				{Note: ZimT takes the voltage in the tVG_file (in parstr) to be the external voltage, Vext}
				{if simtype=1, then Vint=Vext.}
				{if simytype=2,3, then we're simulating Voc, so Vint&Vext will need to be solved, thus the exact value doesn't matter}
				{if simtype=4, then we only know Vext and need to solve for Vint.}
				{therefore, for all simtypes, we simply set Vext=Vint for new:} 		
				Vint:=Vext; 
				{check if V is not too large or small:}
				IF Vext*stv.Vti < -1.95 * LN(Max_Value_myReal) THEN Stop_Prog('V is too small.');
				IF Vext*stv.Vti > 1.95 * LN(Max_Value_myReal) THEN Stop_Prog('V is too large.');
			END;
			parstr:=Copy2SpaceDel(varline); {contains the third parameter}

			astate.Gehp:=StrToFloat(Copy2SpaceDel(parstr));
			foundtVG:=TRUE;
		
			{now determine which kind of simulation we have to do:}
			astate.SimType:=1; {default: we know the internal voltage and it's equal to the external one (Rseries=0)}
			IF SimOC AND (astate.tijd=0) THEN astate.SimType:=2; {steady-state, open-circuit}
			IF SimOC AND (astate.tijd>0) THEN astate.SimType:=3; {transient, open-circuit}
			IF (NOT SimOC) AND (par.Rseries>0) THEN astate.SimType:=4; {steady-state or transient, not open-circuit, Rseries important}
		END;
	EXCEPT {reading didn't work, raise exception}
		WRITELN('Error while reading from file ',par.tVG_file);
		WRITELN('Offending line: ');
		WRITELN(orgline);
		Stop_Prog('See Reference Manual for details.');
	END 
END;

PROCEDURE Open_and_Read_tVG_file(VAR inv : TEXT; VAR new : TState; CONSTREF par : TInputParameters); 
{open tVG file, read header and first time/voltage/Gehp}
VAR foundHeader : BOOLEAN;
BEGIN
	{open input file with times, voltage and generation rate}
	IF NOT FileExists(par.tVG_file) {the file with input par. is not found}
        THEN Stop_Prog('Could not find file '+par.tVG_file);
	ASSIGN(inv, par.tVG_file);
	RESET(inv);
	
	WRITELN('Reading t, Vext, Gehp from file ',par.tVG_file);
	{now try to read from input file until we find 't Vext Gehp'}
	REPEAT
		READLN(inv, dumstr);
		foundHeader:=LeftStr(DelWhite(dumstr),9) ='tVextGehp'; {cut relevant part of string}
	UNTIL foundHeader OR EOF(inv);
	
	IF NOT foundHeader THEN Stop_Prog('Could not find correct header ''t Vext Gehp'' in ' + par.tVG_file + '.');
	
	{DelWhite removes all white spaces (incl. tabs, etc.) from a string, see myUtils}
	Read_tVG(new, 0, inv, foundtVG); {try to read a line of t, Va, Gehp}
	IF NOT(foundtVG) OR (new.tijd<>0) THEN Stop_Prog('tVG_file did not specify steady-state (t=0)');
	IF new.SimType=2 THEN new.SimType:=3; {both 2 and 3 are open-circuit, but 2 = S-S and 3 = transient. that doesn't matter and we take 3}
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
		4 	:  {SS or transient, not open-circuit, Rseries significant}
				BEGIN
					Residual_Current_Voltage:=new.Vext - VextTarget;	{what we should be calculating is the difference between the target Vext and the realised one}
					new.Vext:=VextTarget; {restore Vext in state new}
				END
		ELSE Stop_Prog('Invalid case in function ResCurr');
	END; {case}
END;


PROCEDURE Bracket_device_voltage(VAR Vmin, Vmax, Vguess, ResJmin : myReal; VAR curr, new : TState;
								 CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
{This procedure finds bracketing voltages around Vint, based on a guess (Vguess). Vmin and Vmax bracket Vint, unless unsuccessful}
{See Numerical recipes section 9.1 on Bracketing and Bisection}
CONST ntry = 50; {max. number of iterations}
	  factor = 1.6; {multiplication factor for enlarging interval}

VAR
	ResJmax : myReal;
	it : INTEGER;
	conv, success : BOOLEAN;

BEGIN
	{start with a guess for Vmin,max around Vguess}
	Vmin:=Vguess - 20*par.TolVint;
	Vmax:=Vguess + 20*par.TolVint;
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

	IF NOT success THEN Stop_Prog('Could not find bracketing values for Vint at time ' + FloatToStrF(new.tijd, ffGeneral,10,0));
END;

BEGIN {main program}
	Print_Welcome_Message(ProgName, version); {Welcomes the user, and shows the name and version of the program and the authors}

    {if '-h' or '-H' option is given then display some help and exit:}
    IF hasCLoption('-h') THEN Display_Help_Exit(ProgName);
    IF hasCLoption('-tidy') THEN Tidy_Up_Parameter_File(TRUE); {tidy up file and exit}
    IF NOT Correct_Version_Parameter_File(ProgName, version) THEN Stop_Prog('Version of ZimT and '+parameter_file+' do not match.');
    
{Initialisation:}
    Read_Parameters(MsgStr, par, stv, ProgName); {Read parameters from input file}
    Check_Parameters(stv, par, ProgName); {perform a number of chekcs on the paramters. Note: we need Vt}
    Prepare_Log_File(log, MsgStr, par, version); {open log file}
    IF par.AutoTidy THEN Tidy_Up_Parameter_File(FALSE); {clean up file but don't exit!}
    
    Make_Grid(stv.h, stv.x, stv.i1, stv.i2, par); {Initialize the grid}
    Define_Layers(stv, par); {define layers: Note, stv are not CONSTREF as we need to change them}
	Init_Trap_Distribution(log, stv, par); {Places all types of traps (bulk and interface) in the device at places determined by define_layers.}
    Init_nt0_And_pt0(stv, par); {inits nt0b and pt0 arrays needed for SRH recombination}

	Open_and_Read_tVG_file(inv, new, par); {open tVG file, read header and first time/voltage/Gehp}

	WITH new DO Init_Pot_Dens_Ions_Traps(V, Vgn, Vgp, n, p, nion, pion, f_tb, f_ti, Vint, stv, par); {init. (generalised) potentials and densities}
	new.UpdateIons:=TRUE; {in ZimT this is always true as we don't artificially fix the ions like we can in SimSS}
	{just to make sure these are initialised!}
	curr:=new;
	prev:=curr;

	Prepare_tJV_File(uitv, par.tJ_file, TRUE);   {create the tJV-file}
	IF par.StoreVarFile THEN Prepare_Var_File(stv, par, TRUE); {create a new var_file with appropriate heading}

	WRITELN('The calculation has started, please wait.');

    {Init all parameters and start solving for steady-state (tijd=0)}
	Init_Generation_Profile(stv, log, par); {init. the stv.orgGm array. This is the SHAPE of the profile}
	Update_Generation_Profile(stv.orgGm, new.Gm, new.Gehp, stv, par); {update current Gm array to correct value}
	CountNotConv:=0; {counts the number of times the transient Poisson/cont. eq. solver failed to converge}
	CounttVGPoints:=0; {count number of transient (t,V,G) points}
	countStatic:=0;
	conv:=FALSE;
	keepGoing:=foundtVG; {=true at this point!}

	WHILE keepGoing DO
	BEGIN {main loop over times and t,V,G from input file:}
		Update_Generation_Profile(stv.orgGm, new.Gm, new.Gehp, stv, par); {update current Gm array to correct value}
		prev:=curr;
		
		IF new.SimType = 1 
		THEN Main_Solver(curr, new, MainIt, conv, StatusStr, stv, par) {easy: Vint = Vext, Vext is in input tVG_file}
		ELSE BEGIN {SimType > 1: Voc or Rseries <> 0: in both cases, we don't know Vint}
			{first find Vmn and Vmx, voltages that bracket the device voltage new.Vint}
			Bracket_device_voltage(Vmn, Vmx, curr.Vint, ResVmn, curr, new, stv, par);
			{then use bisection to find Vint}
			MaxItVint:=ROUND(LN((Vmx-Vmn)/par.TolVint)/LN(2)); {max. number of required bisection steps}
			FOR ItVint:=0 TO MaxItVint DO {bisection loop}
			BEGIN
				new.Vint:=0.5*(Vmx + Vmn); {new guess}
				ResVint:=Residual_Current_Voltage(new.Vint, curr, new, conv, stv, par);
				{note, this time in ResCurr we use ResetVnp=false, as we want to keep the new values of V,n,p,nion,pion and Vgn, Vgp}
				IF ResVint*ResVmn < 0 THEN Vmx:=new.Vint ELSE Vmn:=new.Vint; {do bisection}
			END;

			new.Vint:=0.5*(Vmn+Vmx);
		END; {SimType > 1}

		IF conv 
		THEN acceptNewSolution:=TRUE
		ELSE BEGIN {o dear, now what?}
			INC(CountNotConv);
			{put error messages in log file:}
			WRITELN(log);
			WRITELN(log, 'Messages from main solver at time = ',FloatToStrF(new.tijd, ffGeneral,10,0),':');
			WRITELN(log, StatusStr);
			FLUSH(log);
			{now assess whether we accept the new solution, or skip it, or quit:}
			CASE par.FailureMode OF
				0 : Stop_Prog('Convergence failed at time = ' + FloatToStrF(new.tijd, ffGeneral,10,0)+ '. Maybe try smaller time steps.');
				1 : acceptNewSolution:=TRUE;
				2 : acceptNewSolution:=(new.dti=0) {is true if steady-state, false otherwise}
			END;
		END;

		IF acceptNewSolution {OK, new solution is good (even if conv might be FALSE)}
		THEN BEGIN
			curr:=new;
			{output:}
			Write_To_tJV_File(uitv, curr, prev, stv, par, TRUE);
		END
		ELSE BEGIN 
			WRITELN(log, 'Skipping (t,Vext,Gehp) point at time ',FloatToStrF(new.tijd, ffGeneral,10,0)); 
			FLUSH(log) 
		END;

		IF CounttVGPoints MOD par.OutputRatio = 0 THEN {output to screen and var_file every OutputRatio timesteps}
        BEGIN
			IF NOT Conv THEN TextColor(LightRed) ELSE TextColor(LightGray); {Reset to default font colour: it's not white, it's light grey!}
			WITH new DO WRITELN('Time:',tijd:12,' Vext: ',Vext:6:4,' Gehp:',Gehp:11,' Jext:',Jext  :8:3,' convIndex: ',convIndex);
			IF par.StoreVarFile AND acceptNewSolution THEN Write_Variables_To_File(curr, stv, par, TRUE)
        END;
	
		INC(CounttVGPoints);
		
		{is the state (Va, Gehp, ...) still changing significantly?}
		{count the number of points in time where the system hardly changes:}
		IF acceptNewSolution AND (curr.dti*ABS(prev.Vint-curr.Vint)<1) AND (prev.Gehp=curr.Gehp) AND (curr.dti*ABS(prev.Jint-curr.Jint)<1)
			THEN INC(CountStatic)
			ELSE CountStatic:=0; {reset counter}

		staticSystem:=(CountStatic >= MinCountStatic) AND (par.Autostop);
		
		{now see if there's more in the tVG_file:}
		Read_tVG(new, curr.tijd, inv, foundtVG); {try to read a line of t, Va, Gehp}
		IF staticSystem AND par.AutoStop AND foundtVG THEN {we're going to stop prematurely}
		BEGIN
			WRITELN(log, 'Stopping ZimT as the system does not change anymore.');
			TextColor(LightRed); {switch to different colour}
			WRITELN('Stopping ZimT as the system does not change anymore.');
			TextColor(LightGray); {switch back to normal}
		END; 
		keepGoing:=foundtVG AND NOT staticSystem;
    END; {main loop}

	TextColor(LightGray); {whatever happened, switch back to normal colour}
    CLOSE(uitv);

	WRITELN(CounttVGPoints,' (t,V,G)-points');
	IF CountNotConv > 0 THEN WRITELN(CountNotConv,' did not converge.');

    WRITE('Output written in file(s) ', par.tj_file);
	IF par.StoreVarFile THEN WRITE(' and ',par.Var_file);
	WRITELN('.');
	WRITELN(log,CountNotConv,' did not converge.');
	CLOSE(log);

    WRITELN('Finished, press enter to exit');
    IF par.Pause_at_end = 1 THEN READLN {pause at the end of the program}

END.
