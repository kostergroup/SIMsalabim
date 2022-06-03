PROGRAM SimSS;

{
SIMsalabim:a 1D drift-diffusion simulator 
Copyright (c) 2020, 2021, 2022 Dr T.S. Sherkar, Dr V.M. Le Corre, M. Koopmans,
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

{$UNITPATH ../Units/} {first tell compiler where our own units are located}


USES {our own, generic ones:}
	 TypesAndConstants,
     InputOutputUtils, 
     NumericalUtils,
     {our drift-diffusion stuff:}
	 DDTypesAndConstants,
	 DDRoutines,
	 {and normal units:}
	 Math,  
	 StrUtils,
	 SysUtils; 

CONST
    ProgName = TProgram.SimSS;  
    version = '4.35';   {version, 1.00 = 10-03-2004}

{first: check if the compiler is new enough, otherwise we can't check the version of the code}
{$IF FPC_FULLVERSION < 30200} {30200 is 3.2.0}
	{$STOP FPC VERSION SHOULD BE AT LEAST 3.2.0}
{$ENDIF}
{now check to see if the versions of the units match that of this code:}
{$IF (DDRoutinesVersion <> version) OR (DDTypesAndConstantsVersion <> version)} 
	{$STOP Wrong version of the DD units!}
{$ENDIF}


VAR 
	curr, new : TState; {store the previous point in time, the current one and the new}
	{curr: solved, stored, done, 1 time step ago
	new: to be solved, latest time that was read}

	stv : TStaticVars; {all variables that are calculated at the start of the simulation and then remain constant}

	par : TInputParameters; {all input parameters}

	VCount, MainIt : INTEGER; 

	quit_Voc, Conv_Main, acceptNewSolution : BOOLEAN;

	uitv , log : TEXT; {the output files}

	JVSim, JVExp : TJVList; {stores the current-voltage characteristics}
    
	MsgStr : ANSISTRING = ''; {Ansistrings have no length limit, init string to ''}
	
	StatusStr : ANSISTRING; 

FUNCTION Applied_voltage(Vcount : INTEGER; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters) : myReal;
{computes the applied voltage based on Vcount (=i)}
VAR i : INTEGER;

    FUNCTION LogarithmicV(i : INTEGER) : myReal;
    {computes a logarithmic V distribution, stepsize becomes zero at Vacc, so
    Vacc should lie outside [Vmin, Vmax] }
    VAR d : myReal;
    BEGIN
		d:=par.Vacc-par.Vmax;
		LogarithmicV:=par.Vacc - d*EXP((1-i/(stv.NJV-1))*LN((par.Vacc-par.Vmin)/d))
    END;

BEGIN

    CASE Sign(par.Vscan) OF
      -1 : i:=stv.NJV+1-Vcount; {scan down}
      1 : i:=Vcount {scan up}
    END;

    CASE par.Vdistribution of
      1 : Applied_voltage:=par.Vmin + par.Vstep*(i-1);
      2 : Applied_voltage:=LogarithmicV(i-1);
    END;
END;

PROCEDURE Read_Experimental_JV(VAR JVExp : TJVList; VAR log : TEXT; VAR stv : TStaticVars; VAR par : TInputParameters);
{reads the experiment JV curve and determines Vmin, Vmax, Vstep, and Vscan}
{note: This overrides the Vmin,max,step,scan,distribution in the parameter file AND the command line.}
{we could have used a linked list instead...}
{NOTE: stv and par are passed as VAR parameters as we need to set stv.NJV, par.Vmin and par.Vmax}
VAR i : INTEGER;
    inp : TEXT;
    dumstr : STRING;
    V, J : ARRAY[1..maxExpData] OF myReal;
BEGIN
    IF NOT FileExists(par.ExpJV) THEN Stop_Prog('Could not find file '+par.ExpJV, FALSE);
    ASSIGN(inp, par.ExpJV);
    RESET(inp);
    
    {first try to read a header}
    TRY 
		READLN(inp, dumstr);
    EXCEPT
		Stop_Prog('Cannot read a header from '+par.ExpJV+', is it empty?', FALSE);
    END;

    i:=0;
    WHILE (NOT EOF(inp)) AND (i < maxExpData) DO
    BEGIN
        INC(i);
        TRY
            READLN(inp, V[i], J[i]);
        EXCEPT
            Stop_Prog('The experimental JV curve in '+par.ExpJV+' can only contain voltage and current density, etc.', FALSE);
        END
    END;

    CLOSE(inp);

    {now check a number of things:}
    stv.NJV:=i; {this overrides NJV calculated from Vmin,Vmax,Vstep, or from direct input}
    IF stv.NJV < minExpData THEN Stop_Prog('Not enough experimental data points.', FALSE);
    IF stv.NJV = maxExpData THEN Stop_Prog('It looks like the experimental JV file contains more points than I can handle.', FALSE);

    {now we can copy the arrays V and J into ExpJV}
    SetLength(JVExp, stv.NJV+1); {the regular voltages start at 1}
    FOR i:=1 TO stv.NJV DO {note, we only use part of the JVList record}
    BEGIN
        JVExp[i].Vint:=V[i];
        JVExp[i].Vext:=V[i];
        JVExp[i].Jint:=J[i];
        JVExp[i].Jext:=J[i];
        JVExp[i].Use:=TRUE
    END;

    {Vmin, Vmax: we will need these later on}
    par.Vmin:=MIN(JVexp[1].Vint, JVexp[stv.NJV].Vint); {we use min, max function as Vscan might be -1}
    par.Vmax:=MAX(JVexp[1].Vint, JVexp[stv.NJV].Vint);

    {now document what happened and write to log file:}
    WRITELN(log);
    WRITELN(log, 'Read experiment JV curve from ', par.ExpJV,'.');
    WRITELN(log, 'This overrides the voltage distribution that is in ',parameter_file);
    WRITELN(log, 'and any voltage parameters passed via the command line.');
    WRITELN(log, 'Vmin: ',par.Vmin:6:4,' Vmax: ',par.Vmax:6:4);
    WRITELN(log);
    WRITELN('Read experimental JV curve from ', par.ExpJV,'.');
END;

PROCEDURE Init_Voltages_and_Tasks(VAR JVSim, JVExp : TJVList; VAR VCount : INTEGER; VAR log : TEXT;
				VAR stv : TStaticVars; VAR par : TInputParameters);
{determines which voltages need to be simulated. This is either based on 
the normal parameters Vdistribution, Vmin, Vmax, Vstep, Vacc, NJV, or
based on experimental input. 
NOTE: stv and par are passed as VAR parameters as Read_Experimental_JV needs to 
change them!}
VAR i, k : INTEGER;
BEGIN
    IF par.UseExpData THEN Read_Experimental_JV(JVExp, log, stv, par); {calling this overrides NJV calculated from Vmin,Vmax,Vstep, or from direct input}
			
    SetLength(JVSim, stv.NJV+1); {if there is a pre-bias, we need an extra data point, so we make the array one longer than NJV}
    IF par.PreCond THEN {the pre-bias will be stored in point 0}
	WITH JVSim[0] DO {pre-bias point}
	BEGIN
	    Vint:=par.Vpre;
	    UpdateIons:=TRUE; {yes, ions are moving. We checked in Read_Parameters that if PreCond then CIM <>0}
	    Store:=FALSE; {however, don't store the pre-bias point}
	END;
    
    k:=ORD(NOT par.PreCond); {so if PreCond=true then k=0, we will use this below}
    {now for the rest of the voltages:}
    FOR i:=1 TO stv.NJV DO {the regular voltages start at 1}
    BEGIN
		WITH JVSim[i] DO
		BEGIN
			IF par.UseExpData
				THEN Vint:=JVExp[i].Vint
				ELSE Vint:=Applied_voltage(i, stv, par);
			IF par.ion_red_rate=0 THEN {means ions are fixed}
            BEGIN
				UpdateIons:=(i=1) AND (NOT par.PreCond);
				{so if the ions are fixed, but we're doing the first point (i=1)
				then update ions unless we already did a preconditioning}
				Store:=TRUE {by default, store this JV point. If this point does not converge, we will set Store to FALSE}
			END
			ELSE BEGIN {so ions are not always fixed}
				UpdateIons:=((par.CNI<>0) OR (par.CPI<>0)) AND ( (i-k) MOD par.ion_red_rate = 0 );
				{note the k! whether we update the ions also depends on the pre-conditioning}
				Store:=((par.CNI=0) AND (par.CPI=0)) OR ( (i-k+1) MOD par.ion_red_rate = 0)
			END;
		END;
    END;
    
    {VCount is used in the main program as the index of the voltages that should be simulated}
    IF par.PreCond {if we need to pre condition, then use index 0}
		THEN VCount:=0 {Counter of the number of voltages which have been computed}
		ELSE VCount:=1; {index VCount=0 is reserved for the pre-bias voltage (if any)}
    
    {now write to log to show what we are going to do:}
    WRITELN(log);
    WRITELN(log,'The following voltages will be simulated:');
    WRITELN(log,'  i       V    UpdateIons   Store');
    FOR i:=ORD(NOT par.PreCond) TO stv.NJV DO {note: if not PreCond we start at i=1}
	WITH JVSim[i] DO
	    WRITELN(log,i:3,'  ',Vint:8:3,'    ',UpdateIons:5,'      ',Store);
    WRITELN(log);
    FLUSH(log);	
END;

PROCEDURE Init_States(VAR curr, new : TState; Vapp : myReal; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters); 
{inits new (time, Va, Gehp, V, n, p, etc) and sets curr:=new}
BEGIN
    WITH new DO 
    BEGIN
		tijd:=0; {tijd: Dutch for time. SIMsalabim not (yet) time dependent, so set to zero.}
		dti:=0; {dti: inverse of delta t, so if dti=0, the time-step is infinite => steady-state!}
		SimType:=1; {1 means that this is steady-state, not open-circuit}
		Vint:=Vapp; {set applied voltage}
		Gehp:=par.Gehp*par.Gfrac; {set Gehp equal to the value in parameter file.}
		Update_Generation_Profile(stv.orgGm, Gm, par.Gehp*par.Gfrac, stv, par); {Set current Gm array to correct value}
		Init_Pot_Dens_Ions_Traps(V, Vgn, Vgp, n, p, nion, pion, f_tb, f_ti, Vint, stv, par); {init. (generalised) potentials and densities}
	END;
	
	curr:=new;{just to make sure curr is initialised!}
END;


PROCEDURE Find_Solar_Cell_Parameters(JVChar : TJVList; VAR SCPar : TSCPar);
{Finds Jsc, Voc, MPP (and thus FF) by interpolating the V-J data points}
VAR Nusable, i, k, i_start, InterpolationOrder : INTEGER;
    c1, c2, c3 : myReal;
    x_dat, y_dat : Row;
    IntSuccess : BOOLEAN;
BEGIN
    {init the SCPar:}
    WITH SCPar DO {if myReal=extended, then un-init. vars have value 'NaN'.}
    BEGIN
		Jsc:=0; Vmpp:=0; MPP:=0; FF:=0; Voc:=0; 
		ErrJsc:=0; ErrVmpp:=0; ErrMPP:=0; ErrFF:=0; ErrVoc:=0;
		calcSC:=FALSE; calcMPP:=FALSE; calcFF:=FALSE; calcOC:=FALSE;
    END;
    
    {first count how many JV points are usable}
    Nusable:=0;  
    FOR i:=1 TO LENGTH(JVChar)-1 DO 
		IF JVChar[i].Use THEN INC(Nusable);
    {we start at i=1 as the 0th point (if any) corresponds to a pre-bias}
    {note that the length of JVChar includes index 0, so we have to stop at i=LENGHT()-1}
  
    {now create 2 arrays for the J and V points:}
    SetLength(x_dat, Nusable+1); {note: such arrays always start at index 0}
    SetLength(y_dat, Nusable+1); {but we will take them to start at 1.}
    
    {now copy the usable points into x_dat and y_dat:}
    k:=1; {index couter for x_dat, y_dat arrays}
    FOR i:=1 TO LENGTH(JVChar)-1 DO {we don't need i=0 as it corresponds to Vpre (if any)}
		WITH JVChar[i] DO
			IF Use THEN 
			BEGIN
				x_dat[k]:=Vext;
				y_dat[k]:=Jext;
				INC(k)
			END;
	
	InterpolationOrder:=MIN(MaxInterpolationOrder, Nusable-1);

    IF Nusable >= 2 THEN {we need at least 2 points}
    BEGIN
		{first calculate at short-circuit}
		IntSuccess:=InterExtraPolation(x_dat, y_dat, 0, SCPar.Jsc, SCPar.ErrJsc, InterpolationOrder, 2);
		{note: we allow some extrapolation for Jsc (ExtraIndex=1)}
		SCPar.calcSC:=IntSuccess AND (SCPar.ErrJsc < threshold_err*ABS(SCPar.Jsc));

		{now calculate Voc by swapping x_dat and y_dat:}
		IntSuccess:=InterExtraPolation(y_dat, x_dat, 0, SCPar.Voc, SCPar.ErrVoc, InterpolationOrder, 0);
		SCPar.calcOC:=IntSuccess AND (SCPar.ErrVoc < threshold_err*ABS(SCPar.Voc));

		{sometimes (e.g. if until_Voc=1) then it's better to try linear interpolation:}
		IF (NOT SCPar.calcOC) AND (InterpolationOrder>1) THEN
		BEGIN {redo with linear interpolation:}
			IntSuccess:=InterExtraPolation(y_dat, x_dat, 0, SCPar.Voc, SCPar.ErrVoc, 1, 1);
			{now that we use linear interpolation, we allow for some extrapolation as well: ExtraIndex=1}
			SCPar.calcOC:=IntSuccess AND (SCPar.ErrVoc < threshold_err*ABS(SCPar.Voc));
		END;

		{Calculate the Maximum Power Point (MPP) and the fill-factor (FF)}
		{copy the power into y_dat:}
		FOR i:=0 TO Nusable-1 DO y_dat[i]:=-x_dat[i]*y_dat[i];
		i_start:=0;
		FOR i:=1 TO Nusable-3 DO
			IF y_dat[i] > y_dat[i_start] THEN i_start:=i;
		{i_start = i where power is maximum}
		SCPar.calcMPP:=FALSE;
		SCPar.calcFF:=FALSE;
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
			IF (SCPar.Voc<>0) AND (SCPar.Jsc<>0)  AND SCPar.calcOC AND SCPar.calcSC THEN
			BEGIN
				SCPar.FF:=-SCPar.MPP/(SCPar.Voc*SCPar.Jsc); { the fill-factor}
				SCPar.calcFF:=(SCPar.FF<1) AND (SCPar.FF>0); {make sure the value of FF makes sense. If not, set calcFF to FALSE}
				{now calc. the error in FF based on the propagation of errors:}
				SCPar.ErrFF:=SQRT( SQR(SCPar.ErrMPP/(SCPar.Jsc*SCPar.Voc)) + SQR(SCPar.ErrJsc*SCPar.FF/SCPar.Jsc) + SQR(SCPar.ErrVoc*SCPar.FF/SCPar.Voc));
			END
		END;
    END
END;


PROCEDURE Calc_and_Output_Solar_Cell_Parameters(JVExp, JVSim : TJVList; CONSTREF par : TInputParameters);
VAR SCParExp, SCParSim : TSCPar; {Experimental and simulated solar cell parameters}
    digits : INTEGER;
    uitv : TEXT;
    str : STRING;
BEGIN
{note: we use 3 constants (tab1,2,3) to control the layout. They are defined in unit DDTypesAndConstants}
    WRITELN;
    str:='Simulated';
    str:=AddChar(' ',str,tab1+Length(str));
    {add white space: AddCharR adds character ' ' to str until the length is tabpos:}
    str:=AddCharR(' ',str,tab2);
    IF par.UseExpData THEN
    BEGIN
		str:=str + 'Experimental';
		str:=AddCharR(' ',str,tab3); 
		str:=str + 'Deviation';
    END;
    WRITELN(str);
    
    {note: we will use Vext, Jext in calculating the Voc, Jsc, etc. as this is what is
    in the external circuit}
    Find_Solar_Cell_Parameters(JVSim, SCParSim);

    IF par.UseExpData 
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
	IF calcOC AND calcSC AND calcFF AND calcMPP THEN
	BEGIN
	    ASSIGN(uitv, par.scPars_file);  {create the scPars_file}
	    REWRITE(uitv);
	    WRITELN(uitv,'Jsc ErrJsc Voc ErrVoc Vmpp ErrVmpp MPP ErrMPP FF ErrFF');
	    WRITELN(uitv,Jsc:6:4,' ',ErrJsc:6:4,' ',Voc:6:4,' ',ErrVoc:6:4,' ',Vmpp:6:4,' ',ErrVmpp:6:4,' ',MPP:6:4,' ',ErrMPP:6:4,' ',FF:6:4,' ',ErrFF:6:4);
	    CLOSE(uitv);
	END
    ELSE WRITELN('Could not determine (all) solar cell parameters.');
END;


PROCEDURE Compare_Exp_Sim_JV(JVExp, JVSim : TJVList; CONSTREF stv : TStaticVars; CONSTREF par : TInputParameters);
VAR rms, Jmin, Jmax : myReal;
    i, j, count : INTEGER;
    disgardedPoints, bracketed : BOOLEAN;
BEGIN
    {If there is series resistance we need to interpolate the experimental data to match the voltages of the simulation.}
    IF (par.Rseries > 0) THEN
    BEGIN
		j:=1;
		FOR i:=1 TO stv.NJV DO
		BEGIN
			{Some simuation voltages are outside the range of the experiment data, so we can not use these}
			IF (JVSim[i].Vext < par.Vmin) OR (JVSim[i].Vext > par.Vmax) THEN JVExp[i].Use:= FALSE
			ELSE
			BEGIN
				{Set the voltage and interpolate the corresponding experimental current density.}
				JVExp[i].Vext:=JVSim[i].Vext;
				bracketed:=(JVExp[j].Vint < JVExp[i].Vext) AND (JVExp[j+1].Vint > JVExp[i].Vext);
				WHILE (NOT bracketed) AND (j < stv.NJV-1) DO
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
    
    FOR i:=1 TO stv.NJV DO
		{If the simulation did not converge at we can not use it. Similarly if we could not interpolate a point in 
		the experimental data we can not use this voltage}
		IF JVSim[i].Use AND JVExp[i].Use
		THEN 
		BEGIN
			{Calculate the sum of squared residuals}
			CASE par.rms_mode OF
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
    

    IF (par.rms_threshold <= count/stv.NJV) AND (Jmax-Jmin>tolReal) THEN {we need at least 1 data point and a measurable interval [Jmin, Jmax]}
    BEGIN
		{Here we calculate the root mean square error and normalise with respect to the interval [Jmin,Jmax].}
		CASE par.rms_mode OF
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
        IF par.rms_mode = logarithmic THEN WRITELN('-too many voltages were sim. and exp. currents have different signs. Try using rms_mode = lin.');
        WRITELN
    END
END;

BEGIN {main program}
	Print_Welcome_Message(ProgName, version);

    {if '-h' or '-H' option is given then display some help and exit:}
    IF hasCLoption('-h') THEN Display_Help_Exit(ProgName);
    IF hasCLoption('-tidy') THEN Tidy_Up_Parameter_File(TRUE); {tidy up file and exit}
    IF NOT Correct_Version_Parameter_File(ProgName, version) THEN Stop_Prog('Version of SIMsalabim and '+parameter_file+' do not match.');

    {Initialisation:}
    Read_Parameters(MsgStr, par, stv, ProgName); {Read parameters from input file}
    Check_Parameters(stv, par, ProgName); {perform a number of chekcs on the paramters. Note: we need Vt}
    Prepare_Log_File(log, MsgStr, par, version); {open log file}
    IF par.AutoTidy THEN Tidy_Up_Parameter_File(FALSE); {clean up file but don't exit!}

    Make_Grid(stv.h, stv.x, stv.i1, stv.i2, par); {Initialize the grid}
    Define_Layers(stv, par); {define layers: Note, stv are not CONSTREF as we need to change them}
	Init_Trap_Distribution(log, stv, par); {Places all types of traps (bulk and interface) in the device at places determined by define_layers.}
    Init_nt0_And_pt0(stv, par); {inits nt0b and pt0 arrays needed for SRH recombination}
	
    Init_Voltages_and_Tasks(JVSim, JVExp, VCount, log, stv, par);
    Init_Generation_Profile(stv, log, par); {init. the stv.orgGm array. This is the SHAPE of the profile}
	Init_States(curr, new, JVSim[VCount].Vint, stv, par); {inits new (time, Va, Gehp, V, n, p, etc) and sets curr:=new}

    Prepare_tJV_File(uitv, par.JV_file, FALSE);   {create the JV-file}
    IF par.StoreVarFile THEN Prepare_Var_File(stv, par, FALSE); {Create a new var_file with appropriate heading if required}

    quit_Voc:=FALSE;

    WHILE (Vcount<=stv.NJV) AND (NOT quit_Voc) DO  {loop over voltages}
    BEGIN
		new.Vint:=JVSim[VCount].Vint; 
		new.UpdateIons:=JVSim[VCount].UpdateIons;
		
		{now use Main_Solver to iterative solve the Poisson eq and continuity equations:}
		Main_Solver(curr, new, MainIt, Conv_Main, StatusStr, stv, par);

		{First copy the result into JVSim:}		
		IF JVSim[VCount].Store THEN {we store this even if Conv_Main is false}
			WITH JVSim[VCount] DO
			BEGIN
				Jint:=new.Jint; {the internal current}
				Jext:=new.Jext; {internal current}
				Vext:=new.Vext; {external voltage}
				Vint:=new.Vint; {internal voltage}
				Use:=Conv_Main; {store whether convergerence was achieved}
			END;
	
		{now see what else we need to do with the solution}
		IF Conv_Main THEN
		BEGIN
			acceptNewSolution:=TRUE; {new is a keeper, but we'll take care of that later on}
			{output:}
			WRITELN('At Vint=', new.Vint:4:3,' converged in',MainIt:4,' loop(s), Jint=',new.Jint:7:4,' convIndex: ',new.convIndex);

			{we only write the point to the JV_file if Conv_Main and if Store:}
			IF JVSim[VCount].Store THEN Write_To_tJV_File(uitv, new, curr, stv, par, FALSE);						
		END
		ELSE BEGIN {o dear, now what?}
			{tell user of failed convergence:}
			WRITELN('No convergence for Vint=',new.Vint:4:3);
			{put error messages in log file:}
			WRITELN(log);
			WRITELN(log, 'Messages from main solver at voltage= ',new.Vint:4:3,':');
			WRITELN(log, StatusStr);
			FLUSH(log);
			{now assess whether we accept the new solution, or skip it, or quit:}
			CASE par.FailureMode OF
				0 : Stop_Prog('Convergence failed at voltage = ' + FloatToStrF(new.Vint, ffGeneral,5,0)+ '. Maybe try smaller voltage steps?');
				1 : acceptNewSolution:=TRUE;
				2 : acceptNewSolution:=(new.dti=0) {is true if steady-state, false otherwise}
			END;
			
		END;
		
		IF acceptNewSolution {OK, new solution is good (even if conv_main might be FALSE)}
			THEN curr:=new {we keep the newest solution}
			ELSE WRITELN(log, 'Skipping (J,V) point at voltage ',new.Vint); 
		  
		{should we store the internal variables?}
		IF par.StoreVarFile THEN {FIRST check if it should be stored as we cannot compute (Vcount MOD OutputRatio) if OutputRatio=0!}
			IF (Vcount MOD par.OutputRatio = 0) AND acceptNewSolution 
				THEN Write_Variables_To_File(curr, stv, par, FALSE);	
			
        quit_Voc:=(JVSim[VCount].Jext>0) AND par.until_Voc AND (par.Gehp<>0); {only stop past Voc if there is light!}

        VCount:=VCount + 1;
    END; {loop over voltages}

    IF (par.Gehp * par.Gfrac <> 0) AND (stv.V0 <> stv.VL) THEN  {we do have a solar cell}
        Calc_and_Output_Solar_Cell_Parameters(JVExp, JVSim, par);

    IF par.UseExpData THEN Compare_Exp_Sim_JV(JVExp, JVSim, stv, par);
    {note: if Rseries <> 0 then Vext in JVExp and JVSim will be different so we can't do a direct comparison}

    WRITELN('The JV characteristic is written in ', par.JV_file,'.');
	IF par.OutputRatio>0 THEN WRITELN('Stored internal variables in ', par.Var_file,'.');

	WRITELN('Finished, press enter to exit');
    IF par.Pause_at_end = 1 THEN READLN; {pause at the end of the program}
    CLOSE(log);
END.
