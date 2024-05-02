unit NumericalUtils;

{
SIMsalabim: a 1D drift-diffusion simulator 
Copyright (c) 2020, 2021, 2022, 2023, 2024, S. Heester, Dr T.S. Sherkar, V.M. Le Corre, Dr M. Koopmans,
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
email: l.j.a.koster@rug.nl
surface mail: 
L.J.A. Koster
Zernike Institute for Advanced Materials
Nijenborgh 4, 9747 AG Groningen, the Netherlands
}

interface

{$MODE OBJFPC} {force OBJFPC mode}

uses TypesAndConstants,
	 InputOutputUtils, {stop_prog} 
	 Math, {power}
	 SysUtils, {for reading the command line}
	 StrUtils; {for DelSpace}

FUNCTION RombergIntegration(f : MathFunc; a, b, TolRomb : myReal; MaxIt : INTEGER; wait_if_faulty : BOOLEAN) : myReal;
{Romberg integration of function f from a to b. Rom[1,1] and Rom[2,1] are trapezoidal
approximations, Rom[1,j] and Rom[2,j] are Richardson extrapolations.}

FUNCTION RombergIntegrationValues(f : MathFuncValues; vals : Row; a, b, TolRomb : myReal; MaxIt : INTEGER; wait_if_faulty : BOOLEAN) : myReal;
{Romberg integration of function f from a to b. Rom[1,1] and Rom[2,1] are trapezoidal
approximations, Rom[1,j] and Rom[2,j] are Richardson extrapolations.}
{This version can take an array (vals) which is passed on to the MathFunc f}

procedure SolveQuadraticEq(a, b, c : myReal; var nr : integer; var x1, x2 : myReal);
{solves the equation a*x*x + b*x + c = 0, and returns roots x1 and x2}
{nr : number of roots}

FUNCTION Norm_Inf(vec : vector; istart, ifinish : integer) : myReal;
{computes the infinity norm of vector vec}

FUNCTION Norm_Eucl(vec : vector; istart, ifinish : integer) : myReal;
{computes the Euclidean norm of vector vec}

FUNCTION Calc_RMS_Error(y : vector; avy : myReal; istart, ifinish : INTEGER) : myReal;
{computes the rms error of y given average y (avy)}

FUNCTION sf(r : myReal) : INTEGER;
{step function}

FUNCTION Difference(a, b : vector; istart, ifinish : integer) : vector;
{computes the vector a-b}

PROCEDURE Tridiag(VAR x : vector; a, b, c, r : vector; i0, i1 : INTEGER);
{Solves the system A x=r, with A a tridiagonal matrix, a: lower diagonal,
b: main diagonal, and c: upper diagonal.
The total system looks like:
(b[i0] c[i0]                   )  (x[0])      (r[i0])
(a[i0+1] b[i-+1] c[i0+1]       )    .            .
(        ......                )    .       =    .
(               c[i1-1]        )    .            .
(                    a[i1]b[i1]) (x[i1])      (r[i1])    }

PROCEDURE Neville(x0 : myReal; VAR y, err : myReal; n : INTEGER; x_pts, y_pts : Row);
{Neville interpolation. x0: x of desired point. y, err: estimated y and its error, n: order of polynomial}

FUNCTION InterExtraPolation(x, y : Row; x0 : myReal; VAR y_estimate, err_estimate : myReal; Order : INTEGER; 
							ExtraIndex : INTEGER = 0) : BOOLEAN; 
{interpolation / extrapolation routine. This is a wrapper for Neville's algorithm}
{x,y: are the original data, first index 1. Should be either ascending or descending.}
{x0: desired x, y_estimate, err_estimate: estimated y and its error. Order: order of polynomial}
{ExtraIndex: 0=> no extrapolation, 1=> extrapolation limited to average x-spacing, 2=> full extrapolation}
{Returns TRUE if successful, FALSE if not}

PROCEDURE Fit_Parabola(x1, x2, x3, y1, y2, y3 : myReal; VAR a, b, c : myReal);
{Fits a parabola through [(x1, y1), (x2, y2), (x3,y3)], y = ax^2 + bx + c}

FUNCTION BilinearInterpolation(x1, x2 : myReal; VAR x1a, x2a : row; m, n : INTEGER; VAR ya : Table; extra_x2 : BOOLEAN) : myReal;
{Uses bilinear Interpolation for interpolation in 2D. See Numer. Recipes in Pascal 1st Ed. p. 107, paragraph 3.6}
{ data is given as ya[j,k] = y(x1a[j], x2a[k]). We want to know ya at x1, x2}
{ j=1,..,m, k=1,..,n}
{ x1a and x2a need to be ascending!}

FUNCTION Bessel(x : myReal) : myReal;
{calculates the Bessel function of order 1, with some scaling}


implementation

function RandomGaussian : myReal;
{Generates a Gaussianly distributed random number, with mean = 0, SD = 1}
var Phi, Rho, X, dum : myReal ;
begin
     Phi := 2.0*Pi*Random ;
     repeat
		dum:=random
	 until dum <> 0; {guard against dum=0}
     Rho := sqrt(-2.0*Ln(dum)) ;
     X := cos(Phi)*Rho ; {Y := sin(Phi)*Rho ; }
     { A polar to cartesian coordinate transformation is taking place -
     observe that here both X and Y will have normal Gaussian distribution }
     RandomGaussian := X { or Y }
end ;

FUNCTION RombergIntegration(f : MathFunc; a, b, TolRomb : myReal; MaxIt : INTEGER; wait_if_faulty : BOOLEAN) : myReal;
{Romberg integration of function f from a to b. Rom[1,1] and Rom[2,1] are trapezoidal
approximations, Rom[1,j] and Rom[2,j] are Richardson extrapolations.}
VAR i, j, k, m : INTEGER;
    l, h, sum, diffnew, diffold : myReal;
    Rom : ARRAY OF ARRAY OF myReal;
BEGIN
    RombergIntegration:=0; {define return value of qRomb}
    SetLength(Rom, 2, MaxIt);
    {This sets the length of Rom as array[0..1, 0..MaxRombIt-1] of myReal}
    h:=b-a;
    Rom[0,0] := (f(a) + f(b))/2 * h;
    i:=1;
    diffnew:=1; {difference between two successive answers determines convergence}
    diffold:=1;
{   To guard against the possibility that two consecutive row elements agree with
    each other but not with the value of the integral being approximated, it is
    common to genarate approximations until not only |Rn-1,n-1 - Rn,n| (diffnew)
    is within the tolerance, but also |Rn-2,n-2 - Rn-1,n-1| (diffold). }
    WHILE ((diffnew > TolRomb) OR (diffold > TolRomb)) AND (i < MaxIt) DO
    BEGIN
        diffold:=diffnew;
        i:=i+1;
{       approximation from Trapezoidal method                   }
        sum:=0;
        m:=ROUND(POWER(2, i-2));
        FOR k := 1 TO m DO sum:=sum + f(a + (k - 0.5 ) * h);
        Rom[1,0]:=0.5*(Rom[0,0] + h * sum);
{       extrapolation                                           }
        FOR j:= 2 TO i DO
        BEGIN
            l:=POWER(4, j-1);
            Rom[1,j-1] := Rom[1,j-2]+(Rom[1,j-2]-Rom[0,j-2])/(l-1)
        END;
        h:=0.5*h;
{       Only use relative tol. when Rom[1, i-1] <> 0 }
        IF Rom[1, i-1] = 0
            THEN diffnew:=ABS(Rom[1, i-1]-Rom[1, i-2])
            ELSE diffnew:=ABS((Rom[1, i-1]-Rom[0, i-2])/Rom[1, i-1]);
{       since only two rows are kept in storage, this step      }
{       is to prepare for the next row.                         }
{       update row 1 of R                                       }
        FOR j:=0 TO i-1 DO Rom[0,j] := Rom[1,j]
    END;
    IF (diffnew < TolRomb) AND (diffold < TolRomb) {did we converge?}
        THEN RombergIntegration:=Rom[1, i-1]
        ELSE Stop_Prog('Romberg integration did not converge, increase MaxRombIt or decrease TolRomb.',EC_NumericalFailure, wait_if_faulty)
        {no convergence, then stop routine}
END;



FUNCTION RombergIntegrationValues(f : MathFuncValues; vals : Row; a, b, TolRomb : myReal; MaxIt : INTEGER; wait_if_faulty : BOOLEAN) : myReal;
{Romberg integration of function f from a to b. Rom[1,1] and Rom[2,1] are trapezoidal
approximations, Rom[1,j] and Rom[2,j] are Richardson extrapolations.}
{This version can take an array (vals) which is passed on to the MathFunc f}
VAR i, j, k, m : INTEGER;
    l, h, sum, diffnew, diffold : myReal;
    Rom : ARRAY OF ARRAY OF myReal;
BEGIN
    RombergIntegrationValues:=0; {define return value of RombergIntegrationValues}
    SetLength(Rom, 2, MaxIt);
    {This sets the length of Rom as array[0..1, 0..MaxIt-1] of myReal}

    h:=b-a;
    Rom[0,0] := (f(a, vals) + f(b, vals))/2 * h;
    i:=1;
    diffnew:=1; {difference between two successive answers determines convergence}
    diffold:=1;

{   To guard against the possibility that two consecutive row elements agree with
    each other but not with the value of the integral being approximated, it is
    common to genarate approximations until not only |Rn-1,n-1 - Rn,n| (diffnew)
    is within the tolerance, but also |Rn-2,n-2 - Rn-1,n-1| (diffold). }

    WHILE ((diffnew > TolRomb) OR (diffold > TolRomb)) AND (i < MaxIt) DO
    BEGIN
        diffold:=diffnew;
        i:=i+1;
{       approximation from Trapezoidal method                   }
        sum:=0;
        m:=ROUND(POWER(2, i-2));
        FOR k := 1 TO m DO sum:=sum + f(a + (k - 0.5 ) * h, vals);
        Rom[1,0]:=0.5*(Rom[0,0] + h * sum);
{       extrapolation                                           }
        FOR j:= 2 TO i DO
        BEGIN
            l:=POWER(4, j-1);
            Rom[1,j-1] := Rom[1,j-2]+(Rom[1,j-2]-Rom[0,j-2])/(l-1)
        END;
        h:=0.5*h;
{       Only use relative tol. when Rom[1, i-1] <> 0 }
        IF Rom[1, i-1] = 0
            THEN diffnew:=ABS(Rom[1, i-1]-Rom[1, i-2])
            ELSE diffnew:=ABS((Rom[1, i-1]-Rom[0, i-2])/Rom[1, i-1]);
{       since only two rows are kept in storage, this step      }
{       is to prepare for the next row.                         }
{       update row 1 of R                                       }
        FOR j:=0 TO i-1 DO Rom[0,j] := Rom[1,j]
    END;

    IF (diffnew < TolRomb) AND (diffold < TolRomb) {did we converge?}
        THEN RombergIntegrationValues:=Rom[1, i-1]
        ELSE Stop_Prog('Romberg integration did not converge, increase MaxRombIt or decrease TolRomb.', EC_NumericalFailure, wait_if_faulty)
        {no convergence, then stop routine}
        
END;



procedure SolveQuadraticEq(a, b, c : myReal; var nr : integer; var x1, x2 : myReal);
{solves the equation a*x*x + b*x + c = 0, and returns roots x1 and x2}
{note: a should not be zero.}
{nr : number of roots}
var q, det : myReal;
    sgn : shortInt;
begin
     det:=b*b-4*a*c;
     if det < 0 then nr:=0
     else
     begin {det >=0}
          if det = 0 then {1 root}
          begin
               nr:=1;
               x1:=-b/(2*a);
               x2:=x1
          end
          else begin {two roots}
               nr:=2;
               if b > 0 then sgn:=1 else sgn:=-1;
               q:=-0.5*(b + sgn*Sqrt(det));
               x1:=q/a;
               x2:=c/q;
          end; {two roots}
     end;{det >=0}
end;

FUNCTION Norm_Inf(vec : vector; istart, ifinish : integer) : myReal;
VAR i : INTEGER;
    ans : myReal;
{computes the infinity norm of vector vec}
BEGIN
    ans:=0;
    FOR i:=istart TO ifinish DO
        IF ABS(vec[i]) > ans THEN ans:=ABS(vec[i]);
    Norm_Inf:=ans
END;

FUNCTION Norm_Eucl(vec : vector; istart, ifinish : integer) : myReal;
VAR i : INTEGER;
    ans : myReal;
{computes the Euclidean norm of vector vec}
BEGIN
    ans:=0;
    FOR i:=istart TO ifinish DO
        ans:=ans + vec[i]*vec[i];
    Norm_Eucl:=SQRT(ans)
END;

FUNCTION Calc_RMS_Error(y : vector; avy : myReal; istart, ifinish : INTEGER) : myReal;
{computes the rms error of y given average y (avy)}
VAR i : INTEGER;
    sum : myReal;
BEGIN
    IF ifinish <= istart THEN Stop_Prog('Invalid ifinish and istart passed to Calc_RMS_Error.',EC_ProgrammingError);
    sum:=0;
    FOR i:=istart TO ifinish DO
	sum:=sum + SQR(y[i]-avy);
    Calc_RMS_Error:=SQRT(sum/(ifinish-istart));
END;


FUNCTION sf(r : myReal) : INTEGER;
{step function}
BEGIN
	IF r>= 0 THEN sf:=1 ELSE sf:=0;
END;


FUNCTION Difference(a, b : vector; istart, ifinish : integer) : vector;
{computes the vector a-b}
VAR i : INTEGER;
BEGIN
    FOR i:=istart TO ifinish DO Difference[i]:=a[i]-b[i]
END;

PROCEDURE Tridiag(VAR x : vector; a, b, c, r : vector; i0, i1 : INTEGER);
{Solves the system A x=r, with A a tridiagonal matrix, a: lower diagonal,
b: main diagonal, and c: upper diagonal.
The total system looks like:
(b[i0] c[i0]                   )  (x[0])      (r[i0])
(a[i0+1] b[i-+1] c[i0+1]       )    .            .
(        ......                )    .       =    .
(               		c[i1-1])    .            .
(                    a[i1]b[i1]) (x[i1])      (r[i1])    }
{IMPORTANT: THIS FORM OF THE ALGORITM MODIFIES THE ORIGINAL COEFICIENTS, SO 
WE NEED TO PASS THEM AS VALUE (NOT VAR) PARAMETERS}
{This is based on the tridiagonal matrix algorithm, also known as the Thomas algorithm}
VAR w : myReal;
    i : INTEGER ;
BEGIN
	FOR i:=i0+1 TO i1 DO
	BEGIN
		w:=a[i]/b[i-1];
		b[i]:=b[i] - w*c[i-1];
		r[i]:=r[i] - w*r[i-1];
	END;
	x[i1]:=r[i1]/b[i1];
	FOR i:=i1-1 DOWNTO i0 DO
		x[i]:=(r[i]-c[i]*x[i+1])/b[i]
END;


PROCEDURE Neville(x0 : myReal; VAR y, err : myReal; n : INTEGER; x_pts, y_pts : Row);
{Neville interpolation. x0: x of desired point. y, err: estimated y and its error, n: order of polynomial}
VAR i, j : INTEGER;
	p : ARRAY OF ARRAY OF myReal;
BEGIN
	{first check if inputs make sense:}
	IF Length(x_pts) <> Length(y_pts) THEN Stop_Prog('Error in Neville routine: x_pts and y_pts not of same length.', EC_ProgrammingError);
	IF n>=Length(x_pts) THEN Stop_Prog('Error in Neville routine: n cannot be equal or larger than number of points.', EC_ProgrammingError);
	FOR i:=1 TO Length(x_pts)-2 DO
		IF SameValue(x_pts[i],x_pts[i+1]) THEN Stop_Prog('Error in Neville routine: 2 consecutive x-values are equal.', EC_ProgrammingError);

	SetLength(p, n+1, n+1); {our matrix for storing the intermediate results}
	{copy original y-values into the first column of the matrix:}
	FOR i:=1 TO n DO
		p[i,1]:=y_pts[i];

	{now use Neville's formula to calc the polynomial approximations:}
	FOR i:=2 TO n DO
		FOR j:=i TO n DO
			p[j,i]:=((x0 - x_pts[j-i+1]) * p[j,i-1] - (x0 - x_pts[j]) * p[j-1,i-1]) / (x_pts[j] - x_pts[j-i+1]);
			
	y:=p[n,n]; {our final result}
	err:=MIN(ABS(p[n,n-1]-p[n,n]), ABS(p[n-1,n-1]-p[n,n]));	
END;

FUNCTION InterExtraPolation(x, y : Row; x0 : myReal; VAR y_estimate, err_estimate : myReal; Order : INTEGER; 
							ExtraIndex : INTEGER = 0) : BOOLEAN; 
{interpolation / extrapolation routine. This is a wrapper for Neville's algorithm}
{x,y: are the original data, first index 1. Should be either ascending or descending.}
{x0: desired x, y_estimate, err_estimate: estimated y and its error. Order: order of polynomial}
{ExtraIndex: 0=> no extrapolation, 1=> extrapolation limited to average x-spacing, 2=> full extrapolation}
{Returns TRUE if successful, FALSE if not}
VAR i, i1, N, istart, ifin : INTEGER;
	delx : myReal;
	x_selected, y_selected : Row; {selected points close to x0}
	x0Bracketed : BOOLEAN;
	

	FUNCTION Locate(x : Row; x0 : myReal; imin, imax : INTEGER) : INTEGER;
	{simple, slow, (hopefully) robust routine to an index within 
	row x (from imin to imax). x0 is some target we're looking for.
	The resulting index is such that x0 sits between x[i] and x[i+1]
	if x0 is outside the interval determined by x[imin], x[imax] then it 
	returns imin-1 or imax+1.}
	VAR i, r : INTEGER;
	BEGIN
		{check if x is descending or ascending:}
		r:=SIGN(x[imax]-x[imin]); {so r=1 if ascending, r=-1 if descending}
		IF r=0 THEN Stop_Prog('Error in Locate: row x is constant?', EC_ProgrammingError);

		i:=imin; {first init i!}
		{Check if x0 sits beyond the interval, either left or right:}
		IF r*(x0-x[imax]) > 0 THEN i:=imax+1;
		IF r*(x0-x[imin]) < 0 THEN i:=imin-1;
		IF SameValue(x0, x[imax]) THEN i:=imax;
	
		{now: if i is still imin, then we have to look within the interval:}
		IF (i=imin) THEN 
			WHILE NOT ((r*(x[i]-x0)<=0) AND (r*(x[i+1]-x0)>0)) AND (i<imax-1) DO
				INC(i);

		Locate:=i
	END;
	

	
BEGIN
	InterExtraPolation:=FALSE;
	{just make sure these are initialised:}
	y_estimate:=0;
	err_estimate:=0;
	
	{check input}
	N:=LENGTH(x)-1; {max index in our arrays that run from 1...N}
	IF Order>=N THEN Stop_Prog('Error in InterExtraPol: Not enough points for interpolation order.', EC_ProgrammingError);
	IF Order<1 THEN Stop_Prog('Error in InterExtraPol: Order should be at least 1.', EC_ProgrammingError);

	{does the interval bracket x0?}
	x0Bracketed:=(x[1]-x0)*(x[N]-x0)<=0; {note: we include the 0 as x0 might be exactly equal to either point}

	IF (ExtraIndex>0) OR x0Bracketed THEN 
	{only enter this block if either we allow extrapolation or we don't need it:}
	BEGIN
		i1:=Locate(x, x0, 1, N);
		{NOTE: i1 might be 0 or N+1 if x0Bracketed is false}
		
		{by now, i1 might be outside the interval [1,..,N], so we map it back to 1,..,N:}
		i1:=MAX(1, i1);
		i1:=MIN(N, i1);
		
		{istart and ifin define the points that will be passed on to Neville's routine:}
		istart:=MAX(FLOOR(i1 - 0.5*Order), 1); {try to put istart close to i1-Order/2, but >= 1}
		ifin:=Order + istart;
		IF ifin>N THEN 
		BEGIN
			ifin:=N;
			istart:=ifin-Order			
		END;

		{now we should have istart, ifin within [1,..,N] and i1 as close as possible to the middle}
		{copy selected points into local arrays so we can pass them on to the Neville routine:}
		SetLength(x_selected, Order+2); {Order+2? array starts at 0, the last index should be Order+1}
		SetLength(y_selected, Order+2);
		FOR i:=istart TO ifin DO
		BEGIN
			x_selected[i-istart+1]:=x[i];
			y_selected[i-istart+1]:=y[i];
		END;
	
		{do the actual interpolation/extrapolation:}
		Neville(x0, y_estimate, err_estimate, Order+1, x_selected, y_selected);
	
		{OK, now check if extrapolation (if any!) was acceptable}
		delx:=MIN(ABS(x_selected[1]-x0), ABS(x_selected[Order+1]-x0));
		InterExtraPolation:=(ExtraIndex=2) OR x0Bracketed OR (delx <= ABS(x_selected[1]-x_selected[Order+1])/Order) 
		{so IF there was extrapolation, then we accept this if either ExtraIndex=2, or if 
		 the delx ("amount" of extrapolation) is limited to the average x-spacing in the selected points}
	END
END;

PROCEDURE Fit_Parabola(x1, x2, x3, y1, y2, y3 : myReal; VAR a, b, c : myReal);
{Fits a parabola through [(x1, y1), (x2, y2), (x3,y3)], y = ax^2 + bx + c}
BEGIN
    b:=( (SQR(x3)-SQR(x1))*(y2-y1) - (SQR(x2)-SQR(x1))*(y3-y1) )
        /( (x2-x1)*(SQR(x3)-SQR(x1)) -(x3-x1)*(SQR(x2)-SQR(x1)) );
    a:=( y3 - y1 - b*(x3-x1) )/(SQR(x3) - SQR(x1) );
    c:=y1 - a*SQR(x1) - b* x1
END;

FUNCTION BilinearInterpolation(x1, x2 : myReal; VAR x1a, x2a : row; m, n : INTEGER; VAR ya : Table; extra_x2 : BOOLEAN) : myReal;
{Uses bilinear Interpolation for interpolation in 2D.}
{ data is given as ya[j,k] = y(x1a[j], x2a[k]). We want to know ya at x1, x2}
{ j=1,..,m, k=1,..,n}
{ x1a and x2a need to be ascending!}

{For more information see:
Numer. Recipes in Pascal 1st Ed. p. 107, paragraph 3.6
or: 
https://en.wikipedia.org/wiki/Bilinear_interpolation}

VAR j, k : INTEGER;
    y1, y2, y3, y4, t, u : myReal;
BEGIN
    {first find the grid square in which the point (x1, x2) falls}
    {we're looking for j and k such that: x1,2a[j,k] <= x1,2 <= x1,2a[j,k+1]}

    IF (m<2) OR (n<2) THEN Stop_Prog('Not enough points to do BilinearInterpolation', EC_ProgrammingError);

    {in the following two while loops, use the if-statement (and remove the and (j<m-1) in while condition)
    to ensure actual interpolation, the routine as is will extrapolate if necessary!}
    j:=1;
    WHILE NOT ( (x1a[j] <= x1) AND (x1 <= x1a[j+1]) ) AND (j<m-1) DO
		INC(j);

    k:=1;
    WHILE NOT ( (x2a[k] <= x2) AND (x2 <=x2a[k+1]) ) AND (k<n-1) DO
		INC(k);

    {now define the four tabulated points that surround the desired interior point x1,x2}
    {however, there is no guarantee that we have been able to bracket the point.}
    y1:=ya[j,k];
    y2:=ya[j+1,k];
    y3:=ya[j+1,k+1];
    y4:=ya[j,k+1];

    t:=(x1-x1a[j])/(x1a[j+1]-x1a[j]);
    {limit u to 1 if we don't want extrapolation}
    IF NOT extra_x2 AND (x2 > x2a[k+1])
		THEN u:=1
        ELSE u:=(x2-x2a[k])/(x2a[k+1]-x2a[k]);

    BilinearInterpolation:=(1-t)*(1-u)*y1+t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4;
END;

FUNCTION Bessel(x : myReal) : myReal;
{calculates the Bessel function of order 1, with some scaling}
CONST a1 =1/3;
	  a2 = 1/6;
	  a3 = 1/10;
	  a4 = 1/15;
	  a5 = 1/21;
	  a6 = 1/28;
	  a7 = 1/36;
	  a8 = 1/45;
	  a9 = 1/55;
	  a10= 1/66;
BEGIN
   Bessel:=1 + x*(1+a1*x*(1+a2*x*(1+a3*x*(1+a4*x*(1+a5*x*(1+a6*x*(1+a7*x*(1+a8*x*(1+a9*x*(1+a10*x))))))))))
END;


begin

end.
