unit NumericalUtils;

{
SIMsalabim: a 1D drift-diffusion simulator 
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

FUNCTION Norm(vec : vector; istart, ifinish : integer) : myReal;
VAR i : INTEGER;
    ans : myReal;
{computes the infinity norm of vector vec}

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

PROCEDURE Interpolation_rec(x_pts, y_pts : Row; n_pts : INTEGER; x : myReal; VAR y_est, dy : myReal);
    {Given two arrays, x_ptx and y_pts, of length n_pts and a value x at which to interpolate, this
    procedure returns an estimate of y at x. The error estimate is the difference between the y estimate
    using a polynomial of order n_pts - 2 and n_pts - 1.}

PROCEDURE Interpolation(x_array, y_array : Row; x_target : myReal; VAR y_estimate, y_error_estimate : myReal; interp_order : INTEGER);
    {Given two arrays, x_ptx and y_pts, of length n_pts and a value x at which to interpolate, this
    procedure returns an estimate of y at x. The error estimate is the difference between the y estimate
    using a polynomial of order n_pts - 2 and n_pts - 1.}

PROCEDURE Fit_Parabola(x1, x2, x3, y1, y2, y3 : myReal; VAR a, b, c : myReal);
{Fits a parabola through [(x1, y1), (x2, y2), (x3,y3)], y = ax^2 + bx + c}

FUNCTION BilinearInterpolation(x1, x2 : myReal; VAR x1a, x2a : row; m, n : INTEGER; VAR ya : Table; extra_x2 : BOOLEAN) : myReal;
{Uses bilinear Interpolation for interpolation in 2D. See Numer. Recipes in Pascal 1st Ed. p. 107, paragraph 3.6}
{ data is given as ya[j,k] = y(x1a[j], x2a[k]). We want to know ya at x1, x2}
{ j=1,..,m, k=1,..,n}
{ x1a and x2a need to be ascending!}

{use an approximation to the Bernoulli function, much (~ 3x) faster!}
FUNCTION B(x : myReal) : myReal; {The Bernoulli function, an approximation}
{we use several approximations to the Bernoulli function:
 close to 0: 1st order, then 2nd order, then 4th order. Finally, for large (pos.) x B=0 and for
 x << 0 we have B=-x}


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
        ELSE Stop_Prog('Romberg integration did not converge, increase MaxRombIt or decrease TolRomb.', wait_if_faulty)
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
    RombergIntegrationValues:=0; {define return value of qRomb}
    SetLength(Rom, 2, MaxIt);
    {This sets the length of Rom as array[0..1, 0..MaxRombIt-1] of myReal}
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
        ELSE Stop_Prog('Romberg integration did not converge, increase MaxRombIt or decrease TolRomb.', wait_if_faulty)
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


PROCEDURE Interpolation_rec(x_pts, y_pts : Row; n_pts : INTEGER; x : myReal; VAR y_est, dy : myReal);
    {Given two arrays, x_ptx and y_pts, of length n_pts and a value x at which to interpolate, this
    procedure returns an estimate of y at x. The error estimate is the difference between the y estimate
    using a polynomial of order n_pts - 2 and n_pts - 1.}
VAR 
    i_nearest, idx : INTEGER;
    d_c, d_d, CD_prtl, dx_old, dx_new : myReal;
    C_cor, D_cor : ARRAY OF myReal;

    PROCEDURE estimate_y(i_col : INTEGER);
    {This calculates one column of the Neville table, increasing in order of the polynomials
    with increasing column index.}
    VAR 
        i_row : INTEGER;
    BEGIN
        {calculate the corrections to our estimate of y for every row.}
        FOR i_row := 0 TO n_pts-i_col-1 DO 
        BEGIN 
            d_c := x_pts[i_row] - x; 
            d_d := x_pts[i_row+i_col] - x;
            IF d_c = d_d THEN BEGIN 
                Stop_Prog('Error in interpolation routine.');
            END;
            CD_prtl := (C_cor[i_row+1] - D_cor[i_row]) / (d_c - d_d);  
            C_cor[i_row] := d_c * CD_prtl;
            D_cor[i_row] := d_d * CD_prtl;
        END;
        IF 2*(i_nearest+1) < n_pts-i_col THEN 
            dy := C_cor[i_nearest+1]
        ELSE 
        BEGIN 
            dy := D_cor[i_nearest];
            i_nearest := i_nearest - 1
        END;

        y_est := y_est + dy;
        IF i_col < n_pts-1 THEN 
        BEGIN 
            estimate_y(i_col+1);
        END;
    END;
{start of interpolation procedure.}
BEGIN 
    SETLENGTH(C_cor, n_pts);
    SETLENGTH(D_cor, n_pts);
    i_nearest := 0;
    dx_old := ABS(x-x_pts[0]);
    FOR idx := 0 TO n_pts-1 DO 
    BEGIN
        dx_new := ABS(x-x_pts[idx]);
        IF dx_new < dx_old THEN 
        BEGIN 
            i_nearest := idx;
            dx_old := dx_new
        END; 
        C_cor[idx] := y_pts[idx];
        D_cor[idx] := y_pts[idx]
    END;
    y_est := y_pts[i_nearest];
    i_nearest := i_nearest - 1;

    {now estimate y for increasing polynomyal order, starting at order 1.}
    estimate_y(1);
END;


PROCEDURE Interpolation(x_array, y_array : Row; x_target : myReal; VAR y_estimate, y_error_estimate : myReal; interp_order : INTEGER);
VAR 
    closest_sum, sum_pts : myReal;
    x_selected, y_selected : ARRAY OF myReal;
    len_x_arr, len_y_arr, idx, idx_sum, idx_closest : INTEGER;

    PROCEDURE Neville_interpolation(x_pts, y_pts : Row; n_pts : INTEGER; x_target : myReal; VAR y_est, dy : myReal);
        {Given two arrays, x_ptx and y_pts, of length n_pts and a value x at which to interpolate, this
        procedure returns an estimate of y at x. The error estimate is the difference between the y estimate
        using a polynomial of order n_pts - 2 and n_pts - 1.}
    VAR 
        i_nearest, idx : INTEGER;
        d_c, d_d, CD_prtl, dx_old, dx_new : myReal;
        C_cor, D_cor : ARRAY OF myReal;

        PROCEDURE estimate_y(i_col : INTEGER);
        {This calculates one column of the Neville table, increasing in order of the polynomials
        with increasing column index.}
        VAR 
            i_row : INTEGER;
        BEGIN
            {calculate the corrections to our estimate of y for every row.}
            FOR i_row := 0 TO n_pts-i_col-1 DO 
            BEGIN 
                d_c := x_pts[i_row] - x_target; 
                d_d := x_pts[i_row+i_col] - x_target;
                IF d_c = d_d THEN BEGIN 
                    Stop_Prog('Error in interpolation routine, two subsequent x points are of the same value.');
                END;
                CD_prtl := (C_cor[i_row+1] - D_cor[i_row]) / (d_c - d_d);  
                C_cor[i_row] := d_c * CD_prtl;
                D_cor[i_row] := d_d * CD_prtl;
            END;
            IF 2*(i_nearest+1) < n_pts-i_col THEN 
                dy := C_cor[i_nearest+1]
            ELSE 
            BEGIN 
                dy := D_cor[i_nearest];
                i_nearest := i_nearest - 1
            END;

            y_est := y_est + dy; {this is the new estimate of y, the error will be dy, the last step size}
            IF i_col < n_pts-1 THEN 
            BEGIN 
                estimate_y(i_col+1);
            END;
        END;
    
    BEGIN {start of interpolation procedure.}
        SETLENGTH(C_cor, n_pts);
        SETLENGTH(D_cor, n_pts);
        i_nearest := 0;
        dx_old := ABS(x_target-x_pts[0]);

        FOR idx := 0 TO n_pts-1 DO 
        BEGIN
            dx_new := ABS(x_target-x_pts[idx]);
            IF dx_new < dx_old THEN 
            BEGIN 
                i_nearest := idx;
                dx_old := dx_new
            END; 
            C_cor[idx] := y_pts[idx];
            D_cor[idx] := y_pts[idx]
        END;
        y_est := y_pts[i_nearest];
        i_nearest := i_nearest - 1;
        
        {now estimate y for increasing polynomyal order, starting at order 1.}
        estimate_y(1);
    END;

BEGIN {start preprocessing the arrays to interpolate to have the number of points matching the interpolation order}
    SETLENGTH(x_selected, interp_order + 1);
    SETLENGTH(y_selected, interp_order + 1);

    len_x_arr := length(x_array);
    len_y_arr := length(y_array);

    IF len_x_arr <> len_y_arr THEN Stop_Prog('Interpolation: arrays not of equal length.');
    IF len_x_arr <= interp_order THEN Stop_Prog('Interpolation: more points required for requested order of interpolation.');
    closest_sum := 0; {we need to initialize because of a compiler warning, this serves no purpose.}
    FOR idx := 0 TO length(x_array) - interp_order - 1 DO {find closest sequence of points to x_target}
    BEGIN 
        sum_pts := 0;
        FOR idx_sum := 0 TO interp_order DO 
            sum_pts := sum_pts + x_array[idx+idx_sum];

        IF (idx = 0) OR (ABS(sum_pts - interp_order * x_target) < ABS(sum_pts - closest_sum)) THEN 
        BEGIN 
            idx_closest := idx;
            closest_sum := sum_pts
        END;
    END;

    FOR idx := 0 TO interp_order DO {add the selected points from the arrays to the interpolation point arrays.}
    BEGIN
        x_selected[idx] := x_array[idx_closest+idx];
        y_selected[idx] := y_array[idx_closest+idx];
    END;

    {Use Neville interpolation on closest points to x_target}
    Neville_interpolation(x_selected, y_selected, interp_order + 1, x_target, y_estimate, y_error_estimate);
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

    IF (m<2) OR (n<2) THEN Stop_Prog('Not enough points to do BilinearInterpolation');

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


{use an approximation to the Bernoulli function, much (~ 3x) faster!}
FUNCTION B(x : myReal) : myReal; {The Bernoulli function, an approximation}
{we use several approximations to the Bernoulli function:
 close to 0: 1st order, then 2nd order, then 4th order. Finally, for large (pos.) x B=0 and for
 x << 0 we have B=-x}
CONST C1 = 0.083333333333; {1/12}
	  C2 = 1.388888889e-3; {1/720}
	  C3 = 0.25; {cut-off}
	  C4 = 1.2; {another cut-off}
	  C5 = 3.9; {3rd cut-off}
VAR absx : myReal;
BEGIN
    absx:=ABS(x);
    IF absx < C3 {if x close to zero (=most common case!)}
		THEN B:=1 - 0.5*x {then use 1st order Taylor}
		ELSE
		BEGIN {absx >= C3}
			IF absx < C4 THEN B:=1 - x*(0.5 - C1 * x) {use 2nd order}
			ELSE
			BEGIN {absx >= C4}
				IF absx < C5 THEN B:=1 + x*(-0.5 + x*(C1 - x*x*C2)) {use 4th order Taylor}
				ELSE BEGIN
					IF x < -C5 
						THEN B:=-x {x < -C5}
						ELSE B:=x/(EXP(x)-1); {x > C5, will be a very small value}
				END; {absx >= C5}
			END {absx >= C4}
		END; {absx >= C3}
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
