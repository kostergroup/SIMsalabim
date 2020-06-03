unit TypesAndConstants;
{provides the types and constants used by SIMsalabim}

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

{$MODE OBJFPC} {force OBJFPC mode}

interface

const
    q = 1.6022e-19;  	{C} {elementary charge}
    k = 1.3807e-23;     {J/K} {Boltzmann's constant}
    eps_0 = 8.8542e-12; {F/m} {dielectric constant of vacuum}
    Bc = 2.8;           {min. number of bonds for percolation to occur}
    Max_NP = 10000;     {max number of grid points except contacts}
    maxExpData = 10000; {max number of experimental JV points}
	minExpData = 4;     {min number of experimental JV points}
	tolReal = 1e-6;     {tolerance: when are 2 floats equal?}
	
type myReal = DOUBLE; 
	 vector = ARRAY[0..Max_NP + 1] OF myReal;
     MathFunc = FUNCTION(x : myReal) : myReal;
     Row = ARRAY OF myReal; {used in find_solar_cell_parameters and Interpolation, and for n_vals, p_vals, Fn_vals, and Fp_vals}
     Table = ARRAY OF ARRAY OF myReal; {used to store mob_tab, table with elec. mob. as a function of F and n}
     JVList = ARRAY OF RECORD {for storing the JV curve}
                        Vint, Vext : myReal;
                        Jint, Jext : myReal;
                        UpdateIons, Store, Use : BOOLEAN
                      END;
	 TSCPar = RECORD {for storing solar cell parameters}
				Jsc, Vmpp, MPP, FF, Voc, ErrJsc, ErrVmpp, ErrMPP, ErrFF, ErrVoc : myReal;
				calcSC, calcMPP, calcFF, calcOC : BOOLEAN
			  END;
     Tfitmode = (linear, logarithmic);

implementation

begin 

end.
