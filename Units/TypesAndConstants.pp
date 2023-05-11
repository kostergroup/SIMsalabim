unit TypesAndConstants;
{provides the types and constants used by SIMsalabim}

{
SIMsalabim: a 1D drift-diffusion simulator 
Copyright (c) 2020, 2021, 2022, 2023 Dr T.S. Sherkar, V.M. Le Corre, Dr M. Koopmans,
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

{$MODE OBJFPC} {force OBJFPC mode}

interface

uses Ucomplex;

const
    Max_NP	 = 1000;     {max number of grid points except contacts}
	Max_NEtr = 20;	     {max number of trap levels}

	{exit codes:}
	EC_Warning = 3;
	EC_DevParCorrupt = 90; 
	EC_InvalidInput = 91;
	EC_InvalidCLInput = 92;
	EC_NumericalFailure = 93;
	EC_ConverenceFailedHalt = 94;
	EC_ConverenceFailedNotHalt = 95;
	EC_FileNotFound = 96;
	EC_TimeOut = 97;
	EC_ProgrammingError = 99;
	

type myReal = EXTENDED; 
	{note: you can put myReal = single, double, or extended. However, extended may not be available in which case the compiler
	 will simply take double. The size of the real type is saved in the log file, so you can check}
	 vector			= ARRAY[0..Max_NP + 1] OF myReal;
	 TrapArray		= ARRAY[0..Max_NP + 1, 1..Max_NEtr] OF myReal;
	 TrapEnArray	= ARRAY[1..Max_NEtr] OF myReal;
     ShortIntVector	= ARRAY[0..Max_NP + 1] OF ShortInt;
     Row			= ARRAY OF myReal; 
     Table			= ARRAY OF ARRAY OF myReal; {used to store mob_tab, table with elec. mob. as a function of F and n}
     MathFunc		= FUNCTION(x : myReal) : myReal;
     MathFuncValues = FUNCTION(x : myReal; vals : Row) : myReal;
     intArray 		= ARRAY OF ShortInt;
     StringArray 	= ARRAY OF STRING;
     ComplexMatrix 	= ARRAY OF ARRAY OF COMPLEX;

implementation

begin 

end.
