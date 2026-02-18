program Updater;

{
This program updates the simulation setup and layer files of the SIMsalabim project.
It updates such files from minVersion (=oldest version that this program can handle) to newer ones. 
 
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

{$UNITPATH ../../Units/} {first tell compiler where our own units are located}

uses {our own, generic ones:}
	 TypesAndConstants,
     InputOutputUtils, 
     {our drift-diffusion stuff:}
	 DDTypesAndConstants,
	 HelperUnit, {This contains some types and routine that we need in updating from 456 to 500}
	 {and normal units:} 
	 SysUtils, StrUtils, Math;

const minVersion = 419;
	  maxVersion = 526;
	
	  {define shorthand for DirectorySeparator (OS-dependent):}  
	  slash = DirectorySeparator;  

var SetupFile, ProgName : string; 
	version : integer;
	force : boolean = false;
	dumstr : string;


procedure DisplayHelpExit;
begin
	writeln;
	writeln('To use this program, please specify the name of the simulation setup.');
	writeln('For example: ');
{$IFDEF WINDOWS}
	writeln('Updater.exe simulation_setup.txt');
	writeln('If you would like to ignore warnings, etc (NOT recommended), then use -f:');
	writeln('Updater.exe -f simulation_setup.txt');
{$ELSE}
	writeln('./Updater simulation_setup.txt');
	writeln('If you would like to ignore warnings, etc (NOT recommended), then use -f:');
	writeln('./Updater -f simulation_setup.txt');
{$ENDIF}	
	writeln;
	Stop_Prog('''-h''     : displays this help message', 0)
end;

procedure init_from_command_line;
begin
	case ParamCount of
		1 : SetupFile:=ParamStr(1);
		2 : begin {the parameters must look like -f filename.txt}
				if ParamStr(1)='-f' {the first parameter MUST be -f}
					then force:=true
					else DisplayHelpExit; 
				{parameter 1 can only be -f, check for this}
				SetupFile:=ParamStr(2); {this must be the parameter file. We'll check later on to see if it exists}
			end;
		else DisplayHelpExit {input wasn't valid, so show help and exit}
	end	
end;

procedure GetProgNameVersion(var progname : string; var version : integer; fn : string);
{searches the file fn for the name of the progam and the version number.	
If there was a version number, then version>0
If we could find the name of the program, then progname <> ''}
const MaxCount = 1000;
var inp : text;
	count : integer;
	dumStr, line : string;
begin
	assign(inp, fn);
	reset(inp);
   
    version:=0; 
    progname:='';
    
    WHILE NOT EOF(inp) DO
    BEGIN
        READLN(inp, line); {read a line from the file}
        line:=LOWERCASE(line); {convert all characters to lower case to make search easier}
        {first: did we find the version number:}
        IF (POS('version', line) > 0) THEN
        BEGIN {OK, we have a line that contains the word 'version', now let's extract the version number:}
			{now we make use of the fact that version numbers are ALWAYS (assumed) of the form x.xx (where x=0..9)}
			count:=0;
			REPEAT 
				inc(count);
				dumStr:=FormatFloat('0.00', 0.01 * count);
				if pos(dumStr, line) > 0 then version:=count;			
			UNTIL (version>0) or (count > MaxCount);		
		END;
     
        {did we find the correct name of the program:}
        IF POS('simsalabim', line) > 0 THEN progname:='SIMsalabim';
        IF POS('simss', line) > 0 THEN progname:='SimSS';
        IF POS('zimt', line) > 0 THEN progname:='ZimT';
        {POS returns the index of Substr in S, if S contains Substr. In case Substr isn't found, 0 is returned.}
    END;

	close(inp);

end;

procedure UpdateTo420(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=419 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
			{in case we need to add/remove parameters, the code codes here:}
			if pos('ri=', LowerCase(DelSpace(line))) <> 0 then {simply skip ri entry}
				CopyLine:=false; 
			if pos('rf=', LowerCase(DelSpace(line))) <> 0 then {replace rf with accDens}
			begin
				CopyLine:=false;
				buffer:=buffer + 'accDens = 0.9 * accelation parameter for density solver, 0 < accDens < 2' + LineEnding
			end;
			
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo421(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=420 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo422(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=421 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo423(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=422 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo424(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=423 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo425(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=424 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;

			if pos('simsalabim', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				{Replace all occurances of 'SIMsalabim' in line with 'SimSS':}
				buffer:=buffer + ReplaceText(line, 'SIMsalabim', 'SimSS') + LineEnding
			end;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo426(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=425 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo427(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=426 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo428(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=427 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
			{in case we need to add/remove parameters, the code codes here:}
			if pos('etrap=', LowerCase(DelSpace(line))) <> 0 then {add new stuff with traps from file(s)}
			begin
				CopyLine:=false;
				{first: fix notation:}
				line:=ReplaceStr(line, 'Etrap', 'ETrap'); {first replace 'Etrap' with 'ETrap'}
				line:=ReplaceStr(line, 'ETrap', 'ETrapSingle'); {next: replace 'ETrap' with 'ETrapSingle'}
				buffer:=buffer + line + LineEnding;
				{next introduce 2 new parameters: }
				buffer:=buffer + 'BulkTrapFile = none * name of file with bulk trap energy profile (or ''none''). If specified, overrides ETrapSingle' + LineEnding;
				buffer:=buffer + 'IntTrapFile = none * name of file with interface trap energy profile (or ''none''). If specified, overrides ETrapSingle' + LineEnding
			end;
			
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo429(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=428 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo430(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=429 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo431(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=430 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo432(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=431 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo433(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=432 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo434(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=433 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo435(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=434 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
	
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			{in case we need to add/remove parameters, the code codes here:}
			if pos('accdens=', LowerCase(DelSpace(line))) <> 0 then {insert new parameter couplePC}
			begin
				buffer:=buffer + 'couplePC = 4 * >= 0, coupling between Poisson equation and continuity equations' + LineEnding; 
				CopyLine:=true;
			end;						
			
			{now remove TolRomb, MaxRombIt, LowerLimBraun, and UpperLimBraun, and their comments:}
			if pos('tolromb=', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			if pos('*startingpoint', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			if pos('maxrombit=', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			if pos('lowerlimbraun=', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			if pos('*shouldbenon-zero(', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			if pos('upperlimbraun=', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			if pos('*onlyrelevantiftherm', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			if pos('*10seemssufficientfor', LowerCase(DelSpace(line))) <> 0 then CopyLine:=false;
			
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo436(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=435 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo437(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=436 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo438(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=437 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo439(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=438 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo440(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=439 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo441(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=440 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo442(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=441 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo443(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=442 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo444(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=443 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo445(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=444 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo446(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=445 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
				
			{in case we need to add/remove parameters, the code codes here:}
			if pos('outputratio=', LowerCase(DelSpace(line))) <> 0 then {insert new parameter LimitOutput}
			begin
				buffer:=buffer + 'LimitDigits = 1 * if 1, then number of digits in output is limited' + LineEnding; 
				CopyLine:=true;
			end;	
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;	
			

procedure UpdateTo447(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=446 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
						
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo448(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=447 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			{add new 'optics' block:}
			if pos('p_0=', LowerCase(DelSpace(line))) = 1 then
			begin
				buffer:=buffer + line + LineEnding + LineEnding; {first copy line with p_0 and add extra empty line}
				buffer:=buffer + '**Optics****************************************************************************' + LineEnding;
				buffer:=buffer + 'L_TCO = 110E-9 * m, thickness of the TCO. Set to 0 if layer is not used' + LineEnding;
				buffer:=buffer + 'L_BE = 200e-9 * m, thickness of back electrode, must be >0' + LineEnding;
				buffer:=buffer + 'nk_substrate = ..'+slash+'Data'+slash+'nk_SiO2.txt * name of file with n,k values of substrate' + LineEnding;
				buffer:=buffer + 'nk_TCO	= ..'+slash+'Data'+slash+'nk_ITO.txt * name of file with n,k values of TCO' + LineEnding;
				buffer:=buffer + 'nk_active = ..'+slash+'Data'+slash+'nk_P3HTPCBM_BHJ.txt * name of file with n,k values of active layer' + LineEnding;
				buffer:=buffer + 'nk_BE = ..'+slash+'Data'+slash+'nk_Al.txt  * name of file with n,k values of back electrode' + LineEnding;
				buffer:=buffer + 'spectrum = ..'+slash+'Data'+slash+'AM15G.txt * name of file that contains the spectrum' + LineEnding;
				buffer:=buffer + 'lambda_min = 350E-9 * m, lower bound wavelength' + LineEnding;
				buffer:=buffer + 'lambda_max = 800E-9  * m, upper bound wavelength' + LineEnding;
				CopyLine:=false
			end;

			{add two new parameters to Transport layers block:}
			if pos('nc_ltl', LowerCase(DelSpace(line))) <> 0 then
			begin
				buffer:=buffer + 'nk_LTL = ..'+slash+'Data'+slash+'nk_PEDOT.txt * name of file with n,k values of the left TL' + LineEnding;
				buffer:=buffer + 'nk_RTL = ..'+slash+'Data'+slash+'nk_Ca.txt          * name of file with n,k values of the right TL' + LineEnding;
				CopyLine:=true
			end;
		
			
			{change parameter TLsAbsorb to TLsGen:}
			if pos('tlsabsorb=', LowerCase(DelSpace(line))) <> 0 then {insert new parameter LimitOutput}
			begin
				buffer:=buffer + 'TLsGen = 0 * TLs generate electrons/hole yes(1)/no(0), overrides profile!' + LineEnding;
				CopyLine:=false
			end;	
			
			{make sure simss/zimt is correct in parameter autotidy:}
			if pos('autotidy=', LowerCase(DelSpace(line))) <> 0 then 
			begin
				{in older versions, the comment line for autotidy refers to SimSS, also in ZimT!}
				{so, simply replace with 'the program'}
				line:=StringReplace(line, 'SimSS', 'the program', [rfReplaceAll, rfIgnoreCase]);
				line:=StringReplace(line, 'ZimT', 'the program', [rfReplaceAll, rfIgnoreCase]);
				CopyLine:=true
			end;	
	
			{change comment for Gen_profile}
			if pos('gen_profile=', LowerCase(DelSpace(line))) <> 0 then 
			begin
				line:=LeftStr(line, Pos('*',line)); {remove comment}
				line:=line + ' name of file generation profile (or ''none'' or ''calc'')'; {add new comment}
				CopyLine:=true
			end;
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo449(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=448 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			{replace tolJ with tolDens}
			if pos('tolj=', LowerCase(DelSpace(line))) = 1 then
			begin
				buffer:=buffer + 'tolDens = 1e-8 * relative tolerance of density solver' + LineEnding;
				CopyLine:=false
			end;

			{remove MinRelChange}
			if pos('minrelchange=', LowerCase(DelSpace(line))) = 1 then
				CopyLine:=false;
						
			{remove MinAbsJDark}
			if pos('minabsjdark=', LowerCase(DelSpace(line))) = 1 then
				CopyLine:=false;		

			{replace accDens with minAcc and maxAcc}
			if pos('accdens=', LowerCase(DelSpace(line))) = 1 then
			begin
				buffer:=buffer + 'minAcc = 0.05 * >0, min. acceleration parameter' + LineEnding;
				buffer:=buffer + 'maxAcc = 0.95 * <2, max. acceleration parameter' + LineEnding;
				CopyLine:=false
			end;
		
			{now set grad to 4, as this is best!}
			if pos('grad=', LowerCase(DelSpace(line))) = 1 then
			begin
				buffer:=buffer + 'grad = 4 * determines shape of exp. grid, increase grad for smaller h[1]' + LineEnding;
				CopyLine:=false
			end;
					
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo450(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=449 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			{Set CurrDiffInt to 2 as that is the prefered mode as of v4.49:}
			if (pos('currdiffint=', LowerCase(DelSpace(line))) <> 0) then 
			begin
				CopyLine:=false;
				buffer:=buffer + 'CurrDiffInt = 2 * Calc. current from differential (1) or integral (2) expression' + LineEnding;
			end;
			
			{in SimSS, the by-line of Var_file is not quite correct, so we change it:}
			if (pos('var_file=', LowerCase(DelSpace(line))) <> 0) and (ProgName = 'SimSS') then 
			begin
				line:=LeftStr(line, Pos('*', line)); {remove the bit after the *}
				line:=line + ' name of the file with internal variables (x,V,n,p,Jn,Jp,etc.)';
				CopyLine:=true
			end;
					
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo451(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=450 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			{insert new parameter timeout before Pause_at_end:}
			if (pos('pause_at_end=', LowerCase(DelSpace(line))) <> 0) then 
			begin
				buffer:=buffer + 'timeout = -1 * s, max run time, use negative value for unlimited run time.' + LineEnding;
				CopyLine:=true
			end;
					
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo452(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=451 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo453(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=452 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo454(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=453 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo455(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=454 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo456(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=455 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo457(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=456 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo458(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=457 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo459(ProgName : string; var version : integer; fn : string);
{updates the parameter file (fn) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	buffer : ansistring;
	CopyLine : boolean;
begin
	if version=458 then {we only update 1 version number!}
	begin
		assign(inp, fn);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
	
			{in case we need to add/remove parameters, the code codes here:}

			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, fn);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo500(ProgName : string; var version : integer; var fn : string);
{updates the parameter file (fn) to version 5.00. This is different from all previous (and future?) updates as we
need to split the single parameter file into multiple files.}
const setupFilev5 = 'simulation_setup.txt';
	  layerFileName = '_parameters.txt';
	  newVersionString = '5.00';
				
var LTLexists, RTLexists : boolean;
	mainLyr, i : integer;
	
	par : TInputParameters456;
	parnew : TInputParameters500;
	path, msg : AnsiString;

begin

	path:=ExtractFilePath(fn);

	if abs(version - 458)<=1 then {we only update 1 version number!}
	begin
		Read_Parameters_456(fn, msg, par, ProgName); {use 456 reader to read the parameters}
		if msg<>'' then Stop_Prog('Something is wrong.', 3, not force); {is msg <> '', then the reader found things like values from the command line. That's not OK}

		{next: change the filename to setupFile}
		writeln('Will change parameter file name (',fn,') to ',setupFilev5);
		if not force then {warn user}
		begin
			writeln('Please press enter to continue.');
			readln
		end;
		
		if not RenameFile(setupFile, setupFilev5) then {this is the actual statement that renames the file}
			stop_prog('Could not rename '+setupFile+' to '+setupFilev5+'.', 3, not force);
		fn:=setupFilev5;

		{now check a few things to see if the parameters are compatible with the new version:}
		if par.num_GBs <> 0 then stop_prog('Sorry, grain boundaries are not (yet) supported by DevParUpdater.', 3, not force);
		if progName = 'SimSS' then 
		begin
			if par.ion_red_rate > 1 then stop_prog('Sorry, SIMsalabim 5.00 only supports ion_red_rate = 0 or 1.', 3, not force);
			if par.mob_ion_spec <> 0 then stop_prog('Sorry, SIMsalabim 5.00 only supports mob_ion_spec = 0.', 3, not force);
		end;
		
		{Now we need to see how many layers we have in the original file, so look for L_LTL and L_RTL in orgParStr:}
		parnew.NLayers:=1; {there is always at least 1 layer!}
		LTLexists:=par.L_LTL <> 0;
		if LTLexists then begin
			mainLyr:=2;
			inc(parnew.NLayers)
		end
		else 
			mainLyr:=1;
		RTLexists:=par.L_RTL <> 0;
		if RTLexists then inc(parnew.NLayers);

		writeln('Will create ',parnew.NLayers,' layers.');
		SetLength(parnew.lyr, parnew.NLayers+1);

{copy parameters of LTL:}
		if LTLexists then
		with parnew.lyr[mainLyr-1] do begin
			layerFile:='L'+IntToStr(mainLyr-1)+layerFileName;
			L:=par.L_LTL;
			eps_r:=par.eps_r_LTL;
			CB:=par.CB_LTL;
			VB:=par.VB_LTL;
			Nc:=par.Nc_LTL;
			n_0:=-min(0, par.doping_LTL);
			p_0:=max(0, par.doping_LTL);
			mun_0:=par.mob_LTL;
			mup_0:=par.mob_LTL;
			mob_n_dep:=0;
			mob_p_dep:=0;
			gamma_n:=0;
			gamma_p:=0;
			nu_int:=par.nu_int_LTL;
			St:=par.St_L;
			ETrapInt:=par.ETrapSingle;
			IntTrapFile:=par.IntTrapFile;
			Tr_type_Int:=par.Tr_type_L;
			CnInt:=par.Cn;
			CpInt:=par.Cp;
			if par.IonsInTLs then begin
				CNI:=par.CNI;
				IonsMayEnter:=true;
				CPI:=par.CPI
			end
			else begin
				CNI:=0;
				IonsMayEnter:=false;
				CPI:=0			
			end;
			
			if ProgName='SimSS' then begin
				mobnion:=1e-14;
				mobpion:=1e-14
			end
			else begin
				mobnion:=par.mobnion;
				mobpion:=par.mobpion
			end;
			
			layerGen:=par.TLsGen;
			nk_layer:=par.nk_LTL;
			Field_dep_G:=par.Field_dep_G;
			P0:=par.P0;
			a:=par.a;
			ThermLengDist:=par.ThermLengDist;
			kf:=par.kf;
			kdirect:=par.kdirect;
			Lang_pre:=par.Lang_pre;
			UseLangevin:=par.UseLangevin;

			if par.TLsTrap then 
				Bulk_tr:=par.Bulk_tr
			else
				Bulk_tr:=0;
			CnBulk:=par.Cn;
			CpBulk:=par.Cp;
			ETrapBulk:=par.ETrapSingle;
			BulkTrapFile:=par.BulkTrapFile;
			Tr_type_B:=par.Tr_type_B

		end; {copying of LTL parameters is now done!}

{Copy parameters of main absorber layer:}
		with parnew.lyr[mainLyr] do begin
			layerFile:='L'+IntToStr(mainLyr)+layerFileName;
			L:=par.L-par.L_LTL-par.L_RTL;
			eps_r:=par.eps_r;
			CB:=par.CB;
			VB:=par.VB;
			Nc:=par.Nc;
			n_0:=par.n_0;
			p_0:=par.p_0;
			mun_0:=par.mun_0;;
			mup_0:=par.mup_0;
			mob_n_dep:=par.mob_n_dep;
			mob_p_dep:=par.mob_p_dep;
			gamma_n:=par.gamma_n;
			gamma_p:=par.gamma_p;
			nu_int:=par.nu_int_RTL;
			St:=par.St_R;
			ETrapInt:=par.ETrapSingle;
			IntTrapFile:=par.IntTrapFile;
			Tr_type_Int:=par.Tr_type_R;
			CnInt:=par.Cn;
			CpInt:=par.Cp;
			CNI:=par.CNI;
			CPI:=par.CPI;
			IonsMayEnter:=true;
			
			if ProgName='SimSS' then begin
				mobnion:=1e-14;
				mobpion:=1e-14
			end
			else begin
				mobnion:=par.mobnion;
				mobpion:=par.mobpion
			end;
			
			layerGen:=TRUE;
			nk_layer:=par.nk_active;
			Field_dep_G:=par.Field_dep_G;
			P0:=par.P0;
			a:=par.a;
			ThermLengDist:=par.ThermLengDist;
			kf:=par.kf;
			kdirect:=par.kdirect;
			Lang_pre:=par.Lang_pre;
			UseLangevin:=par.UseLangevin;

			Bulk_tr:=par.Bulk_tr;
			CnBulk:=par.Cn;
			CpBulk:=par.Cp;
			ETrapBulk:=par.ETrapSingle;
			BulkTrapFile:=par.BulkTrapFile;
			Tr_type_B:=par.Tr_type_B

		end; {copying of main layer parameters is now done!}

{copy parameters of RTL:}
		if RTLexists then
		with parnew.lyr[mainLyr+1] do begin
			layerFile:='L'+IntToStr(mainLyr+1)+layerFileName;
			L:=par.L_RTL;
			eps_r:=par.eps_r_RTL;
			CB:=par.CB_RTL;
			VB:=par.VB_RTL;
			Nc:=par.Nc_RTL;
			n_0:=-min(0, par.doping_RTL);
			p_0:=max(0, par.doping_RTL);
			mun_0:=par.mob_RTL;
			mup_0:=par.mob_RTL;
			mob_n_dep:=0;
			mob_p_dep:=0;
			gamma_n:=0;
			gamma_p:=0;
			nu_int:=par.nu_int_RTL;
			St:=par.St_R;
			ETrapInt:=par.ETrapSingle;
			IntTrapFile:=par.IntTrapFile;
			Tr_type_Int:=par.Tr_type_R;
			CnInt:=par.Cn;
			CpInt:=par.Cp;
			if par.IonsInTLs then begin
				CNI:=par.CNI;
				IonsMayEnter:=true;
				CPI:=par.CPI
			end
			else begin
				CNI:=0;
				IonsMayEnter:=false;
				CPI:=0			
			end;
			
			if ProgName='SimSS' then begin
				mobnion:=1e-14;
				mobpion:=1e-14
			end
			else begin
				mobnion:=par.mobnion;
				mobpion:=par.mobpion
			end;
			
			layerGen:=par.TLsGen;
			nk_layer:=par.nk_RTL;
			Field_dep_G:=par.Field_dep_G;
			P0:=par.P0;
			a:=par.a;
			ThermLengDist:=par.ThermLengDist;
			kf:=par.kf;
			kdirect:=par.kdirect;
			Lang_pre:=par.Lang_pre;
			UseLangevin:=par.UseLangevin;

			if par.TLsTrap then 
				Bulk_tr:=par.Bulk_tr
			else
				Bulk_tr:=0;
			CnBulk:=par.Cn;
			CpBulk:=par.Cp;
			ETrapBulk:=par.ETrapSingle;
			BulkTrapFile:=par.BulkTrapFile;
			Tr_type_B:=par.Tr_type_B

		end; {copying of RTL parameters is now done!}

{copy general parameters:}
		with parnew do 
		begin
			T:=par.T;
			W_L:=par.W_L;
			W_R:=par.W_R;
			Sn_L:=par.Sn_L;
			Sp_L:=par.Sp_L;
			Sn_R:=par.Sn_R;
			Sp_R:=par.Sp_R;
			Rshunt:=par.Rshunt;
			Rseries:=par.Rseries;
			if ProgName='SimSS' then begin
				Gehp:=par.Gehp;
				Gfrac:=par.Gfrac
			end;
			Gen_profile:=par.Gen_profile;
			L_TCO:=par.L_TCO;
			L_BE:=par.L_BE;
			nk_substrate:=par.nk_substrate;
			nk_TCO:=par.nk_TCO;
			nk_BE:=par.nk_BE;
			spectrum:=par.spectrum;
			lambda_min:=par.lambda_min;
			lambda_max:=par.lambda_max;
			NP:=par.NP;
			tolPois:=par.tolPois;
			maxDelV:=par.maxDelV;
			maxItPois:=par.MaxItPois;
			maxItSS:=par.maxItSS;
			if ProgName='ZimT' then 
			begin
				MaxItTrans:=par.MaxItTrans;
				TolVint:=par.TolVint
			end;
			CurrDiffInt:=par.CurrDiffInt;
			tolDens:=par.tolDens;
			couplePC:=par.couplePC;
			minAcc:=par.minAcc;
			maxAcc:=par.maxAcc;
			IgnoreNegDens:=par.IgnoreNegDens;
			FailureMode:=par.FailureMode;
			grad:=par.grad;
			if ProgName='SimSS' then begin
				Vdistribution:=par.Vdistribution;
				PreCond:=par.PreCond;
				Vpre:=par.Vpre;
				fixIons:=par.ion_red_rate = 0; 
				Vscan:=par.Vscan;
				Vmin:=par.Vmin;
				Vmax:=par.Vmax;
				Vstep:=par.Vstep;
				Vacc:=par.Vacc;
				NJV:=par.NJV;
				until_Voc:=par.until_Voc
			end; {simss bit about voltages}
			timeout:=par.timeout;
			Pause_at_end:=par.Pause_at_end;
			AutoTidy:=par.AutoTidy;
			if ProgName='SimSS' then 
			begin
				UseExpData:=par.UseExpData;
				ExpJV:=par.ExpJV;
				rms_mode:=par.rms_mode;
				rms_threshold:=par.rms_threshold;
				JV_file:=par.JV_file;
				scPars_file:=par.scPars_file 
			end
			else begin {specific for ZimT}
				AutoStop:=par.AutoStop;
				tVG_file:=par.tVG_file;
				tj_file:=par.tj_file
			end;
			
			Var_file:=par.Var_file;
			LimitDigits:=par.LimitDigits;
			OutputRatio:=par.OutputRatio;
			log_file:=par.log_file;
			
		end; {copying of general parameters}
	
		{now we know all the input, so we can write the files:}
		write_simulation_setup_500(parnew, setupFile, path, ProgName); {produce general setup file}
		write_layer_500(parnew.lyr[mainLyr], path); {produce layer file for main absober}
		{new, do the TLs, if any:}
		if LTLExists then 
			write_layer_500(parnew.lyr[mainLyr-1], path);
		if RTLExists then 
			write_layer_500(parnew.lyr[mainLyr+1], path);

		{cleaning:}	
		Tidy_Up_Parameter_Files(setupFile, path, FALSE, parnew);
		
		version:=500;
		writeln('Updated to version ', newVersionString);
		writeln;

		{now tell user which files can be used from now on:}
		writeln('SIMsalabim v5.00 (and later) can now be used with files:');
		writeln(setupFilev5, ' (main setup file)');
		writeln('The parameter files pertaining to each layer are:');
		for i:=1 to parnew.NLayers do
			writeln(parnew.lyr[i].layerFile);
		if not force then 
		begin
			writeln('Please press enter to continue.');
			readln
		end
	end
end;

procedure UpdateTo501(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile501(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=500 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
			{in case we need to add/remove parameters, the code codes here:}
			
			{we add an option to W_L and W_R, so the comment changes:}
			if (pos('w_l=', LowerCase(DelSpace(line))) <> 0) then 
			begin
				line:=LeftStr(line, Pos('*', line)); {remove the bit after the *}
				line:=line + ' eV, work function left electrode (= cathode), or ''sfb''';
				CopyLine:=true
			end;
			if (pos('w_r=', LowerCase(DelSpace(line))) <> 0) then 
			begin
				line:=LeftStr(line, Pos('*', line)); {remove the bit after the *}
				line:=line + ' eV, work function right electrode (= cathode), or ''sfb''';
				CopyLine:=true
			end;
					
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile501(path, LayerFile); 
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo502(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile502(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			{we cannot simply replace 'St' as many words would contain this string, we start from the left}
			if pos('st=', LowerCase(DelSpace(line))) <> 0 then 
				line:=StringReplace(line, 'St', 'N_t_int',[]); {default behaviour: 1) case sensitive, and 2) only first occurrence will be replaced}
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
		buffer:=StringReplace(buffer, LineEnding+'CB', LineEnding+'E_c',[rfReplaceAll]); 		
		buffer:=StringReplace(buffer, LineEnding+'VB', LineEnding+'E_v',[rfReplaceAll]); 		
		buffer:=StringReplace(buffer, LineEnding+'Nc', LineEnding+'N_c',[rfReplaceAll]); 		
		buffer:=StringReplace(buffer, LineEnding+'n_0', LineEnding+'N_D',[rfReplaceAll]); 		
		buffer:=StringReplace(buffer, LineEnding+'p_0', LineEnding+'N_A',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'mun_0', LineEnding+'mu_n',[rfReplaceAll]); 		
		buffer:=StringReplace(buffer, LineEnding+'mup_0', LineEnding+'mu_p',[rfReplaceAll]); 		
		buffer:=StringReplace(buffer, LineEnding+'mob_n_dep', LineEnding+'mobnDep',[rfReplaceAll]); 		
		buffer:=StringReplace(buffer, LineEnding+'mob_p_dep', LineEnding+'mobpDep',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'ETrapInt', LineEnding+'E_t_int',[rfReplaceAll]);
		buffer:=StringReplace(buffer, 'EtrapSingle', 'E_t_int', [rfReplaceAll]); {this corrects a mistake in the comment of variable E_t_int} 
		buffer:=StringReplace(buffer, 'overrides EtrapBulk', 'overrides E_t_bulk', [rfReplaceAll]); {this corrects a mistake in the comment of variable bulkTrapFile} 
		buffer:=StringReplace(buffer, LineEnding+'IntTrapFile', LineEnding+'intTrapFile',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'Tr_type_Int', LineEnding+'intTrapType',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'CnInt', LineEnding+'C_n_int',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'CpInt', LineEnding+'C_p_int',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'CNI', LineEnding+'N_anion',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'CPI ', LineEnding+'N_cation',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'mobnion', LineEnding+'mu_anion',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'mobpion', LineEnding+'mu_cation',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'IonsMayEnter', LineEnding+'ionsMayEnter',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'nk_layer', LineEnding+'nkLayer',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'Field_dep_G', LineEnding+'fieldDepG',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'ThermLengDist', LineEnding+'thermLengDist',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'kf', LineEnding+'k_f',[rfReplaceAll]); 
		buffer:=StringReplace(buffer, LineEnding+'kdirect', LineEnding+'k_direct',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'Lang_pre', LineEnding+'preLangevin',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'UseLangevin', LineEnding+'useLangevin',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'Bulk_tr', LineEnding+'N_t_bulk',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'ETrapBulk', LineEnding+'E_t_bulk',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'BulkTrapFile', LineEnding+'bulkTrapFile',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'Tr_type_B', LineEnding+'bulkTrapType',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'CnBulk', LineEnding+'C_n_bulk',[rfReplaceAll]); 			
		buffer:=StringReplace(buffer, LineEnding+'CpBulk', LineEnding+'C_p_bulk',[rfReplaceAll]); 			
	
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=501 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
			{in case we need to add/remove parameters, the code codes here:}
			
			{search and replace many variable names:}
			{ReplaceText: 1) not case-sensitive, 2) replaces all occurrances}
			line:=ReplaceText(line, 'Sn_L', 'S_n_L'); 		
			line:=ReplaceText(line, 'Sp_L', 'S_p_L'); 		
			line:=ReplaceText(line, 'Sn_R', 'S_n_R'); 		
			line:=ReplaceText(line, 'Sp_R', 'S_p_R'); 		
			line:=ReplaceText(line, 'Rshunt', 'R_shunt'); 		
			line:=ReplaceText(line, 'Rseries', 'R_series'); 		
			line:=ReplaceText(line, 'Gehp', 'G_ehp'); 		
			line:=ReplaceText(line, 'Gfrac', 'G_frac'); 		
			line:=ReplaceText(line, 'Gen_profile', 'genProfile'); 		
			line:=ReplaceText(line, 'nk_substrate', 'nkSubstrate'); 		
			line:=ReplaceText(line, 'nk_TCO', 'nkTCO'); 		
			line:=ReplaceText(line, 'nk_BE', 'nkBE'); 		
			line:=ReplaceText(line, 'MaxItPois', 'maxItPois'); 		
			line:=ReplaceText(line, 'MaxItSS', 'maxItSS'); 		
			line:=ReplaceText(line, 'MaxItTrans', 'maxItTrans'); 		
			line:=ReplaceText(line, 'CurrDiffInt', 'currDiffInt'); 		
			line:=ReplaceText(line, 'IgnoreNegDens', 'ignoreNegDens'); 		
			line:=ReplaceText(line, 'FailureMode', 'failureMode'); 		
			line:=ReplaceText(line, 'TolVint', 'tolVint'); 		
			line:=ReplaceText(line, 'Vdistribution', 'Vdist'); 		
			line:=ReplaceText(line, 'PreCond', 'preCond'); 		
			line:=ReplaceText(line, 'FixIons', 'fixIons');
			line:=ReplaceText(line, 'until_Voc', 'untilVoc'); 		
			line:=ReplaceText(line, 'Pause_at_end', 'pauseAtEnd'); 		
			line:=ReplaceText(line, 'AutoTidy', 'autoTidy'); 		
			line:=ReplaceText(line, 'AutoStop', 'autoStop'); 		
			line:=ReplaceText(line, 'UseExpData', 'useExpData'); 		
			line:=ReplaceText(line, 'ExpJV', 'expJV'); 		
			line:=ReplaceText(line, 'rms_mode', 'rmsMode'); 		
			line:=ReplaceText(line, 'rms_threshold', 'rmsThreshold'); 		
			line:=ReplaceText(line, 'tVG_file', 'tVGFile'); 		
			line:=ReplaceText(line, 'JV_file', 'JVFile'); 		
			line:=ReplaceText(line, 'tj_file', 'tJFile'); 		
			line:=ReplaceText(line, 'Var_file', 'varFile'); 		
			line:=ReplaceText(line, 'LimitDigits', 'limitDigits'); 		
			line:=ReplaceText(line, 'OutputRatio', 'outputRatio'); 		
			line:=ReplaceText(line, 'scPars_file', 'scParsFile'); 		
			line:=ReplaceText(line, 'log_file', 'logFile'); 					
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile502(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo503(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile503(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
		{check if trap files are used => not supported as the user needs to change the files manually}
		dumstr:=LowerCase(DelSpace(buffer)); {remove all whitespace}
		if pos('inttrapfile=none', dumstr)=0 then Warn_User('Please change the intTrapFile contents: trap densities are absolute as of this version!', true);
		if pos('bulktrapfile=none', dumstr)=0 then Warn_User('Please change the bulkTrapFile contents: trap densities are absolute as of this version!', true);
		
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=502 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}


{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile503(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo504(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile504(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=503 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}


{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile504(path, LayerFile); 
				
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo505(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile505(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding;
			{Duplicate the line of the interface transfer velocity to replace later with electrons and holes respectively}
			if (pos('nu_int=', LowerCase(DelSpace(line))) <> 0) THEN 
				buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}

		{replace parameter nu_int with nu_int_n,p. First replace both with nu_int_p, then replace the first occurence again to nu_int_n}
		buffer:=StringReplace(buffer, LineEnding+'nu_int', LineEnding+'nu_int_p',[rfReplaceAll]);
		buffer:=StringReplace(buffer, LineEnding+'nu_int_p', LineEnding+'nu_int_n',[]); {Only the first occurence will be replaced}

		{Update the description with electrons/holes. First replace both with holes, then replace the first occurence again to electrons}
		buffer:=StringReplace(buffer, 'interface transfer velocity', 'interface transfer velocity of holes',[rfReplaceAll]);
		buffer:=StringReplace(buffer, 'interface transfer velocity of holes', 'interface transfer velocity of electrons',[]);{Only the first occurence will be replaced}
		
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=504 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}

{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile505(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo506(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile506(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=505 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}


{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile506(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo507(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile507(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=506 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}


{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile507(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo508(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion, i : integer;
	oldG_ehp, newG_ehp : myReal;
	L, Lgen, Ltot : myReal;
	gen : integer;
	LayerFiles : array of string;
	
	procedure PreReadLayerFile(path, LayerFile : string; var L : myReal; var gen : integer);
	var inp : text;
		line : string;
	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		while not eof(inp) do begin
			readln(inp, line);
		
			{try to find thickness L:}
			if pos('l=',LowerCase(DelWhite(line))) <> 0 then
			begin
				try
					dum1:=ExtractDelimited(2, DelWhite(line), ['=','*']);
					L:=StrToFloat(dum1)
				except 
					stop_prog('Failed to convert L from '+LayerFile+' to a float.', 3, true)
				end;			
			end;

			{try to find layerGen:}
			if pos('layergen=',LowerCase(DelWhite(line))) <> 0 then
			begin
				try
					dum1:=ExtractDelimited(2, DelWhite(line), ['=','*']);
					gen:=StrToInt(dum1)
				except 
					stop_prog('Failed to convert layerGen from '+LayerFile+' to an integer.', 3, true)
				end			
			end
				
		end;
		
		close(inp)

	end;
	
	procedure UpdateLayerFile508(path, LayerFile : string; newG_ehp : myReal); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			{this is where we insert G_ehp:}
			if pos('layergen=', LowerCase(DelWhite(line))) <> 0 then
			begin
				buffer:=buffer + 'G_ehp = ' + FloatToStr(newG_ehp) + '     * m^-3 s^-1, generation rate of electron-hole pairs in this layer ' + LineEnding
			end;
			
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=507 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		Lgen:=0;
		Ltot:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
		
			{in this update, we move G_ehp from the simulation setup to the layers. We also need to recalculate its value...}
			
			if pos('g_ehp=', LowerCase(DelWhite(line))) <> 0 then
			begin
				CopyLine:=false; {make sure we don't copy G_ehp from simulation setup file as we'll have to move it to the layer files}
				{we know the value of G_ehp sits between '=' and '*' plus any white space}
				try
					dum1:=ExtractDelimited(2, DelWhite(line), ['=','*']);
					oldG_ehp:=StrToFloat(dum1)
				except 
					stop_prog('Failed to convert G_ehp from '+SetupFile+' to a float.', 3, true)
				end;
			end;

			{next, we warn the user that user-supplied generation profiles need to be rescaled (not by the updater!!!)}
			if pos('genprofile', LowerCase(DelWhite(line))) <> 0 then 
			begin
				CopyLine:=true;
				dum1:=LowerCase(ExtractDelimited(2, DelWhite(line), ['=', '*']));
				if (dum1 <> 'none') and (dum1 <> 'calc') then warn_user('Please update the generation profile you are using: G_ehp should be absolute, in m-3 s-1.', not force);
				if (dum1 = 'none') and (ProgName = 'ZimT') then warn_user('For ZimT, please manually update G_ehp in the layer files.', not force)
			end;		

{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					SetLength(LayerFiles, LayerCounter+1); {add 1 element to array of file names}
					LayerFiles[LayerCounter]:=Trim(dum1);
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				{first: we read the file and check its thickness and whether it generates}
				PreReadLayerFile(path, LayerFiles[LayerCounter], L, gen);
				Ltot:=Ltot + L; {total thickness of all layers}
				Lgen:=Lgen + L*gen; {sum of thicknesses of layers that generated electron-hole pairs}

				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		
		{re-process the layer files, write the correct G_ehp:}
		if ProgName='ZimT' then
			oldG_ehp:=1; {some random value!}
		
		if Lgen <> 0 then
			newG_ehp:=oldG_ehp*Ltot/Lgen
		else
			newG_ehp:=0;
	
		for i:=1 to LayerCounter do {and update layer file:}
		begin
			GetProgNameVersion(dum1, LayerVersion, path + LayerFiles[i]); 
			if LayerVersion = version then {only update the layer file if we need to!}
				UpdateLayerFile508(path, LayerFiles[i], newG_ehp)
		end;
		
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;




procedure UpdateTo509(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile509(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=508 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
			buffer:=StringReplace(buffer, LineEnding+'rmsMode', LineEnding+'fitMode',[rfReplaceAll]); 	
			buffer:=StringReplace(buffer, LineEnding+'rmsThreshold', LineEnding+'fitThreshold',[rfReplaceAll]); 	
			buffer:=StringReplace(buffer, 'rms error', 'fit error',[rfReplaceAll]); 	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile509(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo510(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile510(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=509 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}

{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile510(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo511(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile511(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=510 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile511(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo512(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile512(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=511 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile512(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo513(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile513(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=512 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile513(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo514(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile514(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=513 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile514(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo515(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile515(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=514 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile515(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo516(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile516(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=515 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile516(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo517(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile517(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=516 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile517(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo518(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile518(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=517 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile518(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo519(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile519(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=518 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile519(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo520(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile520(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=519 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile520(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo521(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile521(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=520 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}

			{insert new parameter leftElec before W_L:}
			if pos('W_L', line) > 0 then
			begin
				CopyLine:=true;
				buffer:=buffer + 'leftElec = -1                 * left electrode is the cathode (-1) or the anode (1)' + LineEnding
			end;

			{remove ' (= cathode)' from the simulation setup file:}
			line:=ReplaceText(line, ' (= cathode)', '');
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile521(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;


procedure UpdateTo522(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile522(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=521 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile522(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo523(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile523(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=522 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile523(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo524(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile524(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=523 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile524(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo525(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile525(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=524 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile525(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo526(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile526(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=525 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile526(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
	end;
end;

procedure UpdateTo527(ProgName : string; var version : integer; SetupFile : string);
{updates the setup file (SetupFile) and layer files (LayerFile) by 1 version number}
var inp, outp : text;
	line, newVersionString, LayerFile : string;
	dum1, dum2 : string;
	buffer, path : ansistring;
	CopyLine : boolean;
	LayerCounter, LayerVersion : integer;
	
	procedure UpdateLayerFile527(path, LayerFile : string); 
	{only works on a single layer file}
	var inp, outp : text;
		buffer : ansistring;
		line, dumstr : string;
		posVersion, posVersionNumber : integer;

	begin
		if FileExists(path + LayerFile) then 
			assign(inp, path + LayerFile)
		else
			Stop_Prog('Could not find file '+LayerFile, 3, not force);
		reset(inp);
		
		{read entire file:}
		buffer:='';
		while not eof(inp) do begin
			readln(inp, line);
			buffer:=buffer + line + LineEnding
		end;
		close(inp);
		
		{change version number:}
		dumstr:=FormatFloat('0.00', 0.01*version);
		posVersion:=pos('version:', LowerCase(buffer));
		posVersionNumber:=pos(dumstr, buffer);
		if (posVersion <> 0) and (posVersion < posVersionNumber) then
			begin {now we simply cut and replace the bit of buffer with the version number:}
				Delete(buffer, posVersionNumber, length(dumstr));
				dumstr:=FormatFloat('0.00', 0.01*(version+1));
				Insert(dumstr, buffer, posVersionNumber)
			end
		else
			Stop_Prog('Could not find correct version number in file '+LayerFile, 3, not force);

{other code that modifies LayerFile goes here:}
	
{end of code that modifies parameters in LayerFile}
		
		{write new layer file:}	
		assign(outp, path + LayerFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp)	
	end;
	
begin
	if version=526 then {we only update 1 version number!}
	begin
		path:=ExtractFilePath(SetupFile); {we need the path to get to the right LayerFile!}
		
		assign(inp, SetupFile);
		reset(inp);
		
		buffer:=''; {this will be our output, empty at first}
		newVersionString:=FormatFloat('0.00', 0.01 * (version+1)); {this converts the version to the next one, so version=420 will yield '4.21'}
		LayerCounter:=0;
		
		while not eof(inp) do
		begin
			readln(inp, line); {read a line from the file}
			{if we find the version number, then we simply replace the version number}
			{otherwise, we simply copy}
			CopyLine:=true;
			if pos('version:', LowerCase(line)) <> 0 then
			begin
				CopyLine:=false;
				buffer:=buffer + '** version: '+newVersionString + LineEnding
			end;
			
{in case we need to add/remove parameters, the code goes here:}

			{insert new parameter tolCurr before tolDens:}
			if pos('tolDens', line) > 0 then
			begin
				CopyLine:=true;
				buffer:=buffer + 'tolCurr = 1E-3                    * relative tolerance on current density' + LineEnding
			end;
	
			{insert new parameter convVar before failureMode:}
			if pos('failureMode', line) > 0 then
			begin
				CopyLine:=true;
				buffer:=buffer + 'convVar = 1                       * integer 1-4, selects which variable to monitor for convergence' + LineEnding
						+ '                                  * 1: densities, 2: current, 3: densities OR current, 4: densities AND current' + LineEnding
			end;
	
{till here}
	
			{now treat the layer(s)}
			dum1:='l'+IntToStr(LayerCounter+1)+'=';
			dum2:=LowerCase(DelSpace(line)); {remove all whitespace}
			if LeftStr(dum2, length(dum1)) = dum1 then {found layer number (LayerCounter+1)}
			begin
				{now try to extract name of layer file:}
				try
					dum1:=ExtractWord(2, line, ['=', '*']); {the file name sits between = and a *}
					inc(LayerCounter);
					LayerFile:=Trim(dum1)
				except
					Stop_Prog('Error reading name of layer '+IntToStr(LayerCounter)+' in file '+SetupFile+'.', 3, not force)
				end;
				
				{and update layer file:}
				GetProgNameVersion(dum1, LayerVersion, path + LayerFile); 
				if LayerVersion = version then {only update the layer file if we need to!}
					UpdateLayerFile527(path, LayerFile); 
					
				CopyLine:=true
			end;	
		
			if CopyLine then
				buffer:=buffer + line + LineEnding {otherwise, simply copy}	
		end;
		
		close(inp);
		{now write the file:}
		assign(outp, SetupFile);
		rewrite(outp);
		write(outp, buffer);
		close(outp);
		
		inc(version);
		writeln('Updated to version ', newVersionString);
		if not force then 
		begin
			writeln('This version introduces two important variables: convVar and tolCurr.');
			writeln('Please check the manual to see how to use them to your advantage!');
			writeln('Press enter to continue...');
			readln
		end
	end
end;

begin
	{init and say hello}
	if hasCLoption('-h') then DisplayHelpExit;
	init_from_command_line; {read parameters from the command line in case there are any}
	writeln('Welcome to Updater.');
    writeln('Copyright (C) 2021, 2022, 2023, 2024, 2025, 2026');
    writeln('Prof. Dr. L.J.A. Koster, University of Groningen.');  
	writeln;

	{first check if file is new enough and exists:}
	if not FileExists(SetupFile) then Stop_Prog('Could not find file '+SetupFile, 3, not force);
	GetProgNameVersion(ProgName, version, SetupFile);
	if ProgName='' then Stop_Prog('Could not find name of the program in file '+SetupFile, 3, not force);
	if version=0 then Stop_Prog('Could not find version number in '+SetupFile, 3, not force);
	if version<minVersion then Stop_Prog('The specified parameter file is too old. I need version '+FormatFloat('0.00', 0.01 * (minVersion))+' or newer.', 3, not force);
	if version>maxVersion then Stop_Prog('The specified parameter file is already up to date.', 3, not force);
	
	writeln('Working on file ',SetupFile);
	writeln('Which is version ',FormatFloat('0.00', 0.01 * (version)),' of ',ProgName);
	writeln;
	
	{warn user we will overwrite the original file(s):}
	writeln('This program will overwrite the original file(s).');
	if not force then 
	begin
		writeln('Please confirm you wish to proceed. [yes/no]');
		readln(dumstr);
		if Trim(LowerCase(dumstr)) <> 'yes' then stop_prog('OK, we will exit.', 0, false);
	end;
	
	{now update step-by-step}
	UpdateTo420(ProgName, version, SetupFile);
	UpdateTo421(ProgName, version, SetupFile);
	UpdateTo422(ProgName, version, SetupFile);
	UpdateTo423(ProgName, version, SetupFile);
	UpdateTo424(ProgName, version, SetupFile);
	UpdateTo425(ProgName, version, SetupFile);
	UpdateTo426(ProgName, version, SetupFile);
	UpdateTo427(ProgName, version, SetupFile);
	UpdateTo428(ProgName, version, SetupFile);
	UpdateTo429(ProgName, version, SetupFile);
	UpdateTo430(ProgName, version, SetupFile);
	UpdateTo431(ProgName, version, SetupFile);
	UpdateTo432(ProgName, version, SetupFile);
	UpdateTo433(ProgName, version, SetupFile);
	UpdateTo434(ProgName, version, SetupFile);
	UpdateTo435(ProgName, version, SetupFile);
	UpdateTo436(ProgName, version, SetupFile);
	UpdateTo437(ProgName, version, SetupFile);
	UpdateTo438(ProgName, version, SetupFile);
	UpdateTo439(ProgName, version, SetupFile);
	UpdateTo440(ProgName, version, SetupFile);
	UpdateTo441(ProgName, version, SetupFile);
	UpdateTo442(ProgName, version, SetupFile);
	UpdateTo443(ProgName, version, SetupFile);
	UpdateTo444(ProgName, version, SetupFile);
	UpdateTo445(ProgName, version, SetupFile);
	UpdateTo446(ProgName, version, SetupFile);
	UpdateTo447(ProgName, version, SetupFile);
	UpdateTo448(ProgName, version, SetupFile);
	UpdateTo449(ProgName, version, SetupFile);
	UpdateTo450(ProgName, version, SetupFile);
	UpdateTo451(ProgName, version, SetupFile);
	UpdateTo452(ProgName, version, SetupFile);
	UpdateTo453(ProgName, version, SetupFile);
	UpdateTo454(ProgName, version, SetupFile);
	UpdateTo455(ProgName, version, SetupFile);
	UpdateTo456(ProgName, version, SetupFile);
	UpdateTo457(ProgName, version, SetupFile);
	UpdateTo458(ProgName, version, SetupFile);
	UpdateTo459(ProgName, version, SetupFile);

	{now we make a bigger step to v5.00}
	{this also changes the name of the file (SetupFile) to the default string}
	UpdateTo500(ProgName, version, SetupFile);

	{the rest are relatively minor updates:}
	UpdateTo501(ProgName, version, SetupFile);
	UpdateTo502(ProgName, version, SetupFile);
	UpdateTo503(ProgName, version, SetupFile);
	UpdateTo504(ProgName, version, SetupFile);
	UpdateTo505(ProgName, version, SetupFile);
	UpdateTo506(ProgName, version, SetupFile);
	UpdateTo507(ProgName, version, SetupFile);
	UpdateTo508(ProgName, version, SetupFile);
	UpdateTo509(ProgName, version, SetupFile);
	UpdateTo510(ProgName, version, SetupFile);
	UpdateTo511(ProgName, version, SetupFile);
	UpdateTo512(ProgName, version, SetupFile);
	UpdateTo513(ProgName, version, SetupFile);
	UpdateTo514(ProgName, version, SetupFile);
	UpdateTo515(ProgName, version, SetupFile);
	UpdateTo516(ProgName, version, SetupFile);
	UpdateTo517(ProgName, version, SetupFile);
	UpdateTo518(ProgName, version, SetupFile);
	UpdateTo519(ProgName, version, SetupFile);
	UpdateTo520(ProgName, version, SetupFile);
	UpdateTo521(ProgName, version, SetupFile);
	UpdateTo522(ProgName, version, SetupFile);
	UpdateTo523(ProgName, version, SetupFile);
	UpdateTo524(ProgName, version, SetupFile);
	UpdateTo525(ProgName, version, SetupFile);
	UpdateTo526(ProgName, version, SetupFile);
	UpdateTo527(ProgName, version, SetupFile);
	
	writeln;
	writeln('Done!');
end.
