unit InputOutputUtils;
{SIMsalabim input and output procedures.}

{
SIMsalabim: a 1D drift-diffusion simulator 
Copyright (c) 2020, 2021, 2023 Dr T.S. Sherkar, V.M. Le Corre, Dr M. Koopmans,
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

{$MODE DELPHI} {force DELPHI mode}

interface

uses TypesAndConstants, 
	 SysUtils, {for reading the command line and FileExists}
	 StrUtils; {for DelSpace}

type TSetChar = Set of Char;

function myBoolStr(val : boolean) : string;
{my implementation of BoolToStr. Returns 'TRUE' or 'FALSE'}

function DelWhite(str : string) : string;
{returns a copy of str with all white spaces (ASCII code 9,..13, and 32) removed from it.}

function DelWhite1(str : string) : string;
{returns a copy of str with all white spaces (ASCII code 9,..13, and 32) reduced to 1 space}

procedure getRealfromCL(ch : String; var foundit : boolean; var r : myReal; CaseSensitive : boolean = false);
{retrieves real number put after option ch in command line}

procedure getStringfromCL(ch : String; var foundit : boolean; var r : string; CaseSensitive : boolean = false);
{retrieves string put after option ch in command line}

function hasCLoption(ch : String; CaseSensitive : boolean = false) : boolean;
{checks if program was started with option ch in command line}

procedure Stop_Prog(msg : STRING; exitCode : INTEGER; wait_for_user : boolean = false);
{displays message (msg), waits for a hard return if wait_for_user and halts the program}

procedure Warn_User(msg : STRING; wait_for_user : boolean = false);
{displays message (msg), waits for a hard return if wait_for_user, does NOT halt the program}

function FindFile(key : string) : string;
{looks for a file with part of its name = key}

PROCEDURE Read_Number(VAR input : TEXT; VAR r : myReal);
{This procedure reads a number from the file input and puts it into the var r}
{The input file may contain comment lines starting with '*' or other characters
preceding the actual value, for instance: 'L = 100e-9 *m' will result in assigning
the value 100e-9 to r. After reading the number the line will be read until eoln}

PROCEDURE Read_Name(VAR inv : TEXT; var_name : STRING; VAR file_name : STRING);
{tries to get a file name for variable var_name from inv file
It discards any spaces, white lines, and *}
OVERLOAD;

PROCEDURE Read_Name(VAR inv : TEXT; var_name : STRING; VAR file_name : SHORTSTRING);
{tries to get a file name for variable var_name from inv file
It discards any spaces, white lines, and *}
OVERLOAD;

PROCEDURE Read_Integer(VAR inxut : TEXT; VAR k : INTEGER);
{reads an integer from file inxut using proc. Read_Number}

PROCEDURE Get_Float(VAR inv, log : TEXT; name_variable : STRING; VAR r : myReal); 
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be written in file log.}
OVERLOAD;

PROCEDURE Get_Float(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : myReal);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
OVERLOAD;

PROCEDURE Get_Integer(VAR inv, log : TEXT; name_variable : STRING; VAR r : INTEGER);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be written in file log.}
OVERLOAD;

PROCEDURE Get_Integer(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : INTEGER);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
OVERLOAD;

PROCEDURE Get_String(VAR inv, log : TEXT; name_variable : STRING; VAR r : STRING);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be written in file log.} 
OVERLOAD;

PROCEDURE Get_String(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : STRING);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
OVERLOAD;

PROCEDURE Get_String(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : SHORTSTRING);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
OVERLOAD;

procedure Read_Table(fileName : string; var data : Table; NumCol : integer; var NumLines : integer;
					 Delims : TSetChar = [#0..' ', ',', ';']; 
					 CommentSym : char = '*'; ErrorHandling : integer = 2;
					 header : string = '');
{Reads data from a file in a table format. Comments (after CommentSym) are ignored. Delims define the delimiters}
{NumCol: The first NumCol columns will be used. A line that contains fewer columns (but is not empty, or comments) will result in an error.}
{NumLines: the number of lines with data that was read}
{Delims: delimiters used to extract the data, so typically a ' ', or a comma}
{CommentSym: comments may appear after this character. Note: should NOT appear in Delims, nor be a decimal separator or a number}
{ErrorHandling: if 0, faulty data is ignored without warning. 1 => warning, but no halt. 2: program halts.}
{Header: Routine looks for this header and it must be there (not case-sensitive), but it can contain comments, etc.}

procedure Read_XY_Table(var x, y : Row; fileName, header : string; var NumLines : integer);
{Wrapper routine for Read_Table. Reads x,y data from a file, delimiter: white space. x, y are set to the correct length (index starting at 0)}

procedure Read_XYZ_Table(var x, y, z : Row; fileName, header : string; var NumLines : integer);
{Wrapper routine for Read_Table. Reads x,y,z data from a file, delimiter: white space. x, y, z are set to the correct length (index starting at 0)}

FUNCTION Count_Substring_In_String(SubStr, Str : STRING) : INTEGER;
{counts the number of occurances of character ch in string Str}
implementation

function myBoolStr(val : boolean) : string;
{my implementation of BoolToStr. Returns 'TRUE' or 'FALSE'}
begin
    if val 
        then myBoolStr:='TRUE'
        else myBoolStr:='FALSE'
end;

function DelWhite(str : string) : string;
{removes all white space from a string, so all characters 9-13 and 32(=space)}
var i : integer;
begin	
	{first convert all white spaces to normal spaces (char. 32):}
	for i:=1 to length(str) do 
		if ( str[i] in [#9..#13] ) then str[i]:=#32;
	{now remove all the spaces from the string:}
	str:=DelSpace(str); 
	{DelSpace returns a copy of S with all spaces (ASCII code 32) removed from it.}
	DelWhite:=str;
end;

function DelWhite1(str : string) : string;
{returns a copy of str with all white spaces (ASCII code 9,..13, and 32) reduced to 1 space}
var i : integer;
begin	
	{first convert all white spaces to normal spaces (char. 32):}
	for i:=1 to length(str) do 
		if ( str[i] in [#9..#13] ) then str[i]:=#32;
	{now reduce all spaces to 1 space}
	str:=DelSpace1(str); 
	DelWhite1:=str;
end;

procedure getRealfromCL(ch : String; var foundit : boolean; var r : myReal; CaseSensitive : boolean = false);
{tries to find ch in command line and returns the real value it preceeds}
var i : integer;
	foundkey : boolean;
begin
     i:=1;
     foundit:=false;
     while (i <= ParamCount) and (not foundit) do
     begin
          {did we find the key?}
          foundkey:= (ParamStr(i)=ch) or ((not CaseSensitive) and (LowerCase(ParamStr(i))=LowerCase(ch)));
          if foundkey and (i < ParamCount) then {we have found ch and there's something coming...}
          begin
               Try
                  r:=StrToFloat(ParamStr(i+1));
                  foundit:=true;
               except {note that if an exception occurs foundit stays false}
                     On val : Exception do
                     Stop_Prog('Exception when reading from command line : '+val.Message, EC_InvalidCLInput);
               end;
          end;
          i:=i+1;
     end;
end;

procedure getStringfromCL(ch : String; var foundit : boolean; var r : string; CaseSensitive : boolean = false);
{tries to find ch in command line and returns the string it preceeds}
var i : integer;
	foundkey : boolean; 
begin
     i:=1;
     foundit:=false;
     while (i <= ParamCount) and (not foundit) do
     begin
          {did we find the key?}
          foundkey:= (ParamStr(i)=ch) or ((not CaseSensitive) and (LowerCase(ParamStr(i))=LowerCase(ch)));
          if foundkey and (i < ParamCount) then {we have found ch and there's something coming...}
          begin
               Try
                  r:=ParamStr(i+1);
                  foundit:=true;
               except {note that if an exception occurs foundit stays false}
                     On val : Exception do
                     Stop_Prog('Exception when reading from command line : '+val.Message, EC_InvalidCLInput);
               end;
          end;
          i:=i+1;
     end;
end;

function hasCLoption(ch : String; CaseSensitive : boolean = false) : boolean;
{tries to find ch in command line}
var i : integer;
    foundit: boolean;
begin
     i:=1;
     foundit:=false;
     while (i <= ParamCount) and (not foundit) do
     begin
          foundit:= (ParamStr(i)=ch) or ((not CaseSensitive) and (LowerCase(ParamStr(i))=LowerCase(ch)));
          i:=i+1;
     end;
     hasCLoption:=foundit;
end;

PROCEDURE Stop_Prog(msg : STRING; exitCode : INTEGER; wait_for_user : boolean);
{Prints a message (msg) and stops the program}
BEGIN
    WRITELN('Program will be terminated.');
    WRITELN(msg);
    IF wait_for_user THEN
    BEGIN
		WRITELN('Please press enter.');
		READLN
	END;
    HALT(exitCode) {This terminates the program and returns the exit code so the outside world knows there was an error.}
END;

procedure Warn_User(msg : STRING; wait_for_user : boolean = false);
{displays message (msg), waits for a hard return if wait_for_user, does NOT halt the program}
BEGIN
    WRITELN('WARNING:');
    WRITELN(msg);
    IF wait_for_user THEN
    BEGIN
		WRITELN('Please press enter.');
		READLN
	END
END;

function FindFile(key : string) : string;
{looks for a file with part of its name = key}
var Info : TSearchRec;
begin
     If FindFirst (key, faAnyFile, Info) = 0
     then FindFile:=Info.Name
     else FindFile:='';
     FindClose(Info);
end;

PROCEDURE Read_Number(VAR input : TEXT; VAR r : myReal);
{This procedure reads a number from the file input and puts it into the var r}
{The input file may contain comment lines starting with '*' or other characters
preceding the actual value, for instance: 'L = 100e-9 *m' will result in assigning
the value 100e-9 to r. After reading the number the line will be read until eoln}
CONST numerical = ['0'..'9', '-', '+'];
VAR ch : CHAR;
    str_num : STRING;
    code : INTEGER;
    got_number, okay : BOOLEAN;
BEGIN
    okay:=FALSE;
    got_number:=FALSE; {we did not yet get the number}
    str_num:='';
    WHILE (NOT got_number) AND (NOT EOF(input)) DO
    BEGIN
        READ(input, ch);
        IF (ch IN numerical) THEN
        BEGIN
            str_num:=str_num + ch;
            got_number:=TRUE; {we did get to the number}
            WHILE (ch<>'*') AND (NOT EOLN(input)) DO
            BEGIN {read the rest of the number until further comment or eoln}
                READ(input, ch);
                IF (ch IN numerical) OR (ch='E') OR (ch='e') OR (ch='.') {e for exponent}
                    THEN str_num:=str_num + ch
            END;
        END
        ELSE IF ch='*' THEN READLN(input) {the line contains only comment}
    END;
    IF got_number
    THEN BEGIN
            READLN(input); {file position is at beginning of next line}
            VAL(str_num, r, code); {this converts the string into it's value}
            IF code=0 THEN okay:=TRUE    {reading of number is successful}
         END
    ELSE r:=0;  {reading of number is unsuccessful}
    IF NOT okay THEN {reading or conversion to real value is unsuccessful}
        Stop_Prog('Parameter file is corrupted.', EC_DevParCorrupt)

END;

PROCEDURE Read_Name(VAR inv : TEXT; var_name : STRING; VAR file_name : STRING);
{tries to get a file name for variable var_name from inv file
It discards any spaces, white lines, and *}
VAR line : STRING;
	i : INTEGER;
BEGIN
    REPEAT
		READLN(inv, line); {read line from input file}
		{FindPart: Search for a substring in a string, using wildcards}
		i:=FindPart(var_name, line); {try to find var_name in line}
		{i = 0 if var_name is not in line}
		{otherwise i indicates position of 1st character of var_name in line}
	UNTIL EOF(inv) OR (i>0);
	{if i=0 var_name not found, terminate program:}
	IF i=0 THEN Stop_Prog('Could not find variable ' + var_name + ' in input file.',EC_DevParCorrupt);

	{now get the part we're interested in:}
	i:=Length(line) - FindPart('=', line); {first get position of '=' sign, counting from the right}
	line:=RightStr(line, i); {remove anything in front of (and including) '=' sign}
	line:=TrimLeft(line); {TrimLeft: Trim whitespace from the beginning of a string}
	i:=FindPart('*', line); {find asterix if there is one}
	IF i > 0 THEN line:=LeftStr(line, i-1); {if there is an asterix, remove it}
	line:=TrimRight(line); {remove any whitespace at the end}

	{we should be done now, check result:}
	IF Length(line) > 0
		THEN file_name:=line
		ELSE Stop_Prog('Could not find correct file name for ' + var_name, EC_DevParCorrupt);

END;

PROCEDURE Read_Name(VAR inv : TEXT; var_name : STRING; VAR file_name : SHORTSTRING);
{tries to get a file name for variable var_name from inv file
It discards any spaces, white lines, and *}
VAR line : STRING;
	i : INTEGER;
BEGIN
    REPEAT
		READLN(inv, line); {read line from input file}
		{FindPart: Search for a substring in a string, using wildcards}
		i:=FindPart(var_name, line); {try to find var_name in line}
		{i = 0 if var_name is not in line}
		{otherwise i indicates position of 1st character of var_name in line}
	UNTIL EOF(inv) OR (i>0);
	{if i=0 var_name not found, terminate program:}
	IF i=0 THEN Stop_Prog('Could not find variable ' + var_name + ' in input file.',EC_DevParCorrupt);

	{now get the part we're interested in:}
	i:=Length(line) - FindPart('=', line); {first get position of '=' sign, counting from the right}
	line:=RightStr(line, i); {remove anything in front of (and including) '=' sign}
	line:=TrimLeft(line); {TrimLeft: Trim whitespace from the beginning of a string}
	i:=FindPart('*', line); {find asterix if there is one}
	IF i > 0 THEN line:=LeftStr(line, i-1); {if there is an asterix, remove it}
	line:=TrimRight(line); {remove any whitespace at the end}

	{we should be done now, check result:}
	IF Length(line) > 0
		THEN file_name:=line
		ELSE Stop_Prog('Could not find correct file name for ' + var_name, EC_DevParCorrupt);

END;

PROCEDURE Read_Integer(VAR inxut : TEXT; VAR k : INTEGER);
{reads an integer from file inxut using proc. Read_Number}
VAR dummy : myReal; {dummy is used for reading integer values, these are read as
                    TmyFloat and then converted to integers}
BEGIN
	Read_Number(inxut, dummy);
	k:=ROUND(dummy);
END;

PROCEDURE Get_Float(VAR inv, log : TEXT; name_variable : STRING; VAR r : myReal);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be written in file log.}
VAR gotit : BOOLEAN;
	dum : myReal;
BEGIN
	Read_Number(inv, r); {FIRST get value from parameter_file}
	{note: if we do not read r from dev.par. file, then the program crashes as it has lost its way in the file}
	getRealfromCL('-'+name_variable, gotit, dum); {try to get it from command line}
	{note: getRealfromCL is not case-sensitive by default}
	IF gotit
		THEN BEGIN
			r:=dum; {if we got it, then assign dum to r}
			WRITELN(log, name_variable,' = ',r); {and write value in log file}
		END;
END;

PROCEDURE Get_Float(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : myReal);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
VAR gotit : BOOLEAN;
	dum : myReal;
BEGIN
	Read_Number(inv, r); {FIRST get value from parameter_file}
	{note: if we do not read r from dev.par. file, then the program crashes as it has lost its way in the file}
	getRealfromCL('-'+name_variable, gotit, dum); {try to get it from command line}
	{note: getRealfromCL is not case-sensitive by default}
	IF gotit
		THEN BEGIN
			r:=dum; {if we got it, then assign dum to r}
			msgstr:=msgstr + name_variable + ' = '+FloatToStr(r) + LineEnding;
			{and store its name and value in the msgstr}
		END;
END;

PROCEDURE Get_Integer(VAR inv, log : TEXT; name_variable : STRING; VAR r : INTEGER);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be written in file log.}
VAR gotit : BOOLEAN;
	dum : myReal;
BEGIN
	Read_Integer(inv, r); {FIRST get value from parameter_file}
	{note: if we do not read r from dev.par. file, then the program crashes as it has lost its way in the file}
	getRealfromCL('-'+name_variable, gotit, dum); {try to get it from command line}
	{note: getRealfromCL is not case-sensitive by default}
		IF gotit
		THEN BEGIN
			r:=ROUND(dum); {if we got it, then assign dum to r}
			WRITELN(log, name_variable,' = ',r); {and write value in log file}
		END;
END;

PROCEDURE Get_Integer(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : INTEGER);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
VAR gotit : BOOLEAN;
	dum : myReal;
BEGIN
	Read_Integer(inv, r); {FIRST get value from parameter_file}
	{note: if we do not read r from dev.par. file, then the program crashes as it has lost its way in the file}
	getRealfromCL('-'+name_variable, gotit, dum); {try to get it from command line}
	{note: getRealfromCL is not case-sensitive by default}
		IF gotit
		THEN BEGIN
			r:=ROUND(dum); {if we got it, then assign dum to r}
			msgstr:=msgstr + name_variable + ' = '+IntToStr(r) + LineEnding
			{and store its name and value in the msgstr}
		END;
END;

PROCEDURE Get_String(VAR inv, log : TEXT; name_variable : STRING; VAR r : STRING);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be written in file log.} 
VAR gotit : BOOLEAN;
	dum : STRING;
BEGIN
	Read_Name(inv, name_variable, r); {FIRST get value from parameter_file}
	{note: if we do not read r from dev.par. file, then the program crashes as it has lost its way in the file}
	getStringfromCL('-'+name_variable, gotit, dum); {try to get it from command line}
	{note: getStringfromCL is not case-sensitive by default}
	IF gotit
		THEN BEGIN
			r:=dum; {if we got it, then assign dum to r}
			WRITELN(log, name_variable,' = ',r); {and write value in log file}
		END;
END;

PROCEDURE Get_String(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : STRING);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
VAR gotit : BOOLEAN;
	dum : STRING;
BEGIN
	Read_Name(inv, name_variable, r); {FIRST get value from parameter_file}
	{note: if we do not read r from dev.par. file, then the program crashes as it has lost its way in the file}
	getStringfromCL('-'+name_variable, gotit, dum); {try to get it from command line}
	{note: getStringfromCL is not case-sensitive by default}
	IF gotit
		THEN BEGIN
			r:=dum; {if we got it, then assign dum to r}
			msgstr:=msgstr + name_variable + ' = '+ r + LineEnding
			{and store its name and value in the msgstr}
		END;
END;

PROCEDURE Get_String(VAR inv : TEXT; VAR msgstr : ANSISTRING; name_variable : STRING; VAR r : SHORTSTRING);
{Reads variable 'name_variable' from input file 'inv' and assigns its value to 'r'
Next, it tries to get this value from the command line. If a command line value is
found, its name and value will be stored in the msgstr.}
VAR gotit : BOOLEAN;
	dum : STRING;
BEGIN
	Read_Name(inv, name_variable, r); {FIRST get value from parameter_file}
	{note: if we do not read r from dev.par. file, then the program crashes as it has lost its way in the file}
	getStringfromCL('-'+name_variable, gotit, dum); {try to get it from command line}
	{note: getStringfromCL is not case-sensitive by default}
	IF gotit
		THEN BEGIN
			r:=dum; {if we got it, then assign dum to r}
			msgstr:=msgstr + name_variable + ' = '+ r + LineEnding
			{and store its name and value in the msgstr}
		END;
END;

procedure Remove_Comments(var str : string; CommentSym : char); 
{remove anything after (and including) CommentSym from string}
var PosComment : integer;
begin
	PosComment:=Pos(CommentSym, str); {this is the position of the comment symbol. If 0, then no such character in str!}
	if PosComment > 0 then
		str:=Copy(str, 1, PosComment-1) {copy str until we get to the CommentSym}
end;

procedure Clean_Line(var str : string; CommentSym : char);
{removes excess white space and anything after CommentSym}
begin
	Remove_Comments(str, CommentSym); {remove anything after CommentSym from string}
	str:=DelWhite1(str); {returns a copy of str with all white spaces (ASCII code 9,..13, and 32) reduced to 1 space}
	str:=Trim(str); {remove any spaces on the left and the right of aline}
end;

procedure Read_Table(fileName : string; var data : Table; NumCol : integer; var NumLines : integer;
					 Delims : TSetChar = [#0..' ', ',', ';']; 
					 CommentSym : char = '*'; ErrorHandling : integer = 2;
					 header : string = '');
{Reads data from a file in a table format. Comments (after CommentSym) are ignored. Delims define the delimiters}
{NumCol: The first NumCol columns will be used. A line that contains fewer columns (but is not empty, or comments) will result in an error.}
{NumLines: the number of lines with data that was read}
{Delims: delimiters used to extract the data, so typically a ' ', or a comma}
{CommentSym: comments may appear after this character. Note: should NOT appear in Delims, nor be a decimal separator or a number}
{ErrorHandling: if 0, faulty data is ignored without warning. 1 => warning, but no halt. 2: program halts.}
{Header: Routine looks for this header and it must be there (not case-sensitive), but it can comments, etc.}

var inp : text;
	i, j, lineCount : integer;
	aline : string; {here we store a string line from the file}
	dataRow : array of myReal; {this is where we store a line of the data}
	FoundHeader : boolean;
	
	procedure Handle_Error(msg : string; ErrorHandling : integer);
	{local routine to handle error messages and warn/halt if need be}
	begin
		case ErrorHandling of
			0 : ; {just for completeness sake. We simply ignore any errors and proceed}
			1 : Warn_User(msg);
			2 : Stop_Prog(msg, EC_InvalidInput);
		otherwise Stop_Prog('Invalid ErrorHandling passed on to Read_Table.', EC_ProgrammingError);
		end;
	end;
	
begin

	if not FileExists(fileName) 
        then Stop_Prog('Cannot find file '+fileName, EC_FileNotFound);
    assign(inp, fileName); {once we get here, we're sure the file exists}
    reset(inp);
  
	i:=0;
	lineCount:=0;
	SetLength(dataRow, NumCol); 
	
	{first: check the header}
	Clean_Line(header, CommentSym); {first, we clean up the header string}
	if header <> '' then begin {OK, we're supposed to look for a header}
		FoundHeader:=false;
		
		{read file and check for header:}
		while (not eof(inp)) and (not FoundHeader) do begin
			readln(inp, aline);
			inc(lineCount);
			Clean_Line(aline, CommentSym);
			FoundHeader:=LowerCase(header) = LowerCase(LeftStr(aline, Length(header)))
		end;
		
		if not FoundHeader then {error, in case we still haven't found the header!}
			Handle_Error('Cannot find correct header in '+fileName+'.', ErrorHandling)
	end;
	
	{now read the data:}
	while not eof(inp) do
	begin
		readln(inp, aline);
		inc(lineCount);
		Clean_Line(aline, CommentSym); {removes excess white space and anything after CommentSym}
	
		{check if number of fields equals number of columns we're supposed to find:}
		if WordCount(aline, Delims) >= NumCol then begin {extract fields/numbers from the line}
			{OK, let's try to extract the data:}
			try
				for j:=0 to NumCol-1 do
					dataRow[j]:=StrToFloat(ExtractWord(j+1, aline, Delims)); {now extract the fields and convert to float:}
			except
				Handle_Error('Cannot process line '+IntToStr(lineCount)+' in '+filename+'.', ErrorHandling);
			end;
			
			{now the data are there, we need to copy it to the main table: data}
			inc(i);
			SetLength(data, i, NumCol); {add a line of data}
			for j:=0 to NumCol-1 do
				data[i-1, j]:=dataRow[j]
				
		end
		else 
			if Length(aline) <> 0 then
				Handle_Error('Cannot process line '+IntToStr(lineCount)+' in '+filename+'.', ErrorHandling);		

	end;
	
	NumLines:=i;
	close(inp)
	
end;

procedure Read_XY_Table(var x, y : Row; fileName, header : string; var NumLines : integer);
{Wrapper routine for Read_Table. Reads x,y data from a file, delimiter: white space. x, y are set to the correct length (index starting at 0)}
var data : Table;
	i : integer;
begin
	Read_Table(fileName, data, 2, NumLines, [#0..' '], '*', 2, header);
	
	{now set the length of x and y to the correct value, such that we can access indeces 0..NumLines-1:}
	SetLength(x, NumLines);
	SetLength(y, NumLines);
	
	{now copy data into x and y:}
	for i:=0 to NumLines-1 do
	begin
		x[i]:=data[i,0];
		y[i]:=data[i,1]
	end
end;

procedure Read_XYZ_Table(var x, y, z : Row; fileName, header : string; var NumLines : integer);
{Wrapper routine for Read_Table. Reads x,y,z data from a file, delimiter: white space. x, y, z are set to the correct length (index starting at 0)}
var data : Table;
	i : integer;
begin
	Read_Table(fileName, data, 3, NumLines, [#0..' '], '*', 2, header);
	
	{now set the length of x, y and z to the correct value, such that we can access indeces 0..NumLines-1:}
	SetLength(x, NumLines);
	SetLength(y, NumLines);
	SetLength(z, NumLines);
	
	{now copy data into x, y and z:}
	for i:=0 to NumLines-1 do
	begin
		x[i]:=data[i,0];
		y[i]:=data[i,1];
		z[i]:=data[i,2]
	end
end;

FUNCTION Count_Substring_In_String(SubStr, Str : STRING) : INTEGER;
{counts the number of occurances of substring in string Str}
VAR n : INTEGER;
BEGIN
    n:=0;
    WHILE NPos(SubStr, Str, n+1) > 0 DO
	INC(n);
    Count_Substring_In_String:=n
END;

begin

end.
