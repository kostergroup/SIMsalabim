unit TransferMatrix;
{Calculates the generation profile for the device}

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
email: l.j.a.koster@rug.nl
surface mail: 
L.J.A. Koster
Zernike Institute for Advanced Materials
Nijenborgh 4, 9747 AG Groningen, the Netherlands
}

{$MODE OBJFPC} 
{$NOTES OFF} {Ignore the inline Notes (Note: Call to subroutine "operator ... marked as inline is not inlined) from the UComplex unit}
INTERFACE

USES 
    TypesAndConstants,
    InputOutputUtils, 
    NumericalUtils,
    DDTypesAndConstants,
	Math,  
	StrUtils,
	SysUtils,
    Ucomplex;

CONST   TransferMatrixVersion = '5.13'; {version of this unit}

        lambda_step = 1E-9; {lambda step size}
        xstep = 1E-9; {grid step size (Only for the TCO and BE layers.)}
TYPE 
    TFMatrix = ARRAY[0..1,0..1] OF COMPLEX;

PROCEDURE Calc_TransferMatrix(VAR stv : TSTaticVars; VAR log : TEXT; CONSTREF par : TInputParameters);
{Main procedure to calculate a generation profile based on input n,k values and spectrum}

IMPLEMENTATION

VAR  E : ComplexMatrix; {Use a global variable because when implementing as a local variable, program crashes with memory issues}

FUNCTION Array_Size(xmin, xmax, xstep : myReal) : INTEGER;
{Calculate the size of an evenly spaced bound array}
BEGIN
    Array_Size := Round((xmax - xmin)/xstep);
END;

FUNCTION Cumulative_Sum(size : INTEGER; arr : Row) : Row;
{Calculate the Cumulative Sum of an Array of length: size}
VAR a: INTEGER;
    retsum: Row;
BEGIN
    SetLength(retsum,size+1); {Initialize array size}
    For a:=0 TO size DO
    BEGIN
        IF (a = 0) THEN {Set the thickness of the first layer manually for indexing}
        retsum[a] := arr[a]
        ELSE 
        retsum[a] := retsum[a-1] + arr[a]; {Rest of the layers}
    END;
    Cumulative_Sum := retsum;
END;

FUNCTION I_Matrix(n1,n2 : COMPLEX) : TFMatrix;
{Calculate the Transfermatrix at an interface between two materials (1,2) using their complex refractive index}
VAR rI,tI: COMPLEX;
    retI: TFMatrix;
BEGIN
    rI := (n1-n2)/(n1+n2); {Reflection coeffient}
    tI := (2*n1)/(n1+n2); {Transmission coefficient}
    {Transfer Matrix of size: [2 x 2]. Fill each element individually}
    retI[0][0] := 1/tI;
    retI[0][1] := rI/tI;
    retI[1][0] := rI/tI;
    retI[1][1] := 1/tI;
    I_Matrix := retI;
END;

FUNCTION L_Matrix(n: COMPLEX; d,lambda: myReal) : TFMatrix;
{Calculate the Transfermatrix for a layer characterized by complex refractive index and thickness (n,d) for a single wavelength (lambda) (1,2)}
VAR xiL: COMPLEX;
    retL: TFMatrix;
BEGIN
    xiL := 2*PI*n/lambda; {layer phase thickness corresponding to the phase change the wave experiences as it traverses the layer}
    {Transfer Matrix of size: [2 x 2]. Fill each element individually. }
    retL[0][0] := CEXP(-1*i*xiL*d);
    retL[0][1] := 0;
    retL[1][0] := 0;
    retL[1][1] := CEXP(1*i*xiL*d);
    L_Matrix := retL;
END;

FUNCTION Matrix_Product(mat1,mat2 : TFMatrix) : TFMatrix;
{Perform standard Matrix multiplication for two matrices. Matrix dimensions: mat1 [a x c], mat2 [c x b]}
VAR a,b,c : INTEGER;
    retmat : TFMatrix;
BEGIN
    FOR a:=0 TO length(mat1)-1 DO 
    BEGIN
        FOR b:=0 TO length(mat2[0])-1 DO
        BEGIN
            retmat[a][b] := 0;
            FOR c:= 0 TO length(mat1[0])-1 DO
                retmat[a][b] := retmat[a][b] + mat1[a][c]*mat2[c][b];
        END;
    END;
    Matrix_Product := retmat;
END;

FUNCTION Interpolate_1D_Linear(xIntSize, NumValues: INTEGER; xInt, xFile, yFile : Row) : Row;
{Make a linear interpolation based on the InterExtraPolation function based on Neville's algorithm for point xFile,yFile. 
Calculate the 'y' values (retInt) for a provided input 'x' array (xInt)}
VAR a,b : INTEGER;
    x,y, retInt: Row;
    IntSuccess: BOOLEAN;
    y_error: myREAL;
BEGIN
    IntSuccess := TRUE; {Boolean to indicate whether interpolation has succeeded}
    SetLength(x, NumValues+1); 
    SetLength(y, NumValues+1); 
    SetLength(retInt, xIntSize); {Init array for interpolated values. Size must match number of x positions}

    FOR a:=0 TO NumValues-1 DO {Create an array with all 'non-zero' values from the file read}
    BEGIN
        x[a+1] := xFile[a];
        y[a+1] := yFile[a];
    END;

    FOR b:= 0 TO xIntSize-1 DO {Calculate an 'interpolated' value for every xInt}
    BEGIN
        IF (IntSuccess) THEN
        begin
            IntSuccess := InterExtraPolation(x,y,xInt[b],retInt[b],y_error,1,2);
        END;
    END;
    Interpolate_1D_Linear := retInt;
END;

FUNCTION Calc_Absorption_Coefficient(sizeLambdas, SizeLayers : INTEGER; lambdas : row; nTotal : ComplexMatrix) : Table;
{Calculate the absorption coefficient for each layer for each lambda}
VAR a,b: INTEGER;
    retAbs : Table;
BEGIN
    SetLength(retAbs, SizeLayers+1, sizeLambdas); {Initialize array size}    
    FOR a:=0 TO SizeLayers DO
    BEGIN
        FOR b:=0 TO sizeLambdas-1 DO
            retAbs[a][b] := (4*PI*nTotal[a][b].im)/(lambdas[b]);
    END;
    Calc_Absorption_Coefficient := retAbs;
END;
 
PROCEDURE Init_Layers_Thicknesses(CONSTREF par : TInputParameters; VAR stv : TSTaticVars; VAR Thicknesses : Row; VAR Layers : StringArray; VAR StartIdx : INTEGER );
{ Create arrays with all layer thicknesses and layer names. Add the TCO layer and BE layer. 
If the thickness of the TCO layer is 0 the StartIdx indictor is set to 1. 
When this is this case, the TCO layer idx/row/column is skipped/ignored for the rest of the Calc_TransferMatrix procedure.}
    
VAR i : INTEGER;
BEGIN
    Thicknesses[0] := par.L_TCO;
    Layers[0] := par.nkTCO;

    IF par.L_TCO > 0 THEN StartIdx := 0;

    FOR i := 1 TO stv.NLayers DO
    BEGIN
        Thicknesses[i] := par.lyr[i].L;
        Layers[i]:= par.lyr[i].nkLayer;
    END;

    Thicknesses[i+1] := par.L_BE;
    Layers[i+1] := par.nkBE;
END;

PROCEDURE Read_nk_Material_From_File(idx, NoOfLambdas: INTEGER; Material : STRING; lambdas : row; VAR nTotal : ComplexMatrix);
{Read lambda,n,k values from a  file. 
Interpolate the read values to assign n,k values for each lambda provided. Store the complex variant for each lambda in nTotal as n + ik }
VAR l,NumValues : INTEGER;
    nInt, kInt, lFile, nFile, kFile : Row;
BEGIN
    writeln('Reading nk file: ',Material);
    Read_XYZ_Table(lFile, nFile, kFile, Material, 'lambda n k', NumValues);

    {Check the whether the spacing between wavelenths is at least the minimal value as defined in the DDTypesAndConstants unit}
    FOR l:=1 to NoOfLambdas-1 DO
        BEGIN
            IF (lFile[l]-lFile[l-1]) < 0.999*minDeltaLambda THEN {Use a tolerance of 0.1% to account for floating point variations}
			    Stop_Prog('Spacing between consecutive wavelengths in the nk-values in file ' + Material + ' is too small, must be >= ' + FloatToStr(minDeltaLambda) + 'm', EC_InvalidInput)
        END;
    
    {Interpolate values to match provided lambdas}
    nInt := Interpolate_1D_Linear(NoOfLambdas, NumValues, lambdas, lFile, nFile); {n value (real part of nTotal)}
    kInt := Interpolate_1D_Linear(NoOfLambdas, NumValues, lambdas, lFile, kFile); { k value (imaginary part of nTotal)}
    FOR l:= 0 TO NoOfLambdas-1 DO {Merge n,k into the complex index of refraction and store in an array}
    BEGIN
        nTotal[idx][l].re := nInt[l]; {Real part}
        nTotal[idx][l].im := kInt[l];  {Imaginary part}
    END;
END;

PROCEDURE Read_AM_From_File(filename: STRING; NoOfLambdas:INTEGER; VAR lambdas, AMValue : Row);
{Read AM value from a file. 
Interpolate the read value to assign AMvalue for each lambda provided }
VAR l, numAMValues : INTEGER;
    amInt,lFile,vFile : Row;
BEGIN
    writeln('Reading spectrum: ',filename);
    Read_XY_Table(lFile, vFile, filename, 'lambda I', NumAMValues);

    {Check the whether the spacing between wavelenths is at least the minimal value as defined in the DDTypesAndConstants unit}
    FOR l:=1 to NumAMValues-1 DO
        BEGIN
            IF (lFile[l]-lFile[l-1]) < 0.999*minDeltaLambda THEN {Use a tolerance of 0.1% to account for floating point variations}
			    Stop_Prog('Spacing between consecutive wavelengths in the used spectrum is too small, must be >= ' + FloatToStr(minDeltaLambda) + 'm', EC_InvalidInput)
        END;

    AMInt := Interpolate_1D_Linear(NoOfLambdas, numAMValues, lambdas, lFile, vFile); {Interpolate AM value}
    FOR l:= 0 TO NoOfLambdas-1 DO 
        AMValue[l] := AMInt[l]; 
END;

PROCEDURE Make_Grid_TM(CONSTREF par : TInputParameters; CONSTREF stv : TSTaticVars; StartIdx : INTEGER; VAR NoOfXPos : INTEGER; VAR xPos : Row; VAR xMat : intArray);
{Create a grid consisting of the existing grid from SIMsalabim, extended with additional layers needed for the TransferMatrix Model. }
VAR f, NoOfXPosTCO, NoOfXPosBE : INTEGER;
BEGIN
    {Determine the number of grid points per layer and sum them to get the total number of grid points}
    NoOfXPosTCO := 0;
    IF par.L_TCO > 0 THEN {Check if TCO layer defined}
    BEGIN
        NoOfXPosTCO := Array_Size(0,par.L_TCO,xstep); {TCO layer}
        IF NoOfXPosTCO*xstep = par.L_TCO THEN {When the last x position of the TCO layer matches the layer boundary, exclude it because it will be taken into account in the next layer }
            NoOfXPosTCO := NoOfXPosTCO - 1;
    END;
    NoOfXPosBE := Array_Size(0,par.L_BE,xstep) +1; {Back electrode, include the very last point}
    NoOfXPos := NoOfXPosTCO + (par.NP + 2) + NoOfXPosBE; {Total number of grid points}
    {Create an Array (grid) with x positions. For the TCO and BE layers, create an unifrom grid.  
    Create matching xMat array (same length and indices) which contains the layer idx for the grid point. Idx 0 is reserved for the TCO, even when it has a thickness of 0.}
    SetLength(xPos,NoOfXPos);
    SetLength(xMat,NoOfXPos);

    IF NoOfXPosTCO <> 0 THEN {TCO layer, uniform grid}
    BEGIN
        FOR f:=0 TO NoOfXPosTCO-1 DO
        BEGIN
            xPos[f] := f*xstep;
            xMat[f] := 0;
        END;
    END;

    FOR f:=0 TO par.NP +1 DO {The existing grid from SIMsalabim}
    BEGIN
        xPos[f+NoOfXPosTCO] := stv.x[f] + par.L_TCO;
        xMat[f+NoOfXPosTCO] := stv.lid[f];
    END;

    FOR f := 1 TO NoOfXPosBE DO {Back electrode, uniform grid}
    BEGIN
        xPos[f+NoOfXPosTCO + par.NP +1] := par.L_TCO + stv.Ltot + (f)*xstep;
        xMat[f+NoOfXPosTCO + par.NP +1] := stv.NLayers +1;
    END;
END;

PROCEDURE Calc_E(sizeLambdas, sizeLayers, sizeX, StartIdx : INTEGER; nTotal, nSubstrate : ComplexMatrix; lambdas, RGlass : Row; xPos, thicknesses, TCumSum : Row; xMat : intArray; VAR E: ComplexMatrix);
{Calculate the optical electric field for an array of lambdas for the entire device based on the Transfer Matrix method}
VAR a,j,k,m, matIdx, posE, xCount : INTEGER;
    S, SPrime, SDPrime : TFMatrix;
    R, T,x : Row;
    xi, numer, denom : COMPLEX;
    dj : myReal;
BEGIN 
    {Initialise array sizes}
    SetLength(R,sizeLambdas);
    SetLength(T,sizeLambdas);
    SetLength(x,sizeX+1);

    FOR k:=0 TO sizeLambdas-1 DO {Calculate each lambda individually}
    BEGIN
        posE:=0; {Reset E Array idx for each layer}

        {First interface}
        S := I_Matrix(nSubstrate[0][k],nTotal[StartIdx][k]); {Initial S matrix for first interface}
        FOR j:=StartIdx TO (sizeLayers-2) DO
            S := Matrix_Product(S,Matrix_Product(L_Matrix(nTotal[j][k],thicknesses[j],lambdas[k]),I_Matrix(nTotal[j][k],nTotal[j+1][k])));
        R[k] := (cmod(S[1][0]/S[0][0]))**2; {Reflection}
        T[k] := cmod((2/(1+nSubstrate[0,k])))/(csqrt(1-(RGlass[k]*R[k]))).re; {Transmission, Returns COMPLEX type but we only need the real part of the number}
        {Other interfaces/layers}
        FOR m:=StartIdx TO sizeLayers - 1 DO {Exclude the first layer, because we already handled it}
        BEGIN
            xi := 2*PI*nTotal[m,k]/Lambdas[k]; 
            dj := thicknesses[m]; {layer thickness}

            xCount:=0; {Reset x position counter, because x refers to the position with respect to the left bound of the layer, not the device}
            FOR a:=0 TO sizeX-1 DO {Clear the array to set correct number of values for new material}
                x[a] :=0;
            FOR a:=0 to sizeX-1 DO {Fill x array with new positions for the selected material}
            BEGIN
                IF (xMat[a] = m) THEN {x position material index matches looped material, thus x position must be used for this run }
                BEGIN
                    IF (m = StartIdx) THEN
                        {The first layer starts at x = 0 anyway, so no need to shift the x position.}
                        x[xCount] := xPos[a]
                    ELSE
                        x[xCount] := xPos[a]-TCumSum[m -1];{Store original pos left shifted to 0 of material in a new array}
                    INC(xCount); {Count the number of x positions for this layer. Used for looping over all x positions in the current layer}
                END;
            END;

            {Calculate S'}
            SPrime := I_Matrix(nSubstrate[0][k],nTotal[StartIdx][k]); {m = 1 (we excluded the first substrate/glass layer already)}
            IF m <> StartIdx THEN
            BEGIN
                FOR matIdx := StartIdx + 1 TO m DO
                    SPrime := Matrix_Product(SPrime,Matrix_Product(L_Matrix(nTotal[matIdx-1][k],thicknesses[matIdx-1],lambdas[k]),I_Matrix(nTotal[matIdx-1][k],nTotal[matIdx][k])));
            END;

            {Init and Calculate S'''}
            SDPrime[0][0] := 1;
            SDPrime[1][0] := 0;
            SDPrime[0][1] := 0;
            SDPrime[1][1] := 1;
            IF m <> sizeLayers-1 THEN {Skip last iteration, because there is no layer beyond the last layer.}
            BEGIN
                FOR matIdx := m TO (sizeLayers - 2) DO
                    SDPrime := Matrix_Product(SDPrime,Matrix_Product(I_Matrix(nTotal[matIdx][k],nTotal[matIdx+1][k]),L_Matrix(nTotal[matIdx+1][k],thicknesses[matIdx+1],lambdas[k])));
            END;

            {Calculate E for every x position (layer specific), but store value per lambda at unshofted x position/index via posE .}
            FOR a:=0 TO xCount-1 DO
            BEGIN
                {Split in numerator/denominator for readability}
                numer := T[k]*(SDPrime[0][0]*CEXP(-1*i*xi*(dj-x[a])) + SDPrime[1][0]*CEXP(1*i*xi*(dj-x[a])));
                denom := SPrime[0][0]*SDPrime[0][0]*CEXP(-1*i*xi*dj) + SPrime[0][1]*SDPrime[1][0]*CEXP(1*i*xi*dj);
                E[posE][k] := numer/denom; {Epos represents the actual x position in the device, not the x position in the layer itself}
                INC(posE);
            END;
        END;
    END;
END;

PROCEDURE Calc_Generation_Rate(CONSTREF stv : TSTaticVars; CONSTREF par : TInputParameters; sizeLambdas, sizeX : INTEGER; lambdastep : myReal; xMat : intArray; lambdas, AMValue : Row; AbsCoef : Table; nTotal : ComplexMatrix; VAR NoOfGenRatePos : INTEGER; VAR genRate : vector);
{Calculate the generation rate based on the Electric Field / Q for all layer. Do not return the TCO and BE layer.}
VAR a,k : INTEGER;
    Q, genRateFull : Row;
BEGIN
    SetLength(Q, sizeX);
    SetLength(genRateFull, sizeX);

    FOR a:=0 TO sizeX - 1 DO
    BEGIN
        Q[a] :=0; {Initialize which is used as input for SimSS or ZimT.Q Array}
        FOR k:=0 TO sizeLambdas -1 DO
            Q[a] := Q[a] + lambdastep*(AbsCoef[xMat[a]][k]*nTotal[xMat[a]][k].re*AMvalue[k]*(cmod(E[a][k])**2)*lambdas[k])/(h*c); {Calculate and Sum Q over lambda for all x positions}
        genRateFull[a]:= Q[a]; {Set the generation rate for a xpos for full device}
        {If the x position is in either the transport layers or active layers, it must be written to the generation profile file}
        IF (xMat[a] <> 0) AND (xMat[a] <> stv.NLayers+1) THEN {When a layer is not present, position indicator is 0, which cannot exist and is thus skipped}
        BEGIN
            IF par.lyr[xMat[a]].layerGen THEN
				genRate[NoOfGenRatePos] := genRateFull[a] {Shift the index of the genRate variable to exclude irrelevant layers}
			ELSE genRate[NoOfGenRatePos] := 0; {layer does not absorb/generate electron-hole pairs}
            INC(NoOfGenRatePos); {Count number of entries}
        END;
    END;
END;

PROCEDURE Calc_TransferMatrix(VAR stv : TSTaticVars; VAR log : TEXT; CONSTREF par : TInputParameters);
{Main procedure to calculate a generation profile based on input n,k values and spectrum}

{Calculate the generation profile for the device defined in the device parametes file, including the substrate, TCO and BE.} 
VAR j,k,NoOfLambdas, NoOfXPos, NoOfGenRatePos,StartIdx : INTEGER;
    Layers : StringArray;
    nTotal, nSubstrate : ComplexMatrix;
    Thicknesses, TCumSum, lambdas, AMValue, xPos, RGlass : Row;
    xMat : intArray;
    AbsCoef : Table;
    genRate : vector;
    lambdaInd, lambdaFormat : myReal;

BEGIN
    WRITELN;
    WRITELN('Started calculation of the generation profile.');

    SetLength(Layers,stv.NLayers+2); {Add 2 to account for the TCO and BE}
    SetLength(Thicknesses, stv.NLayers+2);
    SetLength(TCumSum, stv.NLayers+2);
    StartIdx := 1; {Indicate whether TCO is defined. 1 means the TCO is not defined. Set in Init_Layers_Thicknesses procedure }

    Init_Layers_Thicknesses(par, stv, Thicknesses, Layers, StartIdx); {Create an array with layer thicknesses from input parameters.}

    IF par.lambda_min <> par.lambda_max THEN {Range of lambdas defined, create an evenly spaced array}
    BEGIN
        NoOfLambdas := Array_Size(par.lambda_min, par.lambda_max, lambda_step)+1; {Get the number of elements for lambda/wavelength}
        SetLength(lambdas,NoOfLambdas); {Set size of lambdas array}
        FOR k:=0 TO NoOfLambdas-1 DO {Fill the evenly spaced lambdas array}
        BEGIN
            lambdaInd := par.lambda_min + k*lambda_step;
            lambdaFormat := StrToFloat(Format('%.6e', [lambdaInd])); {Limit the number of decimals to 6 to prevent interpolation issues when processing the nk- and spectrum values}
            lambdas[k] := lambdaFormat;
        END;
    END
    ELSE {One lambda defined (lambda_min = lambda_max). Create an array with just one element}
    BEGIN
        NoOfLambdas := 1;
        SetLength(lambdas,NoOfLambdas); {Size of lambdas array = 1}
        lambdas[0]:=par.lambda_min;
    END;

    SetLength(AbsCoef, stv.NLayers+2, NoOfLambdas); {Initialise Absorption coefficient Matrix}
    SetLength(nTotal, stv.NLayers+2, NoOfLambdas); {Initialise nTotal Matrix}
    SetLength(nSubstrate, 1, NoOfLambdas); {Initialise nSubstrate Matrix}
    SetLength(RGlass, NoOfLambdas); {Initialise RGlass array}

    {The nk vlaues for the substrate are need for the calculation, but are not treated as a layer. They require a seperate parameter.}
    Read_nk_Material_From_File(0, NoOfLambdas, par.nkSubstrate, lambdas, nSubstrate); {Sets complex index of refraction(n+ik) for all lambdas in Lambdas Array}

    FOR j:=StartIdx TO stv.NLayers+1 DO {Calculate complex refracive index based on the material type/name and material data files}
        Read_nk_Material_From_File(j, NoOfLambdas, Layers[j], lambdas, nTotal); {Sets complex index of refraction(n+ik) for all lambdas in Lambdas Array}

    FOR k:=0 TO NoOfLambdas-1 DO {Calculate Reflection/Transmission Coefficient for every lambda}
        RGlass[k] := cmod(((1-nSubstrate[0,k])/(1+nSubstrate[0,k]))**2);

    TCumSum := Cumulative_Sum(stv.NLayers +1, Thicknesses); {Create Array with cumulative sum of Thicknesses}

    Make_Grid_TM(par, stv, StartIdx, NoOfXPos, xPos, xMat); {Create a new temporary grid to include all layers}

    SetLength(E, NoOfXPos, NoOfLambdas); {Initialise size E matrix}

    Calc_E(NoOfLambdas, stv.Nlayers + 2, NoOfXPos, StartIdx, nTotal, nSubstrate, lambdas, RGlass, xPos, Thicknesses, TCumSum, xMat, E); {Calcualte total optical electric field at every position xPos for all lambdas}

    NoOfGenRatePos := 0; {Init number of xpos for generation profile}
    AbsCoef := Calc_Absorption_Coefficient(NoOfLambdas, stv.NLayers + 1, lambdas, nTotal); {Absorption Coefficient per layer}

    SetLength(AMValue,NoOfLambdas); {Initialise AM values Array}
    Read_AM_From_File(par.spectrum, NoOfLambdas, lambdas, AMValue);  {Read AM File}

    Calc_Generation_Rate(stv, par, NoOfLambdas, NoOfXPos, lambda_step, xMat, lambdas, AMValue, AbsCoef, nTotal, NoOfGenRatePos, genRate); {Calculate generation rate (genRate)}
    
    stv.orgGm := genRate; {Fits on the existing grid (par.NP + 1)}

    WRITELN(log);
    WRITELN(log, 'Files used to calculate the generation profile:');
    
    FOR j:=StartIdx TO stv.NLayers+1 DO
        WRITELN(log, '- ' + Layers[j]);
     WRITELN(log, '- ' + par.spectrum);
     WRITELN(log);

    WRITELN('Calculation of the generation profile was successful');
    WRITELN;
END;

BEGIN
END.

