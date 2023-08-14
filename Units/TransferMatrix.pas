unit TransferMatrix;
{Calculates the generation profile for the device}

{
SIMsalabim:a 1D drift-diffusion simulator 
Copyright (c) 2021, 2022, 2023 S. Heester, Dr T.S. Sherkar, Dr V.M. Le Corre, 
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

CONST   TransferMatrixVersion = '4.56'; {version of this unit}

        lambda_step = 1E-9; {lambda step size}
        xstep = 1E-9; {grid step size (not TLs o and active layer)}

        MaxNoOfLayers = 6; {Max Number of Layers in device. Needs to be defined here untill move to N layers}

TYPE 
    TFMatrix = ARRAY[0..1,0..1] OF COMPLEX;

PROCEDURE Calc_TransferMatrix(VAR stv : TSTaticVars; CONSTREF par : TInputParameters);
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
    FOR a:=1 TO SizeLayers DO
    BEGIN
        FOR b:=0 TO sizeLambdas-1 DO
            retAbs[a][b] := (4*PI*nTotal[a][b].im)/(lambdas[b]);
    END;
    Calc_Absorption_Coefficient := retAbs;
END;

PROCEDURE Init_Layers_Thicknesses(CONSTREF par : TInputParameters; VAR Thicknesses : Row; VAR Layers : StringArray; VAR NoOfLayers, posA, posLTL, posRTL, posTCO, posBE : INTEGER );
{Build Arrays with the layer names and thicknesses based on the input parameters. If a layer thickness = 0 (except substrate layer), skip it and do not add it to either array. 
Array size is recorded with NoOfLayers and position indicators if relevant. If such a layer is not present, positon indicator remains 0}
BEGIN
    Thicknesses[NoOfLayers] := 0; {Substrate}
    layers[NoOfLayers] := par.nk_substrate; 
    INC(NoOfLayers);

    IF par.L_TCO > 0 THEN  {TCO}
    BEGIN
        Thicknesses[NoOfLayers] := par.L_TCO;
        layers[NoOfLayers] := par.nk_TCO;
        posTCO := NoOfLayers;
        INC(NoOfLayers);
    END;

    IF par.L_LTL > 0 THEN {LTL}
    BEGIN
        Thicknesses[NoOfLayers] := par.L_LTL;
        layers[NoOfLayers] := par.nk_LTL;
        posLTL := NoOfLayers;
        INC(NoOfLayers);
    END;

    Thicknesses[NoOfLayers] := (par.L - par.L_LTL - par.L_RTL); {Active}
    layers[NoOfLayers] := par.nk_active; 
    posA := NoOfLayers;
    INC(NoOfLayers);

    IF par.L_RTL > 0 THEN {LTL}
    BEGIN
        Thicknesses[NoOfLayers] := par.L_RTL;
        layers[NoOfLayers] := par.nk_RTL;
        posRTL := NoOfLayers;
        INC(NoOfLayers);
    END;

    Thicknesses[NoOfLayers] := par.L_BE; {Back electrode}
    layers[NoOfLayers] := par.nk_BE; 
    posBE := NoOfLayers;
    INC(NoOfLayers);
END;

PROCEDURE Read_nk_Material_From_File(idx, NoOfLambdas: INTEGER; Material : STRING; lambdas : row; VAR nTotal : ComplexMatrix);
{Read lambda,n,k values from a  file. 
Interpolate the read values to assign n,k values for each lambda provided. Store the complex variant for each lambda in nTotal as n + ik }
VAR l,NumValues : INTEGER;
    nInt, kInt, lFile, nFile, kFile : Row;
BEGIN
    Read_XYZ_Table(lFile, nFile, kFile, Material, 'lambda n k', NumValues);

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
    Read_XY_Table(lFile, vFile, filename, 'lambda I', NumAMValues);
    AMInt := Interpolate_1D_Linear(NoOfLambdas, numAMValues, lambdas, lFile, vFile); {Interpolate AM value}
    FOR l:= 0 TO NoOfLambdas-1 DO 
        AMValue[l] := AMInt[l]; 
END;

PROCEDURE Make_Grid_TM(CONSTREF par : TInputParameters; CONSTREF stv : TSTaticVars; posA, posLTL, posRTL, posTCO, posBE : INTEGER; VAR NoOfXPos : INTEGER; VAR xPos : Row; VAR xMat : intArray);
{Create a grid consisting of the existing grid for the TLs and active layer, extended with additional layers needed for the TransferMatrix Model. }
VAR f, NoOfXPosTCO, NoOfXPosBE, count : INTEGER;
BEGIN
    {Determine the number of grid points per layer and sum them to get the total number of grid points}
    NoOfXPosTCO := 0;
    IF par.L_TCO > 0 THEN {Check if TCO layer defined}
    BEGIN
        NoOfXPosTCO := Array_Size(0,par.L_TCO,xstep); {TCO layer}
        IF NoOfXPosTCO*xstep = par.L_TCO THEN {When the last x position of the TCO layer matches the layer boundary, exclude it because it will be taken into account in the next layer }
            NoOfXPosTCO := NoOfXPosTCO - 1;
    END;
    NoOfXPosBE := Array_Size(0,par.L_BE,xstep); {Back electrode}
    NoOfXPos := NoOfXPosTCO + par.NP + 1 + NoOfXPosBE; {Total number of grid points}
    {Create an Array with x positions. Use existing grid for Active layer and TLs. For all other layers, create a uniform grid. 
    Create matching xMat array (same length and indices) which contains the layer idx for the grid point.}
    SetLength(xPos,NoOfXPos);
    SetLength(xMat,NoOfXPos);

    count:=0;
    IF NoOfXPosTCO <> 0 THEN {TCO layer, uniform grid}
    BEGIN
        FOR f:=0 TO NoOfXPosTCO-1 DO
        BEGIN
            xPos[count] := f*xstep;
            xMat[count] := posTCO;
            INC(count);
        END;
    END;

    FOR f := 0 TO par.NP DO {TLs and active layer, exponential symmetric grid}
    BEGIN
        {i1 is the last point in the left insulator (or 0 if there isn't any)
        i2 is the first point in the right insulator (or NP+1 if there is none)}
        xPos[count] := stv.x[f] + par.L_TCO;
        IF f <= stv.i1 THEN 
            xMat[count] := posLTL
        ELSE IF (f > stv.i1) AND (f < stv.i2) THEN
            xMat[count] := posA
        ELSE
            xMat[count] := posRTL;

        INC(count);
    END;

    FOR f := 0 TO NoOfXPosBE -1 DO {Back electrode, uniform grid}
    BEGIN
        xPos[count] := par.L_TCO + par.L + (f)*xstep;
        xMat[count] := posBE;
        INC(count);
    END;
END;

PROCEDURE Calc_E(sizeLambdas, sizeLayers, sizeX : INTEGER; nTotal : ComplexMatrix; lambdas, RGlass : Row; xPos, thicknesses, TCumSum : Row; xMat : intArray; VAR E: ComplexMatrix);
{Calculate the optical electric field for an array of lambdas for the entire device based on the Transfer Matrices method}
VAR a,j,k,m, matIdx, posE, xCount : INTEGER;
    S, SPrime, SDPrime : TFMatrix;
    R, T,x : Row;
    xi, numer, denom : COMPLEX;
    dj : myReal;
BEGIN 
    {Initialise array sizes}
    SetLength(R,sizeLambdas);
    SetLength(T,sizeLambdas);
    SetLength(x,sizeX);

    FOR k:=0 TO sizeLambdas-1 DO {Calculate each lambda individually}
    BEGIN
        posE:=0; {Reset E Array idx for each layer}

        {First interface}
        S := I_Matrix(nTotal[0][k],nTotal[1][k]); {Initial S matrix for first interface}

        FOR j:=1 TO (sizeLayers-1) DO
            S := Matrix_Product(S,Matrix_Product(L_Matrix(nTotal[j][k],thicknesses[j],lambdas[k]),I_Matrix(nTotal[j][k],nTotal[j+1][k])));
        R[k] := (cmod(S[1][0]/S[0][0]))**2; {Reflection}
        T[k] := cmod((2/(1+nTotal[0,k])))/(csqrt(1-(RGlass[k]*R[k]))).re; {Transmission, Returns COMPLEX type but we only need the real part of the number}
        {Other interfaces/layers}
        FOR m:=1 TO sizeLayers DO {Exclude the first layer, because we already handled it}
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
                    x[xCount] := xPos[a]-TCumSum[m -1];{Store original pos left shifted to 0 of material in a new array}
                    INC(xCount); {Count the number of x positions for this layer. Used for looping over all x positions in the current layer}
                END;
            END;

            {Calculate S'}
            SPrime := I_Matrix(nTotal[0][k],nTotal[1][k]); {m = 1 (we excluded the first substrate/glass layer already)}
            IF m <> 1 THEN
            BEGIN
                FOR matIdx := 2 TO m DO
                    SPrime := Matrix_Product(SPrime,Matrix_Product(L_Matrix(nTotal[matIdx-1][k],thicknesses[matIdx-1],lambdas[k]),I_Matrix(nTotal[matIdx-1][k],nTotal[matIdx][k])));
            END;

            {Init and Calculate S'''}
            SDPrime[0][0] := 1;
            SDPrime[1][0] := 0;
            SDPrime[0][1] := 0;
            SDPrime[1][1] := 1;
            IF m <> sizeLayers THEN {Skip last iteration, because there is no layer beyond the last layer.}
            BEGIN
                FOR matIdx := m TO (sizeLayers - 1) DO
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

PROCEDURE Calc_Generation_Rate(sizeLambdas, sizeX, posA, posLTL, posRTL : INTEGER; lambdastep : myReal; xMat : intArray; lambdas, AMValue : Row; AbsCoef : Table; nTotal : ComplexMatrix; VAR NoOfGenRatePos : INTEGER; VAR genRate : vector);
{Calculate the generation rate based on the Electric Field / Q for all layer. Return only the transport and active layers}
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
        IF (xMat[a] = posA) OR (xMat[a] = posLTL) OR (xMat[a] = posRTL) THEN {When a layer is not present, position indicator is 0, which cannot exist and is thus skipped}
        //IF (xMat[a] = posA) THEN {When a layer is not present, position indicator is 0, which cannot exist and is thus skipped}
        BEGIN
            genRate[NoOfGenRatePos] := genRateFull[a]; {Shift the index of the genRate variable to exclude irrelevant layers}
            INC(NoOfGenRatePos); {Count number of values in both the transport and the active layers}
        END;
    END;
END;

PROCEDURE Calc_TransferMatrix(VAR stv : TSTaticVars; CONSTREF par : TInputParameters);
{Main procedure to calculate a generation profile based on input n,k values and spectrum}

{Calculate the generation profile for the device defined in the device parametes file. 
This method takes a 6 layer device as input, but can handle fewer layers. However the substrate, active layer and back electrode must always be defined}
VAR j,k,NoOfLayers,NoOfLambdas, NoOfXPos, NoOfGenRatePos, posA, posLTL, posRTL, posTCO, posBE : INTEGER;
    Layers : StringArray;
    nTotal : ComplexMatrix;
    Thicknesses, TCumSum, lambdas, AMValue, xPos, RGlass : Row;
    xMat : intArray;
    AbsCoef : Table;
    genRate : vector;
    uitv : TEXT;
    ii: INTEGER;
BEGIN
    SetLength(Layers,MaxNoOfLayers);
    SetLength(Thicknesses, MaxNoOfLayers);
    SetLength(TCumSum, MaxNoOfLayers);

    NoOfLayers := 0; {Total number of layers in the device. Minimum 3 (front,active,back)}
    posTCO := 0; {Position of the TCO layer. 0 if not present in the device}
    posBE := 0; {Position of thwhich is used as input for SimSS or ZimT.e back}
    posA := 0; {Position of the active layer}
    posLTL := 0; {Position of the left transport layer. 0 if not present in the device}
    posRTL := 0; {Position of the right transport layer. 0 if not present in the device}

    Init_Layers_Thicknesses(par, Thicknesses, Layers, NoOfLayers, posA, posLTL, posRTL, posTCO, posBE); {Create an array with layer thicknesses from input parameters. Identify the position of active and transport layers}

    IF par.lambda_min <> par.lambda_max THEN {Range of lambdas defined, create an evenly spaced array}
    BEGIN
        NoOfLambdas := Array_Size(par.lambda_min, par.lambda_max, lambda_step); {Get the number of elements for lambda/wavelength}
        SetLength(lambdas,NoOfLambdas); {Set size of lambdas array}
        FOR k:=0 TO NoOfLambdas-1 DO {Fill the evenly spaced lambdas array}
            lambdas[k] := (par.lambda_min + k*lambda_step);
    END
    ELSE {One lambda defined (lambda_min = lambda_max). Create an array with just one element}
    BEGIN
        NoOfLambdas := 1;
        SetLength(lambdas,NoOfLambdas); {Size of lambdas array = 1}
        lambdas[0]:=par.lambda_min;
    END;
    SetLength(AbsCoef, NoOfLayers+1, NoOfLambdas); {Initialise Absorption coefficient Matrix}
    SetLength(nTotal, NoOfLayers+1, NoOfLambdas); {Initialise nTotal Matrix}
    SetLength(RGlass, NoOfLambdas); {Initialise RGlass array}

    FOR j:=0 TO NoOfLayers-1 DO {Calculate complex refracive index based on the material type/name and material data files}
        Read_nk_Material_From_File(j, NoOfLambdas, Layers[j], lambdas, nTotal); {Sets complex index of refraction(n+ik) for all lambdas in Lambdas Array}

    FOR k:=0 TO NoOfLambdas-1 DO {Calculate Reflection/Transmission Coefficient for every lambda}
        RGlass[k] := cmod(((1-nTotal[0,k])/(1+nTotal[0,k]))**2);

    Thicknesses[0] := 0; {Force the first layer size to 0 (Glass)}
    TCumSum := Cumulative_Sum(NoOfLayers - 1, Thicknesses); {Create Array with cumulative sum of Thicknesses}

    Make_Grid_TM(par, stv, posA, posLTL, posRTL, posTCO, posBE, NoOfXPos, xPos, xMat); {Create a new temporary grid to include all layers}
 
    SetLength(E, NoOfXPos, NoOfLambdas); {Initialise size E matrix}
    Calc_E(NoOfLambdas, NoOfLayers - 1, NoOfXPos, nTotal, lambdas, RGlass, xPos, Thicknesses, TCumSum, xMat, E); {Calcualte total optical electric field at every position xPos for all lambdas}

    NoOfGenRatePos := 0; {Init number of xpos for generation profile}
    AbsCoef := Calc_Absorption_Coefficient(NoOfLambdas, NoOfLayers - 1, lambdas, nTotal); {Absorption Coefficient per layer}

    SetLength(AMValue,NoOfLambdas); {Initialise AM values Array}
    Read_AM_From_File(par.spectrum, NoOfLambdas, lambdas, AMValue);  {Read AM File}

    Calc_Generation_Rate(NoOfLambdas, NoOfXPos, posA, posLTL, posRTL, lambda_step, xMat, lambdas, AMValue, AbsCoef, nTotal, NoOfGenRatePos, genRate); {Calculate generation rate (genRate)}
    
    stv.orgGm := genRate; {Fits on the existing grid (par.NP + 1)}

END;
{CHANGELOG
Added a procedure to calculate the generation profile based on the device structure and properties and placed it in the unit TransferMatrix.
The procedure is derived from the TransferMatrix method from Burkhard (2011). The standard device strucure is: Substrate | ITO | LTL | Active | RTL | Back, but not all layers are mandatory, only the Substrate, active and back layer are.
Each layer is characterised by a thickness and material n,k values, defined in the device parametes.
The result of this procedure is a generation profile (with only the LTL, RTL and active layer, if present) with the real value of Gehp, not just the profile.
}
BEGIN
END.

