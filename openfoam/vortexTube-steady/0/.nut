/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{

    stator
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    rotor
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    shroud
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    inlet
    {
        //type            calculated;
        //value           uniform 0;
type cyclic;
    }

    outlet
    {
        //type            calculated;
        //value           uniform 0;
type cyclic;    
}



}


// ************************************************************************* //
