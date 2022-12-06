/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    centralCourantNo

Description
    Calculates the mean and maximum wave speed based Courant Numbers.

\*---------------------------------------------------------------------------*/
        
/*        surfaceScalarField phiv_pos2("phiv_pos2", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg2("phiv_neg2", U_neg & mesh.Sf());
        
        surfaceScalarField cSf_pos2
        (
            "cSf_pos2",
            interpolate(c, pos)*mesh.magSf()
        );

        surfaceScalarField cSf_neg2
        (
            "cSf_neg2",
            interpolate(c, neg)*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos2 + cSf_pos2, phiv_neg2 + cSf_neg2), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos2 - cSf_pos2, phiv_neg2 - cSf_neg2), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);


if (mesh.nInternalFaces())
{
    scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumAmaxSf)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Mean and max Courant Numbers = "
    << meanCoNum << " " << CoNum << endl; */

/*
// Set the reciprocal time-step from the local Courant number
    dt11.ref() = max
    (
       1/dimensionedScalar(dimTime, maxDeltaT),
        fvc::surfaceSum(amaxSf)()()
       /((2*maxCo)*mesh.V())
    );

    forAll(dt11,celli)
{
    Info << "dt[" << celli << "] = " << dt11[celli] << endl;
}*/    volScalarField maxpercell_ucSf1("maxpercell_ucSf1", maxpercell_ucSf);
    
    surfaceScalarField ucSf1( mag( fvc::interpolate(U) & mesh.Sf() )+ fvc::interpolate(mag(c))*mesh.magSf()); 

    forAll(ucSf1.primitiveFieldRef(),fi) 
    {
        ucSf1.primitiveFieldRef()[fi]/= (mesh.V()[mesh.owner()[fi]]/mesh.nGeometricD());
        ucSf1.primitiveFieldRef()[fi]=1/ucSf1.primitiveFieldRef()[fi];
   
        if(ucSf1.primitiveFieldRef()[fi]<maxpercell_ucSf1[mesh.owner()[fi]])
            {
                maxpercell_ucSf1[mesh.owner()[fi]]=ucSf1.primitiveFieldRef()[fi];
            }

        if (ucSf1.primitiveFieldRef()[fi]<maxpercell_ucSf1[mesh.neighbour()[fi]])
        {
                maxpercell_ucSf1[mesh.neighbour()[fi]]=ucSf1.primitiveFieldRef()[fi];
        }
    } 

    maxpercell_ucSf= maxpercell_ucSf1;

Info << gMax(maxpercell_ucSf) << endl;

Info << "Mean and max Courant Numbers = " << meanCoNum << " " << CoNum << endl;










// ************************************************************************* //
