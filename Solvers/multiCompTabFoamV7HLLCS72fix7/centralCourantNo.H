/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

if (mesh.nInternalFaces()){
    scalar i = 0;

    /*surfaceScalarField cSf_pos2("cSf_pos2",interpolate(c, pos)*mesh.magSf() );
        surfaceScalarField cSf_neg2("cSf_neg2",interpolate(c, neg)*mesh.magSf() );
        surfaceScalarField phiv_pos2("phiv_pos2", (U_pos & mesh.Sf()) );//mesh.magSf()
        surfaceScalarField phiv_neg2("phiv_neg2", (U_neg & mesh.Sf()) );//mesh.magSf()
        surfaceScalarField amaxSf("amaxSf", max( max(mag(phiv_pos2+cSf_pos2),mag(phiv_neg2+cSf_neg2)) 
                                                    , max(mag(phiv_pos2-cSf_pos2),mag(phiv_neg2-cSf_neg2)) ) );*/
    
    surfaceScalarField ucSf( mag( fvc::interpolate(U) & mesh.Sf() )+fvc::interpolate(mag(c))*mesh.magSf() ); 
    //surfaceScalarField ucSf2( mag( fvc::interpolate(U) & mesh.Sf() )+fvc::interpolate(mag(c))*mesh.magSf() ); 

        
    forAll(ucSf.primitiveFieldRef(),fi) {
        //ucSf.internalField()[fi]/=(mesh.V()[mesh.owner()[fi]]/mesh.nGeometricD());
          ucSf.primitiveFieldRef()[fi]/= (mesh.V()[mesh.owner()[fi]]/mesh.nGeometricD() );
   } 

   forAll(ucSf, j){
       if( ucSf[j]*runTime.deltaTValue() >= 1 )
            {
                i++;
            }      
   }

    //scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());
    //phi = ucSf;//linearInterpolate(U) & mesh.Sf();

    //CoNumrho = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();
    
    CoNum = gMax( ucSf )*runTime.deltaTValue();
    CoNumMin = gMin( ucSf )*runTime.deltaTValue();
    //CoNumhy = gMax( ucSf2/mesh.V().field() ) * runTime.deltaTValue();
    //CFL.ref() = 0.5*fvc::surfaceSum(sumAmaxSf)/mesh.V().field() * runTime.deltaTValue();
    //CFL.correctBoundaryConditions();
    meanCoNum = gAverage( ucSf )*runTime.deltaTValue();
    //meanCoNumhy = gSum(ucSf2)/gSum(mesh.V().field())*runTime.deltaTValue(); 
    //meanCoNumrho =
    //    0.5*(gSum(sumAmaxSf)/gSum(mesh.V().field()))*runTime.deltaTValue();  
    
    Info << "Cells with Co >= 1 = " << i << endl;
}


Info << "Mean, min and max Courant Numbers = " << meanCoNum << " " << CoNumMin << " " << CoNum << endl;
//Info << "hybrid mean and max Co = " << meanCoNumhy << " " << CoNumhy << endl;
//Info << "rhocentralfoam one: MEand and max = "  << meanCoNumrho << " " << CoNumrho << endl;


// ************************************************************************* //
