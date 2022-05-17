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
    collapseMethod

Description
    Collapse methodology based on the papers of Mihatsch et Al. and Arabnejad
    et Al.

\*---------------------------------------------------------------------------*/

//alphav and alpha_old are the complete lists for alpha at time steps t and t-1 (respectively)
volScalarField alpha_old = alphav.oldTime();

//Divergence of the velocity at time steps t and t-1 (respectively)
volScalarField divU( fvc::div(U) );
volScalarField divUOld( fvc::div( U.oldTime() ) );

int counterCelliCollapse = 0;
Info << "\n";
forAll( alphav, celli ){
    if ( alphav[celli] == 0.0 && Y[celli] == 1 && alpha_old[celli] > 0.0 )
    {
        /*labelList neighCellsi = mesh.cellCells()[celli]; //list of cells adjacent to celli
        
        int candidateCell = 0;
        forAll(neighCellsi, zz)
        {
            if ( alphav[neighCellsi[zz]] == 0.0 && Y[neighCellsi[zz]] == 1.0 ) //non va bene, tutte le 6 facce deveono avere essere ok allo stesso momemtno porco dio
                {
                    candidateCell++;*/
                    if ( divU[celli] > 0.0 && divUOld[celli] < 0.0)
                    {
                        counterCelliCollapse++;
                        pCollapse[celli] = p[celli]; //collapse pressure at celli, other cells have pressure 0 at time t
                        Info << "We have a collapse at Cell ID -> " << celli << "; with a pressure -> " << pCollapse[celli] << " [Pa]" << endl;    
                    }
                //}
        //}
    } 
}

Info << "\nNumber of cells in which we have collapse " << counterCelliCollapse << "\n" << endl; 






    