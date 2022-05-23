
#include "fvCFD.H"
#include "dynamicFvMesh.H"
//#include "turbulentTransportModel.H"
#include "fixedRhoFvPatchScalarField.H"
//#include "incompressibleTwoPhaseMixture.H"
#include "twoPhaseMixtureThermo.H"
#include "turbulentFluidThermoModel.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "tableInfo.H" //thermo class definition
#include "functions.H" //functions implementation
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //#define NO_CONTROL
    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    turbulence->validate();
    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include <fstream>
    #include <iostream>
    #include <string>
    #include <stdio.h>
    #include <stdlib.h>
    //#include <algorithm>
    //#include <vector>

    dimensionedScalar uBounded("uBounded", dimLength/dimTime, 0.01);
  // surfaceScalarField ucSf( mag( fvc::interpolate(U) & mesh.Sf() )+fvc::interpolate(mag(c))*mesh.magSf() );
   OFstream resFile("residual.txt");

    word nameOfTable(" ");
    runTime.controlDict().readIfPresent("nameOfTable", nameOfTable); //read file from controlDict
    tableInfo unif(nameOfTable); //call constructor
    
    #include "calculateThermo.H"

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar CoNumMin = 0.0;

    volScalarField rho_0("rho_0",rho);                                      //NADA
    volScalarField rhoY_0("rhoY_0",rhoY);
	volVectorField rhoU_0("rhoU_0",rhoU);	
	volScalarField rhoE_0("rhoE_0",rhoE);
    
    volScalarField rho_1("rho_1",rho);                                      //NADA
    volScalarField rhoY_1("rhoY_1",rhoY);
	volVectorField rhoU_1("rhoU_1",rhoU);	
	volScalarField rhoE_1("rhoE_1",rhoE);

    volScalarField rho_2("rho_2",rho);                                      //NADA
    volScalarField rhoY_2("rhoY_2",rhoY);
	volVectorField rhoU_2("rhoU_2",rhoU);	
	volScalarField rhoE_2("rhoE_2",rhoE);   
   
    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {

    //EMPIEZA
    volScalarField rhoOld("rhoOld",rho);                                    //NADA
    volScalarField rhoYOld("rhoYOld",rhoY);
	volVectorField rhoUOld("rhoUOld",rhoU);	
	volScalarField rhoEOld("rhoEOld",rhoE);

	rho1.ref() = rho.ref();                                                 //NADA
    rhoU1.ref()= rhoU.ref();
    rhoE1.ref() = rhoE.ref();
    rhoY1.ref() = rhoY.ref();
    p1.ref() = p.ref();
    c1.ref() = c.ref();

    U.correctBoundaryConditions();                                           //NADA
    Y.correctBoundaryConditions();
    rho.correctBoundaryConditions();
    e.correctBoundaryConditions();
    rhoY.boundaryFieldRef() == rho.boundaryField()*Y.boundaryField();
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
    rhoE.boundaryFieldRef() == rho.boundaryField()*( e.boundaryField() + 0.5*magSqr(U.boundaryField()));
    p.correctBoundaryConditions();
    c.correctBoundaryConditions();

    rho1.correctBoundaryConditions();                                         //NADA
    rhoY1.boundaryFieldRef() == rho.boundaryField()*Y.boundaryField();
    rhoU1.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
    rhoE1.boundaryFieldRef() == rho.boundaryField()*( e.boundaryField() + 0.5*magSqr(U.boundaryField()));
    p1.correctBoundaryConditions();
    c1.correctBoundaryConditions();

    #include "createFlux.H"      

    rhok[0]() = fvc::div(rhoPhi);                                            // SUMMATORY FLUXES
    rhoYk[0]() = fvc::div(rhoYPhi);
    rhoUk[0]()= fvc::div(phiUp)-fvc::laplacian(muEff, U) - fvc::div(tauMC);
    rhoEk[0]() = fvc::div(phiEUp) - fvc::div(sigmaDotU)- fvc::laplacian(k_eff,T);

    rho_0 = rho;
    rhoY_0 = rhoY;
    rhoU_0 = rhoU;
    rhoE_0 = rhoE;

    for( int cycle = 0; cycle < 3; cycle++)
    {

	    Info << "starting cycle " << RK4values.size()<< cycle << endl;

        volScalarField& rhokcycle  = rhok[cycle]();                            //NADA
   	    volScalarField& rhoYkcycle = rhoYk[cycle]();
	    volVectorField& rhoUkcycle = rhoUk[cycle]();
	    volScalarField& rhoEkcycle = rhoEk[cycle]();
        
        switch(cycle) 
        {
        case 0:
	    solve(fvm::ddt(rho)  == -rhokcycle);                                     //UPDATE CONSERVATIVE VARIABLES
        solve(fvm::ddt(rhoY) == -rhoYkcycle);
        solve(fvm::ddt(rhoU) == -rhoUkcycle);
        solve(fvm::ddt(rhoE) == -rhoEkcycle);

        rho_1 = rho;
        rhoY_1 = rhoY;
        rhoU_1 = rhoU;
        rhoE_1 = rhoE;

        Info<< "\nZERO\n" << endl;
        break;

        case 1:
	    solve(fvm::ddt(rho)  == -rhokcycle);                                  //UPDATE CONSERVATIVE VARIABLES
        solve(fvm::ddt(rhoY) == -rhoYkcycle);
        solve(fvm::ddt(rhoU) == -rhoUkcycle);
        solve(fvm::ddt(rhoE) == -rhoEkcycle);

        rho =  rho*0.25+0.75*rho_0;
        rhoY = rhoY*0.25+0.75*rhoY_0;
        rhoU = rhoU*0.25+0.75*rhoU_0;
        rhoE = rhoE*0.25+0.75*rhoE_0;

        rho_2  = rho;
        rhoY_2 = rhoY;
        rhoU_2 = rhoU;
        rhoE_2 = rhoE;

        Info<< "\nUNO\n" << endl;
        break;

        case 2:
	    solve(fvm::ddt(rho)  == -rhokcycle);                                  //UPDATE CONSERVATIVE VARIABLES
        solve(fvm::ddt(rhoY) == -rhoYkcycle);
        solve(fvm::ddt(rhoU) == -rhoUkcycle);
        solve(fvm::ddt(rhoE) == -rhoEkcycle);

        rho =  2.0/3.0*rho+1.0/3.0*rho_0;
        rhoY = 2.0/3.0*rhoY+1.0/3.0*rhoY_0;
        rhoU=  2.0/3.0*rhoU+1.0/3.0*rhoU_0;
        rhoE = 2.0/3.0*rhoE+1.0/3.0*rhoE_0;

        Info<< "\nDOS\n" << endl;
        Info<< runTime.deltaTValue() << endl;

        break;
        }

        forAll(rho, celli){                                                     //NADA
            rho.ref()[celli] = max(0.25,rho.ref()[celli]);
            CFL.ref()[celli]=0.0;
        }
        Y.ref() = rhoY.ref() / rho.ref();                                        //NADA
	    U.ref() = rhoU.ref() /rho.ref();
        e = rhoE/rho - 0.5*magSqr(U);
        rhoU.ref() = U.ref()*rho.ref();
        rhoE.ref() = e.ref()*rho.ref() + rho.ref()*0.5*magSqr(U.ref());
        rhoY.ref() = Y.ref()*rho.ref();

        surfaceScalarField ucSf( mag( fvc::interpolate(U) & mesh.Sf() )+fvc::interpolate(mag(c))*mesh.magSf() );

        forAll(ucSf.primitiveFieldRef(),iface)  {
        const label& leftCell = mesh.owner()[iface];
        ucSf.primitiveFieldRef()[iface]/= (mesh.V()[mesh.owner()[iface]]/mesh.nGeometricD());
        ucSf.primitiveFieldRef()[iface]= ucSf.primitiveFieldRef()[iface]*runTime.deltaTValue();
        CFL.ref()[leftCell]=max(CFL.ref()[leftCell],ucSf.primitiveFieldRef()[iface]);
        }

        U.correctBoundaryConditions();                                          //NADA
        Y.correctBoundaryConditions();
	    e.correctBoundaryConditions();
        rho.correctBoundaryConditions();
        rhoY.boundaryFieldRef()== rho.boundaryField()*Y.boundaryField();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
        rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

        #include "calculateThermo.H"                                            //COMPUTE THERMODYNAMIC PROPERTIES
        
        p1.correctBoundaryConditions();                                         //NADA
        c1.correctBoundaryConditions();
        
        rho1.ref() == rho.ref();                                                //NADA
        rhoU1.ref() == rhoU.ref();
        rhoE1.ref() == rhoE.ref();
        rhoY1.ref() == rhoY.ref();
        p1.ref() == p.ref();
        c1.ref() == c.ref();
        
		#include "flux.H"                                                        // COMPUTE FLUJOS EN FACES
		
        rhok[cycle+1]()  = fvc::div(rhoPhi);                                            // SUMMATORY FLUXES
        rhoYk[cycle+1]() = fvc::div(rhoYPhi);
		rhoUk[cycle+1]() = fvc::div(phiUp)-fvc::laplacian(muEff, U) - fvc::div(tauMC);
		rhoEk[cycle+1]() = fvc::div(phiEUp)- fvc::div(sigmaDotU)- fvc::laplacian(k_eff,T);

    }//forLoop

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"                                   //NADA
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    } //while runTime
    Info<< "End\n" << endl;                                                                               //NADA
    return 0;
}
