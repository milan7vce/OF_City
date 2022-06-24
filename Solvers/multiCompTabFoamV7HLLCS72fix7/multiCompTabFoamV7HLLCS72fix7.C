
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentTransportModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "incompressibleTwoPhaseMixture.H"
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

    dimensionedScalar uBounded("uBounded", dimLength/dimTime, 1e-15);
    dimensionedScalar dt("dt", dimensionSet(0,0,1,0,0,0,0),1e-15);  

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
   
    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {

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

    rho_0 = rho;
    rhoY_0 = rhoY;
    rhoU_0 = rhoU;
    rhoE_0 = rhoE;

    for( int cycle = 0; cycle < 3; cycle++)
    {

       Info << "starting cycle " << RK4values.size()<< cycle << endl;
        
       dt=runTime.deltaT();
       drho=dt*(fvc::div(rhoPhi));
       drhoY=dt*(fvc::div(rhoYPhi));
       drhoU=dt*(fvc::div(phiUp)- fvc::laplacian(muEff, U) - fvc::div(tauMC));
       drhoE=dt*(fvc::div(phiEUp) - fvc::div(sigmaDotU)- fvc::laplacian(k_eff,T));
        
        switch(cycle) 
        {
        case 0:
	rho = rho -drho;
        rhoY = rhoY -drhoY;
        rhoU = rhoU -drhoU;
        rhoE = rhoE -drhoE;

        Info<< "\nZERO\n" << endl;
        break;

        case 1:
	rho =  (rho -drho)*0.25+rho_0*0.75;
        rhoY = (rhoY -drhoY)*0.25+rhoY_0*0.75;
        rhoU = (rhoU -drhoU)*0.25+rhoU_0*0.75;
        rhoE = (rhoE -drhoE)*0.25+rhoE_0*0.75;

        Info<< "\nUNO\n" << endl;
        break;

        case 2:
        rho =  2.0/3.0*(rho -drho)+1.0/3.0*rho_0;
        rhoY = 2.0/3.0*(rhoY -drhoY)+1.0/3.0*rhoY_0;
        rhoU=  2.0/3.0*(rhoU -drhoU)+1.0/3.0*rhoU_0;
        rhoE = 2.0/3.0*(rhoE -drhoE)+1.0/3.0*rhoE_0;

        Info<< "\nDOS\n" << endl;
        break;
        }

       
        Y.ref() = rhoY.ref() / rho.ref();                                        //NADA
        U.ref() = rhoU.ref() /rho.ref();
        e = rhoE/rho - 0.5*magSqr(U);

         forAll(rho, celli){                                                     //NADA
            rho.ref()[celli] = max(0.25,rho.ref()[celli]);
        }
       
        rhoU.ref() = U.ref()*rho.ref();
        rhoE.ref() = e.ref()*rho.ref() + rho.ref()*0.5*magSqr(U.ref());
        rhoY.ref() = Y.ref()*rho.ref();

        U.correctBoundaryConditions();                                          //NADA
        Y.correctBoundaryConditions();
	e.correctBoundaryConditions();
        rho.correctBoundaryConditions();
        rhoY.boundaryFieldRef()== rho.boundaryField()*Y.boundaryField();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
        rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

        #include "calculateThermo.H"                                            //COMPUTE THERMODYNAMIC PROPERTIES
        
        
        rho1.ref() == rho.ref();                                                //NADA
        rhoU1.ref() == rhoU.ref();
        rhoE1.ref() == rhoE.ref();
        rhoY1.ref() == rhoY.ref();
        p1.ref() == p.ref();
        c1.ref() == c.ref();

        p1.correctBoundaryConditions();                                         //NADA
        c1.correctBoundaryConditions();

        
	#include "flux.H"                                                        // COMPUTE FLUJOS EN FACE

    }//forLoop

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"                                   //NADA
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    } //while runTime
    Info<< "End\n" << endl;                                                                               //NADA
    return 0;
}
