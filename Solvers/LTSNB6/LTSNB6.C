
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
//#include "fixedRhoFvPatchScalarField.H"
//#include "incompressibleTwoPhaseMixture.H"
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
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0); 
    
    OFstream resFile("residual.txt");

    word nameOfTable(" ");
    runTime.controlDict().readIfPresent("nameOfTable", nameOfTable); //read file from controlDict
    tableInfo unif(nameOfTable); //call constructor
    
    #include "calculateThermo.H"

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    //scalar CoNumMin = 0.0;

    volScalarField rho_0("rho_0",rho);                                      //NADA
    volScalarField rhoY_0("rhoY_0",rhoY);
    volVectorField rhoU_0("rhoU_0",rhoU);
    volScalarField rhoE_0("rhoE_0",rhoE);
    //volScalarField maxpercell_ucSf("maxpercell_ucSf",mag(U));
   
    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {

    rho1.ref() = rho.ref();                                                 //NADA
    rhoU1.ref()= rhoU.ref();
    rhoE1.ref() = rhoE.ref();
    rhoY1.ref() = rhoY.ref();
    p1.ref() = p_r.ref();
    c1.ref() = c.ref();

    U.correctBoundaryConditions();                                           //NADA
    Y.correctBoundaryConditions();
    rho.correctBoundaryConditions();
    e.correctBoundaryConditions();
    rhoY.boundaryFieldRef() == rho.boundaryField()*Y.boundaryField();
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
    rhoE.boundaryFieldRef() == rho.boundaryField()*( e.boundaryField() + 0.5*magSqr(U.boundaryField()));
    p_r.correctBoundaryConditions();
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
        
        switch(cycle) 
        {
        case 0:
        solve(fvm::ddt(rho)  == -fvc::div(rhoPhi));
        solve(fvm::ddt(rhoY) == -fvc::div(rhoYPhi));
        solve(fvm::ddt(rhoU) == -fvc::div(phiUp));
        solve(fvm::ddt(rhoE) == -fvc::div(phiEUp)+ fvc::div(sigmaDotU) + fvc::laplacian(k_eff,T_r));	
	#include "corrections.H"
	U.ref() = rhoU.ref() /rho.ref();
	
	forAll(UBf, pachi)
        {
        if((UBf[pachi].patch().name() == "nozzle1") || (UBf[pachi].patch().name() == "nozzlle2") || (UBf[pachi].patch().name() == "needle"))
        {
        forAll(UBf[pachi], facei)
                  {
                        U.boundaryFieldRef()[pachi][facei] = U.boundaryFieldRef()[pachi][facei]*0.0;
                  }
        }
        }

	solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );

        rhoU = rho*U;

        Info<< "\nZERO\n" << endl;
        break;

        case 1:
        solve(fvm::ddt(rho)  == -fvc::div(rhoPhi));
        solve(fvm::ddt(rhoY) == -fvc::div(rhoYPhi));
        solve(fvm::ddt(rhoU) == -fvc::div(phiUp));
        solve(fvm::ddt(rhoE) == -fvc::div(phiEUp)+ fvc::div(sigmaDotU) + fvc::laplacian(k_eff,T_r));
	#include "corrections.H"
	U.ref() = rhoU.ref() /rho.ref();
	U.correctBoundaryConditions();

	forAll(UBf, pachi)
        {
           if((UBf[pachi].patch().name() == "nozzle1") || (UBf[pachi].patch().name() == "nozzlle2") || (UBf[pachi].patch().name() == "needle"))
         {
              forAll(UBf[pachi], facei)
                  {
                        U.boundaryFieldRef()[pachi][facei] = U.boundaryFieldRef()[pachi][facei]*0.0;
                  }
        }
        }

	
	solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );

        rhoU = rho*U;

	rho =  rho*0.25+rho_0*0.75;
        rhoY = rhoY*0.25+rhoY_0*0.75;
        rhoU = rhoU*0.25+rhoU_0*0.75;
        rhoE = rhoE*0.25+rhoE_0*0.75;

        Info<< "\nUNO\n" << endl;
        break;

        case 2:
	solve(fvm::ddt(rho)  == -fvc::div(rhoPhi));
        solve(fvm::ddt(rhoY) == -fvc::div(rhoYPhi));
        solve(fvm::ddt(rhoU) == -fvc::div(phiUp));
        solve(fvm::ddt(rhoE) == -fvc::div(phiEUp)+ fvc::div(sigmaDotU) + fvc::laplacian(k_eff,T_r));
	#include "corrections.H"
	U.ref() = rhoU.ref() /rho.ref();
	        forAll(UBf, pachi)
        {
           if((UBf[pachi].patch().name() == "nozzle1") || (UBf[pachi].patch().name() == "nozzlle2") || (UBf[pachi].patch().name() == "needle"))
         {
              forAll(UBf[pachi], facei)
                  {
                        U.boundaryFieldRef()[pachi][facei] = U.boundaryFieldRef()[pachi][facei]*0.0;
                  }
        }
        }

	solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
        rhoU = rho*U;

        rho =  2.0/3.0*rho+1.0/3.0*rho_0;
        rhoY = 2.0/3.0*rhoY+1.0/3.0*rhoY_0;
        rhoU=  2.0/3.0*rhoU+1.0/3.0*rhoU_0;
        rhoE = 2.0/3.0*rhoE+1.0/3.0*rhoE_0;
        Info<< "\nDOS\n" << endl;
        break;

	case 3:
	           forAll(UBf, pachi)
        {
           if((UBf[pachi].patch().name() == "nozzle1") || (UBf[pachi].patch().name() == "nozzlle2") || (UBf[pachi].patch().name() == "needle"))
         {
              forAll(UBf[pachi], facei)
                  {
                        U.boundaryFieldRef()[pachi][facei] = U.boundaryFieldRef()[pachi][facei]*0.0;
                  }
        }
        }

        solve(fvm::ddt(rhoU) == fvc::div(tauMC)+ fvc::laplacian(muEff, U));
        solve(fvm::ddt(rhoE) == fvc::div(sigmaDotU) + fvc::laplacian(k_eff,T_r));
	break;
        }

        Y.ref() = rhoY.ref() / rho.ref();                                        //NADA
        U.ref() = rhoU.ref() /rho.ref();
        e = rhoE/rho - 0.5*magSqr(U);

        forAll(rho, celli)
        {                                                     //NADA
            rho.ref()[celli] = max(0.25,rho.ref()[celli]);
        }
       
        rhoU.ref() = U.ref()*rho.ref();
        rhoE.ref() = e.ref()*rho.ref() + rho.ref()*0.5*magSqr(U.ref());
        rhoY.ref() = Y.ref()*rho.ref();

        U.correctBoundaryConditions();
        Y.correctBoundaryConditions();
	    e.correctBoundaryConditions();
        rho.correctBoundaryConditions();
        rhoY.boundaryFieldRef()== rho.boundaryField()*Y.boundaryField();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
        rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));
        
        #include "calculateThermo.H"                //COMPUTE THERMODYNAMIC PROPERTIES
        turbulence->correct();
        
        rho1.ref() == rho.ref();                                                //NADA
        rhoU1.ref() == rhoU.ref();
        rhoE1.ref() == rhoE.ref();
        rhoY1.ref() == rhoY.ref();
        p1.ref() == p_r.ref();
        c1.ref() == c.ref();

        p1.correctBoundaryConditions();                                         //NADA
        c1.correctBoundaryConditions();

	    #include "flux.H"                                                        // COMPUTE FLUJOS EN FACE
    }
    
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"                                   //NADA
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    } //while runTime
    Info<< "End\n" << endl;                                                                               //NADA
    return 0;
}
