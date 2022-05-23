
Info << "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField p1
(
    IOobject
    (
        "p1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*
surfaceScalarField p_limit
(
    IOobject
    (
        "p_limit",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("p_limit", dimPressure , 2500000.0)
);


volScalarField pCollapse
(
    IOobject
    (
        "pCollapse",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p*0
);
*/

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("T", dimTemperature , 0.),
"zeroGradient"
);


volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("c", dimLength/dimTime , 0.),
"zeroGradient"
);

volScalarField c1
(
    IOobject
    (
        "c1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("c1", dimLength/dimTime , 0.),
"zeroGradient"
);

/*surfaceScalarField c_limit
(
    IOobject
    (
        "c_limit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("c_limit", dimLength/dimTime , 400.)
);


volScalarField alphav
(
    IOobject
    (
        "alphav",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/

volScalarField CFL
(
    IOobject
    (
        "CFL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0.0)
);


/*volScalarField alphav_limit
(
    IOobject
    (
        "alphav_limit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
       mesh

);
*/

//#include "rhoBoundaryTypes.H"
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho1
(
    IOobject
    (
        "rho1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/*volScalarField e1
(
    IOobject
    (
        "e1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);*/

volScalarField e
(
    IOobject
    (
        "e",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

volScalarField rhoE1
(
    IOobject
    (
        "rhoE1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);


volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volVectorField rhoU1
(
 IOobject
(
 "rhoU1",
 runTime.timeName(),
 mesh,
 IOobject::NO_READ,
 IOobject::NO_WRITE
 ),
 rho*U
 );

volScalarField Y
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),

    mesh
);

/*volScalarField divide
(
    IOobject
    (
        "divide",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),

    mesh
);*/

surfaceScalarField x_limit
(
    IOobject
    (
        "x_limit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),

    mesh,
    dimensionedScalar(dimless, 0.0)
);


/*surfaceScalarField Y_limit
(
    IOobject
    (
        "Y_limit",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),

    mesh,
    dimensionedScalar(dimless, 0.9999)

);

volScalarField rr
(
    IOobject
    (
        "rr",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/

volScalarField rhoY
(
    IOobject
    (
        "rhoY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Y
);

volScalarField rhoY1
(
    IOobject
    (
        "rhoY1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Y
);


surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless, -1.0)
);


volScalarField Nu //kinematic viscosity calculated from interpolation table
(
    IOobject
    (
        "Nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("Nu", dimLength*dimLength/dimTime , 0.),
"zeroGradient"
);


volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
dimensionedScalar("Cp", (dimLength*dimLength)/(dimTime*dimTime*dimTemperature) , 0.),
"zeroGradient"
);


volScalarField k_eff
(
    IOobject
    (
        "k_eff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
dimensionedScalar("k_eff", (dimMass*dimLength)/(dimTime*dimTime*dimTime*dimTemperature) , 0.),
"zeroGradient"
);

Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);

//incompressibleTwoPhaseMixture mixture(U, phi);
twoPhaseMixtureThermo mixture(U, phi);

volScalarField& alphav(mixture.alpha1());
//alphav.oldTime();
 
volScalarField& alphal(mixture.alpha2());

surfaceScalarField rhoPhi("rhoPhi", mesh.Sf() & fvc::interpolate(rhoU));

surfaceScalarField rhoYPhi("rhoYPhi", mesh.Sf() & fvc::interpolate(rhoU)*fvc::interpolate(Y));

/*surfaceScalarField rhoPhi_MFR
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoPhi
);*/
/*surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //dimensionedScalar("rhoPhi", dimensionSet(1, 0, -1, 0, 0, 0, 0), 1.0)
    //fvc::dotInterpolate(mesh.Sf(), rhoU)
    //fvc::flux(rhoU)
    mesh.Sf() & fvc::interpolate(rhoU))
);*/

// create the k values for the fields

	scalarList RK4values(0);
	scalarList RK4values2(0);

//coefficients of time advancement
	RK4values.append(0);
	RK4values.append(0.5);
	RK4values.append(0.5);
	RK4values.append(1);

//coefficients of intermediate stages
	RK4values2.append(0.16666666666);
	RK4values2.append(0.33333333333);
	RK4values2.append(0.33333333333);
	RK4values2.append(0.16666666666);

		Info << RK4values2 << endl;

		Info << RK4values << endl;


    List<autoPtr<volScalarField> >  rhok(RK4values.size());
    List<autoPtr<volScalarField> >  rhoYk(RK4values.size());
    List<autoPtr<volVectorField> >  rhoUk(RK4values.size());
    List<autoPtr<volScalarField> >  rhoEk(RK4values.size());

forAll(rhok, count)
{

    word namerhokName ("rhok" + name(count));

    rhok[count].reset
    (
        new volScalarField 
	(
            IOobject
            (
                namerhokName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(rhoPhi)
        )
    );

    word namerhoYkName ("rhoYk" + name(count));

    rhoYk[count].reset
    (
        new volScalarField 
	(
            IOobject
            (
                namerhoYkName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(rhoYPhi)
        )
    );

    word namerhoUkName ("rhoUk" + name(count));

    rhoUk[count].reset
    (
        new volVectorField 
	(
            IOobject
            (
                namerhoUkName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::grad(p) // for init
        )
    );

    word namerhoEkName ("rhoEk" + name(count));

    rhoEk[count].reset
    (
        new volScalarField 
    	(
            IOobject
            (
                namerhoEkName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::ddt(rho, e)
        )
    );


}


autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        rhoPhi,
        mixture
    )
);

