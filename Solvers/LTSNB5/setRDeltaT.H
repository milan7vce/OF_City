{
    volScalarField& rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff(100.0);
    /*    runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );*/

    rDeltaT.ref() = max
    (
       1/dimensionedScalar(dimTime, maxDeltaT),
        1/maxpercell_ucSf
    );

    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    
    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    trDeltaT3=trDeltaT.ref();
    tiempopaso=1.0/trDeltaT3;
}
