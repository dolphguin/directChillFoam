//WARNING: solutal convergence generally not guaranteed!

for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    const volScalarField alpha0(alphafvModel);

    #include "EEqn.H"

    fvModels.correct();

    rho = thermo.rho();

    alloy.solve(Us, fvModels, fvConstraints, alphafvModel);

    const scalarField residualField
    (
        mag
        (
            alphafvModel.internalField()
          - alpha0.internalField()
        )
    );

    const scalar residual = gMax(residualField);

    Info<< "Energy iteration " << (Ecorr + 1)
        << " melt fraction residual: " << residual << endl;

    if (residual < 1e-8) break;
}
