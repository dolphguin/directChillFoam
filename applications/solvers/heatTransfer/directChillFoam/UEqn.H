    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

//TODO: this should be a runtime-selectable model instead
{
    const volScalarField& mu (thermo.mu());

    const volScalarField& alpha (alphafvModel);

    const volScalarField gl
    (
        rho2
     / (rho1*(1.0/max(alpha, small) - 1.0) + rho2)
      * pos0(alpha - g_env) // else -> gl = 0.0
    );

    Us = U + gl*(rho1 - rho2)*g*sqr(DAS)/(18*mu);
}

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevTau(U)
     ==
        fvModels.source(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvConstraints.constrain(U);
        K = 0.5*magSqr(U);
    }
