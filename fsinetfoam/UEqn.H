// Solve the Momentum equation
MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
     fvm::ddt(U)
    +fvm::div(phi, U)
    +turbulence->divDevReff(U)
);

UEqn.relax();

Nettings.addResistance(UEqn,nu,mesh);


fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
