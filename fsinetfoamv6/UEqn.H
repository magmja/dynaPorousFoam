// Solve the Momentum equation
const scalar porosity(Nettings.porosity());

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    1/porosity*fvm::ddt(U) + fvm::div(phi, U)
  + 1/porosity*MRF.DDt(U)
  + 1/porosity*turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn.relax();

Nettings.addResistance(UEqn,mesh);

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
