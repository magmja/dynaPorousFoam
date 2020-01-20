// Solve the Momentum equation
const scalar porosity(Nettings.porosity());

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    1/porosity*fvm::ddt(U)
  + 1/porosity*fvm::div(phi, U)
  + MRF.DDt(U)
  + 1/porosity*turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

Nettings.addResistance(UEqn,mesh);

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
