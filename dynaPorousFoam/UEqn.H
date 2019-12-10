// Solve the Momentum equation
const scalar porosity(porousZones.porosity());

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    1/porosity*fvm::ddt(U)
  + 1/porosity*fvm::div(phi, U)
  + 1/porosity*turbulence->divDevReff(U)
);

UEqn.relax();

porousZones.addResistance(UEqn,nu,mesh);


fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
