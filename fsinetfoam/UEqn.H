// Solve the Momentum equation
const scalar porosity(Nettings.porosity());

MRF.correctBoundaryVelocity(U);
// todo: check the to U Eqn difference
 fvVectorMatrix UEqn
 (
     1/porosity*fvm::ddt(U)
   + 1/porosity*fvm::div(phi, U)
   + 1/porosity*turbulence->divDevReff(U)
 );


//fvVectorMatrix UEqn
//(
//   fvm::ddt(U)
//  +fvm::div(phi, U)
//  +turbulence->divDevReff(U)
//);

UEqn.relax();

Nettings.addResistance(UEqn,mesh);


fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
// Nettings.updateVelocity(mesh,UEqn);