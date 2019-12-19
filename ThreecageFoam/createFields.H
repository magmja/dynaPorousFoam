Info << "Reading field p\n"
     << endl;
volScalarField p(
    IOobject(
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence(
    incompressible::turbulenceModel::New(U, phi, laminarTransport));

#include "createMRF.H"
#include "createFvOptions.H"

volScalarField nu(laminarTransport.nu());

Info << "Before read netdict "
     << "\n"
     << endl;
//
//- the net dictionary
IOdictionary netDict(
    IOobject(
        "netDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
netPanel porousZones(netDict);

IOdictionary structuralPositions(
    IOobject(
        "posi",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
//- the structural elements
IOdictionary structuralElements(
    IOobject(
        "surf",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

// constrcut net panel
porousZones.readPosi(structuralPositions);
porousZones.readSur(structuralElements);

Info << "Before read porosity field "
     << "\n"
     << endl;

dimensionedScalar poro(
    "poro",
    dimless,
    1.0);

volScalarField porosityField(
    IOobject(
        "porosityField",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    poro);
// Print information to pring and dubug
// Info << "After read porosity field " << porousZones.posi() << "\n"
//      << endl;
// Info << "After read porosity field " << porousZones.sur() << "\n"
//      << endl;
// Info << "After read porosity field " << porousZones.porosity() << "\n"
//      << endl;
// Info << "After read porosity field " << porousZones.F() << "\n"
//      << endl;
// // update the porosity field

porousZones.updatePoroField(porosityField, mesh);

//- write the porosity field
porosityField.write();
