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

IOdictionary structuralFh(
    IOobject(
        "Fh",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

netPanel Nettings(netDict);

Info << "Solidity of the structure is " << Nettings.Sn() << "\n"
     << endl;
// constrcut net panel
Nettings.readPosi(structuralPositions);
Nettings.readSurf(structuralElements);
Nettings.readForce(structuralFh);

dimensionedScalar poro(
    "poro",
    dimless,
    1.0);
//
volScalarField porosityField(
    IOobject(
        "porosityField",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    poro);

// update the porosity field
Nettings.updatePoroField(porosityField, mesh,U);
//- write the porosity field
porosityField.write();
//