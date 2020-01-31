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

//- the net dictionary
IOdictionary netDict(
    IOobject(
        "netDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

netPanel Nettings(netDict);

IOdictionary structuralPositions(
    IOobject(
        "posi",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
Nettings.readPosi(structuralPositions);

IOdictionary structuralElements(
    IOobject(
        "surf",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
Nettings.readSurf(structuralElements);


Info << "Solidity of the structure is " << Nettings.Sn() << "\n"
     << endl;

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
Nettings.updatePoroField(porosityField, mesh);
// Nettings.updateVelocity(mesh,U);
porosityField.write();
//
// generate the mesh center
scalar numberP = Pstream::nProcs(); //gives the n° of processors

List<pointField> gatheredCentres(numberP);
gatheredCentres[Pstream::myProcNo()] = pointField(mesh.C());
Pstream::gatherList(gatheredCentres);