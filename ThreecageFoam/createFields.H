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
netPanel porousZones0(netDict);
netPanel porousZones1(netDict);
netPanel porousZones2(netDict);


// >>>>>>>>>>>>>>>cage0>>>>>>>>>>>>>>>>>>
IOdictionary structuralPositions(
    IOobject(
        "posi0",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
//- the structural elements
IOdictionary structuralElements(
    IOobject(
        "surf0",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

dimensionedScalar poro(
    "poro",
    dimless,
    1.0);

volScalarField porosityField(
    IOobject(
        "porosityField0",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    poro);

porousZones0.readPosi(structuralPositions);
porousZones0.readSur(structuralElements);
porousZones0.updatePoroField(porosityField, mesh);

// >>>>>>>>>>>>>>> cage1 >>>>>>>>>>>>>>>>>>
IOdictionary structuralPositions(
    IOobject(
        "posi1",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
//- the structural elements
IOdictionary structuralElements(
    IOobject(
        "surf1",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

dimensionedScalar poro(
    "poro",
    dimless,
    1.0);

volScalarField porosityField(
    IOobject(
        "porosityField1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    poro);
porousZones1.readPosi(structuralPositions);
porousZones1.readSur(structuralElements);
porousZones1.updatePoroField(porosityField, mesh);




// >>>>>>>>>>>>>>> cage2 >>>>>>>>>>>>>>>>>>
IOdictionary structuralPositions(
    IOobject(
        "posi2",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
//- the structural elements
IOdictionary structuralElements(
    IOobject(
        "surf2",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

dimensionedScalar poro(
    "poro",
    dimless,
    1.0);

volScalarField porosityField(
    IOobject(
        "porosityField2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    poro);

porousZones2.readPosi(structuralPositions);
porousZones2.readSur(structuralElements);
porousZones2.updatePoroField(porosityField, mesh);

//- write the porosity field
porosityField.write();
