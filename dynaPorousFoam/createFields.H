Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"


volScalarField nu(laminarTransport.nu());

Info<< "Before read netdict " << "\n" << endl;
// 
//- the net dictionary
IOdictionary netDict
(
    IOobject
    (
        "netDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

netPanel porousZones(netDict);

Info<< "the porosity of the net is " << porousZones.porosity()<< "\n"<<endl;

IOdictionary structuralPositions
(
    IOobject
    (
        "posi",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
porousZones.readPosi(structuralPositions);
//- the structural elements
IOdictionary structuralElements
(
    IOobject
    (
        "surc",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
porousZones.readSur(structuralElements);
// constrcut net panel
Info<< "the structuralPositions_memb of the net is " << porousZones.posi()<< "\n"<<endl;
Info<< "the structuralElements_memb of the net is " << porousZones.sur()<< "\n"<<endl;


Info<< "Before read porosity field " << "\n" << endl;
volScalarField porosityField
(
    IOobject
    (
        "porosityField",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    1.0
);
Info<< "After read porosity field " << "\n" << endl;
// read from outside 
// need to confirm ... might be wrong data structure
//- the structural point positions

porousZones.updatePoroField(porosityField, mesh);

//- write the porosity field
porosityField.write();
