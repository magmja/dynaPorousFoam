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
netPanel porousZones(netDict);

Info<< "Porosity of the structure is " << porousZones.porosity()<< "\n"<<endl;
// constrcut net panel
porousZones.readPosi(structuralPositions);
porousZones.readSur(structuralElements);

// Foam::fvMesh mesh
// (
//     Foam::IOobject
//     (
//         Foam::fvMesh::defaultRegion,
//         runTime.timeName(),
//         runTime,
//         IOobject::MUST_READ
//     )
// );

// Info<< "Cell sentres  " << "\n"<<mesh.cellCentres() << endl;
// Info<< "Cell volumes  " << "\n"<<mesh.cellVolumes() << endl;

// volScalarField pp2
// (
//     IOobject
//     (
//         "pp2",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );
// Info<< pp2 << endl;
// Info<< pp2.boundaryField()[0] << endl;


// Info<< "Before read porosity field " << "\n" << endl;
dimensionedScalar poro
(
   "poro",
   dimless,
   1.0
);

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
    poro
);

// Info<< "After read porosity field " << "\n" << endl;

// update the porosity field

porousZones.updatePoroField(porosityField, mesh);

//- write the porosity field
porosityField.write();
