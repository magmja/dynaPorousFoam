/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 
 No.

\*---------------------------------------------------------------------------*/

#include "netPanel.H"
#include "volFields.H"
#include "SortableList.H"
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::netPanel::calcArea(
    const point &pointI,
    const point &pointII,
    const point &pointIII) const
{
    const vector a(pointI - pointII);
    const vector b(pointI - pointIII);
    return 0.5 * mag(a ^ b);
}

bool Foam::netPanel::isInPorousZone(
    const point&            x,
    const List<vector>&     structuralPositions_memb,
    const vector&           structuralElementi) const
{
    bool result(false); // initial value

    const point pointI = structuralPositions_memb[structuralElementi[0]];
    const point pointII = structuralPositions_memb[structuralElementi[1]];
    const point pointIII = structuralPositions_memb[structuralElementi[2]];
    vector panelNorm = calcNorm(pointI, pointII, pointIII);  // a unit vector to indicate the normal 
    scalar dis(mag((x - pointI) & panelNorm));
    // define a const scalar as the distance between point x to net panel
    if (dis <= thickness_memb/2) // distance is less than half thickness
    {
        scalar panelarea(calcArea(pointI, pointII, pointIII));
        vector projectedPoint(0, 0, 0);       // initial the projected point is 0,0,0
        if (((x - pointI) & panelNorm) < 0.0) // on the side of normal vector
        {
            projectedPoint = (x + panelNorm * dis);
        }
        else
        {
            projectedPoint = (x - panelNorm * dis);
        }
        // projectedPiont is the projected point on the net panel
        scalar panelarea3(calcArea(pointI, pointII, projectedPoint) +
                          calcArea(pointI, projectedPoint, pointIII) + 
                          calcArea(projectedPoint, pointII, pointIII)
                         ); //  the area of the three trigular shapes.
        if (panelarea3 <= SMALL + panelarea*1.03)// safe factor
        {
            result = true;
        }
    }

    return result;
}
Foam::vector Foam::netPanel::calcNorm(
    const point &pointI,
    const point &pointII,
    const point &pointIII,
    const vector &fluidVelocity) const
{
    const vector a(pointI - pointII);
    const vector b(pointI - pointIII);
    vector norm(a ^ b / (mag(a ^ b) + SMALL));
    if ((norm & fluidVelocity) < 0) // assume the velocity is x+
    {
        norm = -norm;
    }
    return norm; // can point out or in the panel.
}

Foam::vector Foam::netPanel::calcNorm(
    const point &pointI,
    const point &pointII,
    const point &pointIII) const
{
    const vector a(pointI - pointII);
    const vector b(pointI - pointIII);
    const vector norm(a ^ b / (mag(a ^ b) + SMALL));
    return norm; // can point out or in the panel.
}

Foam::scalar Foam::netPanel::calcDist(
    const point &pointI,
    const point &pointII) const
{
    return mag(pointI - pointII);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::netPanel::netPanel(
    const dictionary &netDict)
    : // initial components
      netDict_memb(netDict),
      Sn_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("Sn"))),
      thickness_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("thickness"))),
      ML_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("meshLength"))),
      rho_fluid(readScalar(netDict_memb.subDict("NetInfo1").lookup("fluidDensity")))
{
    // creat the netpanel object
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::netPanel::~netPanel()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::netPanel::addResistance(
    fvVectorMatrix &UEqn,
    const fvMesh &mesh) const
{
    const vectorField &centres(mesh.C());
    const scalarField V = mesh.V();
    vectorField &Usource = UEqn.source();
    forAll(structuralElements_memb, Elementi)
    {
        point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
        point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
        point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);
        scalar area(calcArea(p0, p1, p2));
        vector sourceforce = structuralForces_memb[Elementi];
        forAll(centres, cellI)
        {
            if (
                isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                Usource[cellI] -= sourceforce /rho_fluid * V[cellI] / (thickness_memb * area + SMALL);
            }
        }
    }
}

void Foam::netPanel::updatePoroField(
    volScalarField &porosityField,
    const fvMesh &mesh) const
{
    // step1 set all the cell as 1
    forAll(mesh.C(), cellI)
    {
        porosityField[cellI] = 1.0; 
    }
    // get the center of all the cells
    const vectorField &centres(mesh.C());
    forAll(structuralElements_memb, Elementi) // loop through all the structural emlements
    {
        forAll(centres, cellI) // loop through all the cell,
        {
            if (isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                porosityField[cellI] = Sn_memb;
            }
        }
    }
}
//get velocit based on net panel element.
// defect: if the net panel is too small, the accuracy can be reduced dramatically.
// void Foam::netPanel::updateVelocity(
//     const fvMesh &mesh,
//     const volVectorField &U)
// {
//     List<vector> fluidVelocity(structuralElements_memb.size(), vector::zero);
//     // get the center of all the cells
//     const vectorField &centres(mesh.C());
//     forAll(structuralElements_memb, Elementi) // loop through all the structural emlements
//     {

//         vector fluidVelocityonElement(vector::zero);
//         scalar num_fvmesh(0);
//         forAll(centres, cellI) // loop through all the cell,
//         {
//             if (isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
//             {
//                 num_fvmesh += 1;
//                 fluidVelocityonElement += U[cellI];
//                 Info << "The velocity is \t " << fluidVelocityonElement << "\n" << endl;
//                 Info << "The number of mesh is\t " << num_fvmesh << "\n" << endl;
//                 // sum the velocity in each cell;
//             }
//         }
//         fluidVelocity[Elementi] = fluidVelocityonElement / (SMALL+num_fvmesh);
//     }
//     fluidVelocity_memb = fluidVelocity;
// }

Foam::List<Foam::vector> Foam::netPanel::updateVelocity(
    const fvMesh &mesh,
    const volVectorField &U)
    // const fvVectorMatrix       &UEqn)

{  // Get the velocity at the nearest cell center. 
    // get the center of all the cells
    List<vector> fluidVelocity_memb(structuralPositions_memb.size(), vector::zero);
    // Info <<"initial velocity is  "<<fluidVelocity_memb<<"\n"<<endl; 
    const vectorField &centres(mesh.C());

    // const vectorField &U = UEqn.psi(); // get the velocity field
    // Info <<"All the/ mesh position are"<< structuralPositions_memb<<endl;
    
    forAll(structuralPositions_memb, Pointi) // loop through all the structural emlements
    {
        scalar maxDistance(1);  //started from 2 m ML_memb
        vector nearestCell(vector::zero);
        scalar loops(0);
        forAll(centres, cellI) // loop through all the cell,
        {
            scalar k1(calcDist(centres[cellI], structuralPositions_memb[Pointi]));
            if (k1<maxDistance)
            {
                
                maxDistance=k1;
                fluidVelocity_memb[Pointi]=U[cellI];        
                nearestCell=centres[cellI];
                loops+=1;
            }
        }
    }
    return fluidVelocity_memb;
}

// * * * * * * * * * * * * * * Communication Functions  * * * * * * * * * * * * * * //
// - the following function is used to communicate with FE solver.
// - initial the structural element
void Foam::netPanel::readSurf(
    const dictionary &structuralElements)
{
    scalar listLength(readScalar(structuralElements.lookup("numOfSurf")));
    List<vector> surf(listLength, vector::zero);
    forAll(surf, i)
    {
        word surfname("e" + Foam::name(i));
        surf[i] = structuralElements.lookup(surfname);
    }
    structuralElements_memb = surf;
}

//- update the position and forces
void Foam::netPanel::readPosi(
    const dictionary &structuralPositions)
{
    scalar listLength(readScalar(structuralPositions.lookup("numOfPoint")));
    List<vector> posi(listLength, vector::zero);
    forAll(posi, i)
    {
        word pointname("p" + Foam::name(i));
        posi[i] = structuralPositions.lookup(pointname);
    }
    structuralPositions_memb = posi;
}

void Foam::netPanel::readForce(
    const dictionary &structuralForces)
{
    scalar listLength(readScalar(structuralForces.lookup("numOfFh")));
    List<vector> Fh(listLength, vector::zero);
    forAll(Fh, i)
    {
        word forcename("fh" + Foam::name(i));
        Fh[i] = structuralForces.lookup(forcename);
    }
    structuralForces_memb = Fh;
}

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// ************************************************************************* //
