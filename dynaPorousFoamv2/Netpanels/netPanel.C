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
//- transform coefficients
Foam::scalar Foam::netPanel::calcArea(
    const point &pointI,
    const point &pointII,
    const point &pointIII) const
{ // the area of the triangular net panel
    scalar panelarea(
        0.5 * mag(
                  (pointI - pointII) ^ (pointI - pointIII)));
    return panelarea;
}

bool Foam::netPanel::isInPorousZone(
    const point x,
    const List<vector> &structuralPositions_memb,
    const vector &structuralElementi) const
{
    bool result(false); // initial value

    const point pointI = structuralPositions_memb[structuralElementi[0]];
    const point pointII = structuralPositions_memb[structuralElementi[1]];
    const point pointIII = structuralPositions_memb[structuralElementi[2]];
    vector panelNorm = calcNorm(pointI, pointII, pointIII);
    scalar panelarea = calcArea(pointI, pointII, pointIII);
    scalar dis(mag((x - pointI) & panelNorm));
    // define a const scalar as the distance between point x to net panel
    if (dis <= thickness_memb / 2) // distance is less than half thickness
    {
        vector projectedPoint(0, 0, 0);       // initial the projected point is 0,0,0
        if (((x - pointI) & panelNorm) < 0.0) // on the side of normal vector
        {
            projectedPoint = (x + panelNorm * dis);
        }
        else
        {
            projectedPoint = (x - panelNorm * dis);
        }
        // projectedP is the projected point on the net panel
        scalar triSma0 = calcArea(projectedPoint, pointII, pointIII);
        scalar triSma1 = calcArea(pointI, projectedPoint, pointIII);
        scalar triSma2 = calcArea(pointI, pointII, projectedPoint);
        //  the area of the three trigular shapes.

        if ((triSma0 + triSma1 + triSma2) <= SMALL + panelarea)
        {
            result = true;
        }
    }

    return result;
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::netPanel::netPanel(
    const dictionary &netDict)
    : // initial components
      netDict_memb(netDict),
      Sn_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("Sn"))),
      thickness_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("thickness")))
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
    const volScalarField &nu,
    const fvMesh &mesh) const
{
    const vectorField &centres(mesh.C());
    const scalarField V = mesh.V();
    vectorField &Usource = UEqn.source();
    scalar volumeE(0.0);
    forAll(structuralElements_memb, Elementi)
    {
        // scalar num_cell = 0;
        volumeE = 0;
        // point pointI = structuralPositions_memb[structuralElements_memb[Elementi],0];
        // point pointII = structuralPositions_memb[structuralElements_memb[Elementi],1];
        // point pointIII = structuralPositions_memb[structuralElements_memb[Elementi],2];
        // // scalar panelarea = calcArea(pointI, pointII, pointIII);
        vector sourceforce = structuralForces_memb[Elementi];

        // const vectorField &U = UEqn.psi(); // get the velocity field
        // Info << "before add the source term, the f is " << f_global << "\n" << endl;
        // Info << "before add the source term, the F_memb is " << F_memb << "\n" << endl;
        forAll(centres, cellI)
        {
            if (
                isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                volumeE += V[cellI];
                // num_cell += 1;
            }

            if (
                isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                Usource[cellI] -= V[cellI] / volumeE * sourceforce;
            }
        }
    }
}

void Foam::netPanel::updatePoroField(
    volScalarField &porosityField,
    fvMesh &mesh) const
{
    // step1 set all the cell as 1
    forAll(mesh.C(), cellI)
    {
        porosityField[cellI] = 1.0;
    }
    // get the center of all the cells
    const vectorField &centres(mesh.C());
    // loop through all the structural emlements
    forAll(structuralElements_memb, Elementi)
    {
        // loop through all the cell,
        vector fluidVelocityonElement(vector::zero);
        scalar num_fvmesh(0);
        forAll(centres, cellI)
        {
            if (isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                porosityField[cellI] = Sn_memb;
                num_fvmesh += 1;
                // sum the velocity in each cell;
            }
        }
        fluidVelocityonElement=fluidVelocityonElement/num_fvmesh;
        fluidVelocity_memb[Elementi]=fluidVelocityonElement;
    }
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
