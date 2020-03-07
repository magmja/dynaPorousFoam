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
    const point &x,
    const List<vector> &structuralPositions_memb,
    const vector &structuralElementi) const
{
    bool result(false); // initial value

    const point pointI = structuralPositions_memb[structuralElementi[0]];
    const point pointII = structuralPositions_memb[structuralElementi[1]];
    const point pointIII = structuralPositions_memb[structuralElementi[2]];
    vector panelNorm = calcNorm(pointI, pointII, pointIII); // a unit vector to indicate the normal
    scalar dis(mag((x - pointI) & panelNorm));
    // define a const scalar as the distance between point x to net panel
    if (dis <= thickness_memb / 2) // distance is less than half thickness
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
                          calcArea(projectedPoint, pointII, pointIII)); //  the area of the three trigular shapes.
        if (panelarea3 <= SMALL + panelarea * 1.0)                      // safe factor
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
    return norm; // must be in the same direction with fulow
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
      thickness_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("PorousMediaThickness"))),
      ML_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("halfMeshSize"))),
      fluidrho_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("fluidDensity"))),
      dw_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("twineDiameter"))),
      updateInterval_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("velocityUpdateInterval")))
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
//    Info << "In addResistance, number of mesh is " << centres.size() << endl;
    // Info << "The structural elements are " << structuralElements_memb << endl;
    // vector sourceforce = structuralForces_memb;
    forAll(structuralForces_memb, Elementi)
    {
        point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
        point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
        point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);
        scalar area(calcArea(p0, p1, p2));

        forAll(centres, cellI)
        {
            if (
                isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                Usource[cellI] -= structuralForces_memb[Elementi]* V[cellI]/fluidrho_memb / (thickness_memb * area + SMALL);
            }
        }
    }
}

void Foam::netPanel::updatePoroField(
    volScalarField &porosityField,
    const fvMesh &mesh) const
{
//    Info << "In updatePoroField, number of mesh is " << (mesh.C()).size() << endl;
    // Info << "The structural elements are " << structuralElements_memb << endl;
    // step1 set all the cell as 1
    forAll(mesh.C(), cellI)
    {
        porosityField[cellI] = 1.0;
    }
    // get the center of all the cells
    const vectorField &centres(mesh.C());
    //- step 2 assign sn to the proper mesh
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

// ge the velocity at the net panel center
void Foam::netPanel::updateVelocity(
    const List<pointField> &gatheredU,
    const List<pointField> &gathered_mesh,
    const scalar &thresholdLength,
    const scalar &time_foam)
{
    List<vector> fluidVelocities(structuralElements_memb.size(), vector::zero);
    int test1=time_foam/updateInterval_memb+SMALL;// It has add SMALL. No reason... otherwise it cannot work.
    float test2=time_foam/updateInterval_memb;

    if(float(test1)==float(test2) or time_foam<0.1)
    {
        Info<< " Update velocity at  = "<<time_foam <<endl;

        scalar maxDistance(1);                 //started from 2 m ML_memb
        forAll(structuralElements_memb, Elemi) // loop through all the structural emlements
        {
            point p0(structuralPositions_memb[structuralElements_memb[Elemi][0]]);
            point p1(structuralPositions_memb[structuralElements_memb[Elemi][1]]);
            point p2(structuralPositions_memb[structuralElements_memb[Elemi][2]]);
            point EP_center = (p0 + p1 + p2) / 3.0;
            maxDistance = thresholdLength*10; //started from 2 m ML_memb
            vector nearestCell(vector::zero);
            scalar loops(0);
            forAll(gathered_mesh, processorI) // loop through all the cell,
            {
                if (maxDistance<thresholdLength)
                {
                    break;
                }
                forAll(gathered_mesh[processorI], PointI)
                {
                    scalar k1(calcDist(gathered_mesh[processorI][PointI], EP_center));
                    if (k1 < maxDistance)
                    {
                        maxDistance = k1;
                        fluidVelocities[Elemi] = gatheredU[processorI][PointI];
                        nearestCell = gathered_mesh[processorI][PointI];
                        loops += 1;
                    }
                    if (maxDistance<thresholdLength)
                    {
                        break;
                    }
                }

            }
//            Info << "After " << loops << " times of loop, the nearest cell is " << nearestCell << "to point " << EP_center <<", and the velocity is "<<fluidVelocities[Elemi]<< "\n"
//                 << endl;
            if (maxDistance >= 0.8)
            {
                Info << "Warnning!!! I cannot find the nearest cell to point " << EP_center << " , because the minimum distance to this point is  " << maxDistance << "\n"
                     << endl;
            }
        }
        fluidVelocity_memb = fluidVelocities; // only assige onece
        // Info << "the velocity on elements are  " << fluidVelocity_memb << endl;
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
        word surf_name("e" + Foam::name(i));
        surf[i] = structuralElements.lookup(surf_name);
    }
    structuralElements_memb = surf;
    // Info << "The structural elements are " << structuralElements_memb << endl;
}

//- update the position and forces
void Foam::netPanel::readPosi(
    const dictionary &structuralPositions)
{
    scalar listLength(readScalar(structuralPositions.lookup("numOfPoint")));
    List<vector> posi(listLength, vector::zero);
    forAll(posi, i)
    {
        word point_name("p" + Foam::name(i));
        posi[i] = structuralPositions.lookup(point_name);
    }
    structuralPositions_memb = posi;
}

void Foam::netPanel::readForce(
    const scalar &time_foam,
    const dictionary &structuralForces)
{
    scalar listLength(readScalar(structuralForces.lookup("numOfFh")));
    scalar time_FE(readScalar(structuralForces.lookup("timeInFE")));
    if (mag(time_FE-time_foam)>10)
    {
        Info<<"Warning!!! The difference of time in FE and FV solvers exceeds 10 s!\n" <<endl;
    }

    List<vector> Fh(listLength, vector::zero);
    forAll(Fh, i)
    {
        word force_name("fh" + Foam::name(i));
        Fh[i] = structuralForces.lookup(force_name);
    }
    structuralForces_memb = Fh;
}

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// ************************************************************************* //
