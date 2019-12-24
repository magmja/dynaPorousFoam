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
        if (panelarea3 <= SMALL + panelarea*1.0)// safe factor
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
    return (a ^ b / (mag(a ^ b) + SMALL));
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

Foam::vector Foam::netPanel::calcLifti(
    const point &pointI,
    const point &pointII,
    const point &pointIII,
    const vector &fluidVelocity) const
{
    const vector eN(calcNorm(pointI, pointII, pointIII, fluidVelocity));
    return (fluidVelocity ^ eN ^ fluidVelocity) / (mag(fluidVelocity ^ eN ^ fluidVelocity) + SMALL);
}

Foam::scalar Foam::netPanel::calcTheta(
    const point &pointI,
    const point &pointII,
    const point &pointIII,
    const vector &fluidVelocity) const
{
    const vector eN(calcNorm(pointI, pointII, pointIII, fluidVelocity));
    return acos((eN & fluidVelocity) / (SMALL + mag(fluidVelocity)));
}
Foam::scalar Foam::netPanel::calcArea(
    const point &pointI,
    const point &pointII,
    const point &pointIII) const
{
    const vector a(pointI - pointII);
    const vector b(pointI - pointIII);
    return 0.5 * mag(a ^ b);
}

void Foam::netPanel::addResistance(
    fvVectorMatrix       &UEqn,
    const fvMesh         &mesh) const
{
    const scalarField V = mesh.V(); // volume of cells
    vectorField &Usource = UEqn.source();
    const vectorField &U = UEqn.psi(); // get the velocity field
    vector resistance_total(vector::zero);
    vector resistanceForce_total(vector::zero); 
    const vectorField &centres(mesh.C());// get the center of all the cells
    vector resistanceForce_Net(vector::zero);
    forAll(structuralElements_memb, Elementi)
    {
        resistanceForce_Net=vector::zero;
        point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
        point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
        point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);
        scalar area(calcArea(p0, p1, p2));
        scalar num_cell(0);
        forAll(centres, cellI)          // loop through all the cell,
        {
            if (isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                vector eL(calcLifti(p0, p1, p2, U[cellI]));
                scalar theta(calcTheta(p0, p1, p2, U[cellI]));
                vector Fd = 0.5 * (0.04 + F_memb.value()[0] * cos(theta)) * mag(U[cellI]) * (U[cellI]);    //* area
                vector Fl = 0.5 * F_memb.value()[1] * sin(2 * theta) * mag(U[cellI]) * mag(U[cellI]) * eL; //* area
                num_cell+=1;
                resistanceForce_Net+=(Fd + Fl)*1000.0;  // 1000 is the density of fluid          
                resistance_total+=(Fd + Fl) * V[cellI] / (thickness_memb * 0.95);
                Usource[cellI] -= (Fd + Fl) * V[cellI] / (thickness_memb * 0.95); //0.6 is safe factor
            }
        }
    resistanceForce_Net*=area/(SMALL+num_cell);
    }
    resistanceForce_total+=resistanceForce_Net;

    Info << "The total resistance source term on netting is  " << resistance_total << "\n" << endl;
    Info << "The total resistance force on netting is  " << resistanceForce_total << "\n" << endl;
}

void Foam::netPanel::updatePoroField(
    volScalarField       &porosityField,
    const fvMesh         &mesh) const
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
        forAll(centres, cellI)
        {
            if (isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                porosityField[cellI] = Sn_memb;
            }
        }
    }
}
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

void Foam::netPanel::readSur(
    const dictionary &structuralElements)
{
    scalar listLength(readScalar(structuralElements.lookup("numOfSurf")));
    List<vector> sur(listLength, vector::zero);
    forAll(sur, i)
    {
        word surname("e" + Foam::name(i));
        sur[i] = structuralElements.lookup(surname);
    }
    structuralElements_memb = sur;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::netPanel::netPanel(
    const dictionary &netDict)
    : // initial components
      netDict_memb(netDict),
      Sn_memb(readScalar(netDict_memb.subDict("porousProperties").lookup("Sn"))),
      thickness_memb(readScalar(netDict_memb.subDict("porousProperties").lookup("thickness"))),
      D_memb(netDict_memb.subDict("porousProperties").lookup("D")),
      F_memb(netDict_memb.subDict("porousProperties").lookup("F"))
{
    // creat the netpanel object
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::netPanel::~netPanel()
{
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// ************************************************************************* //
