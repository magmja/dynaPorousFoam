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
void Foam::netPanel::transformCoeffs(
    const dimensionedVector &D,
    tensor &d,
    const List<vector> &structuralPositions_memb,
    const vector &structuralElementi

    ) const
{
    const point pointI = structuralPositions_memb[structuralElementi[0]];
    const point pointII = structuralPositions_memb[structuralElementi[1]];
    const point pointIII = structuralPositions_memb[structuralElementi[2]];

    const vector a(pointI - pointII);
    const vector b(pointIII - pointII);
    const vector c(pointIII - pointI);
    vector norm(a ^ b / (mag(a ^ b) + SMALL));
    const vector x0(1, 0, 0);
    const vector y0(0, 1, 0);
    const vector z0(0, 0, 1);

    vector e1(vector::zero); //r
    vector e2(vector::zero); //g
    vector e3(vector::zero); //b

    if ((norm & x0) < 0) // assume the velocity is x+
    {
        norm = -norm;
    }
    e1 = norm;

    if (mag(a) > mag(b))
    {
        e3 = b / (mag(b) + SMALL);
        if (mag(b) > mag(c))
        {
            e3 = c / (mag(c) + SMALL);
        }
    }
    else
    {
        e3 = a / (mag(a) + SMALL);
    }

    if ((e3 & x0) < 0) // assume the velocity is x+
    {
        e3 = -e3;
    }

    e2 = (e3 ^ e1 / (mag(e3 ^ e1) + SMALL)); // blue
    if ((e2 & z0) < 0)                       // assume the velocity is x+
    {
        e2 = -e2;
    }

    // Info << "p1 is : " << pointI << "\n" << endl;
    // Info << "p2 is : " << pointII << "\n" << endl;
    // Info << "p3 is : " << pointIII << "\n" << endl;

    // Info << "e1 is : " << e1  << "\n" << endl;
    // Info << "e2 is : " << e2  << "\n" << endl;
    // Info << "e3 is : " << e3  << "\n" << endl;
    // Info << "Tensor e is : " << e.value() << "\n" << endl;
    const tensor e(
        e1 & x0 / (mag(e1 * x0) + SMALL), e1 & y0 / (mag(e1 * y0) + SMALL), e1 & z0 / (mag(e1 * z0) + SMALL),
        e2 & x0 / (mag(e2 * x0) + SMALL), e2 & y0 / (mag(e2 * y0) + SMALL), e2 & z0 / (mag(e2 * z0) + SMALL),
        e3 & x0 / (mag(e3 * x0) + SMALL), e3 & y0 / (mag(e3 * y0) + SMALL), e3 & z0 / (mag(e3 * z0) + SMALL));

    d.xx() = D.value().x();
    d.yy() = D.value().y();
    d.zz() = D.value().z();
    d = e & d & e.T();
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
    vector panelNorm = calcNorm2(pointI, pointII, pointIII);
    scalar dis(mag((x - pointI) & panelNorm));
    // define a const scalar as the distance between point x to net panel
    if (dis <= thickness_memb) // distance is less than half thickness
    {
        scalar panelarea // the area of the triangular net panel
            (
                0.5 * mag(
                          (pointI - pointII) ^ (pointI - pointIII)));
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
        scalar panelarea3(0.5 * (mag((projectedPoint - pointI) ^ (projectedPoint - pointII)) +
                                 mag((projectedPoint - pointI) ^ (projectedPoint - pointIII)) +
                                 mag((projectedPoint - pointII) ^ (projectedPoint - pointIII)))); //  the area of the three trigular shapes.
        if (panelarea3 <= SMALL + panelarea)
        {
            result = true;
        }
    }

    return result;
}
Foam::vector Foam::netPanel::calcNorm2(
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
    return acos(mag(eN & fluidVelocity) / (SMALL + mag(fluidVelocity)));
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
    fvVectorMatrix &UEqn,
    const volScalarField &nu,
    const fvMesh &mesh) const
{
    const vectorField &centres(mesh.C());
    const scalarField V = mesh.V(); // volume of cells
    vectorField &Usource = UEqn.source();
    const vectorField &U = UEqn.psi(); // get the velocity field

    forAll(structuralElements_memb, Elementi)
    {
        // Info << "before add the source term, the f is " << f_global << "\n" << endl;
        // Info << "before add the source term, the F_memb is " << F_memb << "\n" << endl;
        forAll(centres, cellI)
        {
            if (
                isInPorousZone(centres[cellI], structuralPositions_memb, structuralElements_memb[Elementi]))
            {
                point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
                point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
                point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);
                vector eL(calcLifti(p0, p1, p2, U[cellI]));
                scalar theta(calcTheta(p0, p1, p2, U[cellI]));
                vector Fd = 0.5 * (0.04+F_memb.value()[0] * cos(theta)) * mag(U[cellI]) * (U[cellI]); //* area
                vector Fl = 0.5 * F_memb.value()[1] * sin(2 * theta) * mag(U[cellI]) * mag(U[cellI]) * eL;                               //* area
                // tensor dragCoeff = nu[cellI] * d_global + 0.5 * mag(U[cellI]) * f_global;
                // Usource[cellI] -= V[cellI] * dragCoeff & (U[cellI]);
                Usource[cellI] -= (Fd + Fl) / (thickness_memb / (SMALL + V[cellI])); //* area
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
      thickness_memb(readScalar(netDict_memb.subDict("porousProperties").lookup("halfthickness"))),
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
