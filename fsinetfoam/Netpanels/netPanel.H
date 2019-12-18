/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 No.
Class
    Foam::netPanel

Description

SourceFiles
    netPanel.H
    netPanel.C

\*---------------------------------------------------------------------------*/
#ifndef netPanel_H
#define netPanel_H
// class guards that prevents double declaration

#include "dictionary.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "scalar.H"
#include "fvMatrices.H"

// might be used
// #include "pointField.H"
// #include "volFields.H"
// #include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class netPanel Declaration
\*---------------------------------------------------------------------------*/

class netPanel
{
protected:
    //- Private data
    const dictionary netDict_memb;
    const scalar thickness_memb; // thickness of the net panel
    const scalar Sn_memb;
    List<vector> structuralPositions_memb; // positions of each noeds
    List<vector> structuralElements_memb;
    List<vector> structuralForces_memb;
    List<vector> fluidVelocity_memb;

    //-Private member function
    vector calcNorm(
        const point &pointI,
        const point &pointII,
        const point &pointIII) const;
    scalar calcArea(
        const point &pointI,
        const point &pointII,
        const point &pointIII) const;
    bool isInPorousZone(
        const point x,
        const List<vector> &structuralPositions,
        const vector &structuralElementi) const;

public:
    // Constructors
    netPanel(
        const dictionary &netDict);

    //- Destructor
    ~netPanel();

    // public Member Functions
    void addResistance(
        fvVectorMatrix &UEqn,
        const volScalarField &nu,
        const fvMesh &mesh) const;

    void updatePoroField(
        volScalarField &porosityField,
        const fvMesh &mesh) const;

    void updateVelocity(
        const fvMesh &mesh,
        const volVectorField &U);



    scalar Sn() const
    {
        return Sn_memb;
    };

    List<vector> posi() const
    {
        return structuralPositions_memb;
    };

    List<vector> surf() const
    {
        return structuralElements_memb;
    };
    List<vector> FluidU() const
    {
        return fluidVelocity_memb;
    };
    List<vector> Fhout() const
    {
        return structuralForces_memb;
    };

    void readPosi(const dictionary &structuralPositions);
    void readSurf(const dictionary &structuralElements);
    void readForce(const dictionary &structuralForces_memb);

    // Friend Functions

    // Friend Operators
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "netPanelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //r&  t
