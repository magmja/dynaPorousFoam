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
        const dictionary             netDict_;
        const dictionary             porousPropertiesDict_;
        const scalar                 thickness_;  // thickness of the net panel
        const scalar                 porosity_;
        const dimensionedVector      D_;  //- porous resistance coefficient F
        const dimensionedVector      F_;  // private Member fucntion
              List<vector>           structuralPositions_;
              List<vector>           structuralElements_;
    //-Private member function    
        vector calcNorm
        (    
            const point& pointI,
            const point& pointII,
            const point& pointIII
        )const;   

        bool isInPorousZone
        (
            const point             x,
            const List<vector>&     structuralPositions,
            const vector&           structuralElementi
        )const;
        
        void transformCoeffs
        (
            const dimensionedVector&      D,
            dimensionedTensor&            d,
            const List<vector>&           structuralPositions,
            const vector&                 structuralElementi
        ) const;

public:
    // Constructors
       netPanel
       (
           const dictionary&        netDict
       );

    //- Destructor
    ~netPanel();

    // public Member Functions
    void addResistance
        (
            fvVectorMatrix&         UEqn,
            const volScalarField&   nu,
            const fvMesh&           mesh,
            const List<vector>&     structuralPositions,
            const List<vector>&     structuralElements
        )const;

    void updatePoroField
         (
             volScalarField&         porosityField,
             fvMesh&                 mesh,
             const List<vector>&     structuralPositions,
             const List<vector>&     structuralElements
         )const;
    scalar porosity() const
        {
            return porosity_;
        };

    void readSE(
            const dictionary&        structuralPositions,
            const dictionary&        structuralElements
    );

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
