/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 Must get permission from Hui Cheng
 Email: hui.cheng@uis.no
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
        const dictionary             netDict_memb;
        const scalar                 thickness_memb;  // thickness of the net panel
        const scalar                 Sn_memb;
        const dimensionedVector      D_memb;  //- porous resistance coefficient F
        const dimensionedVector      F_memb;  // private Member fucntion
              List<vector>           structuralPositions_memb;
              List<vector>           structuralElements_memb;
    //-Private member function    
        vector calcNorm
        (    
            const point& pointI,
            const point& pointII,
            const point& pointIII,
            const vector& fluidVelocity
        )const;
        vector calcNorm
        (    
            const point& pointI,
            const point& pointII,
            const point& pointIII
        )const;

        vector calcLifti
        (    
            const point& pointI,
            const point& pointII,
            const point& pointIII,
            const vector& fluidVelocity
        )const;   

        scalar calcTheta
        (    
            const point& pointI,
            const point& pointII,
            const point& pointIII,
            const vector& fluidVelocity
        )const;   

        scalar calcArea
        (    
            const point& pointI,
            const point& pointII,
            const point& pointIII
        )const;   

        bool isInPorousZone
        (
            const point&            x,
            const List<vector>&     structuralPositions,
            const vector&           structuralElementi
        )const;
        
        void transformCoeffs
        (
            const dimensionedVector&      D,
            tensor&                       d,
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
            const fvMesh&           mesh
        )const;

    void updatePoroField
         (
             volScalarField&         porosityField,
             const fvMesh&           mesh
         )const;
    scalar porosity() const
        {
            return Sn_memb;
        };

    List<vector>  posi() const
        {
            return structuralPositions_memb;
        };
    List<vector>  sur() const
        {
            return structuralElements_memb;
        };

    dimensionedVector  F() const
        {
            return F_memb;
        };
    
    void readPosi(
            const dictionary&        structuralPositions
    );
    void readSur(
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
