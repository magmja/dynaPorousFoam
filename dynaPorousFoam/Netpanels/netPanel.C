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
void Foam::netPanel::transformCoeffs
(
    const dimensionedVector&      D,
    dimensionedTensor&            d,
    const List<vector>&         structuralPositions,
    const vector&                 structuralElementi

) const
{
    
    const point pointI_=structuralPositions[structuralElementi[0]];
    const point pointII_=structuralPositions[structuralElementi[1]];
    const point pointIII_=structuralPositions[structuralElementi[2]];
    
    const vector a(pointI_-pointII_);
    const vector b(pointIII_-pointII_);
    const vector c(pointIII_-pointI_);
          vector norm(a^b/(mag(a^b)+SMALL));
    const vector x0(1,0,0);
    const vector y0(0,1,0);
    const vector z0(0,0,1);

    vector e1(vector::zero); //r
    vector e2(vector::zero); //g
    vector e3(vector::zero); //b

    if((norm&x0)<0)  // assume the velocity is x+
    {
        norm = -norm;
    }
    e1 = norm;  


    if(mag(a)>mag(b))
    {
        e3 = b/(mag(b)+SMALL);
        if(mag(b)>mag(c))
        {
            e3 = c/(mag(c)+SMALL);
        } 
    }else
    {
        e3 = a/(mag(a)+SMALL);
    }
    
    if((e3&x0)<0)  // assume the velocity is x+
    {
        e3 = -e3;
    }


    e2=(e3^e1/(mag(e3^e1)+SMALL));  // blue
    if((e2&z0)<0)  // assume the velocity is x+
    {
        e2 = -e2;
    }

    Info << "p1 is : " << pointI_ << "\n" << endl;
    Info << "p2 is : " << pointII_ << "\n" << endl;
    Info << "p3 is : " << pointIII_ << "\n" << endl;

    Info << "e1 is : " << e1  << "\n" << endl;
    Info << "e2 is : " << e2  << "\n" << endl;
    Info << "e3 is : " << e3  << "\n" << endl;  
    // Info << "Tensor e is : " << e.value() << "\n" << endl;  
    const tensor e(
        e1&x0/(mag(e1*x0)+SMALL), e1&y0/(mag(e1*y0)+SMALL), e1&z0/(mag(e1*z0)+SMALL), 
        e2&x0/(mag(e2*x0)+SMALL), e2&y0/(mag(e2*y0)+SMALL), e2&z0/(mag(e2*z0)+SMALL), 
        e3&x0/(mag(e3*x0)+SMALL), e3&y0/(mag(e3*y0)+SMALL), e3&z0/(mag(e3*z0)+SMALL)); 

    d.value().xx() = D.value().x();
    d.value().yy() = D.value().y();
    d.value().zz() = D.value().z();
    d.value() = (e & d & e.T()).value();

}


bool Foam::netPanel::isInPorousZone
(
    const point x,
    const List<vector>&     structuralPositions,
    const vector&           structuralElementi
)const
{
    bool result(false); // initial value 

    const point pointI_=structuralPositions[structuralElementi[0]];
    const point pointII_=structuralPositions[structuralElementi[1]];
    const point pointIII_=structuralPositions[structuralElementi[2]];
    vector norm_=calcNorm(pointI_,pointII_,pointIII_);
    scalar dis(mag((x-pointI_)&norm_));
    // define a const scalar as the distance between point x to net panel
    if(dis <= thickness_)// distance is less than half thickness
    {
        scalar panelarea  // the area of the triangular net panel
             (
                 0.5*mag
                     (
                          (pointI_-pointII_)
                         ^(pointI_-pointIII_)
                     )
             );
        vector projectedPoint(0,0,0); // initial the projected point is 0,0,0
                if(((x-pointI_)&norm_)<0.0) // on the side of normal vector
        {
            projectedPoint=(x+norm_*dis);
        }
        else
        {
            projectedPoint=(x-norm_*dis);
        }
        // projectedP is the projected point on the net panel
        scalar panelarea3(0.5*(
                                    mag((projectedPoint-pointI_)^(projectedPoint-pointII_))+
                                    mag((projectedPoint-pointI_)^(projectedPoint-pointIII_))+
                                    mag((projectedPoint-pointII_)^(projectedPoint-pointIII_))
                                    )
                                ); //  the area of the three trigular shapes.
        if(panelarea3<=SMALL+panelarea)
        {
            result = true;
        }
    }
       
    return result;
}

Foam::vector Foam::netPanel::calcNorm
(
    const point& pointI,
    const point& pointII,
    const point& pointIII
)const
{
    const vector a(pointI-pointII);
    const vector b(pointI-pointIII);
    const vector norm(a^b/(mag(a^b)+SMALL));
    return norm; // can point out or in the panel.
}

void Foam::netPanel::addResistance
(
            fvVectorMatrix&         UEqn,
            const volScalarField&   nu,
            const fvMesh&           mesh,
            const List<vector>&   structuralPositions,
            const List<vector>&   structuralElements
)const
{
    const vectorField& centres(mesh.C());
    forAll(structuralElements,Elementi)
    {
        dimensionedTensor d_(tensor::zero);
        dimensionedTensor f_(tensor::zero);
        transformCoeffs(D_,d_,structuralPositions,structuralElements[Elementi]);
        transformCoeffs(F_,f_,structuralPositions,structuralElements[Elementi]);
        tensor& dvalue = d_.value();
        tensor& fvalue = f_.value();
        const scalarField V = mesh.V();
        vectorField& Usource = UEqn.source();
        const vectorField& U = UEqn.psi();

        forAll(centres, cellI)
        {
            if(
                isInPorousZone(centres[cellI],structuralPositions,structuralElements[Elementi])
             )
             {
                tensor dragCoeff = nu[cellI]*dvalue + 0.5*mag(U[cellI])*fvalue;
                Usource[cellI] -=V[cellI]*dragCoeff & (U[cellI] );
             }
        }
    }

}

void Foam::netPanel::updatePoroField
(
             volScalarField&         porosityField,
             fvMesh&                 mesh,
             const List<vector>&     structuralPositions,
             const List<vector>&     structuralElements
)const
{
    // step1 set all the cell as 1
    forAll(mesh.C(),cellI)
    {
        porosityField[cellI] = 1.0;
    }
    // get the center of all the cells
    const vectorField& centres(mesh.C());
    // loop through all the structural emlements
    forAll(structuralElements,Elementi)
    {
        // loop through all the cell, 
        forAll(centres, cellI)
        {
            if(isInPorousZone(centres[cellI],structuralPositions,structuralElements[Elementi]))
            {
                porosityField[cellI] = porosity_;
            }
        }

    }
}
void Foam::netPanel::readSE
(
    const dictionary&        structuralPositions,
    const dictionary&        structuralElements
)
{
    
    // todo 
    //read the dictionary and assign the value to the two member 
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::netPanel::netPanel
(
    const dictionary&        netDict
)
:   // initial components
    netDict_(netDict),
    porousPropertiesDict_(netDict_.subDict("porousProperties")),
    porosity_(readScalar(porousPropertiesDict_.lookup("porosity"))),
    thickness_(readScalar(porousPropertiesDict_.lookup("halfthickness")))
    // D_(readVector(porousPropertiesDict_.lookup("D"))),  // initial as zeros
    // F_(readVector(porousPropertiesDict_.lookup("D")))  // initial as zeros
{
    dimensionedVector D_(porousPropertiesDict_.lookup("D"));
    dimensionedVector F_(porousPropertiesDict_.lookup("F"));

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::netPanel::~netPanel()
{}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
