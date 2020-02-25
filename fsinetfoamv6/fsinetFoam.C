/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "netPanel.H"
#include "OFstream.H"
#include "Pstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n"
         << endl;
    OFstream* myOutFilePtr = NULL;
    if (Pstream::master())
    {
        // Open the file
        myOutFilePtr = new OFstream("velocity_on_elements.txt");
    }
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();
        // read from outside
        // need to confirm ... might be wrong data structure
        IOdictionary structuralPositions(
                IOobject(
                        "posi",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE));
        Nettings.readPosi(structuralPositions);
        Nettings.updatePoroField(porosityField, mesh);

        if (exists("./constant/Fh"))
        {
            Info<< "Reading Fh"<<endl;
            IOdictionary structuralFh(
                    IOobject(
                            "Fh",
                            runTime.constant(),
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::NO_WRITE));
            Nettings.readForce(runTime.value(),structuralFh);
            Info<< "Finish reading Fh"<<endl;
            // Info << "The new force on elements are \n"<<Nettings.Fhout()<<endl;
        }

        List<pointField> gatheredU(numberP);
        gatheredU[Pstream::myProcNo()] = pointField(U);
        Pstream::gatherList(gatheredU);

//        Info << "Number of processors = " << numberP <<endl;
//        Info << "The gatheredpoints is  = " << gatheredcentres.size() <<endl;
//        Info << "The gatheredU is  = " << gatheredU[0].size()<<endl;

        Nettings.updateVelocity(gatheredU,gatheredCentres);

//        Info<<"velocity on the center of net panels are \n"<<Nettings.FluidU()<<endl;
        // write the Nettings.fluidVelocity(); to a extrinal files
        if (Pstream::master())
        {
            OFstream& myOutFile = *myOutFilePtr;
            myOutFile
                    << "The velocities at " << runTime.timeName()<< "s are: "<< Nettings.FluidU()  << endl;
        }
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //