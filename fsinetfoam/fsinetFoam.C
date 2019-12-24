/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "netPanel.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Transient solver for incompressible, turbulent flow,"
        " using the PISO algorithm.");

#include "postProcess.H"

#include "addCheckCaseOptions.H"
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
    fileName output("velocityOnNodes.dat");
    OFstream os(output);
    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

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

        // Info << "The new positions are"<<Nettings.posi()<<endl;
        IOdictionary structuralFh(
            IOobject(
                "Fh",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE));
        Nettings.readForce(structuralFh);
        // Info << "The new Force are"<<Nettings.Fhout()<<endl;
        
        Nettings.updatePoroField(porosityField, mesh);
        
        // Info << "Before update the fluid velocity"<<Nettings.FluidU()<<endl;
        Info << "the velocity on nodes are"<< Nettings.updateVelocity(mesh,U) << endl;
        // Info << "After updateporofield fluid velocity are"<<Nettings.FluidU()<<endl;
        
        os << Nettings.updateVelocity(mesh,U) << endl;
        // write the Nettings.fluidVelocity(); to a extrinal files

        runTime.write();
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
