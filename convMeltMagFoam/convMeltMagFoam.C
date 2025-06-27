/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Reference
    Universität Bayreuth
    Lehrstuhl für Technische Thermodynamik und Trasnportprozesse - LTTT
    Fabian Rösler
    Universitätsstraße 30
    95440 Bayreuth
    Tel.: +49 (921) 55-7163

Application
    convMeltFoam

Description
    Solves a convection dominated solid/liquid phase change process.
    Density is constant and convection is induced by Boussinesq approximation.
    Phase change is described by means of a linear function.
    Enthalpy update with fictious heat source method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }
        
        // --- B-PISO loop
        while (pimple.correct())
        {
            fvVectorMatrix BEqn
            (
                fvm::ddt(B)
              + fvm::div(phi, B)
              - fvc::div(phiB, U)
              - fvm::laplacian(DB, B)
            );

            BEqn.solve();

            volScalarField rAB(1.0/BEqn.A());
            surfaceScalarField rABf("rABf", fvc::interpolate(rAB));

            phiB = fvc::flux(B);

            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix pBEqn
                (
                    fvm::laplacian(rABf, pB) == fvc::div(phiB)
                );

                pBEqn.solve(pB.select(pimple.finalInnerIter()));

                if (pimple.finalNonOrthogonalIter())
                {
                    phiB -= pBEqn.flux();
                }
            }

            #include "magneticFieldErr.H"
		}
        

        Info<< "Liquid fraction [-] = " << alpha1.weightedAverage(mesh.V()).value()
            << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
