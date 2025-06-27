    //Reading thermophysical properties
    Info<< "Reading thermophysical properties\n" << endl;

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    //Reading transport properties
    Info<< "Reading transportProperties\n" << endl;

    #include "readTransportProperties.H"


    //Calculating phase-fraction alpha
    Info<< "Calculating phase-fraction alpha\n" << endl;

    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        max(min((T-Ts)/(Tl-Ts),scalar(1)),scalar(0))
    );
    alpha1.write();

    volScalarField alpha2 ("alpha2", scalar(1) - alpha1);
    

    //Creating fields for thermophysical properties
    Info<< "Creating fields for thermophysical properties\n" << endl;

    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        scalar(1)-max(beta*(T-Tl),scalar(0))
    );

    volScalarField DC
    (
        IOobject
        (
            "DC",
            runTime.timeName(),
            mesh
        ),
        DCl*pow(alpha2,2)/(pow(alpha1,3)+DCs)
    );

    volScalarField lambda
    (
        IOobject
        (
            "lambda",
            runTime.timeName(),
            mesh
        ),
        alpha1*lambdal + alpha2*lambdas
    );

    volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            mesh
        ),
        alpha1*cpl + alpha2*cps
    );
    cp.oldTime();

    surfaceScalarField phiCp
    (
        IOobject
        (
            "phi*cp",
            runTime.timeName(),
            mesh
        ),
        phi*fvc::interpolate(cp)
    );


    //Calculating enthalpy and corrected temperature
    Info<< "Calculating enthalpy and corrected temperature\n" << endl;

    volScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cp*(T-Tmelt) + alpha1*L
    );
    h.write();

    volScalarField Tcorr("Tcorr", (Tl-Ts)*alpha1 + Ts);

    //Creating fields for gravity and pressure
    Info<< "Calculating field g.h\n" << endl;

    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue-getRefCellValue(p, pRefCell)
        );
    }
    
    
    
	Info<< "Reading field pB\n" << endl;
	volScalarField pB
	(
		IOobject
		(
			"pB",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);


	Info<< "Reading field B\n" << endl;
	volVectorField B
	(
		IOobject
		(
			"B",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	#include "createPhiB.H"
	dimensionedScalar DB = 1.0/(mu*sigma);
	DB.name() = "DB";

	dimensionedScalar DBU = 1.0/(2.0*mu*rho);
	DBU.name() = "DBU";

