/*   IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Reading density rho
    dimensionedScalar rho(transportProperties.lookup("rho"));

    // Reading specific heat capacity
    dimensionedScalar cps(transportProperties.lookup("cps"));
    dimensionedScalar cpl(transportProperties.lookup("cpl"));

    // Reading thermal conductivity
    dimensionedScalar lambdas(transportProperties.lookup("lambdas"));
    dimensionedScalar lambdal(transportProperties.lookup("lambdal"));

    // Reading kinematic viscosity
    dimensionedScalar nu(transportProperties.lookup("nu"));

    // Reading latent heat of fusion
    dimensionedScalar L(transportProperties.lookup("L"));

    // Reading melting temperature range
    dimensionedScalar Ts(transportProperties.lookup("Ts"));
    dimensionedScalar Tl(transportProperties.lookup("Tl"));
    
    dimensionedScalar Tmelt
    (
        "Tmelt", 
        (Tl+Ts)/scalar(2.0)
    );

    // Reading volume expansion factor
    dimensionedScalar beta(transportProperties.lookup("beta"));

    // Reading large D'arcy-type source term constant DCl
    dimensionedScalar DCl(transportProperties.lookup("DCl"));
    dimensionedScalar DCs(transportProperties.lookup("DCs"));
*/

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Reading density rho
dimensionedScalar rho("rho", dimDensity, transportProperties.get<scalar>("rho"));

// Reading specific heat capacity
dimensionedScalar cps("cps", dimSpecificHeatCapacity, transportProperties.get<scalar>("cps"));
dimensionedScalar cpl("cpl", dimSpecificHeatCapacity, transportProperties.get<scalar>("cpl"));

// Reading thermal conductivity
dimensionedScalar lambdas("lambdas", dimensionSet(1, 1, -3, -1, 0, 0, 0) , transportProperties.get<scalar>("lambdas"));
dimensionedScalar lambdal("lambdal", dimensionSet(1, 1, -3, -1, 0, 0, 0) , transportProperties.get<scalar>("lambdal"));

// Reading kinematic viscosity
dimensionedScalar nu("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), transportProperties.get<scalar>("nu"));

// Reading Magnetic permeability
dimensionedScalar mu("mu", dimensionSet(1, 1, -2, 0, 0, -2, 0), transportProperties.get<scalar>("mu"));

// Reading latent heat of fusion
dimensionedScalar L("L", dimensionSet(0, 2, -2, 0, 0, 0, 0), transportProperties.get<scalar>("L"));

// Reading melting temperature range
dimensionedScalar Ts("Ts", dimTemperature, transportProperties.get<scalar>("Ts"));
dimensionedScalar Tl("Tl", dimTemperature, transportProperties.get<scalar>("Tl"));

dimensionedScalar Tmelt
(
    "Tmelt", 
    dimTemperature,
    (Tl.value() + Ts.value())/2.0
);

// Reading volume expansion factor
dimensionedScalar beta("beta", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties.get<scalar>("beta"));

// Reading large D'arcy-type source term constant
/*
dimensionedScalar DCl("DCl", dimLength*dimLength/dimTime, transportProperties.get<scalar>("DCl"));
dimensionedScalar DCs("DCs", dimLength*dimLength/dimTime, transportProperties.get<scalar>("DCs"));
*/
dimensionedScalar DCl("DCl", dimensionSet(1, -3, -1, 0, 0, 0, 0), transportProperties.get<scalar>("DCl"));
dimensionedScalar DCs("DCs", dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties.get<scalar>("DCs"));



// Magnetic Properties
dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    transportProperties.get<scalar>("sigma")
);
