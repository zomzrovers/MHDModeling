   /* #include "readTimeControls.H"

    const dictionary& meltingDict = mesh.solutionDict().subDict("MELTING");
    int minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));
    int maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));
    scalar alphaTol(readScalar(meltingDict.lookup("alphaTolerance")));
    scalar alphaRel(readScalar(meltingDict.lookup("alphaRelaxatio")));
*/

/*
#include "readTimeControls.H"

const dictionary& meltingDict = mesh.solutionDict().subDict("MELTING");

// Use dictionary get method with default values
int minTCorr = meltingDict.get<int>("minTempCorrector");
int maxTCorr = meltingDict.get<int>("maxTempCorrector");
scalar alphaTol = meltingDict.get<scalar>("alphaTolerance");
scalar alphaRel = meltingDict.get<scalar>("alphaRelaxation"); // Corrected typo
*/

// In your readControls.H file


// Declare these variables BEFORE including readTimeControls.H
bool adjustTimeStep;
scalar maxCo;
scalar maxDeltaT;

#include "readTimeControls.H"


const dictionary& meltingDict = mesh.solutionDict().subDict("MELTING");

// Rest of your existing code...
int minTCorr = meltingDict.get<int>("minTempCorrector");
int maxTCorr = meltingDict.get<int>("maxTempCorrector");
scalar alphaTol = meltingDict.get<scalar>("alphaTolerance");
scalar alphaRel = meltingDict.get<scalar>("alphaRelaxation");
