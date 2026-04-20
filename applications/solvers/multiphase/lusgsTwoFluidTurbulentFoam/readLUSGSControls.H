dictionary lusgsDict = mesh.solutionDict().subDict("LUSGS");

lusgsIntIters = 
    lusgsDict.lookupOrDefault<label>("internalIterations", 1);

lusgsIters = 
    lusgsDict.lookupOrDefault<label>("iterations", 1);

lusgsTolerance =
    lusgsDict.lookupOrDefault<scalar>("tolerance", 0.0);

lusgsRelTol =
    lusgsDict.lookupOrDefault<scalar>("relTol", 0.0);

// relaxation factor 1<omega<=2
lusgsOmega = 
    lusgsDict.lookupOrDefault<scalar>("omega", 1.2);

lusgsRelax = 
    lusgsDict.lookupOrDefault<scalar>("relaxCoeff", 1.0);

LTS = 
    lusgsDict.lookupOrDefault("localTimestepping", false);

turbOnFinalIterOnly =
    lusgsDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", true);
    