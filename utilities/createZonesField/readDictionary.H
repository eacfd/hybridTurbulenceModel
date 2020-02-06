    algorithm =  word(cZFDict.lookup("Algorithm"));
    Info << "Algorithm for " << algorithm << " mesh selected" << endl << endl;

    if (algorithm == "cellZones")
    {
        const dictionary& cellZoneDict = cZFDict.subDict("cellZonesDict");
        LESCellZoneName = word(cellZoneDict.lookup("LESCellZoneName")); 
        blendedCellZoneName = word(cellZoneDict.lookup("BlendedCellZoneName")); 
    }
    else if (algorithm == "custom")
    {
        const dictionary& customDict = cZFDict.subDict("customDict");
        minEdge = readScalar(customDict.lookup("minEdge")); 
        blendedDist = readScalar(customDict.lookup("blendedDist"));
        clearField = readBool(customDict.lookup("clearField"));
        neiCells = readLabel(customDict.lookup("neighbourCells"));
        scaleFactor = readScalar(customDict.lookup("scaleFactor"));
    }
    else
    {
        FatalError
        << "Algorithm " << algorithm << " not found. " <<
        "Availible algorithms are" << endl << 
        "1. cellZones" << endl << "2. custom" << endl  << 
        exit(FatalError);
    }
