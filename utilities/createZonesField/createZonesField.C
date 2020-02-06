/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    createZonesField

Description
    Create a initial volScalarField zones for hybrid RANS/LES approach. This 
    field is required to define LES/Blended/RANS zones. So value 0 - LESZone, 
    0.5 - BlendedZone, 1 - RANS/URANS region.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a zones field for hybrid approache."
    );
    
    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"
    
    const word dictName("createZonalFieldDict");
    #include "setSystemRunTimeDictionaryIO.H"
    const IOdictionary cZFDict(dictIO);

    word algorithm;
    word LESCellZoneName, blendedCellZoneName;
    scalar minEdge = 0.01;
    scalar blendedDist = 4*minEdge;
    scalar scaleFactor = 1.1;
    label neiCells = 3;
    bool clearField = true;
    
    #include "readDictionary.H"
    
    volScalarField LESZone
    (
        IOobject
        (
            "zones",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zones", dimless, scalar(1.0)),
        "zeroGradient"
    );
    
    forAll(LESZone, I)
    {
        LESZone[I] = 1;
    }

    if (algorithm == "cellZones")
    {
        //forAll(mesh.cellZones(), zoneI) {}
        List<word> zoneNames(2);
        zoneNames[0] = blendedCellZoneName;
        zoneNames[1] = LESCellZoneName;
        
        Info << "LES cellZone name " << LESCellZoneName << endl;
        Info << "Blended cellZone name " << blendedCellZoneName << endl;
          
        forAll(zoneNames, zNI)
        {
            const label cellZoneID = mesh.cellZones().findZoneID(zoneNames[zNI]);
            if (cellZoneID < 0)
            {
                FatalError
                << "Unable to find cell zone " << zoneNames[zNI]  << endl
                << exit(FatalError);
            }

            const cellZone& zone = mesh.cellZones()[cellZoneID];
            forAll(zone, I)
            {
                const label cellID = zone[I];
                if (zoneNames[zNI] == LESCellZoneName)
                {
                    LESZone[cellID] = 0;
                }
                else
                {
                    LESZone[cellID] = 0.5;
                }
            }
        }
    }
    else if (algorithm == "custom")
    {
        label countElement = 0;
        scalarField zeroField;
        volVectorField cellCenters = mesh.C();
        
        Info << "minEdge is " << minEdge << endl;
        
        forAll(LESZone, I)
        {
            if ((mesh.V()[I]) <= (1.1*minEdge*minEdge*minEdge)) 
            {
                LESZone[I] = 0;
                countElement = countElement + 1;
            }
        }
        label allLESelem = countElement;
        reduce(allLESelem, sumOp<label>());
        Info << "LESZone elements = " << allLESelem << endl << endl;
        
        if (clearField == true)
        {
            Info << "Deleting single LESZone elemenets" << endl;
            Info << "Neighbour cells set to " << neiCells << ". ScaleFactor is " << scaleFactor << endl;
            zeroField.setSize(countElement);
        
            label count=0;
            forAll(LESZone, I)
            {
                if (LESZone[I] == 0)
                {
                    zeroField[count] = I;
                    count++;
                }
            }
      
            label iter = 0;
        
            forAll(zeroField, I)
            {
                iter = 0;
                scalar relaxion = scaleFactor*Foam::pow(mesh.V()[zeroField[I]],1.0/3.0);
                forAll(zeroField, J)
                {
                    if ((mag(cellCenters[zeroField[J]] - cellCenters[zeroField[I]]) <= relaxion) && (zeroField[I] != zeroField[J])) 
                    {
                        iter = iter + 1;
                    }
                    if (iter >= neiCells) 
                    { 
                        break;
                    }
                }
                if (iter < neiCells) 
                {
                    LESZone[zeroField[I]] = 1;
                    allLESelem = allLESelem - 1;;
                }
            }
            Info << "LESZone elements = " << allLESelem << endl << endl;
        }
        
        Info << "Creating blendedZone" << endl;
        Info << "Blended zone distance is " << blendedDist << endl;
        forAll(LESZone, I)
        {
            if (LESZone[I] == 0)
            {
                forAll(cellCenters, J)
                {
                    if ((I != J) && (LESZone[J] == 1))
                    {
                        if (mag(cellCenters[J] - cellCenters[I]) < blendedDist)
                        {
                            LESZone[J] = 0.5;
                        }
                    }
                }
            }
        }
    }

    Info << endl << "Writing zones field" << endl << endl;
    LESZone.write();
    Info << "End. Time = " << LESZone.time().elapsedCpuTime() << " s" << endl << endl;
}


// ************************************************************************* //
