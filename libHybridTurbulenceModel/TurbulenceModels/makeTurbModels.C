/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IncompressibleTurbulenceModel.H"
#include "transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

//#include "laminarModel.H"
//#include "RASModel.H"
#include "LESModel.H"
#include "HybridModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#define makeHybridTurbulenceModelTypes(Alpha, Rho, baseModel, BaseModel, Transport)  \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel<Transport> Transport##BaseModel;                      \
        typedef LESModel<Transport##BaseModel> LES##Transport##BaseModel;        \
        typedef HybridModel<BaseModel<Transport>> Hybrid##Transport##BaseModel;  \
    }

    makeHybridTurbulenceModelTypes
    (
        geometricOneField,
        geometricOneField,
        incompressibleTurbulenceModel,
        IncompressibleTurbulenceModel,
        transportModel
    );
    
#define makeHybridBaseTurbulenceModel(Alpha, Rho, baseModel, BaseModel, Transport)   \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef TurbulenceModel                                                \
        <                                                                      \
            Alpha,                                                             \
            Rho,                                                               \
            baseModel,                                                         \
            Transport                                                          \
        > Transport##baseModel;                                                \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            Transport##baseModel,                                              \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(Hybrid##Transport##BaseModel, 0); \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (Hybrid##Transport##BaseModel, dictionary);                           \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##baseModel,                                              \
            Hybrid##Transport##BaseModel,                                     \
            dictionary                                                         \
        );                                                                     \
    }

    makeHybridBaseTurbulenceModel
    (
        geometricOneField,
        geometricOneField,
        incompressibleTurbulenceModel,
        IncompressibleTurbulenceModel,
        transportModel
    );
    
#define makeHybridModel(Type)                                                  \
    makeTemplatedTurbulenceModel                                               \
    (transportModelIncompressibleTurbulenceModel, Hybrid, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelIncompressibleTurbulenceModel, LES, Type)

#include "SmagorinskySgs.H"
makeLESModel(SmagorinskySgs);

#include "kEqnSgs.H"
makeLESModel(kEqnSgs);

#include "zonalHybrid.H"
makeHybridModel(zonalHybrid);


// ************************************************************************* //
