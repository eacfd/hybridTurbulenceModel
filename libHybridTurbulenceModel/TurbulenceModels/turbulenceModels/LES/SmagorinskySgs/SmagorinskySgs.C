/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "SmagorinskySgs.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SmagorinskySgs<BasicTurbulenceModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*Ck_*this->delta()*(dev(D) && D));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("kSgs", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
        )
    );
}


template<class BasicTurbulenceModel>
void SmagorinskySgs<BasicTurbulenceModel>::correctNut()
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ = Ck_*this->delta()*sqrt(k);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SmagorinskySgs<BasicTurbulenceModel>::SmagorinskySgs
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscositySgs<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    
    k_
    (
        IOobject
        (
            IOobject::groupName("kSgs", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kSgs", dimensionSet(0, 2, -2, 0, 0, 0, 0), scalar(0.0))
    ),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SmagorinskySgs<BasicTurbulenceModel>::read()
{
    if (LESeddyViscositySgs<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SmagorinskySgs<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SmagorinskySgs<BasicTurbulenceModel>::omega() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));
    volScalarField epsilon(this->Ce_*k*sqrt(k)/this->delta());

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        epsilon/(0.09*k)
    );
}

template<class BasicTurbulenceModel>
void SmagorinskySgs<BasicTurbulenceModel>::correct()
{
    LESeddyViscositySgs<BasicTurbulenceModel>::correct();
    
    k_ = k(fvc::grad(this->U_));
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
