/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "zonalHybrid.H"
#include "fvOptions.H"

#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HybridModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
void zonalHybrid<BasicTurbulenceModel>::correctNut()
{
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    
    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
zonalHybrid<BasicTurbulenceModel>::zonalHybrid
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
    eddyViscosityHybrid<HybridModel<BasicTurbulenceModel>>
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
    
    rasPtr_(Foam::RASModel<BasicTurbulenceModel>::New(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)),
    lesPtr_(Foam::LESModel<BasicTurbulenceModel>::New(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)),

    zonal_
    (
        IOobject
        (
            "zones",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    mutualK_
    (
        IOobject
        (
            "mutualK",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    
    UBlendingFactor
    (
        IOobject
        (
           "UBlendingFactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(zonal_)
    ),
    
    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> zonalHybrid<BasicTurbulenceModel>::epsilon() const
{
    dimensionedScalar nutSmall ("nutSmall", dimensionSet(0, 2, -1, 0, 0, 0, 0), SMALL);
    return Cmu_*sqr(mutualK_) / (this->nut_ + nutSmall);
}

template<class BasicTurbulenceModel>
bool zonalHybrid<BasicTurbulenceModel>::read()
{
    if (eddyViscosityHybrid<HybridModel<BasicTurbulenceModel>>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField> zonalHybrid<BasicTurbulenceModel>::k() const
{
    return mutualK_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> zonalHybrid<BasicTurbulenceModel>::nut() const
{
    return this->nut_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> zonalHybrid<BasicTurbulenceModel>::LESRegion() const
{
    return zonal_;
}


template<class BasicTurbulenceModel>
void zonalHybrid<BasicTurbulenceModel>::correct()
{
    eddyViscosityHybrid<HybridModel<BasicTurbulenceModel>>::correct();

    //correct RAS model variables
    Info << "RAS Correction" << endl;
    rasPtr_->correct();
    
    const volScalarField& rasNut = rasPtr_->nut();
    const volScalarField& rasK = rasPtr_->k();
    //const volScalarField& nu_ = this->nu();

    //correct LES model variables
    Info << "LES Correction" << endl;
    lesPtr_->correct();

    //calculate Hybrid viscousity
    const volScalarField& lesNuSgs = lesPtr_->nut();
    const volScalarField& lesK = lesPtr_->k();
    
    this->nut_ = rasNut;
    mutualK_ = rasK;

    dimensionedScalar kSmall ("kSmall", dimensionSet(0, 2, -2, 0, 0, 0, 0), SMALL);

    Info << "Hybrid Correction" << endl;

    forAll(zonal_, I)
    {
        if (zonal_[I] == 0)
        {
            this->nut_[I] = lesNuSgs[I];
            mutualK_[I] = lesK[I];
        }
        if (zonal_[I] == 0.5) 
        {
            scalar ff = Foam::tanh(sqr(rasNut[I]/Cmu_.value()/sqrt(mutualK_[I])/5/lesPtr_->delta()[I])); //similar to  Xiao blended functon
            //scalar ff = Foam::tanh( pow(rasNut[I]/rasK[I]*max(500*nu_[I]/sqr(y_[I]),sqrt(rasK[I])/Cmu_.value()/y_[I] ),4)); //similar to Frohlich
            this->nut_[I] = (1-ff)*rasNut[I] + ff*lesNuSgs[I];
            mutualK_[I] = (1-ff)*rasK[I] + ff*lesK[I];
        }

    }

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HybridModels
} // End namespace Foam

// ************************************************************************* //
