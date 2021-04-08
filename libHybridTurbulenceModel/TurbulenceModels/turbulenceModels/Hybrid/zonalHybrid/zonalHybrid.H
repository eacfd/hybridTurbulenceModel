/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

Class
    Foam::HybridModels::zonalHybrid

Description
    zonalHybrid turbulence model for incompressible flows base on combining URANS and LES approaches

    Reference:
    \verbatim
        Epikhin, A. (2019), Numerical Schemes and Hybrid Approach for the Simulation of Unsteady Turbulent Flows. 
        Mathematical Models and Computer Simulations, Vol. 11, No. 6, pp. 1019-1031.
    \endverbatim

Note
    For LES region recomended to use GammaM or filteredLinearM schemes, which implemented in libSchemes 
    https://github.com/unicfdlab/libSchemes.

SourceFiles
    zonalHybrid.C

\*---------------------------------------------------------------------------*/

#ifndef zonalHybrid_H
#define zonalHybrid_H

#include "HybridModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "eddyViscosityHybrid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HybridModels
{

/*---------------------------------------------------------------------------*\
                        Class zonalHybrid Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class zonalHybrid
:
    public eddyViscosityHybrid<HybridModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        zonalHybrid(const zonalHybrid&) = delete;

        //- No copy assignment
        void operator=(const zonalHybrid&) = delete;


protected:

    // Protected data

        autoPtr<Foam::RASModel<BasicTurbulenceModel>> rasPtr_;
        autoPtr<Foam::LESModel<BasicTurbulenceModel>> lesPtr_;
            

        // Fields
        volScalarField zonal_;
            
        //-
        volScalarField mutualK_;
            
        //-
        dimensionedScalar Cmu_;
            
        //-use localBlended scheme for div(phi,U)
        surfaceScalarField UBlendingFactor;

            
        //- Wall distance
        //  Note: different to wall distance in parent RASModel
        //  which is for near-wall cells only
        const volScalarField& y_;


    // Protected Member Functions

        void correctNut(const volScalarField& fv1);
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("zonalHybrid");


    // Constructors

        //- Construct from components
        zonalHybrid
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~zonalHybrid()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();
        
        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const;
        
        //- Return SGS kinetic energy
        virtual tmp<volScalarField> nut() const;
        
        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
        
        virtual tmp<volScalarField> epsilon() const;
        
        //- Return sub-grid specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Correct nuTilda and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HybridModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zonalHybrid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
