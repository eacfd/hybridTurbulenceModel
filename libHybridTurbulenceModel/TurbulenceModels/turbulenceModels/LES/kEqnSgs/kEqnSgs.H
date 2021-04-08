/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::LESModels::kEqnSgs

Group
    grpLESTurbulence

Description
    One equation eddy-viscosity model

    Eddy viscosity SGS model using a modeled balance equation to simulate the
    behaviour of k.

    Reference:
    \verbatim
        Yoshizawa, A. (1986).
        Statistical theory for compressible turbulent shear flows,
        with the application to subgrid modeling.
        Physics of Fluids (1958-1988), 29(7), 2152-2164.
    \endverbatim

    The default model coefficients are
    \verbatim
        kEqnSgsCoeffs
        {
            Ck                  0.094;
            Ce                  1.048;
        }
    \endverbatim

SourceFiles
    kEqnSgs.C

\*---------------------------------------------------------------------------*/

#ifndef kEqnSgs_H
#define kEqnSgs_H

#include "LESeddyViscositySgs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class kEqnSgs Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEqnSgs
:
    public LESeddyViscositySgs<BasicTurbulenceModel>
{
protected:

    // Protected data

        // Fields

            volScalarField k_;


        // Model constants

            dimensionedScalar Ck_;


    // Protected Member Functions

        //- No copy construct
        kEqnSgs(const kEqnSgs&) = delete;

        //- No copy assignment
        void operator=(const kEqnSgs&) = delete;


        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEqnSgs");


    // Constructors

        //- Constructor from components
        kEqnSgs
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
    virtual ~kEqnSgs()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return sub-grid dissipation rate
        virtual tmp<volScalarField> epsilon() const;
        
        //- Return sub-grid specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", this->nut_ + this->nu())
            );
        }

        //- Correct eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEqnSgs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
