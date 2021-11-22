/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "usertwoPhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(usertwoPhaseChangeModel, 0);
    defineRunTimeSelectionTable(usertwoPhaseChangeModel, dictionary);
}

const Foam::word Foam::usertwoPhaseChangeModel::phaseChangePropertiesName
(
    "phaseChangeProperties"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::usertwoPhaseChangeModel::createIOobject
(
    const userimmiscibleIncompressibleTwoPhaseMixture& mixture
) const
{
    IOobject io
    (
        phaseChangePropertiesName,
        mixture.U().mesh().time().constant(),
        mixture.U().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::usertwoPhaseChangeModel::usertwoPhaseChangeModel
(
    const word& type,
    const userimmiscibleIncompressibleTwoPhaseMixture& mixture
)
:
    IOdictionary(createIOobject(mixture)),
    mixture_(mixture),
    twoPhaseChangeModelCoeffs_(optionalSubDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::usertwoPhaseChangeModel::correct()
{}


bool Foam::usertwoPhaseChangeModel::read()
{
    if (regIOobject::read())
    {
        twoPhaseChangeModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
