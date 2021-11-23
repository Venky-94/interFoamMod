#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
userimmiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<usertwoPhaseChangeModel> phaseChangePtr
(
    usertwoPhaseChangeModel::New(mixture)
);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const volScalarField& nuModel1 = mixture.nuM1();
const volScalarField& nuModel2 = mixture.nuM2();

const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1, scalar(0)), scalar(1))
    );

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
	rho2 + (rho1 - rho2)*pow(limitedAlpha1, 1.5*(nuModel1/nuModel2) + 0.75)
	//alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Getting dynamic viscocity for postprocessing
volScalarField muproc
(
	IOobject
	(
		"muproc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	//rho2*nuModel2 + (rho1*nuModel1 - rho2*nuModel2)*pow(alpha1, 1.5*(nuModel1/nuModel2) + 0.75)
	mixture.mu()
);

// Getting kinematic viscocity for postprocessing
volScalarField nuproc
(
	IOobject
	(
		"nuproc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mixture.nu()
);

// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
