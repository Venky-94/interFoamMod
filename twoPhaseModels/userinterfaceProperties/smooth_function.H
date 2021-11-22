/*---------------------------------------------------------------------------*\

    This routine is used to smoothen a volScalarField over a 3D computational 
    domain in OpenFOAM.

    The used smoother was proposed in Onno Ubbink's Ph.D thesis "Numerical 
    prediction of two fluid systems with sharp interfaces".

    This routine is written by Duong Anh Hoang, Ph.D student in Multi-scale 
    Physics Department, Delft University of Technology.

    Email: h.a.duong@tudelft.nl

\*---------------------------------------------------------------------------*/

void Foam::userinterfaceProperties::smoothen
(
    volScalarField& smooth_func
) const
{ 
    const fvMesh& mesh = smooth_func.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    scalar n = 2;

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    for(int iter = 0; iter < n; iter++)
    {
    	scalarField smooth_cal(mesh.nCells(),scalar(0));

    	scalarField sum_area(mesh.nCells(),scalar(0));

	surfaceScalarField smoothF = fvc::interpolate(smooth_func);

	for(int facei = 0; facei < nei.size(); facei++)
	{
	    smooth_cal[own[facei]] += smoothF[facei]*mag(Sf[facei]);
	    sum_area[own[facei]] += mag(Sf[facei]);
	}

	forAll(nei,facei)
	{
	    smooth_cal[nei[facei]] += smoothF[facei]*mag(Sf[facei]);
	    sum_area[nei[facei]] += mag(Sf[facei]);
	}

        forAll(mesh.boundary(), patchi)
        {
            const unallocLabelList& pFaceCells = mesh.boundary()[patchi].faceCells();

            const fvsPatchScalarField& pssf = smoothF.boundaryField()[patchi];

            forAll(mesh.boundary()[patchi], facei)
            {
               smooth_cal[pFaceCells[facei]] += pssf[facei]*mag(Sf[facei]);
               sum_area[pFaceCells[facei]] += mag(Sf[facei]);
            }
        }

	forAll(mesh.cells(),celli)
	{
	    smooth_func[celli] = smooth_cal[celli]/sum_area[celli];
	}

	smooth_func.correctBoundaryConditions();

    }
}

/*---------------------------------------------------------------------------*/
