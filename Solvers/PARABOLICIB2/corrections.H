
forAll(rho,celli)
{

if(rho.ref()[celli]<unif.get_rhoMin() || rho.ref()[celli]>unif.get_rhoMax() || 
	e.ref()[celli]<unif.get_eMin() || e.ref()[celli]>unif.get_eMax() || Y.ref()[celli]<unif.get_YMin()
	|| Y.ref()[celli]>unif.get_YMax() )

	{
	labelList faces = mesh.cells()[celli];
    labelList cells;

  //  Info << faces << endl;

	forAll(faces, facei)
	{
		if(mesh.isInternalFace(faces[facei]))
		{
		if(mesh.owner()[faces[facei]] == celli)
		{
			cells.append(mesh.neighbour()[faces[facei]]);
		}
		else
		{
			cells.append(mesh.owner()[faces[facei]]);
		}
	    }
    }

    rho.ref()[celli] = 0.0;
	rhoU.ref()[celli] = vector::zero;
	rhoE.ref()[celli] = 0.0;
	rhoY.ref()[celli] = 0.0;

	forAll(cells, celli_1)
	{
		rho.ref()[celli] += rho.ref()[cells[celli_1]];
		rhoU.ref()[celli] += rhoU.ref()[cells[celli_1]];
		rhoE.ref()[celli] += rhoE.ref()[cells[celli_1]];
		rhoY.ref()[celli] += rhoY.ref()[cells[celli_1]];
	}

	rho.ref()[celli] /= cells.size();
	rhoU.ref()[celli] /= cells.size();
	rhoE.ref()[celli] /= cells.size();
	rhoU.ref()[celli] /= cells.size();

	}
}
