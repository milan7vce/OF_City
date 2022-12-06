        const volScalarField::Boundary& pbf =  rho.boundaryField(); 
        wordList rhoBoundaryTypes = pbf.types();
       
        forAll(pbf, patchi)
        {        
            if (rhoBoundaryTypes[patchi] == "wall")
            {
                forAll(pbf[patchi], facei)
                {
                    label celli = pbf[patchi].patch().faceCells()[facei];

                    //Density
                    rho_pos.boundaryFieldRef()[patchi][facei] = rho[celli];
                    rho_neg.boundaryFieldRef()[patchi][facei] = rho[celli];
                
                    //Density*y
                    rhoY_pos.boundaryFieldRef()[patchi][facei] = rhoY[celli];
                    rhoY_neg.boundaryFieldRef()[patchi][facei] = rhoY[celli];

                    //DENSITY ENERGY
                    rhoE_pos.boundaryFieldRef()[patchi][facei] = rhoE[celli];
                    rhoE_neg.boundaryFieldRef()[patchi][facei] = rhoE[celli];

                    //Momentum
                    rhoU_pos.boundaryFieldRef()[patchi][facei] = rhoU[celli];
                    rhoU_neg.boundaryFieldRef()[patchi][facei] = -rhoU[celli];
                    //direction of the flow (1.0 = positive)  
                    //a positive value means that the flow is in the direction of the normal vector of the face 
                    // the normal vector is pointing out of the cell
                    
                    //Pressure
                    p_pos.boundaryFieldRef()[patchi][facei] = p[celli];
                    p_neg.boundaryFieldRef()[patchi][facei] = p[celli];

                    //Speed of sound
                    cSf_neg.boundaryFieldRef()[patchi][facei] = c[celli];
                    cSf_pos.boundaryFieldRef()[patchi][facei] = c[celli];
                    
                    c.boundaryFieldRef()[patchi][facei] = c[celli];
                    rho.boundaryFieldRef()[patchi][facei] = rho[celli];
                    Y.boundaryFieldRef()[patchi][facei] = Y[celli];
                    U.boundaryFieldRef()[patchi][facei] = U[celli];
                }
            }else
                 forAll(pbf[patchi], facei)
                {
                    label celli = pbf[patchi].patch().faceCells()[facei];

                    //Density
                    rho_pos.boundaryFieldRef()[patchi][facei] = rho[celli];
                    rho_neg.boundaryFieldRef()[patchi][facei] = rho[celli];
                
                    //Density*y
                    rhoY_pos.boundaryFieldRef()[patchi][facei] = rhoY[celli];
                    rhoY_neg.boundaryFieldRef()[patchi][facei] = rhoY[celli];

                    //DENSITY ENERGY
                    rhoE_pos.boundaryFieldRef()[patchi][facei] = rhoE[celli];
                    rhoE_neg.boundaryFieldRef()[patchi][facei] = rhoE[celli];

                    //Momentum
                    rhoU_pos.boundaryFieldRef()[patchi][facei] = rhoU[celli];
                    rhoU_neg.boundaryFieldRef()[patchi][facei] = rhoU[celli];
                    
                    //Pressure
                    p_pos.boundaryFieldRef()[patchi][facei] = p[celli];
                    p_neg.boundaryFieldRef()[patchi][facei] = p[celli];

                    //Speed of sound
                    cSf_neg.boundaryFieldRef()[patchi][facei] = c[celli];
                    cSf_pos.boundaryFieldRef()[patchi][facei] = c[celli];
                    
                    c.boundaryFieldRef()[patchi][facei] = c[celli];
                    rho.boundaryFieldRef()[patchi][facei] = rho[celli];
                    Y.boundaryFieldRef()[patchi][facei] = Y[celli];
                    U.boundaryFieldRef()[patchi][facei] = U[celli];
                }
        }