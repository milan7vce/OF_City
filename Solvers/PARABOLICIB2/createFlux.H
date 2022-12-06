 	#include "readTimeControls.H"

        if (!LTS)
        {
            #include "setDeltaT.H"
            runTime++;

            // Do any mesh changes
            //mesh.update();
        }

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceScalarField rhoY_pos(interpolate(rhoY, pos));
        surfaceScalarField rhoY_neg(interpolate(rhoY, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos));//, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg));//, U.name()));
	
	volVectorField::Boundary& UBf  = U.boundaryFieldRef();
        forAll(UBf, pachi)
        {
           if((UBf[pachi].patch().name() == "nozzle1") || (UBf[pachi].patch().name() == "nozzlle2")|| (UBf[pachi].patch().name() == "needle"))
         {
             rhoU_neg.boundaryFieldRef()[pachi] = -rhoU_neg.boundaryFieldRef()[pachi];
         
	}
	}

        surfaceScalarField rho_pos1(interpolate(rho1, pos));
        surfaceScalarField rho_neg1(interpolate(rho1, neg));

        surfaceScalarField rhoY_pos1(interpolate(rhoY1, pos));
        surfaceScalarField rhoY_neg1(interpolate(rhoY1, neg));
		
        surfaceVectorField rhoU_pos1(interpolate(rhoU1, pos));//, U1.name()));
        surfaceVectorField rhoU_neg1(interpolate(rhoU1, neg));//, U1.name()));
		
        surfaceScalarField  alpha_pos(interpolate(alphav, pos));//, U1.name()));
	    surfaceScalarField  alpha_neg(interpolate(alphav, neg));//, U1.name()));

        surfaceScalarField rhoE_pos(interpolate(rhoE, pos) );
        surfaceScalarField rhoE_neg(interpolate(rhoE, neg) );

        surfaceScalarField rhoE_pos1(interpolate(rhoE1, pos) );
        surfaceScalarField rhoE_neg1(interpolate(rhoE1, neg) );

        surfaceScalarField cSf_pos("cSf_pos",interpolate(c, pos) );
        surfaceScalarField cSf_neg("cSf_neg",interpolate(c, neg) );

        surfaceScalarField cSf_pos1("cSf_pos1",interpolate(c1, pos) );
        surfaceScalarField cSf_neg1("cSf_neg1",interpolate(c1, neg) );

        surfaceScalarField p_pos(interpolate(p_r, pos));
        surfaceScalarField p_neg(interpolate(p_r, neg));
        surfaceScalarField p_pos1(interpolate(p1, pos));
        surfaceScalarField p_neg1(interpolate(p1, neg));
        //surfaceScalarField pSf( 0.5*(p_pos + p_neg) );


         surfaceScalarField x_dimen1("x_dimen1",mesh.Cf().component(0));

           forAll(cSf_pos.primitiveFieldRef(),iface)  {
    if (x_dimen1.primitiveFieldRef()[iface]>x_limit.primitiveFieldRef()[iface])  {
          rho_pos.primitiveFieldRef()[iface]=rho_pos1.primitiveFieldRef()[iface];
          rho_neg.primitiveFieldRef()[iface]= rho_neg1.primitiveFieldRef()[iface];
          rhoY_pos.primitiveFieldRef()[iface]=rhoY_pos1.primitiveFieldRef()[iface];
          rhoY_neg.primitiveFieldRef()[iface]=rhoY_neg1.primitiveFieldRef()[iface];
          rhoU_pos.primitiveFieldRef()[iface]=rhoU_pos1.primitiveFieldRef()[iface];
          rhoU_neg.primitiveFieldRef()[iface]=rhoU_neg1.primitiveFieldRef()[iface];
          rhoE_pos.primitiveFieldRef()[iface]=rhoE_pos1.primitiveFieldRef()[iface];
          rhoE_neg.primitiveFieldRef()[iface]=rhoE_neg1.primitiveFieldRef()[iface];
          cSf_neg.primitiveFieldRef()[iface]=cSf_neg1.primitiveFieldRef()[iface];
          cSf_pos.primitiveFieldRef()[iface]=cSf_pos1.primitiveFieldRef()[iface];
      }
}
//      #include "correctfaces.h"

        surfaceScalarField Y_pos("Y_pos", rhoY_pos/rho_pos);
        surfaceScalarField Y_neg("Y_neg", rhoY_neg/rho_neg);
 
        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField phiv_pos("phiv_pos", (U_pos & mesh.Sf())/mesh.magSf() );
        surfaceScalarField phiv_neg("phiv_neg", (U_neg & mesh.Sf())/mesh.magSf() );

        surfaceScalarField cSf ( max( max(cSf_pos,cSf_neg), uBounded ) );
          
        surfaceScalarField c_bar ( sqrt(cSf_pos*cSf_neg) );  // If we are to follow this path, we can choose sth similar to Roe - Pike average..
        
        surfaceScalarField s_neg( max(phiv_pos + cSf_pos , phiv_neg + cSf_neg) );  //phiv_neg + cSf_neg);//  phiv_bar + c_bar ;   //
        surfaceScalarField s_pos( min(phiv_pos - cSf_pos , phiv_neg - cSf_neg) );  //phiv_pos - cSf_pos);// phiv_bar - c_bar ; //

        surfaceScalarField Sstar ( "Sstar", (rho_pos*phiv_pos*(s_pos - phiv_pos) - rho_neg*phiv_neg*(s_neg - phiv_neg) + (p_neg - p_pos))
                                            / ( (rho_pos*(s_pos - phiv_pos+ uBounded) - rho_neg*(s_neg - phiv_neg-uBounded) )) );

        surfaceVectorField SstarVec("SstarVec", (rho_pos*U_pos*(s_pos - phiv_pos) - rho_neg*U_neg*(s_neg - phiv_neg) + mesh.Sf()*(p_neg - p_pos)/mesh.magSf() )
                                            / (rho_pos*(s_pos - phiv_pos+ uBounded) - rho_neg*(s_neg - phiv_neg-uBounded)) );

        surfaceScalarField pStar_pos ( rho_pos*(phiv_pos - s_pos)*(phiv_pos - Sstar) + p_pos );
        surfaceScalarField pStar_neg ( rho_neg*(phiv_neg - s_neg)*(phiv_neg - Sstar) + p_neg );
        surfaceScalarField pLR ( 0.5*(pStar_neg+pStar_pos) );

        // define regular L(pos)/R(neg) fluxes
        surfaceScalarField rhoPhi_pos (phiv_pos*rho_pos);
        surfaceScalarField rhoPhi_neg (phiv_neg*rho_neg);

        surfaceScalarField rhoPhiY_pos (phiv_pos*rhoY_pos);
        surfaceScalarField rhoPhiY_neg (phiv_neg*rhoY_neg);

        surfaceVectorField phiUp_pos ( rhoPhi_pos*U_pos + p_pos*mesh.Sf()/mesh.magSf() );
        surfaceVectorField phiUp_neg ( rhoPhi_neg*U_neg + p_neg*mesh.Sf()/mesh.magSf() );

        surfaceScalarField phiEUp_pos ( phiv_pos*(rhoE_pos + p_pos) );
        surfaceScalarField phiEUp_neg ( phiv_neg*(rhoE_neg + p_neg) );


        // use formulation 10.73 / 10.39 from Toro for U*L/R 
        surfaceScalarField u_rho_pos ( rho_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded)  );
        surfaceScalarField u_rho_neg ( rho_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded)  );

        surfaceScalarField u_rhoY_pos ( rhoY_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded) );
        surfaceScalarField u_rhoY_neg ( rhoY_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded) );

        surfaceVectorField u_rhoUp_pos  ( rho_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded) * (U_pos - (phiv_pos - Sstar)*mesh.Sf()/mesh.magSf()  ));
        surfaceVectorField u_rhoUp_neg  ( rho_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded) * (U_neg - (phiv_neg - Sstar)*mesh.Sf()/mesh.magSf()  ));

        surfaceScalarField u_rhoEUp_pos ( rho_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded) 
                                                  *(rhoE_pos/rho_pos + (Sstar - phiv_pos)*(Sstar + p_pos/(rho_pos*(s_pos - phiv_pos+uBounded)) ) )  );
												  
        surfaceScalarField u_rhoEUp_neg ( rho_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded) 
                                                  *(rhoE_neg/rho_neg + (Sstar - phiv_neg)*(Sstar + p_neg/(rho_neg*(s_neg - phiv_neg +uBounded)) ) )  );

        // define STAR L(pos)/R(neg) fluxes, eq. 10.38
        surfaceScalarField rhoPhiS_pos ( rhoPhi_pos + s_pos*(u_rho_pos - rho_pos )  );
        surfaceScalarField rhoPhiS_neg ( rhoPhi_neg + s_neg*(u_rho_neg - rho_neg )  );

        surfaceScalarField rhoPhiYS_pos ( rhoPhiY_pos + s_pos*(u_rhoY_pos - rhoY_pos )  );
        surfaceScalarField rhoPhiYS_neg ( rhoPhiY_neg + s_neg*(u_rhoY_neg - rhoY_neg )  );

        surfaceVectorField phiUpS_pos ( phiUp_pos + s_pos*(u_rhoUp_pos - rhoU_pos)  ) ;
        surfaceVectorField phiUpS_neg ( phiUp_neg + s_neg*(u_rhoUp_neg - rhoU_neg)  ) ;            

        surfaceScalarField phiEUpS_pos ( phiEUp_pos + s_pos*(u_rhoEUp_pos - rhoE_pos)  );
        surfaceScalarField phiEUpS_neg ( phiEUp_neg + s_neg*(u_rhoEUp_neg - rhoE_neg)  );
        

        // define mixed L(pos)/R(neg) fluxes 
        surfaceScalarField rhoPhiM_pos ( (rhoPhi_pos*pos0(s_pos) + rhoPhiS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos)  ));
        surfaceScalarField rhoPhiM_neg ( (rhoPhi_neg*neg0(s_neg) + rhoPhiS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg)  ));

        surfaceScalarField rhoPhiYM_pos ( (rhoPhiY_pos*pos0(s_pos) + rhoPhiYS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos)  ));
        surfaceScalarField rhoPhiYM_neg ( (rhoPhiY_neg*neg0(s_neg) + rhoPhiYS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg)  ));

        surfaceVectorField phiUpM_pos ( (phiUp_pos*pos0(s_pos) + phiUpS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos) )) ;
        surfaceVectorField phiUpM_neg ( (phiUp_neg*neg0(s_neg) + phiUpS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg) )) ;            

        surfaceScalarField phiEUpM_pos ( (phiEUp_pos*pos0(s_pos) + phiEUpS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos)));
        surfaceScalarField phiEUpM_neg ( (phiEUp_neg*neg0(s_neg) + phiEUpS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg) ));
        
        surfaceVectorField phiUp ( (rhoPhi_pos*0.0*U_pos + p_pos*0.0*mesh.Sf()/mesh.magSf())*mesh.magSf() );
        surfaceScalarField phiEUp ( (phiv_pos*0.0*(rhoE_pos + p_pos))*mesh.magSf() );// ( v_zero*rhoE_pos );
	    surfaceScalarField ucSf( mag( fvc::interpolate(U) & mesh.Sf() )+fvc::interpolate(mag(c))*mesh.magSf() );
        
        #include "centralCourantNo.H"

	    if (LTS)
        {
        #include "setRDeltaT.H"
	    runTime.setDeltaT(maxDeltaT);
        runTime++;
        } 

        Info<< "Time = " << runTime.timeName() << nl << endl;
         
        rhoPhi = (rhoPhiM_pos*pos0(Sstar) + rhoPhiM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        rhoYPhi = (rhoPhiYM_pos*pos0(Sstar) + rhoPhiYM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        phiUp = (phiUpM_pos*pos0(Sstar) +  phiUpM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        phiEUp = (phiEUpM_pos*pos0(Sstar) + phiEUpM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();


      volScalarField muEff("muEff", (turbulence->nut()+Nu)*rho) ;
      volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        surfaceScalarField sigmaDotU(
            (
               fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)

              + fvc::dotInterpolate(mesh.Sf(), tauMC)
            )
            & SstarVec
            );

