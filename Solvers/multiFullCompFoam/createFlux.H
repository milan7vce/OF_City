        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceScalarField rhoY_pos(interpolate(rhoY, pos));
        surfaceScalarField rhoY_neg(interpolate(rhoY, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos));//, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg));//, U.name()));

        
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

        surfaceScalarField p_pos(interpolate(p, pos));
        surfaceScalarField p_neg(interpolate(p, neg));
        surfaceScalarField p_pos1(interpolate(p1, pos));
        surfaceScalarField p_neg1(interpolate(p1, neg));
        surfaceScalarField pSf( 0.5*(p_pos + p_neg) );

        surfaceScalarField Y_pos("Y_pos", rhoY_pos/rho_pos);
        surfaceScalarField Y_neg("Y_neg", rhoY_neg/rho_neg);

      //  surfaceScalarField x_pos(interpolate(dimx, pos));
       // surfaceScalarField x_limit(interpolate(x_limit, pos));



        surfaceScalarField x_dimen1("x_dimen1",mesh.Cf().component(0));
        //surfaceScalarField x_pos(interpolate(mesh.C().component(1),pos));

      // if  ((( cSf_pos < c_limit)&&(Y_pos > Y_limit))|| ((cSf_neg < c_limit)&&(Y_neg > Y_limit)))  {
     
   //   if ((x_dimen1<x_limit))  {
     // if ((p_pos<p_limit)||(p_neg<p_limit))  {
    if ((x_dimen1<x_limit))  {
          rho_pos=rho_pos1;
          rho_neg= rho_neg1;
          rhoY_pos=rhoY_pos1;
          rhoY_neg=rhoY_neg1;
          rhoU_pos=rhoU_pos1;
          rhoU_neg=rhoU_neg1;
          rhoE_pos=rhoE_pos1;
          rhoE_neg=rhoE_neg1;
          cSf_neg=cSf_neg1;
          cSf_pos=cSf_pos1;
      }

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
        
        // this is placeholder..
        surfaceVectorField phiUp ( (rhoPhi_pos*0.0*U_pos + p_pos*0.0*mesh.Sf()/mesh.magSf())*mesh.magSf() );
        surfaceScalarField phiEUp ( (phiv_pos*0.0*(rhoE_pos + p_pos))*mesh.magSf() );// ( v_zero*rhoE_pos );
        
        #include "centralCourantNo.H"
        #include "readTimeControls.H"  
        #include "setDeltaT.H"
        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
                     
       
       volScalarField muEff("muEff", (turbulence->nut()+Nu)*rho);  
        //volScalarField muEff("muEff", turbulence->nuEff()*rho);
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        surfaceScalarField sigmaDotU( //viscous heating term
            (
               fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              //+ (mesh.Sf() & fvc::interpolate(tauMC))
              + fvc::dotInterpolate(mesh.Sf(), tauMC)
            )
            & SstarVec 
            );

        rhoPhi = (rhoPhiM_pos*pos0(Sstar) + rhoPhiM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        rhoYPhi = (rhoPhiYM_pos*pos0(Sstar) + rhoPhiYM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        phiUp = (phiUpM_pos*pos0(Sstar) +  phiUpM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        phiEUp = (phiEUpM_pos*pos0(Sstar) + phiEUpM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
