
        rho_pos = (interpolate(rho, pos));
        rho_neg = (interpolate(rho, neg));

        rhoU_pos = (interpolate(rhoU, pos));//, U.name()));
        rhoU_neg = (interpolate(rhoU, neg));//, U.name()));

        forAll(UBf, pachi)
        {
           if((UBf[pachi].patch().name() == "nozzle1") || (UBf[pachi].patch().name() == "nozzlle2") || (UBf[pachi].patch().name() == "needle"))
         {
             rhoU_neg.boundaryFieldRef()[pachi] = -rhoU_neg.boundaryFieldRef()[pachi];
        } 
	}
		
        rhoY_pos = (interpolate(rhoY, pos));
        rhoY_neg = (interpolate(rhoY, neg));
		
	rho_pos1 = (interpolate(rho1, pos));
        rho_neg1 = (interpolate(rho1, neg));

        rhoU_pos1 = (interpolate(rhoU1, pos));
        rhoU_neg1 = (interpolate(rhoU1, neg));
		
	rhoY_pos1 = (interpolate(rhoY1, pos));
        rhoY_neg1 = (interpolate(rhoY1, neg));

        rhoE_pos = interpolate(rhoE, pos);
        rhoE_neg = interpolate(rhoE, neg);

	rhoE_pos1 = interpolate(rhoE1, pos);
        rhoE_neg1 = interpolate(rhoE1, neg);

        cSf_pos = interpolate(c, pos);
        cSf_neg = interpolate(c, neg);

        cSf_pos1 = interpolate(c1, pos);
        cSf_neg1 = interpolate(c1, neg);

        alpha_pos=interpolate(alphav, pos);//, U1.name()));
        alpha_neg=interpolate(alphav, neg);//, U1.name()));

        p_pos = (interpolate(p_r, pos));
        p_neg = (interpolate(p_r, neg));

        p_pos1 = (interpolate(p1, pos));
        p_neg1 = (interpolate(p1, neg));


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
        
        Y_pos = rhoY_pos/rho_pos;
        Y_neg = rhoY_neg/rho_neg;
        
        U_pos = rhoU_pos/rho_pos;
        U_neg = rhoU_neg/rho_neg;

        phiv_pos = (U_pos & mesh.Sf())/mesh.magSf() ;// /mesh.magSf()
        phiv_neg = (U_neg & mesh.Sf())/mesh.magSf() ;// /mesh.magSf()
        cSf = max( max(cSf_pos, cSf_neg), uBounded );

        c_bar = ( sqrt(cSf_pos*cSf_neg) );  // If we are to follow this path, we can choose sth similar to Roe - Pike average..
        
        s_neg = ( max(phiv_pos + cSf_pos , phiv_neg + cSf_neg) );  //phiv_neg + cSf_neg);//  phiv_bar + c_bar ;   //
        s_pos = ( min(phiv_pos - cSf_pos , phiv_neg - cSf_neg) );  //phiv_pos - cSf_pos);// phiv_bar - c_bar ; //
		

        Sstar = ( (rho_pos*phiv_pos*(s_pos - phiv_pos) - rho_neg*phiv_neg*(s_neg - phiv_neg) + (p_neg - p_pos))
                                            / ( (rho_pos*(s_pos - phiv_pos + uBounded) - rho_neg*(s_neg - phiv_neg-uBounded) )) );

        SstarVec = ( (rho_pos*U_pos*(s_pos - phiv_pos) - rho_neg*U_neg*(s_neg - phiv_neg) + mesh.Sf()*(p_neg - p_pos)/mesh.magSf() )
                                            / (rho_pos*(s_pos - phiv_pos+uBounded) - rho_neg*(s_neg - phiv_neg-uBounded)) );
											

        pStar_pos = ( rho_pos*(phiv_pos - s_pos)*(phiv_pos - Sstar) + p_pos );
        pStar_neg = ( rho_neg*(phiv_neg - s_neg)*(phiv_neg - Sstar) + p_neg );
        pLR = ( 0.5*(pStar_neg+pStar_pos) );

        // define regular L(pos)/R(neg) fluxes
        rhoPhi_pos =(phiv_pos*rho_pos);
        rhoPhi_neg =(phiv_neg*rho_neg);

        rhoPhiY_pos =(phiv_pos*rhoY_pos);
        rhoPhiY_neg =(phiv_neg*rhoY_neg);

        phiUp_pos =( rhoPhi_pos*U_pos + p_pos*mesh.Sf()/mesh.magSf() );
        phiUp_neg =( rhoPhi_neg*U_neg + p_neg*mesh.Sf()/mesh.magSf() );

        phiEUp_pos =( phiv_pos*(rhoE_pos + p_pos) );
        phiEUp_neg =( phiv_neg*(rhoE_neg + p_neg) );
		

        // use formulation 10.73 / 10.39 from Toro for U*L/R
        u_rho_pos =( rho_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded)  );
        u_rho_neg =( rho_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded)  );

        u_rhoY_pos =( rhoY_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded) );
        u_rhoY_neg =( rhoY_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded) );
		

        u_rhoUp_pos  =( rho_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded) * (U_pos - (phiv_pos - Sstar)*mesh.Sf()/mesh.magSf()  ));
        u_rhoUp_neg  =( rho_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded) * (U_neg - (phiv_neg - Sstar)*mesh.Sf()/mesh.magSf()  ));

        u_rhoEUp_pos =( rho_pos*(s_pos - phiv_pos )/ (s_pos-Sstar+ uBounded) 
                                                  *(rhoE_pos/rho_pos + (Sstar - phiv_pos)*(Sstar + p_pos/(rho_pos*(s_pos - phiv_pos - uBounded)) ) )  );
        u_rhoEUp_neg =( rho_neg*(s_neg - phiv_neg )/ (s_neg-Sstar+ uBounded) 
                                                  *(rhoE_neg/rho_neg + (Sstar - phiv_neg)*(Sstar + p_neg/(rho_neg*(s_neg - phiv_neg+ uBounded)) ) )  );
												  
  
        // define STAR L(pos)/R(neg) fluxes, eq. 10.38
        rhoPhiS_pos =( rhoPhi_pos + s_pos*(u_rho_pos - rho_pos )  );
        rhoPhiS_neg =( rhoPhi_neg + s_neg*(u_rho_neg - rho_neg )  );

        rhoPhiYS_pos =( rhoPhiY_pos + s_pos*(u_rhoY_pos - rhoY_pos )  );
        rhoPhiYS_neg =( rhoPhiY_neg + s_neg*(u_rhoY_neg - rhoY_neg )  );

        phiUpS_pos =( phiUp_pos + s_pos*(u_rhoUp_pos - rhoU_pos)  ) ;
        phiUpS_neg =( phiUp_neg + s_neg*(u_rhoUp_neg - rhoU_neg)  ) ;            

        phiEUpS_pos =( phiEUp_pos + s_pos*(u_rhoEUp_pos - rhoE_pos)  );
        phiEUpS_neg =( phiEUp_neg + s_neg*(u_rhoEUp_neg - rhoE_neg)  );
		

        // define mixed L(pos)/R(neg) fluxes 
        rhoPhiM_pos =( rhoPhi_pos*pos0(s_pos) + rhoPhiS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos) );
        rhoPhiM_neg =( rhoPhi_neg*neg0(s_neg) + rhoPhiS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg) );

        rhoPhiYM_pos =( rhoPhiY_pos*pos0(s_pos) + rhoPhiYS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos) );
        rhoPhiYM_neg =( rhoPhiY_neg*neg0(s_neg) + rhoPhiYS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg) );

        phiUpM_pos =(  phiUp_pos*pos0(s_pos) + phiUpS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos)  ) ;
        phiUpM_neg =(  phiUp_neg*neg0(s_neg) + phiUpS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg)) ;            

        phiEUpM_pos =(  phiEUp_pos*pos0(s_pos) + phiEUpS_pos*neg0(s_pos) )/(pos0(s_pos) + neg0(s_pos) );
        phiEUpM_neg =(  phiEUp_neg*neg0(s_neg) + phiEUpS_neg*pos0(s_neg) )/(neg0(s_neg) + pos0(s_neg));

        rhoPhi = (rhoPhiM_pos*pos0(Sstar) + rhoPhiM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar))*mesh.magSf();
        rhoYPhi = (rhoPhiYM_pos*pos0(Sstar) + rhoPhiYM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        phiUp = (phiUpM_pos*pos0(Sstar) +  phiUpM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();
        phiEUp = (phiEUpM_pos*pos0(Sstar) + phiEUpM_neg*neg0(Sstar) )/(pos0(Sstar) + neg0(Sstar) )*mesh.magSf();

        muEff = (turbulence->nut()+Nu)*rho;
        tauMC = muEff*dev2(Foam::T(fvc::grad(U)));

        sigmaDotU =
                    (
                (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)

                + fvc::dotInterpolate(mesh.Sf(), tauMC)
                )
                & SstarVec
                    );

