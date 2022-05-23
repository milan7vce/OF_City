//interpolation from uniform structured thermodynamic table
//in OF7, bisogna sostituire .ref() con .dimensionedInternalField() definito in OF24x

forAll(rho, celli){//loop 

// bound Y.ref()[celli]

Y.ref()[celli] = max(0.0,min(Y.ref()[celli],1.0));
//rho.ref()[celli] = max(0.2,rho.ref()[celli]);

	if(rho.ref()[celli]<unif.get_rhoMin() || rho.ref()[celli]>unif.get_rhoMax() || e.ref()[celli]<unif.get_eMin() || e.ref()[celli]>unif.get_eMax() || Y.ref()[celli]<unif.get_YMin() || Y.ref()[celli]>unif.get_YMax() ){
	Info<<"out of table!\n"<<endl; 


// simple diagnostic information ..
	Info<<"Cell :" <<  celli <<endl; 
	Info<<"Rho :" <<  rho.ref()[celli] <<endl; 
	Info<<"e :" <<  e.ref()[celli] <<endl; 
	Info<<"Y :" <<  Y.ref()[celli] <<endl;
        Info<<"u :" <<  U.ref()[celli] <<endl;
        Info<<"c :" <<  c.ref()[celli] <<endl;
        Info<<"p :" <<  p.ref()[celli] <<endl;
        Info<<"x :" <<  mesh.C()[celli].component(0) <<endl;
                  

	Info<<"rhomin :" <<  unif.get_rhoMin() <<endl; 
  	Info<<"rhomax:" <<  unif.get_rhoMax() <<endl;
	Info<<"emin :" <<  unif.get_eMin() <<endl; 
  	Info<<"emax:" <<  unif.get_eMax() <<endl;
	Info<<"Ymin :" <<  unif.get_YMin() <<endl; 
  	Info<<"Ymax:" <<  unif.get_YMax() <<endl; 

 	 exit(-1);
	}

	unif.initAux(rho,e,Y,celli);

	/*if ( approximate(unif, "T") < 330)
		{
			//Info<<"Cell " <<  celli << " Fails. Its T is " << approximate(unif, "T")<<endl; 
			p.ref()[celli]=p.ref()[celli];
			T.ref()[celli]=T.ref()[celli];
			c.ref()[celli]=c.ref()[celli];
			alphav.ref()[celli]=alphav.ref()[celli];
			Nu.ref()[celli]=Nu.ref()[celli];
			Cp.ref()[celli]=Cp.ref()[celli];
			k_eff.ref()[celli]=k_eff.ref()[celli];

			alphav.ref()[celli] =
				max
				(
					min(
					alphav.ref()[celli], scalar(1)
					),
				scalar(0)
				);

			alphal.ref()[celli] = 1.0 - alphav.ref()[celli];
		}
		else
		{*/
			p.ref()[celli]=max(approximate(unif, "p"),10.0);
			T.ref()[celli]=approximate(unif, "T");
			c.ref()[celli]= max(approximate(unif, "c"),1.0);
			alphav.ref()[celli]=approximate(unif, "alphav");
			Nu.ref()[celli]=approximate(unif, "nu");
			Cp.ref()[celli]=approximate(unif, "Cp");
			k_eff.ref()[celli]=approximate(unif, "k");

			alphav.ref()[celli] =
				max
				(
					min(
					alphav.ref()[celli], scalar(1)
					),
				scalar(0)
				);

			alphal.ref()[celli] = 1.0 - alphav.ref()[celli];

		//}


}//end loop 

   p.correctBoundaryConditions();
   c.correctBoundaryConditions();
   T.correctBoundaryConditions();
   Nu.correctBoundaryConditions();
   Cp.correctBoundaryConditions();
   k_eff.correctBoundaryConditions();
   alphav.correctBoundaryConditions();
   alphal.correctBoundaryConditions();
