//interpolation from uniform structured thermodynamic table
//in OF7, bisogna sostituire .ref() con .dimensionedInternalField() definito in OF24x

forAll(rho, celli){//loop 

// bound Y.ref()[celli]

Y.ref()[celli] = max(0.0,min(Y.ref()[celli],1.0));
//rho.ref()[celli] = max(0.2,rho.ref()[celli]);

	if(rho.ref()[celli]<unif.get_rhoMin() || rho.ref()[celli]>unif.get_rhoMax() || e.ref()[celli]<unif.get_eMin() || e.ref()[celli]>unif.get_eMax() || Y.ref()[celli]<unif.get_YMin() || Y.ref()[celli]>unif.get_YMax() ){
	Info<<"out of table!\n"<<endl; 


	        p_r.ref()[celli]=1000.0;
			T_r.ref()[celli]=330.0;
			c.ref()[celli]= 1.0;
			alphav.ref()[celli]=1;
			Nu.ref()[celli]=1.0;  //change
			Cp.ref()[celli]=1.0;	//change
			k_eff.ref()[celli]=1.0;	//change
			alphav.ref()[celli] =1.0; 	//change


 	// exit(-1);
	}else{

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
			p_r.ref()[celli]=max(approximate(unif, "p_r"),10.0);
			T_r.ref()[celli]=max(330,approximate(unif, "T_r"));
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

			//alphal.ref()[celli] = 1.0 - alphav.ref()[celli];

		//}

}
}//end loop 

   p_r.correctBoundaryConditions();
   c.correctBoundaryConditions();
   T_r.correctBoundaryConditions();
   Nu.correctBoundaryConditions();
   Cp.correctBoundaryConditions();
   k_eff.correctBoundaryConditions();
   alphav.correctBoundaryConditions();
   //alphal.correctBoundaryConditions();
