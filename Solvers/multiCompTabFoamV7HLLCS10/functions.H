//
//in OF7, bisogna sostituire .ref() con .dimensionedInternalField() definito in OF240

scalar approximate(const tableInfo &tbi, string quantType){

return tbi.thermoInter( tbi.selectQuant.find(quantType)->second);
}



//function for interpolating from the thermodynamic table
scalar tableInfo::thermoInter(scalar *** quantity)const{

/*scalar prop_e1, prop_e2, interpolatedValue;
prop_e1= ( rhoSolver-rhoTable[xPos][yPos] )/( rhoTable[xPos+1][yPos+1]-rhoTable[xPos][yPos] )*( quantity[xPos+1][yPos]-quantity[xPos][yPos] )+quantity[xPos][yPos];
prop_e2= ( rhoSolver-rhoTable[xPos][yPos] )/( rhoTable[xPos+1][yPos+1]-rhoTable[xPos][yPos] )*( quantity[xPos+1][yPos+1]-quantity[xPos][yPos+1] )+quantity[xPos][yPos+1];
interpolatedValue=( eSolver-e_table[xPos][yPos] )/( e_table[xPos+1][yPos+1]-e_table[xPos][yPos] )*( prop_e2-prop_e1 ) + prop_e1;*/

//Info<<"xPos :" <<  xPos <<endl; 
//Info<<"yPos :" <<  yPos <<endl; 
//Info<<"zPos :" <<  zPos <<endl; 

//Info<<"Rho :" <<  rhoSolver << rhoTable[xPos][yPos][zPos]  << rhoTable[xPos+1][yPos+1][zPos] <<endl; 
//Info<<"e :" <<  eSolver << eTable[xPos][yPos][zPos]  << eTable[xPos+1][yPos+1][zPos] <<endl; 
//Info<<"y :" <<  YSolver << YTable[xPos][yPos][zPos]  << YTable[xPos+1][yPos+1][zPos] <<endl; 

//Info<<"Rho :" <<  rhoSolver << rhoTable[xPos][yPos][zPos+1]  << rhoTable[xPos+1][yPos+1][zPos+1] <<endl; 
//Info<<"e :" <<  eSolver << eTable[xPos][yPos][zPos+1]  << eTable[xPos+1][yPos+1][zPos+1] <<endl; 
//Info<<"y :" <<  YSolver << YTable[xPos][yPos][zPos+1]  << YTable[xPos+1][yPos+1][zPos+1] <<endl; 
//cin.get();


	scalar prop_1, prop_2, prop_3, prop_4, var_1, var_2, interpolatedValue;

    prop_1 = ( rhoSolver - rhoTable[xPos][yPos][zPos] )/( rhoTable[xPos+1][yPos+1][zPos] - rhoTable[xPos][yPos][zPos] )*
			( quantity[xPos+1][yPos][zPos] - quantity[xPos][yPos][zPos] ) + quantity[xPos][yPos][zPos];
	prop_2 = ( rhoSolver - rhoTable[xPos][yPos][zPos] )/( rhoTable[xPos+1][yPos+1][zPos] - rhoTable[xPos][yPos][zPos] )*
			( quantity[xPos+1][yPos+1][zPos] - quantity[xPos][yPos+1][zPos] ) + quantity[xPos][yPos+1][zPos];
	/*prop_1 = ( rhoSolver - rhoTable[xPos][yPos][zPos] )/( rhoTable[xPos+1][yPos][zPos] - rhoTable[xPos][yPos][zPos] )*
			( quantity[xPos+1][yPos][zPos] - quantity[xPos][yPos][zPos] ) + quantity[xPos][yPos][zPos];
	prop_2 = ( rhoSolver - rhoTable[xPos][yPos+1][zPos] )/( rhoTable[xPos+1][yPos+1][zPos] - rhoTable[xPos][yPos+1][zPos] )*
			( quantity[xPos+1][yPos+1][zPos] - quantity[xPos][yPos+1][zPos] ) + quantity[xPos][yPos+1][zPos];*/

	var_1 = ( eSolver - eTable[xPos][yPos][zPos] )/( eTable[xPos+1][yPos+1][zPos] - eTable[xPos][yPos][zPos] )*( prop_2 - prop_1) + prop_1;

    prop_3 = ( rhoSolver - rhoTable[xPos][yPos][zPos+1] )/( rhoTable[xPos+1][yPos+1][zPos+1] - rhoTable[xPos][yPos][zPos+1] )*
			( quantity[xPos+1][yPos][zPos+1] - quantity[xPos][yPos][zPos+1] ) + quantity[xPos][yPos][zPos+1];
	prop_4 = ( rhoSolver - rhoTable[xPos][yPos][zPos+1] )/( rhoTable[xPos+1][yPos+1][zPos+1] - rhoTable[xPos][yPos][zPos+1] )*
			( quantity[xPos+1][yPos+1][zPos+1] - quantity[xPos][yPos+1][zPos+1] ) + quantity[xPos][yPos+1][zPos+1];
	/*prop_3 = ( rhoSolver - rhoTable[xPos][yPos][zPos+1] )/( rhoTable[xPos+1][yPos][zPos+1] - rhoTable[xPos][yPos][zPos+1] )*
			( quantity[xPos+1][yPos][zPos+1] - quantity[xPos][yPos][zPos+1] ) + quantity[xPos][yPos][zPos+1];
	prop_4 = ( rhoSolver - rhoTable[xPos][yPos+1][zPos+1] )/( rhoTable[xPos+1][yPos+1][zPos+1] - rhoTable[xPos][yPos+1][zPos+1] )*
			( quantity[xPos+1][yPos+1][zPos+1] - quantity[xPos][yPos+1][zPos+1] ) + quantity[xPos][yPos+1][zPos+1];*/

	var_2 = ( eSolver - eTable[xPos][yPos][zPos+1] )/( eTable[xPos+1][yPos+1][zPos+1] - eTable[xPos][yPos][zPos+1] )*( prop_4 - prop_3 ) + prop_3;

	interpolatedValue = ( YSolver - YTable[xPos][yPos][zPos] )/( YTable[xPos+1][yPos+1][zPos+1] - YTable[xPos][yPos][zPos] )*( var_2 - var_1 ) + var_1;


//Info<<"Val = " <<  interpolatedValue <<endl; 
//cin.get();


	return interpolatedValue;
}



//calculate manipulation parameters based on the CFD solution
void tableInfo::initAux(volScalarField& rho, volScalarField& e, volScalarField& Y, label celli){//initialize auxialiary members
eSolver = e.ref()[celli];
YSolver = Y.ref()[celli];


switch(format){//switch

case 0:
rhoTable=rho_table;
rhoSolver=rho.ref()[celli];
xPos=static_cast<int> (  ( rhoSolver-rhoMin )/dRho  );
break;

case 1:
rhoTable=rhoLog_table;
rhoSolver=Foam::log10( rho.ref()[celli] );
xPos=static_cast<int> (  ( rhoSolver- Foam::log10(rhoMin) )/dRho  );
break;

default:
Info<<"wrong format of table in x-axis \n"<<endl;
exit(-1);
break;

}//end switch

eTable = e_table;
YTable = Y_table;

//xPos=static_cast<int> (  ( rhoSolver-rhoMin )/dRho  );
yPos=static_cast<int> (  ( eSolver-eMin )/de  );
zPos=static_cast<int> (  ( YSolver-YMin )/dY  );

// limit x,y,zPos -- note max value is *Node - 1 !!! 
zPos = max(0,min(YNodes-2,zPos));
yPos = max(0,min(eNodes-2,yPos));
xPos = max(0,min(rhoNodes-2,xPos));


//Info<<"xPos :" <<  xPos <<endl; 
//Info<<"yPos :" <<  yPos <<endl; 
//Info<<"zPos :" <<  zPos <<endl; 
//cin.get();

}//end of function



//constructor: reads from structured equally spaced thermodynamic grid
tableInfo::tableInfo(word fileName){ //constructor
rhoMin=1.e10; rhoMax=-1.e10; eMin=1.e10; eMax=-1.e10; YMin=1.e10; YMax=-1.e10;

table.open(fileName.c_str());

	if (!table) {
	Info << "\nProblem reading the thermodynamic table "<< fileName <<" :(\n" <<endl;
	exit(-1);
	}

table>> rhoNodes >> eNodes >> YNodes >> dRho >> de >> dY >> MW >> format; //read first line

p_table = new scalar **[rhoNodes] ;
T_table = new scalar **[rhoNodes] ;
rho_table = new scalar **[rhoNodes] ;
e_table = new scalar **[rhoNodes] ;
h_table = new scalar **[rhoNodes] ;
s_table = new scalar **[rhoNodes] ;
c_table = new scalar **[rhoNodes] ;
alphav_table = new scalar **[rhoNodes] ;
rhoLog_table = new scalar **[rhoNodes] ;
Cp_table = new scalar **[rhoNodes] ;
k_table = new scalar **[rhoNodes] ;
mu_table = new scalar **[rhoNodes] ;
Y_table = new scalar **[rhoNodes] ;
nu_table = new scalar **[rhoNodes] ;

for(label i=0; i<rhoNodes; i++){
	p_table[i]=new scalar *[eNodes];
	T_table[i]=new scalar *[eNodes];
	rho_table[i]=new scalar *[eNodes];
	e_table[i]=new scalar *[eNodes];
	h_table[i]=new scalar *[eNodes];
	s_table[i]=new scalar *[eNodes];
	c_table[i]=new scalar *[eNodes];
	alphav_table[i]=new scalar *[eNodes];
	rhoLog_table[i]=new scalar *[eNodes];
	Cp_table[i]=new scalar *[eNodes];
	k_table[i]=new scalar *[eNodes];
	mu_table[i]=new scalar *[eNodes];
	Y_table[i]=new scalar *[eNodes];
	nu_table[i]=new scalar *[eNodes];
}

for(label i=0; i<rhoNodes; i++){
	for(label j=0; j<eNodes; j++){
		p_table[i][j]=new scalar[YNodes];
		T_table[i][j]=new scalar[YNodes];
		rho_table[i][j]=new scalar[YNodes];
		e_table[i][j]=new scalar[YNodes];
		h_table[i][j]=new scalar[YNodes];
		s_table[i][j]=new scalar[YNodes];
		c_table[i][j]=new scalar[YNodes];
		alphav_table[i][j]=new scalar[YNodes];
		rhoLog_table[i][j]=new scalar[YNodes];
		Cp_table[i][j]=new scalar[YNodes];
		k_table[i][j]=new scalar[YNodes];
		mu_table[i][j]=new scalar[YNodes];
		Y_table[i][j]=new scalar[YNodes];
		nu_table[i][j]=new scalar[YNodes];
	}
}



for(label k=0; k<YNodes; k++){
	for(label j=0; j<eNodes; j++){
		for(label i=0; i<rhoNodes; i++){
			
		table>> (rhoLog_table[i][j][k]) >> (rho_table[i][j][k]) >> (e_table[i][j][k]) >> (Y_table[i][j][k]) >> (T_table[i][j][k]) >> (p_table[i][j][k]) >> (h_table[i][j][k]) >> (s_table[i][j][k]) >> (Cp_table[i][j][k]) >> (c_table[i][j][k]) >> (k_table[i][j][k]) >> (mu_table[i][j][k]) >> (alphav_table[i][j][k]) ; //read the rest the thermodynamies properties
		nu_table[i][j][k]=mu_table[i][j][k]/rho_table[i][j][k];
		rhoMin=min(rho_table[i][j][k],rhoMin);
		rhoMax=max(rho_table[i][j][k],rhoMax);
		eMin=min(e_table[i][j][k],eMin);
		eMax=max(e_table[i][j][k],eMax);
		YMin=min(Y_table[i][j][k],YMin);
		YMax=max(Y_table[i][j][k],YMax);
		}
	}
}

table.close();
Info<< "\nThermodynamic table "<< fileName <<" has been loaded" << endl;

Info<< "\n Displaying limits: " << endl;
Info<< "\n Rho: " << rhoMin << ".." << rhoMax << endl;
Info<< "\n e: " << eMin << ".." << eMax << endl;
Info<< "\n y: " << YMin << ".." << YMax << endl;


//fill in the container in order to distinct which thermodynamic variable is interpolated each time based on rho,e
selectQuant.insert (it, std::pair<string,scalar ***>("p",p_table));
selectQuant.insert (it, std::pair<string,scalar ***>("T",T_table));
selectQuant.insert (it, std::pair<string,scalar ***>("c",c_table));
selectQuant.insert (it, std::pair<string,scalar ***>("alphav",alphav_table));
selectQuant.insert (it, std::pair<string,scalar ***>("nu",nu_table));
selectQuant.insert (it, std::pair<string,scalar ***>("Cp",Cp_table));
selectQuant.insert (it, std::pair<string,scalar ***>("k",k_table));

} //end of constructor



//release memory used for the thermodynamic table
tableInfo::~tableInfo(){ //destructor

	for(label i=0; i<rhoNodes; i++){
		for(label j=0; j<eNodes; j++){
			delete [] T_table[i][j];
			delete [] rho_table[i][j];
			delete [] e_table[i][j];
			delete [] h_table[i][j];
			delete [] s_table[i][j];
			delete [] c_table[i][j];
			delete [] alphav_table[i][j];
			delete [] rhoLog_table[i][j];
			delete [] Cp_table[i][j];
			delete [] k_table[i][j];
			delete [] mu_table[i][j];
			delete [] Y_table[i][j];
			delete [] nu_table[i][j];
		}
	}
	
	for(label i=0; i<rhoNodes; i++){
	delete [] T_table[i];
	delete [] rho_table[i];
	delete [] e_table[i];
	delete [] h_table[i];
	delete [] s_table[i];
	delete [] c_table[i];
	delete [] alphav_table[i];
	delete [] rhoLog_table[i];
	delete [] Cp_table[i];
	delete [] k_table[i];
	delete [] mu_table[i];
	delete [] Y_table[i];
	delete [] nu_table[i];
	}
delete [] T_table;
delete [] rho_table;
delete [] e_table;
delete [] h_table;
delete [] s_table;
delete [] c_table;
delete [] alphav_table;
delete [] rhoLog_table;
delete [] Cp_table;
delete [] k_table;
delete [] mu_table;
delete [] Y_table;
delete [] nu_table;

//Info<<"\nMEMORY RELEASED\n"<<endl;

} //end of destructor
