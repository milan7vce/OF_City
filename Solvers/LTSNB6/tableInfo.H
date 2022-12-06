#include <map>

class tableInfo{

private:
//table parameters
ifstream table;
scalar ***T_table, ***p_table, ***rho_table, ***e_table, ***h_table, ***s_table, ***c_table, ***alphav_table, ***rhoLog_table, ***Cp_table, ***k_table, ***mu_table, ***nu_table, ***Y_table; //table values at the nodes
label rhoNodes, eNodes, YNodes, format; //format=0: delta(rho)=ct; format=1: delta( log10(rho) )=ct
scalar dRho, de, dY, MW;
scalar rhoMin, rhoMax, eMin, eMax, YMin, YMax;

scalar thermoInter(scalar *** )const; //function for interpolating from the thermodynamic table


public:
//manipulation parameters
scalar rhoSolver, eSolver, YSolver, ***rhoTable, ***eTable, ***YTable;
label xPos, yPos, zPos;

std::map<string,scalar ***> selectQuant; //container for specifying which variable is interpolated
std::map<string,scalar ***>::iterator it = selectQuant.begin();  //iterator for the above container

tableInfo(word); //constructor
~tableInfo();  //destructor

void initAux(volScalarField&, volScalarField&, volScalarField&, label);  //calculate manipulation parameters based on the CFD solution

//access functions
scalar get_rhoMin(){ return rhoMin; }
scalar get_eMin(){ return eMin; }
scalar get_YMin(){ return YMin; }
scalar get_rhoMax(){ return rhoMax; }
scalar get_eMax(){ return eMax; }
scalar get_YMax(){ return YMax; }

friend scalar approximate(const tableInfo&, string);  //friend function; switch between different thermoInter calls

};

