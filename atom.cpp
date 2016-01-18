#include "atom.h"
#include "random_park.h"
#include "constants.h"
#include "mainheader.h"

enum{NONE, BERENDSEN};

Atom::Atom() {


	//defaults
	thermo_created=0;
	thermo_type=NONE;
	initial_temp = 1.0;
	final_temp = 1.0;

	PE=KE=p_PE=p_KE=0;
	

}

Atom::~Atom() {

}

void Atom::set_num_atoms(int n) {
	natoms = n;
	
	for (int i=0; i<natoms; i++) {
		type.push_back(1);
		vx.push_back(0.0);
		vy.push_back(0.0);
		vz.push_back(0.0);
		rx.push_back(0.0);
		ry.push_back(0.0);
		rz.push_back(0.0);
		fx.push_back(0.0);
		fy.push_back(0.0);
		fz.push_back(0.0);
		
	}

}

void Atom::initialize_velocities() {

	RanPark *random = new RanPark(18918);

	double sqrt_T = sqrt(initial_temp);
	
	for(int i=0; i<natoms; i++) {
		vx[i]=sqrt_T*random->gaussian();
		vy[i]=sqrt_T*random->gaussian();
		vz[i]=sqrt_T*random->gaussian();
	}
	delete random;

}


void Atom::calc_force() {


	PE=0;
	p_PE=0;

	double RadiusSQ=rcut*rcut;
	double Xr, Yr, Zr, Xr2, Yr2, Zr2;
	double RijSQ;
	double r2i, r6i;
	double Fij, Fxij, Fyij, Fzij;
	double force;

	for (int i=0; i<natoms; i++) {
		fx[i]=fy[i]=fz[i]= 0.0;
	}

	for(int i=0; i<(natoms-1); i++) {

		for(int j=i+1; j<natoms; j++) {

			//Calculate distance between atoms
			Xr = rx[i] - rx[j];
			Yr = ry[i] - ry[j];
			Zr = rz[i] - rz[j];
			Xr = Xr - box_x*round(Xr/box_x); //periodic boundary condition
			Yr = Yr - box_y*round(Yr/box_y);
			Zr = Zr - box_z*round(Zr/box_z);
 
			Xr2 = Xr*Xr; //square of distance
			Yr2 = Yr*Yr;
			Zr2 = Zr*Zr;

			RijSQ = Xr2 + Yr2 + Zr2;

			if(RijSQ<=RadiusSQ) {
				
				r2i = 1.0/RijSQ; //Intermediate Calcs for the LJ (6-12) potential
				r6i = r2i*r2i*r2i;

				//Fij = 48.0*r2i*r6i*(r6i-0.5); //LJ Force in reduced units

				force = 48.0*pow(1.0/RijSQ,6)-24*pow(1.0/RijSQ,3);
				Fij = force*r2i;

				PE=PE+r6i*(r6i-1.0);

				//Component Forces
				Fxij = Xr*Fij;
				Fyij = Yr*Fij;
				Fzij = Zr*Fij;

				//Sum forces on atom i
				fx[i] +=Fxij;
				fy[i] +=Fyij;
				fz[i] +=Fzij;
				
				//Sum forces on atom j
				fx[j] -= Fxij;
				fy[j] -= Fyij;
				fz[j] -= Fzij;
				
				p_PE -= Fxij*Xr-Fyij*Yr-Fzij*Zr;



			}

		}

	}

	PE*=4;


}



void Atom::thermo_update_berendsen() {

	double lamda;
	double d_natoms=static_cast<double>(natoms);
	sum_KE();
	lamda = sqrt((3.0*d_natoms*final_temp)/(2.0*KE)); //for TauT=tstep,
	
	for (int i=0; i<natoms; i++) {
		vx[i] *= lamda;
		vy[i] *= lamda;
		vz[i] *= lamda;
	}
	
}


void Atom::thermo_create_none() {

	thermo_type = NONE;
	thermo_created = 1;
	cout<<"created no thermo"<<endl;
}

void Atom::thermo_create_berendsen() {
	thermo_type = BERENDSEN;
	thermo_created = 1;
	cout<<"created Berendsen thermostat"<<endl;
}
void Atom::thermo_set_temps(double i_temp, double f_temp) {
	initial_temp = i_temp;
	final_temp = f_temp;
}

void Atom::thermo_update() {
	if(thermo_type==BERENDSEN) {
		thermo_update_berendsen();
	}
}


void Atom::sum_KE_and_P() {
	sum_KE();
	calc_pressure();
}

void Atom::sum_KE() {

	KE=0;

	for (int i=0; i<natoms; i++) {
		KE+= vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
	}

	KE=0.5*KE;
}

void Atom::calc_pressure() {

	p_KE = -KE*2.0;
	p_KE *= (-1.0/(3.0*box_x*box_y*box_z));
	p_PE *= (-1.0/(3.0*box_x*box_y*box_z));
	
	tot_press = p_KE + p_PE;
	
}

	
	
