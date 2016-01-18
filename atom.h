#ifndef ATOM_H
#define ATOM_H
#include <vector>

class Atom
{

 public:
	Atom();
	~Atom();

	int natoms;
	std::vector<double> vx, vy, vz;  //velocity vectors
	std::vector<double> rx, ry, rz;  //position vectors
	std::vector<int> type;           //atom type
	std::vector<double> fx, fy, fz;  //force vectors

	void set_num_atoms(int);
	void initialize_velocities();

	void calc_force();
	
	void set_box (double bx, double by, double bz) {
		box_x=bx; box_y=by; box_z=bz;
	}

	//note all units non-dimensionalized
	double T;           //temperature
	double PE, KE;      //potential and kinetic energy
	double p_PE, p_KE;  //potential (configurational) and kinetic components of pressure
	double tot_press;   //total pressure (scalar)

	//thermostat
	int thermo_step;    //how often to do thermostatting
	int thermo_type;    //type of thermostat (none, berendsen)
  
	double initial_temp;  //initial temperature
	double final_temp;    //desired final temperature
	int thermo_created;   //=1 if thermostat has been created

	void thermo_create_none();
	void thermo_create_berendsen();
	void thermo_set_temps(double, double);

	void thermo_update();

	void sum_KE_and_P();

 private:
	double box_x, box_y, box_z;    //x,y,z components of simulation box

	void thermo_update_berendsen();
	void calc_pressure();
	void sum_KE();

};

#endif

