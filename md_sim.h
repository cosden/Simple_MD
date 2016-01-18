#ifndef MD_SIM_H
#define MD_SIM_H
#include <stdio.h>
#include <string.h>
class Md_sim 
{

 public:
	
	Md_sim();
	~Md_sim();

	class Atom *atom;
	class Verlet *verlet;
	
	void initialize();
	void run();
	void finalize();


	void set_num_t_steps(int n) {
		num_t_steps = n;
	}
	void set_num_unit_cells(int n) {
		num_unit_cells = n;
	}
	void set_rho(double r) {
		rho = r;
	}
	void set_dt(double t) {
		dt = t;
	}
	void set_write_press_freq(int wpf) {
		write_press_freq = wpf;
	}
	void set_write_freq(int wf) {
		write_freq = wf;
	}
	void set_initial_temp(double it) {
		initial_temp = it;
	}
	void set_final_temp(double ft) {
		final_temp = ft;
	}
	void set_runname(char *rn) {
		strcpy(runname,rn);
	}


 private:

	void form_lattice();
	void write_init_to_file();
	void vel_A();
	void vel_B();
	void setup_write();
	void write();

	//note all units non-dimensionalized

	int num_t_steps;        //number of timesteps to run
	int num_unit_cells;     //number of fcc unit cells in each direction (i.e. cubic)
	double rho;             //desired density
	double dt;              //timestep
	double time;            //current simulation time
	int write_press_freq;   //number of timesteps to write pressure info
	int write_freq;         //number of timesteps to write position/velocity info
	double initial_temp;    //initial temperature setpoint
	double final_temp;      //final thermostated temperature setpoint
	char *runname;          //runname to put in output files

	double box_x, box_y, box_z; //dimensions of simulation box
	double dt2;                 //dt/2.0
	double dtsq2;               //1/2 dt^2

	FILE *pos_and_v_file;           //position/velo file
	FILE *P_KE_Press_file;          //KE and pressure file
	char pos_and_v_filename[50];    //position/velo filename
	char P_KE_Press_filename[50];   //KE and pressure filename
	
};

#endif

