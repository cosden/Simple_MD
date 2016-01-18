#include "mainheader.h"
#include "atom.h"
#include "constants.h"
#include "md_sim.h"


Md_sim::Md_sim() {
	atom = new Atom; 
	runname = new char[20];
	
	//defaults
	num_t_steps = 0;
	num_unit_cells = 5;
	rho = 0.6;
	dt = 0.002;
	write_press_freq = 10;
	write_freq = 10;
	initial_temp = 1.0;
	final_temp = 1.0;
	strcpy(runname,"md_out");
}

Md_sim::~Md_sim() {
	delete [] runname;
	delete atom;
	

}
//intialize the simulation
//sets temperatures, creates atoms on a lattice, writes data to file
//sets up files for output
void Md_sim::initialize() {
	
	
	atom->thermo_set_temps(initial_temp,final_temp);
	form_lattice();
	atom->initialize_velocities();
	write_init_to_file();

	//store 1/2*dt and 1/2*dt*dt for later calculations
	dt2 = dt/2.0;
	dtsq2 = dt*dt2;

	setup_write();

}

//write initial positions and velocities to a file
void Md_sim::write_init_to_file(){

	FILE *pos_file;
	char pos_filename [50];
	strcpy(pos_filename, runname);
	strcat(pos_filename, "_initial_r_and_v.out");
	pos_file = fopen(pos_filename,"w");
	
	for (int i=0; i<atom->natoms; i++) {

		fprintf(pos_file, "%d %g  %g  %g  %g  %g %g\n", i, atom->rx[i], atom->ry[i], atom->rz[i], 
		        atom->vx[i], atom->vy[i], atom->vz[i]);
	}

	fclose(pos_file);
 
}

//main time loop
void Md_sim::run() {

	atom->calc_force(); //calculate for first timestep
	
	time = 0;

	cout<<endl<<"time  T*"<<endl;

	for(int t_step=0; t_step<num_t_steps; t_step++) {
		time+=dt;

		vel_A(); 
		atom->calc_force();
		vel_B();
		atom->thermo_update();
		atom->sum_KE_and_P();

		cout<<time<<" "<<((2.0/3.0)*(atom->KE)/atom->natoms)<<endl;
		
		write();

	}

}


void Md_sim::vel_A() {
	
	for(int i=0; i<atom->natoms; i++) {

		//advance positions a full time step
		atom->rx[i] += dt*atom->vx[i] + dtsq2*atom->fx[i];
		atom->ry[i] += dt*atom->vy[i] + dtsq2*atom->fy[i];
		atom->rz[i] += dt*atom->vz[i] + dtsq2*atom->fz[i];
		
		//advance velocities a half time step
		atom->vx[i] = atom->vx[i] + dt2*atom->fx[i];
		atom->vy[i] = atom->vy[i] + dt2*atom->fy[i];
		atom->vz[i] = atom->vz[i] + dt2*atom->fz[i];
		

	}
}

void Md_sim::vel_B() {

	for(int i=0; i<atom->natoms; i++) {
		
		//advance velocities a half time step
		atom->vx[i] = atom->vx[i] + dt2*atom->fx[i];
		atom->vy[i] = atom->vy[i] + dt2*atom->fy[i];
		atom->vz[i] = atom->vz[i] + dt2*atom->fz[i];

	}
	

}


void Md_sim::form_lattice() {
	
	int x_unit_cells, y_unit_cells, z_unit_cells;
	x_unit_cells=y_unit_cells=z_unit_cells=num_unit_cells;

	int num_atoms = num_unit_cells*num_unit_cells*num_unit_cells*4;
	atom->set_num_atoms(num_atoms);
	
	double lattice_const = pow((4.0/rho),(1.0/3.0));
	double a = lattice_const;
	cout<<"num_atoms="<<num_atoms<<" lattice_const="<<lattice_const<<endl;
	box_x = lattice_const*static_cast<double>(num_unit_cells);
	box_y = box_x;
	box_z = box_x;

	atom->set_box(box_x, box_y, box_z);

	//intialize first unit cell
	int i = 0;
	atom->rx[i]=0;
	atom->ry[i]=0;
	atom->rz[i]=0;

	i=1;
	atom->rx[i]=0;
	atom->ry[i]=a/2.0;
	atom->rz[i]=a/2.0;
	
	i=2;
	atom->rx[i]=a/2.0;
	atom->ry[i]=0;
	atom->rz[i]=a/2.0;

	
	i=3;
	atom->rx[i]=a/2.0;
	atom->ry[i]=a/2.0;
	atom->rz[i]=0;


	//Build up first "row" of unit cells along x axis

	for(i=4; i<4*x_unit_cells; i++) {
	 atom->rx[i]=atom->rx[i-4]+a;
	 atom->ry[i]=atom->ry[i-4];
	 atom->rz[i]=atom->rz[i-4];
	}

	int index=4*x_unit_cells;

	for (i=index; i<4*x_unit_cells*y_unit_cells; i++) {
		atom->rx[i]=atom->rx[i-index];
		atom->ry[i]=atom->ry[i-index]+a;
		atom->rz[i]=atom->rz[i-index];
	}

	//Expland x-y plane along Z axis
	index=4*x_unit_cells*y_unit_cells;

	for(i=index; i<4*x_unit_cells*y_unit_cells*z_unit_cells; i++) {
		atom->rx[i]=atom->rx[i-index];
		atom->ry[i]=atom->ry[i-index];
		atom->rz[i]=atom->rz[i-index]+a;
	}


   //Moves origin to center of box
	for (i=0; i<num_atoms; i++) {
		atom->rx[i]-=box_x/2.0;
		atom->ry[i]-=box_y/2.0;
		atom->rz[i]-=box_z/2.0;
	}
	


}


void Md_sim::setup_write() {
	
	strcpy(pos_and_v_filename, runname);
	strcat(pos_and_v_filename, "_pos_v_time.out");
	pos_and_v_file = fopen(pos_and_v_filename,"w");
	fprintf(pos_and_v_file, "index rx  ry  rz  vx  vy vz\n");
	fclose(pos_and_v_file);

	strcpy(P_KE_Press_filename, runname);
	strcat(P_KE_Press_filename, "_P_KE_time.out");
	P_KE_Press_file = fopen(P_KE_Press_filename,"w");
	fprintf(P_KE_Press_file,"Time KE tot_press\n");

	fclose(P_KE_Press_file);
	
 
}

void Md_sim::write() {

	pos_and_v_file = fopen(pos_and_v_filename, "a");

	fprintf(pos_and_v_file,"Time: %g\n", time);
	for (int i=0; i<atom->natoms; i++) {

		fprintf(pos_and_v_file, "%d %g  %g  %g  %g  %g %g\n", i, atom->rx[i], atom->ry[i], atom->rz[i], 
		        atom->vx[i], atom->vy[i], atom->vz[i]);
	}
	fclose(pos_and_v_file);

	P_KE_Press_file = fopen(P_KE_Press_filename, "a");
	fprintf(P_KE_Press_file,"%g %g %g\n", time, atom->KE, atom->tot_press);
	fclose(P_KE_Press_file);

}


void Md_sim::finalize() {


}
