#include "mainheader.h"
#include "input_file.h"
#include "md_sim.h"
#include "atom.h"


void Input_file::read_input_file(int argc, char **argv, Md_sim *md_sim) {


	//set max length of a line in the input file
	const int maxline=400;

	char line[maxline];
	char *command, *arg;
	char copy[maxline]="xx";
	char *nameoffile;

	if (argc==2){
		nameoffile=argv[1];
	}
	else if(argc<3){
		cout<<"ERROR: input file not specified in command line"<<endl;
		exit(1);
	}
	else {
		cout<<"ERROR: too many input commands"<<endl;
		exit(1);
	}


	FILE *infile;
	infile=fopen(nameoffile,"r");
	if (infile==NULL) {
		cout<<"ERROR: cannot open input file: "<<nameoffile<<endl;
		exit(1);
	}



	int z=1;

	//  finput input;

	int counter=0;

	while(z) {
		//get sequential line one by one
		fgets(line,maxline,infile);

		//make a copy to split apart
		strcpy(copy,line);

		//split string into tokens, command is first token
		command=strtok(copy," \t\n\r\f");

		if(command!=NULL) {

			if (!strcmp(command,"End"))
				{
					//cout<<"end of input file"<<endl<<endl;
					return;
				}


			if (!strcmp(command,"num_t_steps")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_num_t_steps(atoi(arg));

			}
			else if (!strcmp(command,"num_unit_cells")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_num_unit_cells(atoi(arg));
			}
			else if (!strcmp(command,"rho")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_rho(atof(arg));
			}
			else if (!strcmp(command,"dt")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_dt(atof(arg));
			}
			else if (!strcmp(command,"write_press_freq")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_write_press_freq(atoi(arg));
			}
 			else if (!strcmp(command,"write_freq")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_write_freq(atoi(arg));
			}
			else if (!strcmp(command,"thermostat")) {
				arg=strtok(NULL, " \t\n\r\f");
				if(!strcmp(arg,"berendsen")) {
					md_sim->atom->thermo_create_berendsen();
				}
				else if (!strcmp(arg,"none")) {
					md_sim->atom->thermo_create_none();
				}
				else {
					cout<<"Unknown Thermostat type: "<<arg<<"  using NO THERMOSTAT"<<endl;
					md_sim->atom->thermo_create_none();
				}
			}		
			else if (!strcmp(command,"initial_temp")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_initial_temp(atof(arg));
			}
			else if (!strcmp(command,"final_temp")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_final_temp(atof(arg));
			}
			else if (!strcmp(command,"runname")) {
				arg=strtok(NULL, " \t\n\r\f");
				md_sim->set_runname(arg);
			}



		}
		counter++;
	}



	return;
}
