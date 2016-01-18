#ifndef INPUT_FILE
#define INPUT_FILE
#include "md_sim.h"


class Input_file {

 public:
	char input_filename[400];


	void read_input_file(int, char **, Md_sim*);


};

#endif
