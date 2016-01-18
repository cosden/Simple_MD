#include"md_sim.h"
#include"input_file.h"

int main(int argc, char**argv) {

	Md_sim md_sim;

	Input_file input_file;
	input_file.read_input_file(argc, argv, &md_sim);

	md_sim.initialize();
	md_sim.run();
	md_sim.finalize();

	

	return 1;
}
