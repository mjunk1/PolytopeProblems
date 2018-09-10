#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>
#include <stabiliser.h>

#include "GLPKConvexSeparation.h"

using namespace std;


int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

string infile;
string outfile;
bool Quiet,mode;


try {

	TCLAP::CmdLine cmd("Program for eliminating redundant points from the convex hull of a set of points.", ' ', "0.1");

	// arguments

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices ", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", false, "out.dat", "string");
	cmd.add(output_arg);

	TCLAP::SwitchArg Quiet_arg ("Q","Quiet","Suppress detailed output to files", cmd, false);
	TCLAP::SwitchArg mode_arg ("m","mode","Use naive method to eliminate", cmd, false);

	cmd.parse(argc, argv);

	infile = input_arg.getValue();
	outfile = output_arg.getValue();
	Quiet = Quiet_arg.getValue();
	mode = mode_arg.getValue();


} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}

cout << "#--------------------------------------------" << endl;
cout << "# Elimination of non-extremal points " << endl;
cout << "#--------------------------------------------" << endl;


// write program call to stdout
cout << endl;
cout << "# Program call" << endl;
cout << "#-------------" << endl;
cout << "   ";
for(unsigned i=0; i<argc; i++) {
	cout << argv[i] << " ";
}
cout << endl << endl;


// read 
vector<LabelledState> states;

if(get_states(states, infile+".mat", infile+".lab") == 0) {

	cout << "   Read " << states.size() << " points." << endl;
}
else {
	cout << "   Failed. Will now quit." << endl;
	return -1;
}

// ----------------------------------
// ----- elimination
// ----------------------------------


// ---- now perform elimination on the remaining points
unsigned dim = states.at(0).object.size();
sort(states.begin(), states.end());

GLPKConvexSeparation lp (dim);
lp.set_verbosity(1);


lp.add_vertices(states);
cout << "   Found " << lp.get_nvertices() << " extremal points." << endl;
lp.print_parameters();

if(mode == true) {
	lp.delete_redundant_points();
	cout << "   Found " << lp.get_nvertices() << " extremal points." << endl;
	lp.print_parameters();
}

states = lp.iget_labelled_vertices();

// add a "c" to the end
for(unsigned i=0; i<states.size(); i++) {
	states.at(i).label += " c";
}

write_states(states, outfile+".mat", outfile+".lab");

}