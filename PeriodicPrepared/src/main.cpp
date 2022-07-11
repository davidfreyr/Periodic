// main.cpp : Defines the entry point for the application.
//

#include "EntSys.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main(void) {
	int NrThreads = 4;
	wall_clock timer;
	timer.tic();

/*
	Lorenz sys;
	umat r_a = MulInds(2, sys.dim);
	uvec resol = { 500, 50, 100 };
	EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
	int N = 4;			   // no. iterations
	double ta = 2.0, tb = 0.0; // t_j=ta/(j+tb)
	ofstream fout1("./results/LorenzEnt.txt");
	ofstream fout2("./results/LorenzRes.txt");
*/

	//VanDerPol sys;
 //   umat r_a = MulInds(2, sys.dim);
 //   uvec resol = { 100,100 };
 //   EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
 //   int N = 10;	// no. iterations
 //   double ta = 1.0, tb = 0.0;  // t_j=ta/(j+tb)
 //   ofstream fout1("./results/VanDerPolMaxVal.txt");
 //   ofstream fout2("./results/VanDerPolRes.txt");
 //   ofstream fout3("./results/VanDerPolResM.txt");

/*
    SpeedControl sys;
    umat r_a = MulInds(0, sys.dim);
    uvec resol = { 1000,1000 };
    EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
    int N = 4000;	// no. iterations
    double ta = 1.0, tb = 0.0;  // t_j=ta/(j+tb)
    ofstream fout1("./results/SpeedControlMaxVal.txt");
    ofstream fout2("./results/SpeedControlRes.txt");
    ofstream fout3("./results/SpeedControlResM.txt");
*/
/*
    JetEngine sys;
    umat r_a = MulInds(4, sys.dim);
    uvec resol = { 1000,1000 };
    EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
    int N = 1000;	// no. iterations
    double ta = 1.0, tb = 0.0;  // t_j=ta/(j+tb)
    ofstream fout1("./results/JetEngineMaxVal.txt");
    ofstream fout2("./results/JetEngineRes.txt");
    ofstream fout3("./results/JetEngineResM.txt");
    */

    //three_d_example sys;
    //umat r_a = MulInds(4, sys.dim);
    //uvec resol = { 500, 50, 100 };
    //EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
    //int N = 4000;			   // no. iterations
    //double ta = 1.0, tb = 0.0; // t_j=ta/(j+tb)
    //ofstream fout1("./results/3Dmaxval.txt");
    //ofstream fout2("./results/3DRes.txt");
    //ofstream fout3("./results/3DResM.txt");


	//UnitCircle sys;
 //   umat r_a = MulInds(4, sys.dim);
 //   uvec resol = { 500, 500 };
 //   EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
 //   int N = 15;			   // no. iterations
 //   double ta = 1.0, tb = 0.0; // t_j=ta/(j+tb)
 //   ofstream fout1("./results/Unitcirclemaxval.txt");
 //   ofstream fout2("./results/Unitcircle.txt");
 //   ofstream fout3("./results/UnitcircleResM.txt");


	//Orbit3D sys;
	//umat r_a = MulInds(3, sys.dim);
	//uvec resol = { 500, 50, 100 };
	//EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
	//int N = 50;				// no. iterations
	//double ta = 1.0, tb = 0.0;	//t_j=ta/(j+tb)
	//ofstream fout1("./results/Orbit3DMaxVal.txt");
	//ofstream fout2("./results/Orbit3DRes.txt");
	//ofstream fout3("./results/Orbit3DResM.txt");

/*
	BouncingBall sys;
	umat r_a = MulInds(3, sys.dim);
	uvec resol = { 1000,1000 };
	EntSysDISC Esys(&sys, r_a, resol, true, NrThreads);
	int N = 200;	// no. iterations
	double ta = 1.0, tb = 0.0;  // t_j=ta/(j+tb)
	ofstream fout1("./results/BBEnt.txt");
	ofstream fout2("./results/BBRes.txt");
*/

/*
	Henon sys;
	umat r_a = MulInds(3, sys.dim);
	uvec resol = { 1000,1000 };
	EntSysDISC Esys(&sys, r_a, resol, true, NrThreads);
	int N = 1000;	// no. iterations
	double ta = 4.0, tb = 0.0;  // t_j=ta/(j+tb)
	ofstream fout1("./results/HenonEnt.txt");
	ofstream fout2("./results/HenonRes.txt");
*/

	Repressilator sys;
	umat r_a = MulInds(2, sys.dim);
	uvec resol = { 50, 50, 50 };
	EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
	int N = 10;			// no. iterations
	double ta = 1.0, tb = 0.0;	//t_j=ta/(j+tb)
	ofstream fout1("./results/RepressilatorMaxVal.txt");
	ofstream fout2("./results/RepressilatorRes.txt");
	ofstream fout3("./results/RepressilatorResM.txt");
	if (sys.s <= 2) {
		cout << "s is too small, it must be > 2";
		exit(0);
	}
	else if (sys.mu <= pow(2 / (sys.s - 2), 1 / sys.s) * (sys.s / (sys.s - 2))) {
		cout << "mu is too small, it must be > " << pow(2 / (sys.s - 2), 1 / sys.s) * (sys.s / (sys.s - 2));
		exit(0);
	}


	int k0, best_k;
	double MaxVal, EstEnt, bestMaxVal = Infinity;
	vec cx, s1, best_a;
	mat s2, best_p;
	double t;
	cout.precision(16);
	fout1.precision(16);
	fout2.precision(16);
	timer.tic();
	for (int j = 0; j < N; j++) {
		tie(MaxVal, cx) = Esys.FindMaximum();
		tie(s1, s2, k0, EstEnt) = Esys.riem_subg(cx);
		//cout << "s1 " << s1 << "s2 " << s2 << endl;
		cout << "itr " << j << ": MaxVal " << MaxVal << " punktur " << cx.t() << endl;
		fout1 << MaxVal << endl;
		if (MaxVal < bestMaxVal) {
			bestMaxVal = MaxVal;
			best_p = Esys.p;
			best_a = Esys.a;
			best_k = j;
		}
		t = ta / (j + 1 + tb);
		Esys.Norms1s2(s1, s2);
		Esys.StepForward(t, s1, s2); // step j+1
        //cout << "a " << Esys.a << endl;
        //cout << "P " << Esys.p << endl;
	}
	cout << N << " iterations computed in " << timer.toc() << " s" << endl;
	fout2 << N << " iterations computed in " << timer.toc() << " sec. : t_j = " << ta << "/(j+" << tb << ")" << endl;
	fout2 << "Best estimate of Restoration Entropy " << endl
    	<< "obtained in iteration " << best_k << " with " << endl
		<< "a : (pow of x_1, pow of x_2, etc. )  : coefficient " << endl;
	for (int i = 0; i < Esys.PolyDim; i++) {
		fout2 << "( ";
		for (int j = 0; j < Esys.dim; j++) {
			fout2 << Esys.r_a(i, j) << " ";
		}
		fout2 << ") : " << best_a(i) << endl;
	}
	fout2 << endl
		<< "p :" << endl;
	best_p.raw_print(fout2);

    for (int i = 0; i < Esys.PolyDim; i++) {
        for (int j = 0; j < Esys.dim; j++) {
            fout3 << Esys.r_a(i, j) << " ";
        }
        fout3 << " " << best_a(i) << endl;
    }
    fout3 << endl
          << "p :" << endl;
    best_p.raw_print(fout3);
}