// main.cpp : Defines the entry point for the application.
//

#include <armadillo>
#include "EntSys.h"
using namespace std;
using namespace arma;

int main(void) {
	int NrThreads = 4;
	wall_clock timer;
	timer.tic();

/*
    UnitCircle sys;
    umat r_a = MulInds(6, sys.dim);
    uvec resol = { 1000,1000 };
    EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
    int N = 4000;			   // no. iterations
    double ta = 1.0, tb = 0.0; // t_j=ta/(j+tb)
    ofstream fout1("./results/Unitcirclemaxval.txt");
    ofstream fout2("./results/Unitcircle.txt");
    ofstream fout3("./results/UnitcircleResM.txt");

*/

    Orbit3D sys;
    umat r_a = MulInds(5, sys.dim);
    uvec resol = {  100, 100, 100 };
    EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
    int N = 2000;			   // no. iterations
    double ta = 3.0, tb = 0.0; // t_j=ta/(j+tb)
    ofstream fout1("./results/3Dmaxval.txt");
    ofstream fout2("./results/3DRes.txt");
    ofstream fout3("./results/3DResM.txt");

/*
 *
    SprottSystem sys;
    umat r_a = MulInds(3, sys.dim);
    uvec resol = { 50, 5, 10 };
    EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
    int N = 1;			   // no. iterations
    double ta = 3.0, tb = 0.0; // t_j=ta/(j+tb)
    ofstream fout1("./results/Sprottmaxval.txt");
    ofstream fout2("./results/SprottRes.txt");
    ofstream fout3("./results/SprottResM.txt");
*/
/*
    repressilator sys;
    umat r_a = MulInds(4, sys.dim);
    uvec resol = { 500, 50, 100 };
    EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
    int N = 500;			   // no. iterations
    double ta = 3.0, tb = 0.0; // t_j=ta/(j+tb)
    ofstream fout1("./results/repressilatormaxval.txt");
    ofstream fout2("./results/repressilatorRes.txt");
    ofstream fout3("./results/repressilatorResM.txt");
*/


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
			bestMaxVal= MaxVal;
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