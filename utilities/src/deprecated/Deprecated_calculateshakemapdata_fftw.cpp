#define _USE_MATH_DEFINES
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <gdal_priv.h>
#include <filesystem>
#include <vector>
#include <fftw3.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "ThreadPool.h"

using namespace std;
struct SingleSpecData {
    string& stationName;
    int type;
    vector<double>& spec;
};

void elastic_spectrum_fft(double sourceDt, vector<double>& ag, int n, 
                          double zeta, double dT, double Tmax) 
{
    int i, agSize;
    agSize = ag.size();

    while (n < agSize) n <<= 1;

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    memset(in, 0, n * sizeof(fftw_complex));
    for(i = 0; i < agSize; i++) in[i][0] = -ag[i];

    fftw_execute(p);

    int length = round(Tmax / dT);
    auto Tn = new double[length];
    auto freq = new double[n];
    auto cf = new double[n];
    gsl_complex H, Utemp;
    double cfn, max_sd;

    for(i = 0; i < n; i++) {
        freq[i] = (i < (n / 2)) ? ((double)i / (sourceDt * n)) : ((double)(i - n) / (sourceDt * n));
    }

    for(i = 0; i < n; i++) {
        cf[i] = 2 * M_PI * freq[i];
    }

    fftw_complex* din = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex* dout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fftw_plan dp = fftw_plan_dft_1d(n, din, dout, FFTW_BACKWARD, FFTW_ESTIMATE);

    for(i = 0; i < length; i++) 
    {
        Tn[i] = dT * (i + 1);
        if (!(Tn[i] == 0.3 || Tn[i] == 1.0 || Tn[i] == 3.0)) continue;
        cfn = 2 * M_PI / Tn[i];

        for (int j = 0; j < n; j++) {
            H = gsl_complex_inverse(gsl_complex_rect(-cf[j] * cf[j] + cfn * cfn, 2 * zeta * cfn * cf[j]));
            Utemp = gsl_complex_mul(gsl_complex_rect(out[j][0], out[j][1]), H);
            din[j][0] = Utemp.dat[0];
            din[j][1] = Utemp.dat[1];
        }
        fftw_execute(dp);
        max_sd = abs(dout[0][0]);
        for(int j = 0; j < n; j++) {
            if(max_sd < abs(dout[j][0])) max_sd = abs(dout[j][0]);
        }
        cout << max_sd / n << endl;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_destroy_plan(dp);
    fftw_free(din);
    fftw_free(dout);
    delete[] Tn;
    delete[] freq;
    delete[] cf;
}

void processOneSingleFile(const filesystem::path& p) {
    ifstream afile;
    afile.open(p, ios::in);
    double t, a;
    vector<double> T;
    vector<double> A;
    T.reserve(sizeof(double) * 20000);
    A.reserve(sizeof(double) * 20000);
    afile >> t;
    while(afile >> t >> a) {
        T.emplace_back(t);
        A.emplace_back(a);
    }
    afile.close();
    
    elastic_spectrum_fft(T[1] - T[0], A, 1 << 15, 0.05, 0.01, 6);
}

void processPath(const filesystem::path& p){
    if(!exists(p)){         //目录不存在直接返回
        return;
    }
    ThreadPool pool(8);
    vector<future<void> > results;
    for(auto& directory: filesystem::directory_iterator(p)){
		if(filesystem::is_directory(directory)){
			string directory_string{directory.path().filename().u8string()};

			for(auto & path: filesystem::directory_iterator(directory.path())) {
				if(filesystem::is_regular_file(path)){
					string path_string{path.path().filename().u8string()};
					if(!path_string.compare("EW.txt")) {
                        processOneSingleFile(path.path());
					} else if(!path_string.compare("NS.txt")) {
						processOneSingleFile(path.path());
					} else if(!path_string.compare("UD.txt")) {
						processOneSingleFile(path.path());
					} else {
						continue;
					}
				}
			}
		}
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    setvbuf(stdin, new char[1 << 20], _IOFBF, 1 << 20);
    setvbuf(stdout, new char[1 << 20], _IOFBF, 1 << 20);
    filesystem::path p("../sample_data");
    processPath(p);
    return 0;
}