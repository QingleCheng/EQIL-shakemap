#define _USE_MATH_DEFINES
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <gdal_priv.h>
#include <filesystem>
#include <vector>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex_math.h>

#include "ThreadPool.h"

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

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
    auto AArray = new double[2 * n];
    memset(AArray, 0, 2 * n * sizeof(double));

    for(i = 0; i < agSize; i++) REAL(AArray, i) = -ag[i];

    gsl_fft_complex_radix2_forward(AArray, 1, n);

    int length = round(Tmax / dT);
    auto Tn = new double[length];
    auto freq = new double[n];
    auto cf = new double[n];
    gsl_complex H, Utemp;
    auto U = new double[2 * n];
    double cfn, max_sd;

    for(i = 0; i < n; i++) {
        freq[i] = (i < (n / 2)) ? ((double)i / (sourceDt * n)) : ((double)(i - n) / (sourceDt * n));
    }

    for(i = 0; i < n; i++) {
        cf[i] = 2 * M_PI * freq[i];
    }

    for(i = 0; i < length; i++) 
    {
        Tn[i] = dT * (i + 1);
        if (!(Tn[i] == 0.3 || Tn[i] == 1.0 || Tn[i] == 3.0)) continue;
        cfn = 2 * M_PI / Tn[i];

        for (int j = 0; j < n; j++) {
            H = gsl_complex_inverse(gsl_complex_rect(-cf[j] * cf[j] + cfn * cfn, 2 * zeta * cfn * cf[j]));
            Utemp = gsl_complex_mul(gsl_complex_rect(REAL(AArray, j), IMAG(AArray, j)), H);
            REAL(U, j) = Utemp.dat[0];
            IMAG(U, j) = Utemp.dat[1];
        }
        gsl_fft_complex_radix2_inverse(U, 1, n);
        max_sd = REAL(U, 0);
        for(int j = 0; j < n; j++) {
            if(max_sd < abs(REAL(U, j))) max_sd = abs(REAL(U, j));
        }
        cout << max_sd << endl;
    }

    delete[] AArray;
    delete[] Tn;
    delete[] freq;
    delete[] cf;
    delete[] U;
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
                        results.emplace_back(pool.enqueue(processOneSingleFile, path.path()));
					} else if(!path_string.compare("NS.txt")) {
						results.emplace_back(pool.enqueue(processOneSingleFile, path.path()));
					} else if(!path_string.compare("UD.txt")) {
						results.emplace_back(pool.enqueue(processOneSingleFile, path.path()));
					} else {
						continue;
					}
				}
			}
		}
    }
    for (auto&& result : results) {
        result.get();
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