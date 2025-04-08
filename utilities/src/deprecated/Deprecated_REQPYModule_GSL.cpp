#define _USE_MATH_DEFINES
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <gdal_priv.h>
#include <filesystem>
#include <vector>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex_math.h>
#include <fftw3.h>
#include <tuple>
#include <algorithm>

#include "convolution_fftw.h"
#include "ThreadPool.h"

using namespace std;

double median(vector<double>& inputV) {
    int size = inputV.size();
    return ((size % 2 == 0) ? ((inputV[size / 2 - 1] + inputV[size / 2]) / 2) : (inputV[size / 2]));
}

void cwtzm(
    vector<double>& inputS, 
    double inputFs,
    vector<double>& inputScales,
    double omega,
    double zeta,
    vector<vector<double> >& outputCoefs
) {
    int nf = inputScales.size();
    double dt = 1 / inputFs;
    int n = inputS.size();

    for(auto& c: outputCoefs) c.clear();
    outputCoefs.clear();
    outputCoefs.reserve(nf);

    vector<double> t, wv;
    double ttemp;
    t.reserve(n);
    wv.reserve(n);

    for(int i = 0; i < n; i++) t.emplace_back(i * dt);

    double centertime = median(t);

    for(int i = 0; i < n; i++) t[i] -= centertime;
	
    FFTW_Convolution::Workspace ws;
    FFTW_Convolution::init_workspace(ws, FFTW_Convolution::LINEAR_SAME, n, 1, n, 1);
    for(int i = 0; i < nf; i++) {
        wv.clear();
		ttemp = t[t.size() - 1] / inputScales[i];
		wv.emplace_back((exp(-zeta * omega * abs(ttemp))*sin(omega * ttemp)) / sqrt(inputScales[i]));
        for(int j = 0; j < t.size() - 1; j++) {
            ttemp = t[j] / inputScales[i];
            wv.emplace_back((exp(-zeta * omega * abs(ttemp))*sin(omega * ttemp)) / sqrt(inputScales[i]));
        }
        vector<double> coef;
        coef.reserve(n);
        FFTW_Convolution::convolve(ws, &inputS[0], &wv[0]);
        for(int j = 0; j < n; j++) {
            coef.emplace_back(ws.dst[j]);
        }
        outputCoefs.emplace_back(move(coef));
    }
    FFTW_Convolution::clear_workspace(ws);
}

void getdetails(
    vector<double>& inputT, 
    vector<double>& inputS, 
    vector<vector<double> >& inputC, 
    vector<double>& inputScales, 
    double omega, 
    double zeta,
    vector<vector<double> >& outputD,
    vector<double>& outputSr
) {
	int NS = inputScales.size();
	int n = inputS.size();

    for(auto& d: outputD) d.clear();
    outputD.clear();
    outputD.reserve(NS);

	double centertime = median(inputT);
	double ttemp;
	
	vector<double> wv;
	wv.reserve(n);
	
	FFTW_Convolution::Workspace ws;
    FFTW_Convolution::init_workspace(ws, FFTW_Convolution::LINEAR_SAME, n, 1, n, 1);

	for(int i = 0; i < NS; i++) {
		wv.clear();
		ttemp = (inputT[inputT.size() - 1] - centertime) / inputScales[i];
		wv.emplace_back((exp(-zeta * omega * abs(ttemp))*sin(omega * ttemp)));
        for(int j = 0; j < inputT.size() - 1; j++) {
            ttemp = (inputT[j] - centertime)/ inputScales[i];
            wv.emplace_back((exp(-zeta * omega * abs(ttemp))*sin(omega * ttemp)));
        }
		vector<double> Dtemp;
		Dtemp.reserve(n);
		FFTW_Convolution::convolve(ws, &inputC[i][0], &wv[0]);
		for(int j = 0; j < n; j++) {
            Dtemp.emplace_back(- ws.dst[j] / pow(inputScales[i], 5.0/2.0));
        }
        outputD.emplace_back(move(Dtemp));
	}

    outputSr.clear();
    outputSr.reserve(n);

	for (int i = 0; i < n; i++) {
		ttemp = 0;
		for (int j = 0; j < NS - 1; j++) {
			ttemp += (inputScales[j + 1] - inputScales[j]) * (outputD[j + 1][i] + outputD[j][i]) / 2;
		}
        outputSr.emplace_back(ttemp);
	}
    double s_max = abs(inputS[0]);
    double sr_max = abs(outputSr[0]);
    for (auto s: inputS) if (s_max < abs(s)) s_max = s;
    for (auto s: outputSr) if (sr_max < abs(s)) sr_max = s;
    auto ff = s_max /sr_max;
    for (auto& d: outputSr) d *= ff;
    for (auto& d: outputD) {
        for (auto& dd: d) {
            dd *= ff;
        }
    }
	FFTW_Convolution::clear_workspace(ws);
}

/*
    Only for sorted data
*/
void interp(
    vector<double>& inputX,
    vector<double>& inputXk,
    vector<double>& inputYk,
    double defaultValue,
    vector<double>& outputAns
) {
    int kIndex = 0;
    int kSize = inputXk.size();
    outputAns.clear();
    for (auto x: inputX) {
        if (x < inputXk[kIndex]) {
            outputAns.emplace_back(defaultValue);
            continue;
        }
        while (kIndex < kSize && inputXk[kIndex] < x) {
            kIndex++;
        }
        if (kIndex == kSize) {
            outputAns.emplace_back(defaultValue);
            continue;
        }
        outputAns.emplace_back(
            inputYk[kIndex - 1] + 
            (x - inputXk[kIndex - 1]) * 
            (inputYk[kIndex] - inputYk[kIndex - 1]) / 
            (inputXk[kIndex] - inputXk[kIndex - 1])
        );
    }
}

void ResponseSpectrumTheta(
    vector<double>& inputT,
    vector<double>& inputS1,
    vector<double>& inputS2,
    double z,
    double dt,
    vector<double>& inputTheta,
    vector<vector<double> >& outputPSA
) {
    for(auto& psa: outputPSA) psa.clear();
    outputPSA.clear();

    vector<double> theta;
    int ntheta = inputTheta.size();
    theta.reserve(ntheta);
    for(auto t: inputTheta) theta.emplace_back(t * M_PI / 180.0);

    int npo = max(inputS1.size(), inputS2.size());
    int nT = inputT.size();

    int nor = npo;

    size_t n = (1 << (int)ceil(log2(npo + 10.0 * inputT[nT - 1] / dt)));

    double fs = 1 / dt;
    double fres = fs / n;
    int nfrs = (int)ceil(n / 2);
    vector<double> freqs, ww;
    freqs.reserve(nfrs + 1);
    ww.reserve(nfrs + 1);
    for(int i = 0; i < nfrs + 1; i++) freqs.emplace_back(fres * i);
    for(auto f: freqs) ww.emplace_back(2 * M_PI * f);

    vector<gsl_complex> ffts1, ffts2;
    ffts1.reserve(n);
    ffts2.reserve(n);

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    memset(in, 0, n * sizeof(fftw_complex));
    for(int i = 0; i < inputS1.size(); i++) {
        in[i][0] = inputS1[i];
    }

    fftw_execute(p);

    for(int i = 0; i < n; i++) {
        ffts1.emplace_back(gsl_complex_rect(out[i][0], out[i][1]));
    }

    memset(in, 0, n * sizeof(fftw_complex));
    for(int i = 0; i < inputS2.size(); i++) {
        in[i][0] = inputS2[i];
    }

    fftw_execute(p);

    for(int i = 0; i < n; i++) {
        ffts2.emplace_back(gsl_complex_rect(out[i][0], out[i][1]));
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    
    outputPSA.reserve(ntheta);
    for(int i = 0; i < ntheta; i++) {
        vector<double> sd;
        sd.resize(nT);
        outputPSA.emplace_back(move(sd));
    }
            
    #pragma omp parallel for   
    for(int kk = 0; kk < nT; kk++) {
        vector<gsl_complex> H1;
        vector<gsl_complex> d1, d2;
        vector<double> gsl_in;
        H1.reserve(n);
        d1.reserve(nor);
        d2.reserve(nor);
        gsl_in.reserve(2 * n);
        int m = 1;
        gsl_complex temp;
        double w ,k ,c;
        
        w = 2 * M_PI / inputT[kk];
        k = m * w * w;
        c = 2 * z * m * w;
        for(int i = 0; i < nfrs + 1; i++) {
            H1.emplace_back(gsl_complex_inverse(gsl_complex_rect(-m * ww[i] * ww[i] + k, c * ww[i])));
        }
        for(int i = nfrs + 1; i < n; i++) {
            H1.emplace_back(gsl_complex_conjugate(H1[n - i - 1]));
        }
        H1[n / 2].dat[1] = 0;

        for(int i = 0; i < n; i++) {
            temp = gsl_complex_mul(H1[i], ffts1[i]);
            gsl_in.emplace_back(temp.dat[0]);
            gsl_in.emplace_back(temp.dat[1]);
        }
        gsl_fft_complex_radix2_inverse(&gsl_in[0], 1, n);

        for(int i = 0; i < nor; i++) d1.emplace_back(gsl_complex_rect(gsl_in[i * 2], gsl_in[i * 2 + 1]));

        gsl_in.clear();
        for(int i = 0; i < n; i++) {
            temp = gsl_complex_mul(H1[i], ffts2[i]);
            gsl_in.emplace_back(temp.dat[0]);
            gsl_in.emplace_back(temp.dat[1]);
        }
        gsl_fft_complex_radix2_inverse(&gsl_in[0], 1, n);

        for(int i = 0; i < nor; i++) d2.emplace_back(gsl_complex_rect(gsl_in[i * 2], gsl_in[i * 2 + 1]));

        for(int i = 0; i < ntheta; i++) {
            double max_drot = 0;
            double max_temp;
            for(int j = 0; j < nor; j++) {
                max_temp = gsl_complex_abs(
                    gsl_complex_add(
                        gsl_complex_mul_real(d1[j], cos(theta[i])), 
                        gsl_complex_mul_real(d2[j], sin(theta[i]))
                    )
                );
                if(max_drot < max_temp) max_drot = max_temp;
            }
            outputPSA[i][kk] = max_drot * (2 * M_PI / inputT[kk]) * (2 * M_PI / inputT[kk]);
        }
    }
}

void REQPYrotdnn(
    vector<double>& inputS1, 
    vector<double>& inputS2, 
    double inputFs, 
	vector<double>& inputDso, 
	vector<double>& inputTo, 
	int nn, 
	double T1, 
	double T2, 
	double zi = 0.05,
	int nit = 15,
	int NS = 100
) {
    size_t n = min(inputS1.size(), inputS2.size());
    double dt = 1 / inputFs;

	vector<double> theta;
    theta.reserve(180);
    for(int i = 0; i < 180; i++) theta.emplace_back(i);

	vector<double> t;
	t.reserve(n);
	for(int i = 0; i < n; i++)
		t.emplace_back(i * dt);
	
    double FF1 = 0.1, FF2 = 1.0 / (2.0 * dt);
    double omega = M_PI;
    double zeta = 0.05;
    if (T1 < (1/FF2)) T1 = 1/FF2;
    if (T2 > (1/FF1)) FF1 = 1 / T2;
    vector<double> freqs, T, scales;
    freqs.reserve(NS);
    T.reserve(NS);
    scales.reserve(NS);
    for (int i = 0; i < NS; i++) {
        freqs.emplace_back(FF2 * exp(i * log(FF1 / FF2) / (NS - 1)));
    }
    for (int i = 0; i < NS; i++) {
        T.emplace_back(1.0 / freqs[i]);
    }
    for (int i = 0; i < NS; i++) {
        scales.emplace_back(1.0 / (2.0 * freqs[i]));
    }
    vector<vector<double> > C1, C2;
    cwtzm(inputS1, inputFs, scales, omega, zeta, C1);
    cwtzm(inputS2, inputFs, scales, omega, zeta, C2);
	
	cout << "Wavelet decomposition performed" << endl;
	
	vector<vector<double> > D1, D2;
	vector<double> sr1, sr2;
	
	getdetails(t, inputS1, C1, scales, omega, zeta, D1, sr1);
	getdetails(t, inputS2, C2, scales, omega, zeta, D2, sr2);

	cout << "Detail functions generated" << endl;
	
    vector<double> ds;
    interp(T, inputTo, inputDso, -1, ds);

    vector<vector<double> > PSA180or;
    ResponseSpectrumTheta(T, inputS1, inputS2, zi, dt, theta, PSA180or);

	for(int loopcount = 0; loopcount < nit; loopcount++) {
		cout << "Now performing iteration " << loopcount + 1 << " of " << nit << endl;
        ResponseSpectrumTheta(T, inputS1, inputS2, zi, dt, theta, PSA180or);
	}
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    setvbuf(stdin, new char[1 << 20], _IOFBF, 1 << 20);
    setvbuf(stdout, new char[1 << 20], _IOFBF, 1 << 20);
    vector<double> s1;
    vector<double> s2;
    vector<double> dso;
    vector<double> To;
    s1.reserve(30000);
    s2.reserve(30000);
    dso.reserve(10000);
    To.reserve(10000);
    fstream inputfile;
    double temp;
    inputfile.open("../REQPY_test/s1.txt");
    while(inputfile >> temp) s1.emplace_back(temp);
    inputfile.close();
    inputfile.open("../REQPY_test/s2.txt");
    while(inputfile >> temp) s2.emplace_back(temp);
    inputfile.close();
    inputfile.open("../REQPY_test/Dso.txt");
    while(inputfile >> temp) dso.emplace_back(temp);
    inputfile.close();
    inputfile.open("../REQPY_test/To.txt");
    while(inputfile >> temp) To.emplace_back(temp);
    inputfile.close();
    REQPYrotdnn(s1, s2, 100, dso, To, 100, 0.05, 6.0, 0.05, 15, 100);
    return 0;
}