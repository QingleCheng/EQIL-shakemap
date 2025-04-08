#define BOOST_PYTHON_STATIC_LIB
#define BOOST_NUMPY_STATIC_LIB
#include <Python.h>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <gsl/gsl_fft_complex.h>

using namespace std;
using namespace boost::python;
namespace np = boost::python::numpy;

np::ndarray REQPYrotdnn(np::ndarray& s1, 
						np::ndarray& s2, 
						double fs, 
						np::ndarray& dso, 
						np::ndarray& To, 
						int nn, 
						double T1, 
						double T2, 
						double zi,
						int nit,
						double NS);

BOOST_PYTHON_MODULE(REQPY) {
	Py_Initialize();
	np::initialize();
	def("REQPYrotdnn", &REQPYrotdnn);
}

np::ndarray REQPYrotdnn(np::ndarray& nds1, 
						np::ndarray& nds2, 
						double fs, 
						np::ndarray& nddso, 
						np::ndarray& ndTo, 
						int nn, 
						double T1, 
						double T2, 
						double zi,
						int nit,
						double NS)
{
	double* s1 = reinterpret_cast<double*>(nds1.get_data());
	double* s2 = reinterpret_cast<double*>(nds2.get_data());
	double* dso = reinterpret_cast<double*>(nddso.get_data());
	double* To = reinterpret_cast<double*>(ndTo.get_data());
	int s1_length = nds1.get_shape()[0];
	int dso_length = nddso.get_shape()[0];
	Py_intptr_t shape[2] = {2, 10};
	np::ndarray result = np::zeros(2, shape, np::dtype::get_builtin<double>());
	auto p_result = reinterpret_cast<double *>(result.get_data());
	p_result[0] = (double)s1_length;
	p_result[19] = (double)dso_length;
	return result;
}