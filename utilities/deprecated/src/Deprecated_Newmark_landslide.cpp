/**
 * @file Newmark_landslide.cpp
 * @author MPCB_Bishop
 * @version 1.0.1
 * @date 2021-05-02
 * @brief Newmark刚体滑块法计算程序，python接口
 */

#define BOOST_PYTHON_STATIC_LIB
#define BOOST_NUMPY_STATIC_LIB
#include <Python.h>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

using namespace std;
using namespace boost::python;
namespace np = boost::python::numpy;

/**
 * @brief Newmark刚体滑块法计算永久位移
 * @param time numpy数组
 * @param acc numpy数组
 * @param size 前两个参数的长度
 * @param ay 屈服加速度大小
 * @return 返回永久位移的大小，单位为m
 */
double NewmarkDisp(np::ndarray& time, np::ndarray& acc, long int size, double ay);

/// 导出至Python
BOOST_PYTHON_MODULE(Newmark_landslide) {
	Py_Initialize();
	np::initialize();
	def("NewmarkDisp", &NewmarkDisp);
}

double NewmarkDisp(np::ndarray& ndtime, np::ndarray& ndacc, long int size, double ay)
{
	double Displacement = 0.0;
	//if (ay<0) {
	//	Displacement = 10.0;
	//}
	//else {
	double* time = reinterpret_cast<double*>(ndtime.get_data());
	double* acc = reinterpret_cast<double*>(ndacc.get_data());
	double v0, v1, d0, d1;
	d0 = d1 = v1 = v0 = 0.0;
	double dt = time[1] - time[0];
	double* abs_acc;
	abs_acc = new double[size];
	double* rel_acc;
	rel_acc = new double[size];
	double* rel_vel;
	rel_vel = new double[size];
	double* rel_dis;
	rel_dis = new double[size];

	abs_acc[0] = 0;		//% total acceleration
	rel_acc[0] = 0;		//% relative acceleration
	rel_vel[0] = 0;		//% relative velocity
	rel_dis[0] = 0;		//% relative displacement

	//for (unsigned int i=0;i<time.size()-1;i++)
	//{
	//	if (acc[i]<ay && v0>0)
	//	{
	//		//Сٶ
	//		//v1 = v0 + acc[i] * dt;
	//		v1 = v0 + (acc[i] + acc[i + 1] - 2 * ay) / 2.0*dt;
	//		d1 = d0 + (v0 + v1) / 2.0*dt;
	//		v0 = v1;
	//		d0 = d1;
	//	} 
	//	else if (acc[i] < ay && v0 < 0)
	//	{
	//	}else
	//	{
	//		v1 = v0+(acc[i] + acc[i + 1]-2*ay) / 2.0*dt;				
	//		d1=d0+ (v0 + v1) / 2.0*dt;
	//		v0 = v1;
	//		d0 = d1;
	//	}
	//}

	for (unsigned int i = 1; i < size; i++)
	{
		abs_acc[i] = ay;
		rel_acc[i] = acc[i] - ay;
		rel_vel[i] = rel_vel[i - 1] + 0.5 * (rel_acc[i - 1] + rel_acc[i]) * dt;

		if (rel_vel[i] < 0)
		{
			abs_acc[i] = 0;
			rel_vel[i] = 0;
			rel_acc[i] = 0;
		}
		rel_dis[i] = rel_dis[i - 1] + rel_vel[i - 1] * dt + (2 * rel_acc[i - 1] + rel_acc[i]) * dt * dt / 6;
	}
	Displacement = rel_dis[size - 1];
	delete[]abs_acc;
	delete[]rel_acc;
	delete[]rel_vel;
	delete[]rel_dis;
	//}
	return Displacement;
}