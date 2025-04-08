/**
 * @file NewmarkDisp.hpp
 * @author MPCB_Bishop
 * @version 1.0.1
 * @date 2021-05-02
 * @brief Newmark刚体滑块法计算
 */
#ifndef NEWMARK_DISP_H
#define NEWMARK_DISP_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

/// Newmark刚体滑块法相关功能函数
namespace NewmarkDisp {
    /**
     * @brief Newmark 刚体滑块法计算永久位移
     * @param time 时间序列
     * @param acc 加速度序列
     * @param ay 屈服加速度的值
     * @return 永久位移量
     */
    double NewmarkDisp(vector<double>& time, vector<double>& acc, double ay) {
        // 获取时间序列长度与时间间隔
        size_t ntime = time.size();
        double dt = time[1] - time[0];

        double* abs_acc;
        abs_acc = new double[ntime];
        double* rel_acc;
        rel_acc = new double[ntime];
        double* rel_vel;
        rel_vel = new double[ntime];
        double* rel_dis;
        rel_dis = new double[ntime];

        abs_acc[0] = 0;		// total acceleration
        rel_acc[0] = 0;		// relative acceleration
        rel_vel[0] = 0;		// relative velocity
        rel_dis[0] = 0;		// relative displacement

        for (size_t i = 1; i < ntime; i++)
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
            rel_dis[i] = rel_dis[i - 1] + rel_vel[i - 1] * dt + 
                         (2 * rel_acc[i - 1] + rel_acc[i]) * dt * dt / 6;
        }

        double displacement = rel_dis[ntime - 1];
        delete[]abs_acc;
        delete[]rel_acc;
        delete[]rel_vel;
        delete[]rel_dis;

        return displacement;
    }

    /**
     * @brief 计算滑坡可能性
     * @note 位移不大于 0.05m 的认为滑坡可能非常低 \n
     * 位移大于 0.05m 且不大于 0.15m 的认为滑坡可能性低 \n
     * 位移大于 0.15m 且不大于 0.3m 的认为滑坡可能性中等 \n
     * 位移大于 0.3m 认为滑坡可能性高
     * @param disp 蒙特卡洛后得到的位移量, 即不同的 C, Phi 值对应的位移量
     * @param Ratio 输出滑坡可能性大小, 分别为非常低、低、中等、高的占总体的比例
     * @return 无返回值, 使用Ratio引用返回
     */
    void GetLandslideRatio(vector<double>& disp, vector<double>& Ratio) {
        Ratio.clear();
        int countVeryLow, countLow, countModerate, countHigh;
        countVeryLow = countLow = countModerate = countHigh = 0;
        const double veryLow = 0.05;
        const double low = 0.15;
        const double moderate = 0.3;
        for (auto& d:disp)
        {
            if (d <= veryLow)
            {
                countVeryLow++;
            }
            else if (d > veryLow && d <= low)
            {
                countLow++;
            }
            else if (d > low && d <= moderate)
            {
                countModerate++;
            }
            else {
                countHigh++;
            }
        }
        double size_disp = disp.size();
        Ratio.push_back(countVeryLow / size_disp);
        Ratio.push_back(countLow / size_disp);
        Ratio.push_back(countModerate / size_disp);
        Ratio.push_back(countHigh / size_disp);
    }

    /// 获取向量的中位数
    double getMedian(vector<double> x) {
        int xsize = x.size();
        sort(x.begin(), x.end());
        return (xsize % 2 == 0) ? ((x[xsize / 2 - 1] + x[xsize / 2]) / 2) : (x[xsize / 2]);
    }
}

#endif // NEWMARK_DISP_H