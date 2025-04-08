/**
 * @file NewmarkDisp.cpp
 * @author MPCB_Bishop
 * @version 1.0.3
 * @date 2021-07-15
 * @brief Newmark刚体滑块法计算程序
 * @details Shell Usage:
 * $ ./NewmarkDisp.exe <num_threads>
 *  </path/to/acc/dir>
 *  </path/to/slope>
 *  </path/to/lithologic>
 *  <calculation_step>
 *  <is_pga_am>
 */

/**
 * @mainpage Landslide 滑坡计算程序
 * 
 * Landslide滑坡计算程序, 本文档仅用于API展示, 
 * 具体使用详见项目
 *  <a href="https://gitee.com/guyi2000/land-slide">README</a>.
 */

#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip>
#include <omp.h>

#include "GeoTiffReader.hpp"
#include "NewmarkDisp.hpp"
#include "StationLithologic.hpp"
#include "StationAcc.hpp"
using namespace std;
namespace fs = std::filesystem;

/// PGA 文件读取类
class PGAReader {
private:
    int __nlons;
    int __nlats;
    int __step;
    int __x;
    int __y;
    vector<double> __pga;
public:
    /// 使用 PGA 文件与间隔大小初始化类
    PGAReader(const char * pgaFilePath, int step): __x(0), __y(0) {
        ifstream ipt;
        ipt.open(pgaFilePath);

        string ig;
        ipt >> ig >> __nlons >>__nlats;

        double ignore, pga_temp;
        while(ipt >> ignore >> ignore >> pga_temp) {
            __pga.emplace_back(pga_temp);
        }

        __step = step;
    }

    /**
     * @brief 读取一个值, 并将指针移动到下个值
     * @note 本函数不是内存安全的, 请确保可以调用
     */
    double iterator() {
        double pga = __pga[__y * __nlons + __x];
        __x += __step;
        if (__x >= __nlons) {
            __y += __step;
            __x = 0;
        }
        return pga;
    }

    /// 指针归零, 可以重新遍历
    void clear() {
        __x = 0; __y = 0;
    }
};

/**
 * @brief 获取滑坡概率,  **并行** 计算
 * @note 蒙特卡洛方法模拟, 默认模拟100次, 
 * 取中位数为滑坡永久移动量, 以此计算概率
 * @param lon 计算点的经度
 * @param lat 计算点的纬度
 * @param slope 此处的坡度
 * @param m 滑坡体饱和度, 0为完全不饱和, 1为完全饱和
 * @param sl 岩性台站数据
 * @param sa 加速度台站数据
 * @return 滑坡发生的概率, 取值范围为 [0, 1]
 * @retval 1 一定会发生滑坡
 * @retval 0 此地不会发生滑坡
 * @retval -nan/nan 计算中屈服加速度未定义
 * @todo 可以增加预判断, 如果明显不会滑坡则不计算
 */
double getDispProb(double lon,
                   double lat,
                   double slope,
                   double m,
                   LithoStation& sl, 
                   AccStation& sa) {
    /// 蒙特卡洛模拟次数, 默认为100
    int n = 100;
    vector<double> C_seq, phi_seq;
    C_seq.reserve(n); phi_seq.reserve(n);

    sl.searchNearestStation(lon, lat); // 获取最近台站的岩性数据
    sl.getSeqC_Phi(n, C_seq, phi_seq);

    auto [time, acc] = sa.searchNearestStation(lon, lat); // 获取最近台站的地震动数据

    vector<double> displacement;
    displacement.resize(n); 
#pragma omp parallel for
    for (int j = 0; j < n; j++) {
        displacement[j] = NewmarkDisp::NewmarkDisp(
            time, acc, sl.getAcBy(
                C_seq[j], phi_seq[j], 20, slope, m, 2.4
            )
        );
    }
    return 0.335 * (1 - exp(-0.048 * pow(NewmarkDisp::getMedian(displacement) * 100, 1.565)));
}

/**
 * @brief 获取滑坡概率,  **并行** 计算, 进行加速度调幅
 * @note 蒙特卡洛方法模拟, 默认模拟100次, 
 * 取中位数为滑坡永久移动量, 以此计算概率
 * @param lon 计算点的经度
 * @param lat 计算点的纬度
 * @param slope 此处的坡度
 * @param m 滑坡体饱和度, 0为完全不饱和, 1为完全饱和
 * @param sl 岩性台站数据
 * @param sa 加速度台站数据
 * @param pgar 各点PGA迭代器
 * @return 滑坡发生的概率, 取值范围为 [0, 1]
 * @retval 1 一定会发生滑坡
 * @retval 0 此地不会发生滑坡
 * @retval -nan/nan 计算中屈服加速度未定义
 * @todo 可以增加预判断, 如果明显不会滑坡则不计算
 */
double getDispProbAndAM(double lon,
                        double lat,
                        double slope,
                        double m,
                        LithoStation& sl, 
                        AccStation& sa,
                        PGAReader& pgar) {
    /// 蒙特卡洛模拟次数, 默认为100
    int n = 100;
    vector<double> C_seq, phi_seq;
    C_seq.reserve(n); phi_seq.reserve(n);

    sl.searchNearestStation(lon, lat); // 获取最近台站的岩性数据
    sl.getSeqC_Phi(n, C_seq, phi_seq);
    auto [time, acc] = sa.searchNearestStationAndAM(lon, lat, pgar.iterator()); // 获取最近台站的地震动数据

    vector<double> displacement;
    displacement.resize(n); 
#pragma omp parallel for
    for (int j = 0; j < n; j++) {
        displacement[j] = NewmarkDisp::NewmarkDisp(
            time, acc, sl.getAcBy(
                C_seq[j], phi_seq[j], 20, slope, m, 2.4
            )
        );
    }
    return 0.335 * (1 - exp(-0.048 * pow(NewmarkDisp::getMedian(displacement) * 100, 1.565)));
}

/**
 * @brief 计算滑坡概率, 并生成计算结果, 存储为`result0.txt`、
 * `result5.txt`、`result9.txt`, 分别代表滑坡体饱和比例为0\%、50\%、90\%
 * @param argc 输入参数数量, 必须等于7, 即输入六个参数
 * @param argv 输入的参数值, 第一个参数为计算并行核心数, 
 * 第二个参数为时程数据目录, 第三个参数为坡度数据, 
 * 第四个参数为岩性数据, 第五个参数为台站分辨率, 
 * 第六个参数为是否进行加速度调幅
 * @retval 0 正常
 * @retval -1 输入参数错误, 程序异常退出
 */
int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    setvbuf(stdin, new char[1 << 20], _IOFBF, 1 << 20);
    setvbuf(stdout, new char[1 << 20], _IOFBF, 1 << 20);

    if(argc != 7) exit(-1);

    omp_set_num_threads(atoi(argv[1])); // 第一个参数为计算并行核心数

    fs::path working_dir(argv[2]); // 第二个参数为时程数据目录

    GeoTiff slope(argv[3]); // 第三个参数为坡度数据
    
    LithoStation sl(argv[4]); // 第四个参数为岩性数据

    AccStation sa(working_dir.c_str()); 

    int station_step = atoi(argv[5]); // 第五个参数为台站分辨率

    if(!atoi(argv[6])) { // 第六个参数为是否进行加速度调幅
        ofstream opt;
        opt.open(working_dir / "result0.txt");
        slope.generator_step(station_step,
            [&](double x, double y, float s) -> int{
                opt << fixed << setprecision(8) << x << " " << y << " " 
                    << getDispProb(x, y, (double)s, 0, sl, sa) << "\n";
                return 0;
            }
        );
        opt.close();

        opt.open(working_dir / "result5.txt");
        slope.generator_step(station_step,
            [&](double x, double y, float s) -> int{
                opt << fixed << setprecision(8) << x << " " << y << " " 
                    << getDispProb(x, y, (double)s, 0.5, sl, sa) << "\n";
                return 0;
            }
        );
        opt.close();

        opt.open(working_dir / "result9.txt");
        slope.generator_step(station_step,
            [&](double x, double y, float s) -> int{
                opt << fixed << setprecision(8) << x << " " << y << " " 
                    << getDispProb(x, y, (double)s, 0.9, sl, sa) << "\n";
                return 0;
            }
        );
        opt.close();
    } else {
        PGAReader pgar((working_dir / "pga.txt").c_str(), station_step);
        ofstream opt;
        opt.open(working_dir / "result0.txt");
        slope.generator_step(station_step,
            [&](double x, double y, float s) -> int{
                opt << fixed << setprecision(8) << x << " " << y << " " 
                    << getDispProbAndAM(x, y, (double)s, 0, sl, sa, pgar) << "\n";
                return 0;
            }
        );
        opt.close();
        pgar.clear();

        opt.open(working_dir / "result5.txt");
        slope.generator_step(station_step,
            [&](double x, double y, float s) -> int{
                opt << fixed << setprecision(8) << x << " " << y << " " 
                    << getDispProbAndAM(x, y, (double)s, 0.5, sl, sa, pgar) << "\n";
                return 0;
            }
        );
        opt.close();
        pgar.clear();

        opt.open(working_dir / "result9.txt");
        slope.generator_step(station_step,
            [&](double x, double y, float s) -> int{
                opt << fixed << setprecision(8) << x << " " << y << " " 
                    << getDispProbAndAM(x, y, (double)s, 0.9, sl, sa, pgar) << "\n";
                return 0;
            }
        );
        opt.close();
    }

    return 0;
}