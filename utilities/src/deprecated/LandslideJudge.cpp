/**
 * @file LandslideJudge.cpp
 * @author MPCB_Bishop
 * @version 1.0.1
 * @date 2021-05-02
 * @brief 滑坡判断程序
 * @details Shell Usage:
 * $ ./NewmarkDisp.exe </path/to/acc/dir>
 *  </path/to/slope>
 *  </path/to/lithologic>
 */

#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip>

#include "GeoTiffReader.hpp"
#include "StationLithologic.hpp"
#include "StationSlope.hpp"
using namespace std;
namespace fs = std::filesystem;

/**
 * @brief 判断滑坡是否发生，并生成计算结果，存储为`result.txt`
 * @param argc 输入参数数量，必须等于4，即输入三个参数
 * @param argv 输入的参数值，第一个参数为NewmarkDisp程序运行结果目录，
 * 第二个参数为输入的坡度数据路径，
 * 第三个参数为输入的岩性数据路径
 * @retval 0 正常
 * @retval -1 输入参数错误，程序异常退出
 */
int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    setvbuf(stdin, new char[1 << 20], _IOFBF, 1 << 20);
    setvbuf(stdout, new char[1 << 20], _IOFBF, 1 << 20);

    if(argc != 4) exit(-1);
    fs::path working_dir(argv[1]); // 第一个参数为NewmarkDisp程序运行结果目录

    GeoTiff slope(argv[2]); // 第二个参数为输入的坡度数据路径
    LithoStation sl(argv[3]); // 第三个参数为输入的岩性数据路径
    SlopeStation ss0((working_dir / "result0.txt").c_str());
    SlopeStation ss5((working_dir / "result5.txt").c_str());
    SlopeStation ss9((working_dir / "result9.txt").c_str());

    auto actual_slope = slope.get_tiff_data();
    ofstream opt;

    // 判断滑坡体饱和度为0%时的滑坡情况
    auto critical_slope0 = slope.generator(
        [&](double x, double y) -> float {
            return ss0.searchNearestStation(
                x, y, sl.searchNearestStation(
                    x, y
                )
            );
        }
    );
    opt.open(working_dir / "judgeresult0.txt");
    for(int i = 0; i < actual_slope.size(); i++) {
        opt << fixed << setprecision(8) << critical_slope0[i].lon << " "
            << critical_slope0[i].lat << " " << critical_slope0[i].data << " "
            << actual_slope[i] << " " << (actual_slope[i] >= critical_slope0[i].data) << endl;
    }
    opt.close();

    // 判断滑坡体饱和度为50%时的滑坡情况
    auto critical_slope5 = slope.generator(
        [&](double x, double y) -> float {
            return ss5.searchNearestStation(
                x, y, sl.searchNearestStation(
                    x, y
                )
            );
        }
    );
    opt.open(working_dir / "judgeresult5.txt");
    for(int i = 0; i < actual_slope.size(); i++) {
        opt << fixed << setprecision(8) << critical_slope5[i].lon << " "
            << critical_slope5[i].lat << " " << critical_slope5[i].data << " "
            << actual_slope[i] << " " << (actual_slope[i] >= critical_slope5[i].data) << endl;
    }
    opt.close();

    // 判断滑坡体饱和度为90%时的滑坡情况
    auto critical_slope9 = slope.generator(
        [&](double x, double y) -> float {
            return ss9.searchNearestStation(
                x, y, sl.searchNearestStation(
                    x, y
                )
            );
        }
    );
    opt.open(working_dir / "judgeresult9.txt");
    for(int i = 0; i < actual_slope.size(); i++) {
        opt << fixed << setprecision(8) << critical_slope9[i].lon << " "
            << critical_slope9[i].lat << " " << critical_slope9[i].data << " "
            << actual_slope[i] << " " << (actual_slope[i] >= critical_slope9[i].data) << endl;
    }
    opt.close();

    return 0;
}