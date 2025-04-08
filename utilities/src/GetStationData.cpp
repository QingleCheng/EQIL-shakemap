/**
 * @file GetStationData.cpp
 * @author MPCB_Bishop
 * @version 1.0.1
 * @date 2021-05-02
 * @brief StationData生成程序
 * @details Shell Usage:
 * $ ./GetStationData.exe </path/to/slope>
 *  </path/to/vs30>
 *  </path/to/output/shakemap/point>
 */

#include <iostream>
#include <fstream>

#include "GeoTiffReader.hpp"

using namespace std;

/**
 * @brief 生成shakemap计算点
 * @param argc 输入参数数量, 必须等于4, 即输入三个参数
 * @param argv 输入的参数值, 第一个参数为坡度数据路径, 
 * 第二个参数为Vs30数据文件路径, 
 * 第三个参数为输出的shakemap数据路径
 * @retval 0 正常
 * @retval -1 输入参数错误, 程序异常退出
 */
int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    setvbuf(stdin, new char[1 << 20], _IOFBF, 1 << 20);
    setvbuf(stdout, new char[1 << 20], _IOFBF, 1 << 20);

    if (argc != 4) exit(-1);

    GeoTiff slope(argv[1]); // 第一个参数为坡度数据路径
    slope.write_shakemap_points(argv[3], argv[2]); // 第三个参数为输出的shakemap数据路径
                                                   // 第二个参数为Vs30数据文件路径
    return 0;
}