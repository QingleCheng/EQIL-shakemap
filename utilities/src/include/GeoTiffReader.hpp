/**
 * @file GeoTiffReader.hpp
 * @author MPCB_Bishop
 * @version 1.0.1
 * @date 2021-05-02
 * @brief 读取 `GeoTiff` 文件内容
 * @details `GeoTiff` 文件读取API \n
 * 实现了一个简单的遍历器用于遍历网格数据
 */
#ifndef GEO_TIFF_READER_H
#define GEO_TIFF_READER_H

#include <iostream>
#include <gdal_priv.h>  // 读取 GeoTiff 文件
#include <streambuf>
#include <fstream>
#include <ctime>
#include <string>
#include <fstream>
#include <iomanip>      // 改变输出的精度
using namespace std;

/**
 * @brief 用于 `generator()` 函数的返回值打包, 
 * 格式为经度、维度、GeoTiff文件数据、函数返回数据
 * @tparam F `GeoTiff` 文件中存储数据的格式
 * @tparam T `generator()` 函数的返回值类型
 */
template <typename F, typename T>
struct packedRasData {
    double lon;         ///< 文件中存储的经度
    double lat;         ///< 文件中存储的纬度
    F rasData;          ///< `GeoTiff` 文件中存储数据
    T data;             ///< 函数返回值数据
};

/// 地面以下 30m 剪切波速处理类, 用于读取相关数据
class Vs30 {
private:
    GDALRasterBand* __band;
    int __nX;
    int __nY;
    int __offX;
    int __offY;
    float* __pafScanData;
    double* __transform;
public:
    /// 使用 Vs30 文件路径初始化类
    Vs30(const char* vs30_file_path):__pafScanData(nullptr), 
                                     __nX(0), __nY(0), 
                                     __offX(0), __offY(0) {
        __transform = new double[6];
        GDALAllRegister();
        GDALDataset* vs30Data = (GDALDataset* )GDALOpen(vs30_file_path, GA_ReadOnly);
        vs30Data->GetGeoTransform(__transform);
        __band = vs30Data->GetRasterBand(1);
    }
    ~Vs30() {
        delete[] __transform;
        delete[] __pafScanData;
    }

    /// 根据经纬度范围读取文件中数据并存入内存
    void get_data(const double lon_min, 
                  const double lon_max, 
                  const double lat_min, 
                  const double lat_max) {
        if(!__pafScanData) {
            delete[] __pafScanData;
        }
        __offX = floor((lon_min + 180) * 120);
        __nX = ceil((lon_max + 180) * 120) - __offX + 1;
        __offY = floor((84 - lat_max) * 120);
        __nY = ceil((84 - lat_min) * 120) - __offY + 1;
        __pafScanData = new float[__nX * __nY];
        if (CE_None != __band->RasterIO(
            GF_Read, __offX, __offY, __nX, __nY,
            __pafScanData, __nX, __nY,
            __band->GetRasterDataType(), 0, 0
        )) throw "Can not read data correct!";
    }

    /**
     * @brief 根据内存数据及经纬度坐标读取并插值得出此处的剪切波波速, 插值方法为 **双线性**
     * @note 请在运行本函数前先运行 `get_data()` 函数, 否则内存中无数据将抛出异常
     * @return 此处的地下 30m 剪切波速大小
     */
    double get_vs30(const double lon, const double lat) {
        int lon_index = floor((lon + 180) * 120);
        int lat_index = floor((84 - lat) * 120);
        double delta_lon = (lon + 180) * 120 - lon_index;
        double delta_lat = (84 - lat) * 120 - lat_index;
        lon_index -= __offX;
        lat_index -= __offY;
        double vs00 = __pafScanData[lat_index * __nX + lon_index];
        double vs01 = __pafScanData[lat_index * __nX + lon_index + 1];
        double vs10 = __pafScanData[(lat_index + 1) * __nX + lon_index];
        double vs11 = __pafScanData[(lat_index + 1) * __nX + lon_index + 1];
        double vs50 = (vs10 - vs00) * delta_lat + vs00;
        double vs51 = (vs11 - vs01) * delta_lat + vs01;
        return (vs51 - vs50) * delta_lon + vs50;
    }
};

/// `GeoTiff` 文件处理类, 用于读取此类文件
class GeoTiff {
private:
    GDALDataset* __geoTiffData;
    GDALRasterBand* __band;
    int __nX;
    int __nY;
    float* __pafScanData;
    double* __transform;
public:
    double lon_min; ///< `GeoTiff` 文件的范围
    double lon_max; ///< `GeoTiff` 文件的范围
    double lat_min; ///< `GeoTiff` 文件的范围
    double lat_max; ///< `GeoTiff` 文件的范围

    /// 使用 `GeoTiff` 文件路径初始化类
    GeoTiff(const char* tiff_file_path) {
        __transform = new double[6];
        GDALAllRegister();
        __geoTiffData = (GDALDataset* )GDALOpen(tiff_file_path, GA_ReadOnly);

        // 获取坐标变换
        __geoTiffData->GetGeoTransform(__transform);
        __nX = __geoTiffData->GetRasterXSize();
        __nY = __geoTiffData->GetRasterYSize();
        __band = __geoTiffData->GetRasterBand(1);
        __pafScanData = new float[__nX * __nY];
        if (CE_None != __band->RasterIO(
            GF_Read, 0, 0, __nX, __nY,
            __pafScanData, __nX, __nY,
            __band->GetRasterDataType(), 0, 0
        )) throw "Can not read data correct!";

        // 计算GeoTiff文件的坐标范围
        lon_min = __transform[0];
        lon_max = __transform[0] + (__nX - 1) * __transform[1];
        lat_min = __transform[3] + (__nY - 1) * __transform[5];
        lat_max = __transform[3];
    }
    ~GeoTiff() {
        delete[] __transform;
        delete[] __pafScanData;
    }

    /// 输出 ShakeMap 的点云文件, 共4列, 分别为经度、纬度、波速、编码（自动递增）
    void write_shakemap_points(const char* outPutFilePath, 
                               const char* vs30FilePath) {
        Vs30 vs30(vs30FilePath);
        vs30.get_data(lon_min, lon_max, lat_min, lat_max);
        ofstream outFile;
        outFile.open(outPutFilePath, ios::out | ios::trunc);
        double x, y;
        for(int j = 0; j < __nY; j++) {
            for(int i = 0; i < __nX; i++) {
                x = __transform[0] + i * __transform[1];
                y = __transform[3] + j * __transform[5];
                outFile << fixed << setprecision(8)  << x << " " << y << " " 
                        << vs30.get_vs30(x, y) << " " << j * __nX + i << "\n";
            }
        }
        outFile.close();
    }

    /**
     * @brief `GeoTiff` 文件迭代器
     * @note 如果对迭代顺序敏感, 请注意本函数使用的迭代顺序为从左到右, 从上到下
     * @tparam F 函数模版, 至少有三个参数, 分别代表经度、维度、坡度, 
     * 例如: `F(lon, lat, slope, ...)`
     * @tparam Args... 函数 `F` 的其他参数（除经纬坐标及坡度）
     * @param f 需要遍历网格数据的函数, 要求至少有三个参数, 
     * 分别代表经度、维度、坡度, 例如: `F(lon, lat, slope, ...)`
     * @param args 函数 `f` 所需要的其他参数
     * @return `packedRasData` 的向量, 每个可以被解压为 `lon`, `lat`, `rasData`, `data`, 
     * 并且 `data` 为函数 `f` 的返回值
     */
    template<class F, class... Args>
    auto generator(F&& f, Args&&... args) -> 
    vector<packedRasData<float, typename result_of<F(double, double, float, Args...)>::type> > {
        using return_type = typename result_of<F(double, double, float, Args...)>::type;

        vector<packedRasData<float, return_type> > ans;
        ans.reserve(__nX * __nY);
        double x, y;
        for(int j = 0; j < __nY; j++) {
            for(int i = 0; i < __nX; i++) {
                x = __transform[0] + i * __transform[1];
                y = __transform[3] + j * __transform[5];
                ans.emplace_back(move(packedRasData<float, return_type>{ 
                    x, y, __pafScanData[j * __nX + i], 
                    forward<F>(f)(x, y, __pafScanData[j * __nX + i], forward<Args>(args)...) 
                }));
            }
        }
        return ans;
    }

    /**
     * @brief `GeoTiff` 文件迭代器
     * @note 如果对迭代顺序敏感, 请注意本函数使用的迭代顺序为从左到右, 从上到下
     * @tparam F 函数模版, 至少有三个参数, 分别代表经度、维度、坡度, 
     * 例如: `F(lon, lat, slope, ...)`
     * @tparam Args... 函数 `F` 的其他参数（除经纬坐标及坡度）
     * @param step 经纬度的间隔（用于节省计算时间）
     * @param f 需要遍历网格数据的函数, 要求至少有三个参数, 
     * 分别代表经度、维度、坡度, 例如: `F(lon, lat, slope, ...)`
     * @param args 函数 `f` 所需要的其他参数
     * @return `packedRasData` 的向量, 每个可以被解压为`lon`, `lat`, `rasData`, `data`, 
     * 并且 `data` 为函数 `f` 的返回值
     */
    template<class F, class... Args>
    auto generator_step(int step, F&& f, Args&&... args) -> 
    vector<packedRasData<float, typename result_of<F(double, double, float, Args...)>::type> > {
        using return_type = typename result_of<F(double, double, float, Args...)>::type;

        vector<packedRasData<float, return_type> > ans;
        double x, y;
        for(int j = 0; j < __nY; j++) {
            for(int i = 0; i < __nX; i++) {
                if (i % step == 0 && j % step == 0) {
                    x = __transform[0] + i * __transform[1];
                    y = __transform[3] + j * __transform[5];
                    ans.emplace_back(move(packedRasData<float, return_type>{ 
                        x, y, __pafScanData[j * __nX + i],
                        forward<F>(f)(x, y, __pafScanData[j * __nX + i], forward<Args>(args)...) 
                    }));
                }
            }
        }
        return ans;
    }
};

#endif // GEO_TIFF_READER_H