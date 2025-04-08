/**
 * @file StationAcc.hpp
 * @author MPCB_Bishop
 * @version 1.0.1
 * @date 2021-07-08
 * @brief 地震动数据获取
 * @details 根据经纬度获取台站地震动数据 \n
 * 使用 `KD-Tree` 加速运算
 */
#ifndef STATION_SLOPE_H
#define STATION_SLOPE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include <filesystem>

#include "nanoflann.hpp" // kdtree

using namespace std;
using namespace nanoflann;
namespace fs = std::filesystem;

/**
 * @brief 坐标点结构体
 * @tparam T 坐标数据类型, 为浮点数类型, 例如: `float`, `double`
 */
template <typename T>
struct AccPoint
{
    vector<double> time;    ///< 台站时间序列
    vector<double> acc;     ///< 台站加速度时程
    T x; 			        ///< 点的坐标
    T y; 			        ///< 点的坐标
};

/**
 * @brief 点云结构体, 存储由 `Point` 构成的向量
 * @tparam T 坐标数据类型, 为浮点数类型, 例如: `float`, `double`
 */
template <typename T>
struct AccPointCloud
{
    /// 由 `Point` 构成的向量
    std::vector<AccPoint<T> >  pts;

    /// Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    /// 欧氏距离的平方, 不影响排序
    inline T kdtree_distance(const T* p1, const size_t idx_p2, size_t /*size*/) const
    {
        const T d0 = p1[0] - pts[idx_p2].x;
        const T d1 = p1[1] - pts[idx_p2].y;
        return d0 * d0 + d1 * d1;
    }

    /**
     * @brief Returns the dim'th component of the idx'th point in the class
     * @details Since this is inlined and the "dim" argument is typically 
     * an immediate value, the "if/else's" are actually solved at compile time.
     */
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0) return pts[idx].x;
        else return pts[idx].y;
    }

    /** 
     * @brief Optional bounding-box computation: return false to default to a standard bbox computation loop.
     * @details Return true if the BBOX was already computed by the class and returned 
     * in "bb" so it can be avoided to redo it again. Look at bb.size() to find out 
     * the expected dimensionality (e.g. 2 or 3 for point clouds)
     */
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

/// 地震动台站类, 内部实现了各类搜索与计算算法
class AccStation {
private:
    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, AccPointCloud<float> >,
        AccPointCloud<float>,
        2
    > my_kd_tree_t;
    // kd-tree
    my_kd_tree_t* __index;

    AccPointCloud<float> __cloud;

    // 仅搜索最近的一个台站
    const size_t __num_results = 1;
    size_t __ret_index;
    float __out_dist_sqr;

    KNNResultSet<float> __resultSet;

    SearchParams __searchParams;

public:

    /// 初始化 Station, `kd-tree`
    AccStation(): __index(nullptr), __resultSet(__num_results) { }

    /**
     * @brief 初始化 Station, 将文件中台站数据读取并保存至`kd-tree`中
     * @param acc_dir 地震动目录
     */
    AccStation(const char * acc_dir): __index(nullptr), 
                                      __resultSet(__num_results) {
        readLithologicData(acc_dir);
    }

    ~AccStation() { if(__index != nullptr) delete __index; }

    /**
     * @brief 将文件中台站数据读取并保存至 `kd-tree` 中
     * @note 数据文件格式一定要正确, 否则将无法读取正确数据
     * @param acc_dir 地震动目录
     */
    void readLithologicData(const char * acc_dir) {
        // 判断kd-tree中是否已有数据, 防止重复读取
        if(__index != nullptr) return;
        __cloud.pts.reserve(1844032);

        string filename; int delim_pos;
        float lon, lat;
        double time_temp, acc_temp;
        ifstream ipt;

        for(auto& directory_item: fs::directory_iterator(acc_dir)) {
            if(directory_item.is_directory()) {
                filename = directory_item.path().filename().u8string();
                delim_pos = filename.find('&');

                lon = stof(filename.substr(0, delim_pos));
                lat = stof(filename.substr(delim_pos + 1));

                vector<double> time, acc;

                ipt.open(directory_item.path() / "Acc.txt");
                while(ipt >> time_temp >> acc_temp) {
                    time.emplace_back(time_temp);
                    acc.emplace_back(acc_temp);
                }
                ipt.close();

                __cloud.pts.emplace_back(AccPoint<float>{
                    time, acc, lon, lat
                });
            }
        }

        // 建树
        __index = new my_kd_tree_t(2, __cloud, KDTreeSingleIndexAdaptorParams(10));
        __index->buildIndex();
    }

    /**
     * @brief 搜索最近的台站, 并将结果保存在 `ret_index` 中, 同时返回此处时间序列与加速度时程
     * @param lon 查询地点的经度
     * @param lat 查询地点的纬度
     * @return 返回查询地点时间序列及其加速度时程, 以 tuple 形式返回
     */
    auto searchNearestStation(float lon, float lat) 
    -> tuple<vector<double>, vector<double> > {
        float search_point[] = {lon, lat};

        __searchParams = SearchParams(10);
        __resultSet.init(&__ret_index, &__out_dist_sqr);
        __index->findNeighbors(__resultSet, search_point, __searchParams);

        return {__cloud.pts[__ret_index].time, __cloud.pts[__ret_index].acc};
    }

    /**
     * @brief 搜索最近的台站, 并将结果保存在 `ret_index` 中, 带 **调幅**
     * @param lon 查询地点的经度
     * @param lat 查询地点的纬度
     * @param pga 查询地点的PGA
     * @return 返回查询地点时间序列及其加速度时程, 以tuple形式返回
     * @details 本函数同时返回此处时间序列与 **调幅后** 的加速度时程
     */
    auto searchNearestStationAndAM(float lon, float lat, double pga) 
    -> tuple<vector<double>, vector<double> > {
        float search_point[] = {lon, lat};

        __searchParams = SearchParams(10);
        __resultSet.init(&__ret_index, &__out_dist_sqr);
        __index->findNeighbors(__resultSet, search_point, __searchParams);

        vector<double> accAM;
        accAM.reserve(__cloud.pts[__ret_index].acc.size());
        double s = *max_element(__cloud.pts[__ret_index].acc.begin(), __cloud.pts[__ret_index].acc.end());
        for(auto& p : __cloud.pts[__ret_index].acc) {
            accAM.emplace_back(p / s * pga);
        }

        return {__cloud.pts[__ret_index].time, accAM};
    }
};

#endif // STATION_SLOPE_H