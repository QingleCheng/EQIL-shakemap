/**
 * @file StationLithologic.hpp
 * @author MPCB_Bishop
 * @version 1.0.1
 * @date 2021-05-02
 * @brief 地质数据获取
 * @details 根据经纬度获取地质数据 \n
 * 使用 `KD-Tree` 加速运算
 */
#ifndef STATION_LITHOLOGIC_H
#define STATION_LITHOLOGIC_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>

#include "nanoflann.hpp" // kdtree

using namespace std;
using namespace nanoflann;

/**
 * @brief 坐标点结构体
 * @tparam T 坐标数据类型, 为浮点数类型, 例如: `float`, `double`
 */
template <typename T>
struct LithoPoint
{
    int LithoGroup;	///< 场地分类
    T x; 			///< 点的坐标
    T y; 			///< 点的坐标
};

/**
 * @brief 点云结构体, 存储由 `Point` 构成的向量
 * @tparam T 坐标数据类型, 为浮点数类型, 例如: `float`, `double`
 */
template <typename T>
struct LithoPointCloud
{
    /// 由 `Point` 构成的向量
    std::vector<LithoPoint<T> >  pts;

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

/// 岩土台站类, 内部实现了各类搜索与计算算法
class LithoStation {
private:
    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, LithoPointCloud<float> >,
        LithoPointCloud<float>,
        2
    > my_kd_tree_t;
    // kd-tree
    my_kd_tree_t* __index;

    LithoPointCloud<float> __cloud;

    // 仅搜索最近的一个台站
    const size_t __num_results = 1;
    size_t __ret_index;
    float __out_dist_sqr;

    KNNResultSet<float> __resultSet;

    SearchParams __searchParams;

    size_t __seed;

    default_random_engine __gen;

public:
    double C;             ///< 土的粘聚力
    double C_std;         ///< 土的粘聚力标准差
    double phi;           ///< 土的内摩擦角
    double phi_std;       ///< 土的内摩擦角标准差
    double gamma;         ///< 土的平均重度, 单位kN/m<sup>3</sup>

    /// 初始化Station, `kd-tree`
    LithoStation(): __index(nullptr), __resultSet(__num_results), 
                    __seed(time(nullptr)), __gen(__seed) { }

    /**
     * @brief 初始化 Station, 将文件中台站数据读取并保存至 `kd-tree` 中
     * @param lithologic_filepath 场地文件, 数据文件为第一行为表头, 
     * 第二行为编号、场地类别、经度、纬度
     */
    LithoStation(const char * lithologic_filepath): __index(nullptr), 
                                                    __resultSet(__num_results), 
                                                    __seed(time(nullptr)), 
                                                    __gen(__seed) {
        readLithologicData(lithologic_filepath);
    }

    ~LithoStation() { if(__index != nullptr) delete __index; }

    /**
     * @brief 将文件中台站数据读取并保存至 `kd-tree` 中
     * @note 数据文件格式一定要正确, 否则将无法读取正确数据
     * @param lithologic_filepath 场地文件, 数据文件为第一行为表头, 
     * 第二行为编号、场地类别、经度、纬度
     */
    void readLithologicData(const char * lithologic_filepath) {
        // 判断 kd-tree 中是否已有数据, 防止重复读取
        if(__index != nullptr) return;
        __cloud.pts.reserve(1844032);

        // 跳过第一行, 表头
        string ig;
        ifstream ipt(lithologic_filepath);
        ipt >> ig >> ig >> ig >> ig;

        int ind, litho_class;
        float lon, lat;

        while(ipt >> ind >> litho_class >> lon >> lat)
            __cloud.pts.emplace_back(LithoPoint<float>{litho_class, lon, lat});
        ipt.close();

        // 建树
        __index = new my_kd_tree_t(2, __cloud, KDTreeSingleIndexAdaptorParams(10));
        __index->buildIndex();
    }

    /**
     * @brief 搜索最近的台站, 并将结果保存在 `ret_index` 中, 同时获取土层岩性参数
     * @param lon 查询地点的经度
     * @param lat 查询地点的纬度
     * @return 返回查询地点的场地分类, 结果为 0, 1, 2, 3, 4
     */
    int searchNearestStation(float lon, float lat) {
        float search_point[] = {lon, lat};

        __searchParams = SearchParams(10);
        __resultSet.init(&__ret_index, &__out_dist_sqr);
        __index->findNeighbors(__resultSet, search_point, __searchParams);

        getC_Phi();

        return __cloud.pts[__ret_index].LithoGroup;
    }

    /**
     * @brief 根据搜索的台站确定岩性参数
     * @note 必须 `ret_index` 已定义后才能运行, 
     * 即如要外部调用, 必须运行完 `searchNearestStation()` 函数
     */
    void getC_Phi() {
        switch (__cloud.pts[__ret_index].LithoGroup) {
        case 0:
        case 1:
            C = 35;
            C_std = 10;
            phi = 40;
            phi_std = 4;
            gamma = 25;
            break;
        case 2:
            C = 30;
            C_std = 9;
            phi = 35;
            phi_std = 3;
            gamma = 23;
            break;
        case 3:
            C = 25;
            C_std = 7;
            phi = 30;
            phi_std = 3;
            gamma = 20;
            break;
        case 4:
            C = 15;
            C_std = 5;
            phi = 20;
            phi_std = 2;
            gamma = 16;
            break;
        }
    }

    /**
     * @brief 根据地质分类确定岩性参数, 本函数没有用到内部变量
     * @param litho_class 地质分类类别
     */
    void getC_Phi(int litho_class) {
        switch (litho_class) {
        case 0:
        case 1:
            C = 35;
            C_std = 10;
            phi = 40;
            phi_std = 4;
            gamma = 25;
            break;
        case 2:
            C = 30;
            C_std = 9;
            phi = 35;
            phi_std = 3;
            gamma = 23;
            break;
        case 3:
            C = 25;
            C_std = 7;
            phi = 30;
            phi_std = 3;
            gamma = 20;
            break;
        case 4:
            C = 15;
            C_std = 5;
            phi = 20;
            phi_std = 2;
            gamma = 16;
            break;
        }
    }

    /**
     * @brief 计算屈服加速度
     * @note 本函数中没有使用`Station`类中的任何成员变量, 假定重力加速度为 9.8m/s<sup>2</sup>
     * @param c_input 土的 c 值
     * @param phi_input 土的内摩擦角
     * @param gamma_input 土的平均重度
     * @param slope_input 坡度, 角度制
     * @param m_input 滑坡体饱和度, 完全饱和为 1, 完全不饱和为 0
     * @param t_input 土条宽度
     * @return 返回计算后的屈服加速度, 单位 m/s<sup>2</sup>
     */
    inline double getAcBy(double c_input, double phi_input, double gamma_input, double slope_input, double m_input, double t_input) {
        double Fs = c_input / gamma_input / sin(slope_input / 180 * M_PI) / t_input
             + tan(phi_input / 180 * M_PI) / tan(slope_input / 180 * M_PI)
             - m_input * 10 * tan(phi_input / 180 * M_PI) / gamma_input / tan(slope_input / 180 * M_PI);

        return (Fs - 1) * sin(slope_input / 180 * M_PI) * 9.8;
    }

    /**
     * @brief 计算台站处屈服加速度
     * @note 本函数使用`Station`类中的成员变量, 假定重力加速度为 9.8m/s<sup>2</sup>, 
     * 滑坡体完全非饱和, 土条宽度 2.7
     * @param slope_input 坡度, 角度制
     * @return 返回计算后的屈服加速度, 单位 m/s<sup>2</sup>
     */
    inline double getAc(double slope_input) {
        return getAcBy(C, phi, gamma, slope_input, 0, 2.7);
    }

    /**
     * @brief 根据标准差生成台站处岩性数据随机序列, 便于蒙特卡洛模拟
     * @note 本函数截断小于 0 的值
     * @param n 生成的向量长度, 即模拟次数
     * @param OutputCSeq 输出的 C 值序列, 会清空原先数据
     * @param OutputPhiSeq 输出的 Phi 值序列, 会清空原先数据
     * @return 无返回值, 通过引用返回
     */
    void getSeqC_Phi(int n, vector<double>& OutputCSeq, vector<double>& OutputPhiSeq) {
        OutputCSeq.clear(); OutputPhiSeq.clear();
        OutputCSeq.reserve(n); OutputPhiSeq.reserve(n);
        double C_temp, phi_temp;
        // 正态分布模拟
        normal_distribution<double> C_dis(C, C_std);
        normal_distribution<double> phi_dis(phi, phi_std);

        for (int i = 0; i < n; i++) {
            C_temp = -1; phi_temp = -1;
            do {
                C_temp = C_dis(__gen);
            } while (C_temp < 0);
            do {
                phi_temp = phi_dis(__gen);
            } while (phi_temp < 0);
            OutputCSeq.emplace_back(C_temp); OutputPhiSeq.emplace_back(phi_temp);
        }
    }

    /// Station输出, 输出为Station的场地分类
    friend ostream& operator<< (ostream& os, LithoStation& s) {
        os << s.__cloud.pts[s.__ret_index].LithoGroup;
        return os;
    }
};

#endif // STATION_LITHOLOGIC_H