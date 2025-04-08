"""shakemap数据后处理

这一模块用于 ShakeMap 计算后的数据处理.

Examples
--------
不更改默认值情况下, 可以如下使用本包::
    
    $ python shakemap2acc.py

"""
import os, re
import math
from typing import List, Tuple, Union
import numpy as np
import h5py
import time
from multiprocessing import Pool

from config import *
from REQPY_Module import REQPY_single
from shakemapdatagenerator import ShakeMapDataGenerator


class ShakeMap2Acc(object):
    """生成各网格点加速度时程
    
    Examples
    --------
    >>> sa = ShakeMap2Acc(
    ...     ShakeMapDataGenerator,
    ...     ROOT,
    ...     DESTINATION_ROOT
    ... )
    >>> sa.get_shakemap_data()
    >>> sa.gen_calculate_station(NUM_SIDE_STATION)
    >>> sa.start_calculate_share()
    """
    
    def __init__(self, 
                 shake_map_data: ShakeMapDataGenerator, 
                 root: str = "output_data", 
                 destination_root: str = "output_data"):
        """该类构造函数

        Parameters
        ----------
        shake_map_data : ShakeMapDataGenerator
            ShakeMap数据生成器类的一个实例
        root : str, optional
            输入数据文件夹（ShakeMap生成的数据文件）, 默认值为 "output_data"
        destination_root : str, optional
            输出数据文件夹, 默认值为 "output_data"
        """
        ## shakemapgenerator实例, 用于获取台站数据及部分常量的值
        self.shake_map_data = shake_map_data
        ## 输入数据目录, 即Shakemap输出数据目录
        self.root = root
        ## 输出数据目录, 即存储结果目录
        self.destination_root = destination_root
        ## 台站数据字典
        self.station_data = None
        ## 点列的经度
        self.lons = None
        ## 点列的纬度
        self.lats = None
        ## 点列的经度个数
        self.nlon = None
        ## 点列的纬度个数
        self.nlat = None
        ## 点列的烈度
        self.mmi  = None
        ## 点列的PGA
        self.pga  = None
        ## 点列的PGV
        self.pgv  = None
        ## 点列的SA(0.1)
        self.sa01 = None
        ## 点列的SA(0.3)
        self.sa03 = None
        ## 点列的SA(1.0)
        self.sa10 = None
        ## 点列的SA(3.0)
        self.sa30 = None
        ## 点列的SA(6.0)
        self.sa60 = None
        ## 需要计算加速度时程的点列序号的列表
        self.calculate_list = []
    
    def get_shakemap_data(self, shakemap_result_file_name: str='shake_result.hdf'):
        """从 ``shake_result.hdf`` 文件中读取计算结果

        Notes
        -----
        ShakeMap数据文件中的加速度单位是 ``ln(g)``, 本函数将其转化为 ``g``
        
        Parameters
        ----------
        shakemap_result_file_name : str, optional
            shakemap计算结果文件的名称, 相对于 ``self.root``, 默认值为 'shake_result.hdf'
        """
        with h5py.File(os.path.join(
            self.root, shakemap_result_file_name
        ), 'r') as shake_result_data:
            self.lons = shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['MMI']['lons'][:]
            self.lats = shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['MMI']['lats'][:]
            self.mmi  = shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['MMI']['mean'][:]
            self.pga  = np.exp(
                shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['PGA']['mean'][:]
            )
            self.pgv  = np.exp(
                shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['PGV']['mean'][:]
            )
            self.sa01 = np.exp(
                shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['SA(0.1)']['mean'][:]
            )
            self.sa03 = np.exp(
                shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['SA(0.3)']['mean'][:]
            )
            self.sa10 = np.exp(
                shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['SA(1.0)']['mean'][:]
            )
            self.sa30 = np.exp(
                shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['SA(3.0)']['mean'][:]
            )
            self.sa60 = np.exp(
                shake_result_data['arrays']['imts']['GREATER_OF_TWO_HORIZONTAL']['SA(6.0)']['mean'][:]
            )
            self.nlat = int(np.unique(self.lons, return_counts=True)[1][0])
            self.nlon = int(np.unique(self.lats, return_counts=True)[1][0])
    
    def gen_calculate_station(self, nums_station_by_side: int = 2):
        """生成计算的台站序号, 生成方法为等距取样
        
        Notes
        -----
        这将严重影响计算速度, 建议取值小于 40
        
        Parameters
        ----------
        nums_station_by_side : int, optional
            每边取样的站点数, 所以总站点数即为其平方, 默认值为 2
        """
        self.calculate_list = []
        for i in range(nums_station_by_side):
            lon_index = math.floor(self.nlon * (i + 1) / (nums_station_by_side + 1))
            for j in range(nums_station_by_side):
                lat_index = math.floor(self.nlat * (j + 1) / (nums_station_by_side + 1))
                self.calculate_list.append(lat_index * self.nlon + lon_index)
        
        if IS_PGA_AM == 1:
            np.savetxt(
                os.path.join(self.destination_root, "pga.txt"), 
                np.vstack((self.lons, self.lats, self.pga * self.shake_map_data.g)).T, 
                fmt="%.8f",
                header="%d %d" % (self.nlon, self.nlat)
            )
        
        print("Calculate following station:")
        for i in self.calculate_list:
            print("Index: " + str(i) + " Lontitude:" + str(self.lons[i]) + " Latitude:" + str(self.lats[i]))
            
        if self.station_data is None:
            self.station_data = self.shake_map_data.gen_station_data()
            
        for station_name, station_loc in self.station_data.items():
            t, acc = self.get_station_acc_by_name(station_name)
            dir_name = os.path.join(self.destination_root, str(station_loc[0]) + "&" + str(station_loc[1]))
            
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
                ShakeMap2Acc.save_station_acc(dir_name, t, acc)
    
    @staticmethod
    def save_station_acc(dir_name: str, dt: np.ndarray, acc: np.ndarray):
        """保存台站的加速度时程

        Parameters
        ----------
        dir_name : str
            保存的路径
        dt : np.ndarray
            时间向量
        acc : np.ndarray
            一个方向的加速度时程
        """
        np.savetxt(os.path.join(dir_name, 'Acc.txt'), np.vstack((dt, acc)).T, '%.8f')
    
    def start_calculate_share(self):
        """开始计算,  **并行** 计算"""  
        pool = Pool(NUM_THREADS)
        for idx in self.calculate_list:
            t, acc = self.get_nearest_station(self.lons[idx], self.lats[idx])
            To, dso = self.interpolation_one_spec(index=idx)
            pool.apply_async(
                calculate_single_acc,
                args=(
                    self.destination_root,
                    self.lons[idx],
                    self.lats[idx],
                    self.shake_map_data.g,
                    To,
                    dso,
                    t[1] - t[0],
                    acc
                )
            )
        pool.close()
        pool.join()
    
    def start_calculate(self):
        """开始计算,  **单线程** 计算"""
        for idx in self.calculate_list:
            t, acc = self.get_nearest_station(self.lons[idx], self.lats[idx])
            To, dso = self.interpolation_one_spec(index=idx)
            calculate_single_acc(
                self.destination_root,
                self.lons[idx],
                self.lats[idx],
                self.shake_map_data.g,
                To,
                dso,
                t[1] - t[0],
                acc
            )
    
    def interpolation_one_spec(self, dt: float=0.001, index: int=None, pga: float=None,
                               sa01: float=None, sa03: float=None, sa10: float=None, 
                               sa30: float=None, sa60: float=None) -> Union[Tuple[np.ndarray, np.ndarray], None]:
        """反应谱插值（线性插值）, 可以输入index, 也可以输入所有其他数据进行插值

        Parameters
        ----------
        dt : float, optional
            插值得到的反应谱分辨率, 默认值为 0.001
        index : int, optional
            需要插值的反应谱的序号, 默认值为 None
        pga : float, optional
            需要插值的反应谱的PGA, 默认值为 None
        sa01 : float, optional
            需要插值的反应谱的SA(0.1), 默认值为 None
        sa03 : float, optional
            需要插值的反应谱的SA(0.3), 默认值为 None
        sa10 : float, optional
            需要插值的反应谱的SA(1.0), 默认值为 None
        sa30 : float, optional
            需要插值的反应谱的SA(3.0), 默认值为 None
        sa60 : float, optional
            需要插值的反应谱的SA(6.0), 默认值为 None

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            返回时间序列和插值后反应谱曲线, 单位和输入单位相同
        None
            插值失败
        """
        if index is not None:
            pga = self.pga[index]
            sa01 = self.sa01[index]
            sa03 = self.sa03[index]
            sa10 = self.sa10[index]
            sa30 = self.sa30[index]
            sa60 = self.sa60[index]
        if pga is None or sa01 is None or sa03 is None or sa10 is None or sa30 is None or sa60 is None:
            return None
        x = np.array([0.0, 0.1, 0.3, 1.0, 3.0, 6.0])
        y = np.array([pga, sa01, sa03, sa10, sa30, sa60])
        t = np.linspace(0.05, 6.0, num=int((6 - 0.05) / dt + 1))
        return (t, np.interp(t, x, y))
    
    def get_nearest_station(self, lon: float, lat: float) -> Tuple[List[float], np.ndarray]:
        """根据坐标得到最近的台站的地震动数据, 方法为遍历, 速度较慢, 但因为台站较少, 因此可以采用

        Parameters
        ----------
        lon : float
            当前点经度
        lat : float
            当前点纬度

        Returns
        -------
        Tuple[List[float], np.ndarray]
            返回时间序列及加速度时程向量
        
        Examples
        --------
        >>> self.get_nearest_station()
        ([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, ...], array([ 5.72400e-05,...1194e-03]))
        """
        if self.station_data is None:
            self.station_data = self.shake_map_data.gen_station_data()
        nearest_dis = float("inf")
        nearest_station = ""
        for station_name, station_loc in self.station_data.items():
            dis = (station_loc[0] - lon) ** 2 + (station_loc[1] - lat) ** 2
            if (dis < nearest_dis):
                nearest_dis = dis
                nearest_station = station_name
        return self.get_station_acc_by_name(nearest_station)

    def get_station_acc_by_name(self, station_name:str) -> Tuple[List[float], np.ndarray]:
        """根据台站名称获取台站的加速度时程信息

        获取最大的加速度时程信息

        Parameters
        ----------
        station_name : str
            台站名称

        Returns
        -------
        Tuple[List[float], np.ndarray]
            返回时间序列及加速度时程向量
        
        Examples
        --------
        >>> self.get_nearest_station("NIG019")
        ([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, ...], array([ 1.7160000e-0...4346e-01]))
        """
        t = []
        with open(os.path.join(ROOT,station_name,"EW.txt"), 'r', encoding='gbk') as f:
            t = []
            acc1 = []
            i = 0
            for line in f.readlines():
                if i >= 1:  # 跳一行
                    line=re.sub(' +', ' ', line).replace(' ','\n')
                    line=line.replace('\t','\n')
                    temp = line.split('\n')
                    if(temp[0]==''):
                        t.append(float(temp[1]))
                        acc1.append(float(temp[2]))
                    else:
                        t.append(float(temp[0]))
                        acc1.append(float(temp[1]))
                i += 1
            acc1 = np.array(acc1)  # 注意单位
            acc1 -= np.average(acc1[0:5])  # 纠偏
        
        with open(os.path.join(ROOT,station_name,"NS.txt"), 'r', encoding='gbk') as f:
            t = []
            acc2 = []
            i = 0
            for line in f.readlines():
                if i >= 1:  # 跳一行
                    line=re.sub(' +', ' ', line).replace(' ','\n')
                    line=line.replace('\t','\n')
                    temp = line.split('\n')
                    if(temp[0]==''):
                        t.append(float(temp[1]))
                        acc2.append(float(temp[2]))
                    else:
                        t.append(float(temp[0]))
                        acc2.append(float(temp[1]))
                i += 1
            acc2 = np.array(acc2)  # 注意单位
            acc2 -= np.average(acc2[0:5])  # 纠偏  
        if (max(np.max(acc1), -np.max(acc1)) > max(np.max(acc2), -np.max(acc2))):
            return(t, acc1)
        return (t, acc2)
    

def calculate_single_acc(destination_root: str,
                         lon: float,
                         lat: float,
                         g: float,
                         To: np.ndarray,
                         dso: np.ndarray,
                         DTs: float,
                         acc: np.ndarray):
    """根据输入数据计算单个点加速度时程, 并保存

    Parameters
    ----------
    destination_root : str
        数据保存路径
    lon : float
        待计算台站经度
    lat : float
        待计算台站纬度
    g : float
        重力加速度, 单位 :math:`m/s^2`
    To : np.ndarray
        反应谱时间序列
    dso : np.ndarray
        反应谱
    DTs : float
        加速度时间分辨率
    acc : np.ndarray
        加速度时程, 单位 :math:`m/s^2`
    """
    dir_name = os.path.join(destination_root, str(lon) + "&" + str(lat))
    if not os.path.exists(dir_name):
        acc /= g
        ccs, *_ = REQPY_single(acc, 1/DTs, dso, To, T1=0.05, nit=15, T2=6, plots=0)
        t = np.linspace(0.0, DTs * (len(ccs) - 1), num=len(ccs))
        ccs *= g
        os.mkdir(dir_name)
        ShakeMap2Acc.save_station_acc(dir_name, t, ccs)

        
if __name__ == "__main__": 
    if not os.path.exists(DESTINATION_ROOT):
        os.mkdir(DESTINATION_ROOT)
        
    output_dir = os.path.join(DESTINATION_ROOT, time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()))
    os.mkdir(output_dir)
    shake_map2acc = ShakeMap2Acc(ShakeMapDataGenerator(
        ROOT, os.path.join(SHAKEMAP_DATA_HOME, ID, "current")
    ), os.path.join(SHAKEMAP_DATA_HOME, ID, "current", "products")
     , output_dir)
    shake_map2acc.get_shakemap_data()
    shake_map2acc.gen_calculate_station(NUM_SIDE_STATION)
    shake_map2acc.start_calculate_share()
    
    # 使用判断程序计算滑坡概率
    start_time = time.perf_counter()
    os.system("./utilities/NewmarkDisp.exe %d %s %s %s %d %d" % (
        NUM_THREADS, output_dir, SLOPE_DATA_FILE, LITHOLOGIC_DATA_FILE, 
        STATION_STEP, IS_PGA_AM
    ))
    end_time = time.perf_counter()
    print("Landslide probability calculating finished. Time elapsed: %.2fs" % (end_time - start_time))