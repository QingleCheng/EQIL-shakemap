"""Vs30 文件读取包

这一模块用于读取 Vs30 栅格数据文件, 由于速度原因, 仅用于 `shakemapgenerator` 包.

Notes
-----
    该类需要 **Vs30** 文件才能正常工作. 经测试, 计算速度约为 C++ 版本的百分之一.

"""
import os
import math
from osgeo import gdal

class Vs30(object):
    """Vs30 栅格数据读取类
    
    Attributes
    ----------
    vs30_data : gdal.Dataset
        Vs30 数据集
    """
    
    def __init__(self, vs30_path: str="global_vs30.grd"):
        """根据 Vs30 数据集路径构造.

        Parameters
        ----------
        vs30_path : str
            Vs30 数据集路径, 默认值为 "global_vs30.grd"
        """
        gdal.AllRegister()
        self.vs30_data = gdal.Open(vs30_path)
    
    def get_point_vs30(self, lon: float, lat: float) -> float:
        """根据经纬度坐标获取对应的剪切波速, 采用 **双线性** 插值

        Parameters
        ----------
        lon : float
            待求点经度
        lat : float
            待求点纬度

        Returns
        -------
        float
            地下 30m 剪切波速的值
        
        Examples
        --------
        >>> self.get_point_vs30(138.83890315, 37.45423190)
        218.09299576077512
        """
        lat_index = math.floor((84 - lat) * 120)
        lon_index = math.floor((lon + 180) * 120)
        delta_lat = (84 - lat) * 120 - lat_index
        delta_lon = (lon + 180) * 120 - lon_index
        k = self.vs30_data.ReadAsArray(lon_index, lat_index, 2, 2)
        vs50 = (k[1][0] - k[0][0]) * delta_lat + k[0][0]
        vs51 = (k[1][1] - k[0][1]) * delta_lat + k[0][1]
        return (vs51 - vs50) * delta_lon + vs50

if __name__ == "__main__":
    vs30 = Vs30(os.path.join("vs30_data", "global_vs30.grd"))
    print(vs30.get_point_vs30(138.83890315, 37.45423190))