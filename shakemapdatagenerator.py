"""ShakeMap 数据生成包.

这一模块用于生成 ShakeMap 需要的数据文件.

Examples
--------
不更改默认值情况下, 可以如下使用本包::
    
    $ python shakemapdatagenerator.py

"""
import os
import re
import json
import time
from typing import Any, Dict, Tuple, Iterator
import numpy as np
from multiprocessing import Pool
import xml.etree.ElementTree as ET

from vs30 import Vs30
from config import *


class ShakeMapDataGenerator(object):
    """ShakeMap 数据生成器
    
    Examples
    --------
    >>> shake = ShakeMapDataGenerator(ROOT, 
    ...                               DESTINATION_ROOT, 
    ...                               STATION_LOCATION_DATA)
    >>> shake.gen_instrumented_data()
    >>> shake.gen_event_xml(
    ...     id=ID, lat="37.2260", lon="138.7790", depth="16.0", 
    ...     mag="6.6", earthquake_time="2004-10-23T08:56:00Z", 
    ...     locstring="Honshu, Japan"
    ... )
    >>> shake.gen_model_conf(SLOPE_DATA_FILE)
    """
    
    def __init__(self, root:str, 
                 destination_root:str = "output_data", 
                 station_file:str = "StationforHtml.txt"):
        """该类的构造函数.

        根据文件路径、预计输出路径、台站数据文件构建生成器.

        Parameters
        ----------
        root : str
            台站地震动根目录
        destination_root : str, optional
            输出文件根目录, 默认值为 "output_data"
        station_file : str, optional
            台站文件路径, 相对于 ``root`` 目录, 默认值为 "StationforHtml.txt"
        """
        ## 傅里叶采样, 默认为2<sup>15</sup>
        self.n = 2 ** 15
        ## 阻尼比, 默认为0.05
        self.zeta = 0.05
        ## 反应谱时间间隔, 默认为0.01秒
        self.dt = 0.01
        ## 反应谱计算最大时间, 默认为6秒
        self.t_max = 6
        ## 重力加速度, 默认为9.81m/s<sup>2</sup>
        self.g = 9.81
        ## 工作目录, 即需要处理的数据目录. 内部应该为多个文件夹, 一个文件夹为一个台站, 内有多组地震动数据
        self.root = root
        ## 输出目录
        self.destination_root = destination_root
        ## 台站数据路径, 相对于 ``root`` 的路径, 存有台站的经纬度信息
        self.station_file = station_file
        ## 结果列表
        self.res_list = None
        ## 存储台站信息
        self.station_dict = None
        ## 存储事件信息
        self.event_data = None
        ## 存储shakemap数据格式
        self.data_dict = None
    
    def set_fft_args(self, n:int, 
                     zeta:float, dt:float, 
                     t_max:float) -> Tuple[int, float, float, float]:
        """设定快速傅里叶变换参数.

        Parameters
        ----------
        n : int
            傅里叶采样值.
        zeta : float
            阻尼比.
        dt : float
            输出反应谱时间间隔.
        t_max : float
            反应谱最大时间.

        Returns
        -------
        Tuple[int, float, float, float]
            返回四个参数本身, 便于后续级联调用.
        """
        self.n = n
        self.zeta = zeta
        self.dt = dt
        self.t_max = t_max
        return n, zeta, dt, t_max
    
    def gen_station_data(self, 
                         destination_file_name:str = "station.json") -> Dict[str, Tuple[float, float]]:
        """台站数据生成.

        生成台站数据文件 ``station.json`` 便于内部查找调用, 将信息存储到 ``self.station_dict`` 中

        Parameters
        ----------
        destination_file_name : str, optional
            处理后的台站数据文件导出路径, 相对 ``root`` 目录, 默认值为 "station.json"

        Returns
        -------
        Dict[str, Tuple[float, float]]
            返回台站名为键, 台站坐标为值的字典.
        """
        if self.station_dict is None:
            with open(os.path.join(self.root, self.station_file), "rb") as r, \
                open(os.path.join(self.destination_root, destination_file_name), "wb") as w:
                w.write("{".encode())
                w.write(r.read())
                w.seek(-1, os.SEEK_END)
                w.write("}".encode())
            
            with open(os.path.join(self.destination_root, destination_file_name), "r") as f:
                station_dict = json.load(f)
                self.station_dict = station_dict
                
        return self.station_dict
    
    def gen_sa_dict(self, sa: float, tn: float) -> Dict[str, Any]:
        """生成 ShakeMap 字段.

        Parameters
        ----------
        sa : float
            反应谱值.
        tn : float
            对应周期.

        Returns
        -------
        Dict[str, Any]
            ShakeMap 字典格式.
        """
        return {
            "name": "sa(%.1f)" % tn,
            "value": round(sa * 100 / self.g, 6),
            "units": "%g",
            "flag": "0",
            "ln_sigma": 0.0
        }
    
    def gen_instrumented_data(self, 
                              destination_file_name:str = "instrumented_dat.json") -> Dict[str, Any]:
        """生成ShakeMap所需数据.

        Parameters
        ----------
        destination_file_name : str, optional
            输出文件名, 相对 ``destination_root`` 目录, 默认值为 "instrumented_dat.json"

        Returns
        -------
        Dict[str, Any]
            返回数据字典
        """
        ## 新建进程池
        pool = Pool(NUM_THREADS)
        res_list = []
        for attr, stationname, filepath in gen_file_path(self.root):
            ## 计算单个反应谱
            res = pool.apply_async(
                func=gen_one_spec, 
                args=(filepath, self.n, self.zeta, self.dt, self.t_max), 
                callback=''
            )
            print("Start calculate " + stationname + " " + attr + " accelaration spectrum")
            res_list.append((attr, stationname, res))
        pool.close()
        pool.join()
        print("========== All calculations done successfully ==========")
        
        print("==========         Start data process         ==========")
        self.res_list = res_list
        
        if self.station_dict is None:
            self.gen_station_data()
            
        vs30 = Vs30(VS30_DATA_FILE)
        station_dict = {k: {
            "point": v, "channel": {}, 
            "vs30": round(vs30.get_point_vs30(v[0], v[1]), 3)
        } for k, v in self.station_dict.items()}
        
        for attr, stationname, res in res_list:
            res_dict = res.get()
            station_dict[stationname]["channel"][attr] = {
                "name": attr,
                "amplitudes": [
                    {
                        "name": "pga",
                        "value": round(res_dict['pga'] * 100 / self.g, 6),
                        "units": "%g",
                        "flag": "0",
                        "ln_sigma": 0.0
                    },
                    {
                        "name": "pgv",
                        "value": round(res_dict['pgv'] * 100, 6),
                        "units": "cm/s",
                        "flag": "0",
                        "ln_sigma": 0.0
                    }
                ]
            }

            station_dict[stationname]["channel"][attr]["amplitudes"].append(
                self.gen_sa_dict(res_dict['sa'][res_dict['tn'] == 0.3][0], 0.3)
            )
            station_dict[stationname]["channel"][attr]["amplitudes"].append(
                self.gen_sa_dict(res_dict['sa'][res_dict['tn'] == 1.0][0], 1.0)
            )
            station_dict[stationname]["channel"][attr]["amplitudes"].append(
                self.gen_sa_dict(res_dict['sa'][res_dict['tn'] == 3.0][0], 3.0)
            )
        
        data_dict = {"type": "FeatureCollection"}
        feature_list = []
        for station_name, station_data in station_dict.items():
            feature_list.append({ 
                "geometry": {"type": "Point", "coordinates": station_data["point"]},
                "type": "Feature",
                "id": "BO." + station_name,
                "properties": {
                    "code": station_name,
                    "name": "",
                    "instrumentType": "UNK",
                    "source": "BO",
                    "network": "BO",
                    "commType": "UNK",
                    "location": "",
                    "intensity": None,
                    "intensity_flag": "",
                    "intensity_stddev": None,
                    "pga": None,
                    "pgv": None,
                    "distance": None,
                    "elev": None,
                    "vs30": station_data["vs30"],
                    "channels": [v for _,v in station_data["channel"].items()]
                }
            })
            
        data_dict["features"] = feature_list
        self.data_dict = data_dict
        with open(os.path.join(self.destination_root, destination_file_name), "w") as f:
            json.dump(data_dict, f)
        print("==========   Data process done successfully   ==========")
        return data_dict
    
    def gen_event_xml(self, id: str, lat: str, lon: str, depth: str, 
                      mag: str, earthquake_time: str, locstring: str, 
                      event_file_name: str = "event.xml") -> ET.Element:
        """创建 ``event.xml`` 数据文件.

        Parameters
        ----------
        id : str
            事件ID
        lat : str
            震源纬度
        lon : str
            震源经度
        depth : str
            震源深度
        mag : str
            震级
        earthquake_time : str
            地震时间
        locstring : str
            地震地点
        event_file_name : str, optional
            输出文件名称, 相对 ``destination_root`` 目录, 默认值为 "event.xml"

        Returns
        -------
        ET.Element
            返回事件ET树
        """
        earthquake_event = ET.Element("earthquake")
        earthquake_event.set("id", id)
        earthquake_event.set("netid", "BO")
        earthquake_event.set("network", "BO")
        earthquake_event.set("lat", lat)
        earthquake_event.set("lon", lon)
        earthquake_event.set("depth", depth)
        earthquake_event.set("mag", mag)
        earthquake_event.set("time", earthquake_time)
        earthquake_event.set("locstring", locstring)
        earthquake_event.set("event_type", "ACTUAL")
        self.event_data = earthquake_event
        
        tree = ET.ElementTree(earthquake_event)
        tree.write(os.path.join(self.destination_root, event_file_name))
        print("==========         event.xml generated        ==========")
        
        return earthquake_event
    
    def gen_model_conf(self, 
                       slope_data_file_path: str, 
                       point_file_name: str = 'sample.shakemap'):
        """创建 ``model.conf`` 文件及所有点数据

        Parameters
        ----------
        slope_data_file_path : str
            坡度数据文件路径
        point_file_name : str, optional
            输出点数据文件名, 默认值为 "sample.shakemap"
        """
        print("==========       Generating model config      ==========")
        os.system("./utilities/GetStationData.exe %s %s %s" % (
            slope_data_file_path,
            VS30_DATA_FILE,
            os.path.join(self.destination_root, point_file_name)
        ))
        with open(os.path.join(self.destination_root, 'model.conf'), 'w') as w:
            w.write('[interp]\n'
                    '    imt_list = '
                    'PGA, PGV, MMI, '
                    'SA(0.1), SA(0.3), SA(1.0), SA(3.0), SA(6.0)\n'
                    '    [[prediction_location]]\n'
                    '        file = %s\n' %
                    os.path.join(self.destination_root, point_file_name))
        print("==========       Model config generated       ==========")
    

def gen_file_path(dirPath:str="./") -> Iterator[Tuple[str, str, str]]:
    """地震动文件迭代器, 使用迭代以节省内存.

    Parameters
    ----------
    dirPath : str, optional
        需要迭代的目录, 默认值为 "./"

    Yields
    -------
    Iterator[Tuple[str, str, str]]
        返回三元组, 分别为地震动名称、台站名称、地震动相对路径
    """
    for path in os.listdir(dirPath):
        subDirPath = os.path.join(dirPath, path)
        if os.path.isdir(subDirPath):
            for filename in os.listdir(subDirPath):
                if str.lower(filename[-3:]) == 'txt':
                    if 'EW' in filename:
                        yield 'EW', path, os.path.join(subDirPath, filename)
                    elif 'NS' in filename:
                        yield 'NS', path, os.path.join(subDirPath, filename)
                    elif 'UD' in filename:
                        yield 'UD', path, os.path.join(subDirPath, filename)

def gen_one_spec(file_path:str, 
                 n:int, 
                 zeta:float, 
                 dT:float, 
                 Tmax:float) -> Dict[str, Any]:
    """生成一个反应谱

    Parameters
    ----------
    file_path : str
        加速度时程路径
    n : int
        傅里叶采样
    zeta : float
        阻尼比
    dT : float
        时间分辨率
    Tmax : float
        反应谱最大时间

    Returns
    -------
    Dict[str, Any]
        返回反应谱字典
    """
    with open(file_path, 'r', encoding='gbk') as f:
        t = []
        acc = []
        i = 0
        for line in f.readlines():
            if i >= 1:  # 跳一行
                line=re.sub(' +', ' ', line).replace(' ','\n')
                line=line.replace('\t','\n')
                temp = line.split('\n')
                if(temp[0]==''):
                    t.append(float(temp[1]))
                    acc.append(float(temp[2]))
                else:
                    t.append(float(temp[0]))
                    acc.append(float(temp[1]))
            i += 1
        dt = t[-1] / (len(t) - 1)
        acc = np.array(acc)  # 注意单位
        acc -= np.average(acc[0:5])  # 纠偏
    v = [0]
    for i in range(len(acc) - 1):
        v.append( v[i] + dt * (acc[i] + acc[i + 1]) / 2)
    v = np.array(v)
    v -= np.average(v[0:5])
    Tn, Sa = elastic_spectrum_fft_onlyneeded(dt, acc, n, zeta, dT, Tmax)
    print("==========    Caculation done successfully    ==========")
    return {
        "pgv": max(np.max(v), -np.min(v)),
        "pga": max(np.max(acc), -np.min(acc)),
        "sa": Sa,
        "tn": Tn
    }

def elastic_spectrum_fft_onlyneeded(t: np.ndarray, 
                                    ag: np.ndarray, 
                                    n: float, 
                                    zeta: float, 
                                    dT: float, 
                                    Tmax: float) -> Tuple[np.ndarray, np.ndarray]:
    """只计算 ShakeMap 有用的反应谱, 即 ``Tn = 0.3, 1.0, 3.0``

    Parameters
    ----------
    t : np.ndarray
        时间序列
    ag : np.ndarray
        加速度时程
    n : float
        傅里叶采样
    zeta : float
        阻尼比
    dT : float
        时间分辨率
    Tmax : float
        反应谱最大时间

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        返回反应谱的时间及反应谱值数组
    """
    sp = np.fft.fft(-ag, n)
    freq = np.fft.fftfreq(n, d=t)
    cf = 2 * np.pi * freq

    Tn = np.linspace(dT, Tmax, int(Tmax / dT))
    cfn = 2 * np.pi / Tn
    Sa = np.array([])
    TnA = np.array([])
    for ind, cfn1 in enumerate(cfn):
        if not (Tn[ind] == 0.3 or Tn[ind] == 1.0 or Tn[ind] == 3.0): continue
        H = 1 / (-cf ** 2 + (1j) * 2 * zeta * cfn1 * cf + cfn1 ** 2)
        U = sp * H
        u = np.fft.ifft(U, n)

        Sd1 = np.max(np.abs(u.real))
        Sa1 = (cfn1 ** 2) * Sd1
        Sa = np.append(Sa, Sa1)
        TnA = np.append(TnA, Tn[ind])

    # add initial point
    TnA = np.append(np.array([0]), TnA)
    ag_max = np.max(ag)
    Sa = np.append(np.array([ag_max]), Sa)

    return (TnA, Sa)

def elastic_spectrum_fft(t: np.ndarray, 
                         ag: np.ndarray, 
                         n: float, 
                         zeta: float, 
                         dT: float, 
                         Tmax: float):
    """计算一条地震动的反应谱

    Parameters
    ----------
    t : np.ndarray
        时间序列
    ag : np.ndarray
        加速度时程
    n : float
        傅里叶采样
    zeta : float
        阻尼比
    dT : float
        时间分辨率
    Tmax : float
        反应谱最大时间

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        返回反应谱的时间及反应谱值数组
    """
    sp = np.fft.fft(-ag, n)
    freq = np.fft.fftfreq(n, d=t)
    cf = 2 * np.pi * freq

    Tn = np.linspace(dT, Tmax, int(Tmax / dT))
    cfn = 2 * np.pi / Tn
    Sa = np.array([])
    TnA = np.array([])
    for ind, cfn1 in enumerate(cfn):
        H = 1 / (-cf ** 2 + (1j) * 2 * zeta * cfn1 * cf + cfn1 ** 2)
        U = sp * H
        u = np.fft.ifft(U, n)

        Sd1 = np.max(np.abs(u.real))
        Sa1 = (cfn1 ** 2) * Sd1
        Sa = np.append(Sa, Sa1)
        TnA = np.append(TnA, Tn[ind])

    # add initial point
    TnA = np.append(np.array([0]), TnA)
    ag_max = np.max(ag)
    Sa = np.append(np.array([ag_max]), Sa)

    return (TnA, Sa)

if __name__ == '__main__':
    # 备份原先数据, 并创建新文件夹
    if not os.path.exists(os.path.join(SHAKEMAP_DATA_HOME, ID)):
        os.mkdir(os.path.join(SHAKEMAP_DATA_HOME, ID))
    if os.path.exists(os.path.join(SHAKEMAP_DATA_HOME, ID, "current")):
        os.rename(os.path.join(SHAKEMAP_DATA_HOME, ID, "current"), 
                  os.path.join(SHAKEMAP_DATA_HOME, ID, time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())))
    os.mkdir(os.path.join(SHAKEMAP_DATA_HOME, ID, "current"))

    shake_map_data_generator = ShakeMapDataGenerator(ROOT, os.path.join(SHAKEMAP_DATA_HOME, ID, "current"), STATION_LOCATION_DATA)
    
    start_time = time.perf_counter()
    shake_map_data_generator.gen_instrumented_data()
    shake_map_data_generator.gen_event_xml(
        id=ID, lat=LAT, lon=LON, depth=DEPTH, 
        mag=MAG, earthquake_time=EARTHQUAKE_TIME, 
        locstring=LOCSTRING
    )
    shake_map_data_generator.gen_model_conf(SLOPE_DATA_FILE)
    end_time = time.perf_counter()
    print("Time elapsed: %.2fs" % (end_time - start_time))