"""ShakeMap 常量定义.

此处定义了 ShakeMap 程序计算过程中的所有常量.

Notes
-----
    请先修改相关配置, 如 CPU 核数, 再运行.

Attributes
----------
ROOT : str
    源数据目录
IS_PGA_AM : int
    是否进行PGA调幅
DESTINATION_ROOT : str
    输出结果目录
ID : str
    生成的事件ID
LON : str
    震源经度
LAT : str
    震源纬度
DEPTH : str
    震源深度
MAG : str
    震级
EARTHQUAKE_TIME : str
    地震发生时间, 格式为 "YYYY-mm-ddTHH:MM:SSZ"
LOCSTRING : str
    地震发生地点
NUM_THREADS : int
    多线程计算核心数
NUM_SIDE_STATION : int
    每边取样台站数量
STATION_STEP : int
    滑坡计算分辨率, 1为最高
SHAKEMAP_DATA_HOME : str
    ShakeMap程序计算时的数据目录, 即generator的输出目录
STATION_LOCATION_DATA : str
    Station位置信息文件
VS30_DATA_FILE : str
    Vs30数据文件路径
SLOPE_DATA_FILE : str
    坡度数据文件路径
LITHOLOGIC_DATA_FILE : str
    岩性数据文件路径
"""
ROOT = "sample_data"                        # 源数据目录
IS_PGA_AM = 1                               # 是否进行PGA调幅
DESTINATION_ROOT = "output_data"            # 输出结果目录
ID = "usp000d6vkkkkk"                       # 生成的事件ID
LON = "138.7790"                            # 震源经度
LAT = "37.2260"                             # 震源纬度
DEPTH = "16.0"                              # 震源深度
MAG = "6.6"                                 # 震级
EARTHQUAKE_TIME = "2004-10-23T08:56:00Z"    # 地震发生时间
LOCSTRING = "Honshu, Japan"                 # 地震发生地点
NUM_THREADS = 8                             # 多线程计算核心数
NUM_SIDE_STATION = 2                        # 每边取样台站数量
STATION_STEP = 10                           # 滑坡计算分辨率, 1为最高

# ShakeMap程序计算时的数据目录, 即generator的输出目录
SHAKEMAP_DATA_HOME = "/home/guyi/shakemap_profiles/default/data"

# Station位置信息文件
STATION_LOCATION_DATA = "StationforHtml.txt"

# Vs30数据文件路径
VS30_DATA_FILE = "./vs30_data/global_vs30.grd"

# 坡度数据文件路径
SLOPE_DATA_FILE = "./tiff_data/2004_slope.tif"

# 岩性数据文件路径
LITHOLOGIC_DATA_FILE = "./lithologic_data/NationalPoint.txt"