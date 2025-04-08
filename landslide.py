"""主程序入口.

Landslide滑坡程序入口, 用于计算各点滑坡概率.

Notes
-----
使用前先阅读使用说明, 见 README_

Examples
--------
请先修改配置文件 ``config.py``, 此后直接在终端运行::

    $ python landslide.py

.. _README: https://gitee.com/guyi2000/land-slide
"""
import os, time

from shakemapdatagenerator import ShakeMapDataGenerator
from shakemapdataprocessor import shakemap_data_process
from shakemap2acc import ShakeMap2Acc
from plot import Plot
from config import *


def main():
    """主程序入口"""
    
    # 备份原先数据, 并创建新文件夹
    if not os.path.exists(os.path.join(SHAKEMAP_DATA_HOME, ID)):
        os.mkdir(os.path.join(SHAKEMAP_DATA_HOME, ID))
    if os.path.exists(os.path.join(SHAKEMAP_DATA_HOME, ID, "current")):
        os.rename(os.path.join(SHAKEMAP_DATA_HOME, ID, "current"), 
                  os.path.join(SHAKEMAP_DATA_HOME, ID, time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())))
    os.mkdir(os.path.join(SHAKEMAP_DATA_HOME, ID, "current"))
    
    # 使用 ShakeMapGenerator 生成 ShakeMap 需要的文件
    shake_map_data_generator = ShakeMapDataGenerator(ROOT, os.path.join(SHAKEMAP_DATA_HOME, ID, "current"), STATION_LOCATION_DATA)

    start_time = time.perf_counter()
    total_start_time = start_time
    shake_map_data_generator.gen_instrumented_data()
    shake_map_data_generator.gen_event_xml(
        id=ID, lat=LAT, lon=LON, depth=DEPTH, 
        mag=MAG, earthquake_time=EARTHQUAKE_TIME, 
        locstring=LOCSTRING
    )
    shake_map_data_generator.gen_model_conf(SLOPE_DATA_FILE)
    end_time = time.perf_counter()
    print("Generate data finished. Time elapsed: %.2fs" % (end_time - start_time))
    
    # 调用 ShakeMap 程序计算
    start_time = time.perf_counter()
    shakemap_data_process(ID)
    end_time = time.perf_counter()
    print("Shakemap calculating finished. Time elapsed: %.2fs" % (end_time - start_time))

    # 插值加速度时程
    start_time = time.perf_counter()
    if not os.path.exists(DESTINATION_ROOT):
        os.mkdir(DESTINATION_ROOT)
    
    output_dir = os.path.join(DESTINATION_ROOT, time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()))
    os.mkdir(output_dir)
    shake_map2acc = ShakeMap2Acc(
        shake_map_data_generator, 
        os.path.join(SHAKEMAP_DATA_HOME, ID, "current", "products"), 
        output_dir
    )
    shake_map2acc.get_shakemap_data()
    shake_map2acc.gen_calculate_station(NUM_SIDE_STATION)
    shake_map2acc.start_calculate_share()
    end_time = time.perf_counter()
    print("Acc calculating finished. Time elapsed: %.2fs" % (end_time - start_time))
    
    # 计算滑坡概率
    start_time = time.perf_counter()
    os.system("./utilities/NewmarkDisp.exe %d %s %s %s %d %d" % (
        NUM_THREADS, output_dir, SLOPE_DATA_FILE, LITHOLOGIC_DATA_FILE, 
        STATION_STEP, IS_PGA_AM
    ))
    end_time = time.perf_counter()
    print("Landslide probability calculating finished. Time elapsed: %.2fs" % (end_time - start_time))
    
    # 绘图
    start_time = time.perf_counter()
    Plot(output_dir).plot()
    end_time = time.perf_counter()
    print("Ploting finished. Time elapsed: %.2fs" % (end_time - start_time))
    print("All calculation done. Total time elapsed: %.2fs" % (end_time - total_start_time))
    
if __name__ == "__main__":
    main()