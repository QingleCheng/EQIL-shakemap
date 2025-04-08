from osgeo import gdal
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Landslide 程序前处理模块, 裁剪 GeoTiff 坡度文件."
    )
    parser.add_argument(
        "src_tiff", type=str, 
        help="输入文件路径"
    )
    parser.add_argument(
        "-o", "--output", type=str, dest="dst_tiff", 
        default="slope.tif", help="输出文件路径"
    )
    args = parser.parse_args()
    
    print("请输入需要计算的中心点坐标 (一般可取震源为计算中心点): ")
    
    lon = float(input("经度值 (格式: xxx.xxxxx, 化为度): "))
    lat = float(input("纬度值 (格式: xxx.xxxxx, 化为度): "))
    
    side = input("请输入每边测点数 (默认值 720, 回车以保持默认): ")
    if "" == side:
        side = 720
    side = int(side)
    
    gdal.AllRegister()
    src_dataset = gdal.Open(args.src_tiff)
    transform = src_dataset.GetGeoTransform()
    
    c_lon_index = round((lon - transform[0]) / transform[1])
    c_lat_index = round((lat - transform[3]) / transform[5])
    
    ul_lon_index = c_lon_index - side // 2
    ul_lat_index = c_lat_index - side // 2
    
    gdal.Translate(
        args.dst_tiff, 
        src_dataset,
        options=gdal.TranslateOptions(
            srcWin=[ul_lon_index, ul_lat_index, side, side]
        )
    )