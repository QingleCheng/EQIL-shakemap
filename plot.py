"""程序输出数据绘图

Examples
--------
修改需要绘图的文件夹名称, 然后运行::

    $ python plot.py

"""
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class Plot(object):
    """数据绘图类

    Attributes
    ----------
    root : str
        绘图文件夹
    
    Examples
    --------
    >>> Plot(ROOT).plot()
    >>> Plot(ROOT).plot().plot_single_file(
    ...     SOME_FILE, OUTPUT_FILE
    ... )
    """
    
    def __init__(self, root: str):
        """构造 Plot 类

        Parameters
        ----------
        root : str
            绘图文件夹
        """
        self.root = root
        
    def plot(self):
        """绘制滑坡概率图

        Returns
        -------
        Plot
            返回Plot类本身, 方便级联调用
        """
        return self.plot_single_file("result0.txt", "plot_result0.png") \
                   .plot_single_file("result5.txt", "plot_result5.png") \
                   .plot_single_file("result9.txt", "plot_result9.png")
        
    def plot_single_file(self, file_name: str, save_file_name: str):
        """绘制热度图

        只要输入文件满足为网格数据, 第一列为经度, 第二列为纬度, 第三列为绘图数据, 
        即可使用此函数绘图.

        Parameters
        ----------
        file_name : str
            需要绘制的文件名
        save_file_name : str
            保存的图片名

        Returns
        -------
        Plot
            返回Plot类本身, 方便级联调用
        """
        plt.figure(figsize=(10.8, 10.8), dpi=100)
        raw_prob_data = np.loadtxt(os.path.join(self.root, file_name))
        nX = np.unique(raw_prob_data[:, 1], return_counts=True)[1][0]
        nY = np.unique(raw_prob_data[:, 0], return_counts=True)[1][0]
        prob_data = raw_prob_data[:, 2]
        column = []
        for y_i in range(nY):
            row = []
            for x_i in range(nX):
                row.append(prob_data[y_i * nX + x_i])
            column.append(row)
        sns.heatmap(pd.DataFrame(column), cmap="YlGnBu")
        plt.savefig(os.path.join(self.root, save_file_name))
        return self
    
    
if __name__ == "__main__":
    start_time = time.perf_counter()
    Plot("./output_data/2021_07_14_12_13_45").plot()
    end_time = time.perf_counter()
    print("Ploting finished. Time elapsed: %.2fs" % (end_time - start_time))
