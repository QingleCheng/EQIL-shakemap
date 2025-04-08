"""调用shakemap程序进行计算

Examples
--------
如果已经配置完成, 可以直接运行::

    $ python shakemapdataprocessor.py

"""
from shakemap.coremods import assemble, select, model
from config import ID

def shakemap_data_process(identity: str, comment: str="Testing"):
    """调用 ShakeMap 程序计算

    Parameters
    ----------
    identity : str
        事件ID
    comment : str, optional
        ShakeMap 计算的标签, 默认值为 "Testing"
    """
    select.SelectModule(identity).execute()
    assemble.AssembleModule(identity, comment).execute()
    print("==========      Running shakemap program      ==========")
    model.ModelModule(identity).execute()
    print("==========     Shakemap calculation done      ==========")

if __name__ == "__main__":    
    shakemap_data_process(ID)