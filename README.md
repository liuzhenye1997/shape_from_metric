# shape_from_metric
这是作者提供的第一个例子的matlab版本的代码，即BunnyFromMetric.目前该代码写的十分粗糙，变量名基本抄的作者的，注释几乎没有，唯一的作用是跑通。
作者提供的这个例子是所有例子中最基础的。例如后面对于存在边界的模型又要进行其它处理。不过核心算法还是有的，所以主要是用来对照论文的,如果要实现其它例子也只需以此为框架。
要运行此代码只需打开main.m运行即可，唯一可调的参数是迭代次数，预设为60次。最后代码所在文件夹会有一个obj文件输出，即为结果。
至于运行速度方面，感觉在不更改算法的情况下很难有较大提升。以迭代60次为例，耗时约为385s，其中极分解，解方程和稀疏矩阵相乘这三部分耗时247s左右。而这三部分我暂时都没想到如何提高速度。事实上平均7s迭代一次的速度与作者给的houdini的例子的运行速度基本一样。
