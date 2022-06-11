# Graduation-Design
## Hybrid mechanism
通过一下命令编译算法对比compare.cpp,测试算法test_FDA (需要在目录下放置数据集 Retail.csv census.csv)
```bash
mkdir build
cd build
cmake ..
make
ctest
./compare
```
会生成不同配置下不同算法产生的均方误差csv类似
```csv
TM,548.236,659.463,783.056,941.657,1127.37,1342.03
BM,2233.94,2697.51,3151.72,3652.29,4125.93,4673.51
Hybrid Mechanism,702.453,1673.39,2377.36,3152.19,3932.33,4751.85
FDA,1444.06,1734.8,1925.78,2143.36,2456.66,3028.23
```
通过绘制图表代码即可绘制论文中的算法对比图表
```
python draw.py
```
## Matrix binary mechanism

通过运行 pyomo_FDA.py pyomo_LRM.py LRM.cpp FDA.cpp 需要安装 Pyomo 以及求解器Ipopt

即可求解出矩阵B，M