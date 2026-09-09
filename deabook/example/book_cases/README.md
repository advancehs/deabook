# 能源与环境案例：每章一个主样本

本目录为2026-09-09修订的正文配套案例。旧chap*/notebook保留历史版本，不作为当前正文的复算入口。

| 章节 | 主样本 | 主题 |
|---|---|---|
| 1—3、11 | 美国1978年19家电厂 | 能源生产与经典效率、编程 |
| 4 | 中国2015—2017年省级火电面板 | 发电生产率 |
| 5、12—13 | 中国省级能源/CO2面板 | 能源效率、程序及对偶 |
| 6—7 | OECD能源/CO2数据 | 环境效率与环境生产率 |
| 8 | 美国1978年19家电厂 | 普通SBM与投入、发电松弛 |
| 9 | 中国省级火电面板 | 物质平衡 |
| 10 | 芬兰89家配电企业 | 能源行业成本监管 |

每章只使用一个数据集，章内全部实证模型比较、图表与配套代码围绕该数据集展开。不同章节允许共用同一数据集。数据来源须核实为已发表期刊论文使用的能源或环境样本；未完成来源核实的案例不得标记为定稿。第4章与第9章、部分编程章节可保留跨章共用安排。菲律宾稻农和混合动力车型均已退出正文主案例。

运行：在项目根目录 `python deabook/example/book_cases/run.py`；可在末尾指定章节号，如 `4 8`。依赖numpy、pandas、scipy、matplotlib、openpyxl、tabulate；第10章另需cvxpy和CLARABEL。环境接口交叉检查另需pyomo、highspy。ch01.ipynb—ch13.ipynb为各章入口。

完整观测保存在results/chNN_full.csv，正文摘要保存在chNN_summary.csv。每章正文表6—9行，按章节采用数据记录、分布统计或时期分解；图片同时保存300 dpi PNG和矢量PDF至正文images。缩短表格不缩减DEA参照集。第8章使用全部19家电厂。

## 数据来源

美国1978年电厂来自deaR Electric_plants，关联Färe、Grosskopf和Kokkelenberg（1989）及Simar和Wilson（1998）；芬兰企业与pyStoNED公开文件一致，关联Kuosmanen（2012）。中国/OECD数据沿用项目原始来源。美国2022年煤电厂候选因未核实期刊论文出处而撤出正式案例；相关文件仅作历史资料。外部数据不因放入本项目即改为MIT授权，需保留原包及官方来源声明。

第8章没有排放变量，只计算普通SBM。第10章为加法成本、半正态非效率教学模型，有2个观测发生效率截断，诊断另存。第13章输出标准化边际替代率，未换算为元/吨碳价。

## 数据存放位置

当前正文案例统一从本目录的 `data/` 读取数据，无需访问旧 `deabook/data/` 目录。正式使用的文件及对应章节见 `data/README.md`。`deabook.data_path()` 同样指向此目录。旧数据归档于 `research/legacy-data-20260909/data/`。

## 安装及输出目录

需要 Python 3.9 或更高版本。在源码根目录执行 `python -m pip install -e ".[book,test]"` 可安装案例与测试依赖；正式计算不再自动加载用户目录中的依赖缓存。

运行全部章节：`python -m deabook.example.book_cases.run`。指定章节和输出目录：

```bash
python -m deabook.example.book_cases.run 4 7 --output-dir ./book-output
python -m pytest tests -q
```

源码环境默认沿用正文图表位置；安装包环境默认写入当前目录的 `book-case-output/results` 和 `book-case-output/images`。仅导入模块不会创建目录。Notebook可先 `from deabook.example.book_cases import run`，调用 `run.configure_output("./book-output")`，再执行 `run.FUN[章节号]()`，通过 `run.OUT` 读取结果。

动态计算支持任意行索引，拒绝重复个体—年份；新进入且缺少前期观测的个体保留缺失指数。同一计算中的距离按观测和参照期缓存，不改变样本或模型定义。未采用的候选脚本和数据归档于 `research/legacy-case-candidates-20260910/`。
