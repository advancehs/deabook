# chap4：能源效率（静态）与能源生产率（动态）

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap4`
- 主要代码文件：
  - `chap 4静态.ipynb`
  - `chap 4静态_province.ipynb`
  - `chap 4动态.ipynb`
  - `chap 4动态 province.ipynb`

## 2. 本章主要模型 / 函数

- DEA.DDF2、MQDEA.MQDEA

## 3. 本章使用的数据

仍使用 `china data.xlsx` 与 `oecd data.xlsx`。

## 4. 本章 import 的逐项含义

- 静态 notebook 主要导入 `DEA`。
- 动态 notebook 主要导入 `MQDEA`。
- `CONTEMPORARY`：动态部分最常用的技术集。
- `CategoricalDtype`：控制动态指标展示顺序。

## 5. 本章参数穷尽说明

### 静态部分：`DEA.DDF2(...)`

与 chap2 的 `DDF2` 参数集合完全相同，区别只在方向向量：

- `sent="K L E=Y"`
- `gx=[0,0,1]`
- `gy=[0]`
- `rts=RTS_VRS1`

这代表：**只在第 3 个投入变量 `E`（能源）上做方向调整**。

### 动态部分：`MQDEA.MQDEA(...)`

参数集合与 chap3 一致，但实际取值是：

- `id=dmuname`
- `year='year'`
- `sent="K L E=Y"`
- `gx=[0,0,1]`
- `gy=[0]`
- `rts=RTS_VRS1`
- `tech=CONTEMPORARY`
- `solver="mosek"`

#### 这组参数会触发的分支

- `gx` 非零、`gy` 为零 → 输入方向动态指数
- `tech=CONTEMPORARY` → 走同期技术集比较
- 返回主列：`MQ`, `MEFFCH`, `MTECHCH`

## 6. 本章输出文件

- 静态：`table4.1.md`、`DDF_VRS.md`、`chap41分...png`
- 动态：`MQDEA_VRS_.md`、`chap42分...png`

## 7. 阅读与运行提醒

- `gx=[0,0,1]` 是本章最关键的代码信号。它决定这章是在算“能源效率”，不是资本效率。
- 静态和动态虽然都围绕能源，但底层类完全不同：静态用 `DDF2`，动态用 `MQDEA`。
