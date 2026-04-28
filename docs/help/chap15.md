# 第十五章：DEA 调试与排错

**主题**：通过在 dea8 函数中故意引入 bug，学习 DEA 模型的调试技术。

**数据文件**：
- `Ex3.dta`

**说明**：本章不使用 deabook 包，使用第十二章开发的 `dea8` 函数及其变体。

---

## 导入

```python
import pandas as pd
from pyomo.environ import ConcreteModel, Set, Var, Objective, Constraint
```

---

## 调试场景 1：变量名错误

### Bug 描述

在 LP 约束中使用了错误的循环变量名：

```python
# 正确代码
for i in model.I:
    ...

# Bug 代码
for i2 in model.I:
    ...  # 但后续使用 i（未定义）
```

### 错误现象

```
NameError: name 'i' is not defined
```

### 诊断方法

1. 检查约束函数内的变量名是否与循环变量一致
2. 确认 Pyomo 约束的 rule 函数签名正确

---

## 调试场景 2：缺失参数导致约束缺失

### Bug 描述

`dea8` 的 `orient` 参数缺失，导致没有匹配的 if/elif 分支：

```python
# dea8 缺少 orient 参数
def dea8(formula, data, rts="crs"):
    # 没有 orient 参数 → 无产出/投入导向分支
    # Constraint 函数返回 None
```

### 错误现象

- 约束未正确构建
- LP 模型可能无解或结果异常

### 诊断方法

1. 使用 `model.pprint()` 查看 LP 模型的完整结构
2. 检查约束是否被正确添加
3. 确认 `orient` 参数传递正确

---

## 调试工具和技术

### 1. model.pprint()

打印完整的 LP 模型结构（变量、约束、目标函数）：

```python
model.pprint()
```

输出包括：
- 所有变量及其 bounds
- 所有约束及表达式
- 目标函数

### 2. print 语句检查数据

在构建约束前打印参考集数据：

```python
print(xref)  # 检查参考集投入数据
print(yref)  # 检查参考集产出数据
```

### 3. 使用 glpk 求解器调试

```python
res = dea8("Y~K L", data, rts="crs", orient="oo", solver="glpk")
```

- `glpk` 是开源求解器，适合调试
- 不依赖商业求解器许可证
- 速度较慢但错误信息更详细

### 4. 逐步构建

从最简版本（dea1）开始，逐步增加功能：
1. 先确保基础 CRS + 产出导向正确
2. 添加 VRS 约束
3. 添加投入导向
4. 添加 evaquery/refquery

---

## 常见错误类型

| 错误类型 | 错误信息 | 原因 | 解决方案 |
|---------|---------|------|---------|
| 变量名错误 | `NameError: name 'i' is not defined` | 循环变量与使用变量不一致 | 检查所有变量名 |
| 约束缺失 | 结果异常或无解 | if/elif 分支不完整 | 添加所有必要的分支 |
| 索引越界 | `IndexError` | 方向向量长度与变量数不匹配 | 检查 gx/gy 长度 |
| 求解失败 | `optimization_status != ok` | 模型不可行或无界 | 检查约束符号和 bounds |

---

## dea8 正确实现要点

```python
def dea8(formula, data, evaquery=None, refquery=None, rts="crs", orient="oo"):
```

关键检查点：
1. `formula` 解析正确（`~` 分隔产出和投入）
2. `orient` 为 `"oo"` 或 `"io"`
3. `rts` 为 `"crs"` 或 `"vrs"`
4. 产出导向：目标函数为 maximize，约束符号正确
5. 投入导向：目标函数为 minimize，约束符号正确
6. VRS：添加 `sum(lamda) == 1` 约束
7. evaquery/refquery 正确筛选 DMU
