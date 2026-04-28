# 第四章：静态与动态能源效率

**主题**：使用 DDF 测量静态能源效率，使用 MQDEA 测量动态能源效率变动。

**数据文件**：
- `oecd data.xlsx`（33个OECD国家，2016-2019）
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP)

---

## 第一部分：静态能源效率（DDF）

### 导入

```python
from deabook import DEA
from deabook.constant import RTS_VRS1, RTS_CRS, OPT_LOCAL
```

### 模型：DDF 投入导向（能源效率）— VRS

```python
model = DEA.DDF2(data, sent="K L E=Y", gy=[0], gx=[0,0,1], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `sent` | `"K L E=Y"` | 3个投入(K,L,E)，1个产出(Y) |
| `gy` | `[0]` | 不调整产出（投入导向） |
| `gx` | `[0,0,1]` | 只沿第三个投入（E，能源）方向缩减 |
| `rts` | `RTS_VRS1` | 可变规模报酬 |

### 返回值

- `optimization_status`：求解状态
- `rho`：距离函数值
- `tei`：能源技术效率（`tei = 1 - rho`，范围 0-1）

### gx 方向向量的含义

`gx=[0,0,1]` 表示只对能源 E 做效率评估：
- `[0,0,1]`：能源效率（EO）
- `[1,0,0]`：资本效率（PEO）
- `[0,1,0]`：劳动效率
- `[1,1,1]`：全要素效率

---

## 第二部分：动态能源效率（MQDEA Malmquist）

### 导入

```python
from deabook import MQDEA
from deabook.constant import RTS_VRS1, RTS_CRS, OPT_LOCAL, TOTAL, CONTEMPORARY
```

### 模型：MQDEA 投入导向（能源）— VRS + 当期技术

```python
model = MQDEA.MQDEA(
    data, id=dmuname, year='year',
    sent="K L E=Y",
    gx=[0,0,1], gy=[0],
    rts=RTS_VRS1,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `gx` | `[0,0,1]` | 只沿能源方向缩减（能源 Malmquist） |
| `gy` | `[0]` | 不调整产出（投入导向） |
| `tech` | `CONTEMPORARY` | 当期生产技术 |

### 返回值

- `D11`：当期能源距离值
- `MQ`：Malmquist 指数
- `MEFFCH`：效率变化分量
- `MTECHCH`：技术变化分量

### 静态 vs 动态对比

| 特征 | 静态（DDF2） | 动态（MQDEA） |
|------|-------------|--------------|
| 时间维度 | 单截面 | 面板跨期 |
| 输入要求 | 仅需数据+sent | 需 data + id + year |
| 输出指标 | tei（单期效率） | MQ, MEFFCH, MTECHCH |
| 效率含义 | 相对当期前沿 | 跨期生产率变动 |
