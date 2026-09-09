"""Book tables with sufficient detail for explanation without full sample listings."""
import numpy as np
import pandas as pd

def expand_table(ch, df, old):
    if ch in (1,11):
        # Fixed original order avoids selecting observations for their outcomes.
        out=df.iloc[:8][['Plant','Labor','Fuel','Capital','Output']].copy()
        out.columns=['电厂','劳动','燃料','装机','发电量']
        if ch==1:out['发电燃料比']=df['发电燃料比'].iloc[:8].values
        if ch==11:out['VRS产出效率']=df['效率'].iloc[:8].values
        return out
    if ch in (4,5,7):
        rows=[]
        for year,g in df.groupby('year'):
            for col,label in [('MQ','生产率指数'),('EC','效率变化'),('TC','技术变化')]:
                v=g[col].dropna();v=v[v>0]
                if v.empty:continue
                rows.append({'期末年份':int(year),'指标':label,'有效数':len(v),
                             '几何均值':np.exp(np.log(v).mean()),'中位数':v.median(),
                             'P10':v.quantile(.1),'P90':v.quantile(.9)})
        return pd.DataFrame(rows)
    columns={2:['Labor','Fuel','Capital','Output'],3:['crs_input','vrs_input','crs_output','vrs_output','scale'],
             6:['弱处置CRS','弱处置VRS共同因子','强处置CRS'],
             8:['SBM-CRS','SBM-VRS','劳动松弛比例','燃料松弛比例','装机松弛比例','发电松弛比例'],
             9:['TE','PIAE','NPIAE','NPOAE','保产出效率'],
             10:['DEA成本效率','StoNED成本效率'],
             12:['普通SBM','含污染SBM','减排方向效率'],
             13:['效率','产出约束乘数','排放约束乘数','标准化边际替代率']}[ch]
    labels={'Labor':'劳动','Fuel':'燃料','Capital':'装机','Output':'发电量','crs_input':'CRS投入','vrs_input':'VRS投入','crs_output':'CRS产出','vrs_output':'VRS产出','scale':'规模效率'}
    stat=[('有效样本数',lambda v:len(v)),('均值',lambda v:v.mean()),('标准差',lambda v:v.std()),
          ('最小值',lambda v:v.min()),('P10',lambda v:v.quantile(.1)),('中位数',lambda v:v.median()),
          ('P90',lambda v:v.quantile(.9)),('最大值',lambda v:v.max())]
    rows=[]
    for name,fn in stat:
        row={'统计量':name}
        for col in columns:
            v=df[col].replace([np.inf,-np.inf],np.nan).dropna()
            row[labels.get(col,col)]=fn(v)
        rows.append(row)
    return pd.DataFrame(rows)
