"""Recalculate book exhibits: python run.py [chapter ...]."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[3]
# Support direct execution from a source checkout without changing import priority.
if __package__ in (None, ''):
    sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy import stats
from deabook.example.book_cases.models import radial,bad_distance,sbm,material_balance
from deabook.example.book_cases.loaders import read_raw
BASE=Path(__file__).resolve().parent
# Importing the module must not write to the source or installation directory.
_SOURCE_CHECKOUT = (ROOT/'docs/deabook文档').is_dir()
OUT = BASE/'results' if _SOURCE_CHECKOUT else Path.cwd()/'book-case-output/results'
IMG = ROOT/'docs/deabook文档/images' if _SOURCE_CHECKOUT else Path.cwd()/'book-case-output/images'


def configure_output(output_dir):
    """Select writable output directories; directories are created only on save."""
    global OUT, IMG
    root = Path(output_dir).expanduser().resolve()
    OUT, IMG = root/'results', root/'images'
    return OUT, IMG

for font in [Path('C:/Windows/Fonts/msyh.ttc'),Path('/mnt/c/Windows/Fonts/msyh.ttc')]:
    if font.exists():
        font_manager.fontManager.addfont(str(font));plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font)).get_name();break
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'figure.dpi':130,'savefig.dpi':300,'pdf.fonttype':42})
BLUE='#24557A';ORANGE='#B65A2B';GREEN='#41785B'

def data(name):
    if name=='plants':return read_raw('plants')
    if name=='firms':return read_raw('firms')
    if name in ['china','oecd']:
        d=read_raw(name).query('2016<=year<=2019').reset_index(drop=True)
        d['id']=d['省份' if name=='china' else '国家'];return d
    if name=='power':
        d=read_raw('power');d.columns=['province','K','L','F','Y','CO2','year','id'];return d.query('2015<=year<=2017').reset_index(drop=True)
    raise ValueError(name)

def arrays(d,x,y,b=None):
    cols=x+y+(b or []);v=d[cols].to_numpy(float)
    if not np.isfinite(v).all() or (v<=0).any():raise ValueError('Nonpositive or missing production data')
    sc=d[cols].mean();z=d[cols]/sc
    return z[x].to_numpy(),z[y].to_numpy(),None if b is None else z[b].to_numpy(),sc

def save(ch,df,summary,fig,caption):
    from deabook.example.book_cases.exhibits import expand_table
    summary = expand_table(ch, df, summary)
    OUT.mkdir(parents=True, exist_ok=True)
    IMG.mkdir(parents=True, exist_ok=True)
    axes = fig.axes
    fig.set_size_inches(7.2, 3.8 if len(axes)>1 else 4.2)
    for i, ax in enumerate(axes):
        if not ax.get_title():
            ax.set_title(f"（{chr(97+i)}）" + ax.get_ylabel(), fontsize=10, pad=12)
        ax.tick_params(labelsize=9)
        ax.grid(axis='y', alpha=.16)
    fig.tight_layout(pad=1.4, w_pad=2.4)

    df.to_csv(OUT/f'ch{ch:02}_full.csv',index=False)
    summary.to_csv(OUT/f'ch{ch:02}_summary.csv',index=False)
    (OUT/f'ch{ch:02}_summary.md').write_text(summary.to_markdown(index=False,floatfmt='.3f'),encoding='utf-8')
    fig.savefig(IMG/f'case_ch{ch:02}.png',bbox_inches='tight')
    fig.savefig(IMG/f'case_ch{ch:02}.pdf',bbox_inches='tight')
    plt.close(fig)
    (OUT/f'ch{ch:02}_meta.json').write_text(json.dumps({'caption':caption,'rows':len(df),'table_rows':len(summary),'columns':list(df.columns)},ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'CH{ch:02} DONE rows={len(df)} table={len(summary)}',flush=True)

def summ(df,cols,labels=None):
    return pd.DataFrame([{'指标':(labels or {}).get(c,c),'均值':df[c].mean(),'中位数':df[c].median(),'P10':df[c].quantile(.1),'P90':df[c].quantile(.9)} for c in cols])

def boxes(ax,df,cols,labels=None):
    ax.boxplot([df[c].dropna() for c in cols],tick_labels=labels or cols,showfliers=True,patch_artist=True,boxprops={'facecolor':'#DCE6EE'},medianprops={'color':ORANGE})
    ax.set_ylabel('效率');ax.set_ylim(0,1.05);ax.grid(axis='y',alpha=.2)

def ch3():
    d=data('plants');X,Y,_,_=arrays(d,['Labor','Fuel','Capital'],['Output']);out=d.copy()
    for rt in ['crs','vrs']:
        for o in ['input','output']:out[f'{rt}_{o}']=[radial(X,Y,x,y,rt,o) for x,y in zip(X,Y)]
    out['scale']=out.crs_input/out.vrs_input
    labels={'crs_input':'CRS投入效率','vrs_input':'VRS投入效率','crs_output':'CRS产出效率','vrs_output':'VRS产出效率','scale':'规模效率'}
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.7));boxes(ax[0],out,['crs_input','vrs_input','scale'],['CRS投入','VRS投入','规模效率'])
    ax[1].scatter(out.crs_input,out.vrs_input,color=BLUE,s=24);ax[1].plot([.75,1],[.75,1],'--',color='gray',lw=.8);ax[1].set(xlabel='CRS投入效率',ylabel='VRS投入效率');fig.tight_layout()
    save(3,out,summ(out,list(labels),labels),fig,'美国19家电厂的效率分布与规模报酬比较')

def dynamic(d, xcols, ycols, mode='output', bcols=None, sequential=False):
    """Adjacent-period indices, invariant to DataFrame labels and row order.

    Unmatched DMUs retain a row with missing distances and indices. Duplicate
    DMU-years are rejected rather than choosing one silently.
    """
    d = d.reset_index(drop=True)
    if d[['id', 'year']].isna().any().any():
        raise ValueError('Panel id and year must not be missing')
    if d.duplicated(['id', 'year']).any():
        raise ValueError('Duplicate id-year observations')
    if mode not in ('output', 'energy'):
        raise ValueError("mode must be 'output' or 'energy'")
    X, Y, B, _ = arrays(d, xcols, ycols, bcols)
    years = sorted(d.year.unique())
    refs = {t: (d.year <= t).to_numpy() if sequential else (d.year == t).to_numpy()
            for t in years}
    cache = {}

    def dist(i, year):
        key = (i, year)
        if key not in cache:
            ref = refs[year]
            if bcols:
                value = bad_distance(X[ref], Y[ref], B[ref], X[i], Y[i], B[i],
                                     allow_infeasible=True)[0]
            else:
                value = radial(X[ref], Y[ref], X[i], Y[i],
                               orient='input' if mode == 'energy' else 'output',
                               energy=0 if mode == 'energy' else None,
                               allow_infeasible=True)
            cache[key] = value
        return cache[key]

    rows = []
    for t0, t1 in zip(years[:-1], years[1:]):
        prev = {v: i for i, v in d.loc[d.year == t0, 'id'].items()}
        for i in d.index[d.year == t1]:
            p = prev.get(d.loc[i, 'id'])
            a, bb, c, e = ((dist(p, t0), dist(i, t0), dist(p, t1), dist(i, t1))
                           if p is not None else (np.nan,) * 4)
            if np.isfinite([a, bb, c, e]).all() and min(a, bb, c, e) > 0:
                mq = np.sqrt(bb / a * e / c)
                ec = e / a
                tc = mq / ec
            else:
                mq = ec = tc = np.nan
            rows.append({'id': d.loc[i, 'id'], 'year': t1,
                         'D00': a, 'D10': bb, 'D01': c, 'D11': e,
                         'MQ': mq, 'EC': ec, 'TC': tc})
    return pd.DataFrame(rows, columns=['id', 'year', 'D00', 'D10', 'D01', 'D11',
                                       'MQ', 'EC', 'TC'])

def geo(s):
    a=np.asarray(s);a=a[np.isfinite(a)&(a>0)];return np.exp(np.log(a).mean()) if len(a) else np.nan

def dynplot(ch,r,title):
    valid=r.dropna(subset=['MQ','EC','TC']);g=valid.groupby('year')[['MQ','EC','TC']].agg(geo)
    table=g.reset_index();table['有效样本']=valid.groupby('year').size().values
    if len(table)>6:
        table=pd.DataFrame([{'时期':f'{int(g.index.min())}—{int(g.index.max())}','MQ':geo(valid.MQ),'EC':geo(valid.EC),'TC':geo(valid.TC),'有效观测':len(valid)}])
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.6));g.rename(columns={'MQ':'生产率','EC':'效率变化','TC':'技术变化'}).plot(ax=ax[0],marker='o',color=[BLUE,ORANGE,GREEN]);ax[0].axhline(1,c='gray',lw=.7);ax[0].set(xlabel='年份',ylabel='几何平均指数');ax[0].legend(fontsize=7)
    ax[1].hist(valid.MQ,bins=14,color=BLUE,alpha=.85);ax[1].axvline(1,c=ORANGE,lw=1);ax[1].set(xlabel='单个观测的生产率指数',ylabel='观测数');fig.tight_layout();save(ch,r,table,fig,title)

def ch4():dynplot(4,dynamic(data('power'),['K','L','F'],['Y']),'中国火电行业发电生产率的年度分解与观测分布')

def ch5():
    d=data('china');r=dynamic(d,['E','K','L'],['Y'],mode='energy');X,Y,_,_=arrays(d,['E','K','L'],['Y']);static=np.array([radial(X,Y,x,y,energy=0) for x,y in zip(X,Y)])
    d['EE']=static;r=r.merge(d[['id','year','EE']],on=['id','year'],how='right');g=r.dropna(subset=['MQ','EC','TC']).groupby('year')[['MQ','EC','TC']].agg(geo)
    tab=d.groupby('year').EE.mean().rename('能源效率').to_frame().join(g).reset_index();fig,ax=plt.subplots(1,2,figsize=(7.2,2.6));ax[0].plot(tab.year,tab['能源效率'],marker='o',c=BLUE);ax[0].set(xlabel='年份',ylabel='相对混合前沿的平均能源效率');g.plot(ax=ax[1],marker='o',color=[BLUE,ORANGE,GREEN]);ax[1].axhline(1,c='gray',lw=.7);ax[1].set(xlabel='年份',ylabel='几何平均指数');fig.tight_layout();save(5,r,tab,fig,'中国省级能源效率与相邻期能源生产率')

def ch6():
    d=data('oecd').query('year==2019').reset_index(drop=True);X,Y,B,_=arrays(d,['K','L','E'],['Y'],['CO2']);out=d[['id','year']].copy()
    for key,rt,strong in [('弱处置CRS','crs',False),('弱处置VRS共同因子','vrs',False),('强处置CRS','crs',True)]:out[key]=[bad_distance(X,Y,B,x,y,b,rt,strong)[0] for x,y,b in zip(X,Y,B)]
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.7));boxes(ax[0],out,list(out.columns[2:]),['弱CRS','弱VRS','强CRS']);ax[1].scatter(out.iloc[:,2],out.iloc[:,3],c=BLUE,s=24);ax[1].plot([0,1],[0,1],'--',c='gray',lw=.7);ax[1].set(xlabel='弱处置CRS效率',ylabel='弱处置VRS共同因子效率');fig.tight_layout();save(6,out,summ(out,list(out.columns[2:])),fig,'2019年OECD碳排放方向效率的技术假设比较')

def ch7():dynplot(7,dynamic(data('oecd'),['K','L','E'],['Y'],bcols=['CO2'],sequential=True),'OECD动态环境绩效的年度分解与分布（时序CRS技术）')

def ch8():
    d=data('plants');X,Y,_,_=arrays(d,['Labor','Fuel','Capital'],['Output']);out=d.copy()
    for rt in ['crs','vrs']:
        res=[sbm(X,Y,x,y,rt) for x,y in zip(X,Y)]
        out['SBM-'+rt.upper()]=[v[0] for v in res]
        if rt=='vrs':
            sl=np.array([v[1] for v in res])
            for k,name in enumerate(['劳动','燃料','装机']):out[name+'松弛比例']=sl[:,k]/X[:,k]
            out['发电松弛比例']=sl[:,3]/Y[:,0]
    cols=['SBM-CRS','SBM-VRS'];fig,ax=plt.subplots(1,2,figsize=(7.2,2.7));boxes(ax[0],out,cols)
    ax[1].bar(['劳动','燃料','装机','发电'],out[['劳动松弛比例','燃料松弛比例','装机松弛比例','发电松弛比例']].median()*100,color=[BLUE,BLUE,BLUE,ORANGE]);ax[1].set(ylabel='相对松弛中位数（%）');fig.tight_layout()
    save(8,out,summ(out,cols),fig,'美国19家电厂的SBM效率与分项松弛')

def ch9():
    d=data('power');cols=['K','L','F','Y','CO2'];sc=d[cols].mean();z=d[cols]/sc;X=z[['K','L','F']].to_numpy();Y=z[['Y']].to_numpy();B=z[['CO2']].to_numpy();coef=(d.CO2.sum()/d.F.sum())*sc.F/sc.CO2
    out=d[['id','year']].copy()
    for level in range(1,5):out[f'B{level}']=[material_balance(X,Y,B,x,y,b,coef,level) for x,y,b in zip(X,Y,B)]
    out['B4_keep']=[material_balance(X,Y,B,x,y,b,coef,4,True) for x,y,b in zip(X,Y,B)]
    out['TE']=out.B1/B[:,0];out['PIAE']=out.B2/out.B1;out['NPIAE']=out.B3/out.B2;out['NPOAE']=out.B4/out.B3;out['保产出效率']=out.B4_keep/B[:,0]
    assert (out.B1+1e-7>=out.B2).all() and (out.B2+1e-7>=out.B3).all()
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.6));boxes(ax[0],out,['TE','NPIAE','NPOAE'],['TE','NPIAE','允许停产配置项']);ax[1].scatter(out.B4/B[:,0],out['保产出效率'],s=14,c=BLUE,alpha=.7);ax[1].set(xlabel='允许产出归零的整体比率',ylabel='增加产出下限后的比率');fig.tight_layout();save(9,out,summ(out,['TE','PIAE','NPIAE','NPOAE','保产出效率']),fig,'物质平衡分解与产出下限敏感性（2015—2017年）')

def ch10():
    import cvxpy as cp
    from scipy.optimize import minimize
    d=data('firms');n=len(d);X=d[['Energy','Length','Customers']].to_numpy();X=X/X.mean(axis=0);c=d.TOTEX.to_numpy()/d.TOTEX.mean()
    # Convex monotone cost CNLS: f_i >= f_j + beta_j'(x_i-x_j).
    f=cp.Variable(n);beta=cp.Variable((n,3),nonneg=True)
    cons=[f>=1e-8]
    for j in range(n):cons.append(f>=f[j]+(X-X[j])@beta[j,:])
    prob=cp.Problem(cp.Minimize(cp.sum_squares(c-f)),cons);prob.solve(solver='CLARABEL',tol_gap_abs=1e-8,tol_feas=1e-8,max_iter=300)
    if prob.status!='optimal':raise RuntimeError(prob.status)
    fitted=f.value;res=c-fitted
    # Centered normal/half-normal cost residual density; mu=E[u].
    def nll(par):
        su,sv=np.exp(par);sig=np.hypot(su,sv);e=res+su*np.sqrt(2/np.pi)
        return -np.sum(np.log(2)-np.log(sig)+stats.norm.logpdf(e/sig)+stats.norm.logcdf(e*su/(sig*sv)))
    sd=res.std();fits=[minimize(nll,np.log([sd*a,sd*(1-a)]),bounds=[(-14,4),(-14,4)],method='L-BFGS-B') for a in [.2,.5,.8]];fit=min(fits,key=lambda q:q.fun)
    if not fit.success:raise RuntimeError(fit.message)
    su,sv=np.exp(fit.x);e=res+su*np.sqrt(2/np.pi);sig2=su*su+sv*sv;mu=e*su*su/sig2;ss=su*sv/np.sqrt(sig2);zz=mu/ss;uhat=mu+ss*np.exp(stats.norm.logpdf(zz)-stats.norm.logcdf(zz))
    out=d.copy();out['row_id']=np.arange(1,n+1);out['CNLS_fitted']=fitted;out['residual']=res;out['JLMS_u']=uhat
    # Additive cost efficiency conditional on the observation; report clipping.
    raw=(c-uhat)/c;out['StoNED成本效率']=np.clip(raw,0,1);out['DEA成本效率']=[radial(c[:,None],X,[ci],xi,'vrs') for ci,xi in zip(c,X)]
    meta={'solver_status':prob.status,'su':su,'sv':sv,'residual_skew':float(stats.skew(res)),'cost_efficiency_clipped':int(((raw<0)|(raw>1)).sum()),'afriat_max_violation':float(max(np.max(fitted[j]+(X-X[j])@beta.value[j]-fitted) for j in range(n)))}
    (OUT/'ch10_diagnostics.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.6));ax[0].scatter(c,fitted,s=16,c=BLUE);lim=max(c);ax[0].plot([0,lim],[0,lim],'--',c='gray',lw=.7);ax[0].set(xlabel='观测总成本／样本均值',ylabel='CNLS拟合成本／样本均值');ax[1].scatter(out['DEA成本效率'],out['StoNED成本效率'],s=16,c=ORANGE);ax[1].plot([0,1],[0,1],'--',c='gray',lw=.7);ax[1].set(xlabel='DEA成本效率',ylabel='StoNED成本效率');fig.tight_layout();save(10,out,summ(out,['DEA成本效率','StoNED成本效率']),fig,'芬兰89家配电企业的成本拟合与效率比较')

def ch12():
    d=data('china').query('year==2019').reset_index(drop=True);X,Y,B,_=arrays(d,['K','L','E'],['Y'],['CO2']);out=d[['id','year']].copy()
    out['普通SBM']=[sbm(X,Y,x,y,'crs')[0] for x,y in zip(X,Y)]
    out['含污染SBM']=[sbm(X,Y,x,y,'crs',B,b)[0] for x,y,b in zip(X,Y,B)]
    out['减排方向效率']=[bad_distance(X,Y,B,x,y,b)[0] for x,y,b in zip(X,Y,B)]
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.6));boxes(ax[0],out,list(out.columns[2:]));ax[1].scatter(out['普通SBM'],out['含污染SBM'],c=BLUE,s=20);ax[1].plot([0,1],[0,1],'--',c='gray',lw=.7);ax[1].set(xlabel='普通SBM',ylabel='含污染SBM');fig.tight_layout();save(12,out,summ(out,list(out.columns[2:])),fig,'同一中国省级样本中不同目标函数的结果比较（2019年）')

def ch13():
    d=data('china').query('year==2019').reset_index(drop=True);X,Y,B,sc=arrays(d,['K','L','E'],['Y'],['CO2']);out=d[['id','year']].copy();rows=[]
    for x,y,b in zip(X,Y,B):
        eff,r=bad_distance(X,Y,B,x,y,b);cost_y=-r.ineqlin.marginals[3];cost_b=-r.eqlin.marginals[0]
        # Magnitudes from minimized -DDF; report only a normalized tradeoff.
        trade=abs(cost_b/cost_y) if abs(cost_y)>1e-8 else np.nan
        # Finite differences check value-function shadow prices on the scaled LP.
        h=1e-5;b2=b.copy();b2[0]+=h
        # Direction depends on b as well, so hold the direction fixed in the audit.
        E=np.column_stack([B.T,b]);A=np.vstack([np.column_stack([X.T,np.zeros(3)]),np.column_stack([-Y.T,np.zeros(1)])])
        from deabook.example.book_cases.models import solve
        rr=solve(np.r_[np.zeros(len(X)),-1],A,np.r_[x,-y],E,b2,bounds=[(0,None)]*len(X)+[(None,1)])
        rows.append({'效率':eff,'产出约束乘数':cost_y,'排放约束乘数':cost_b,'标准化边际替代率':trade,'对偶有限差分误差':abs((rr.fun-r.fun)/h-r.eqlin.marginals[0])})
    out=pd.concat([out,pd.DataFrame(rows)],axis=1);v=out['标准化边际替代率'].dropna()
    tab=pd.DataFrame([{'统计量':k,'数值':val} for k,val in [('总样本',len(out)),('可归一化样本',len(v)),('边际替代率中位数',v.median()),('边际替代率P90',v.quantile(.9)),('最大有限差分误差',out['对偶有限差分误差'].max())]])
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.6));ax[0].hist(v,bins=np.linspace(float(v.min())-0.05,float(v.max())+0.05,11),color=BLUE);ax[0].set(xlabel='标准化产出／排放边际替代率',ylabel='省份数');ax[1].scatter(out['效率'],out['标准化边际替代率'],c=ORANGE,s=20);ax[1].set(xlabel='减排方向效率',ylabel='标准化边际替代率');fig.tight_layout();save(13,out,tab,fig,'中国省级样本的排放约束边际价值与效率（2019年）')

def introductory(ch):
    d=data('plants');cols=['Labor','Fuel','Capital','Output'];out=d.copy()
    if ch==1:
        out['发电燃料比']=d.Output/d.Fuel;tab=summ(out,['发电燃料比']);fig,axs=plt.subplots(1,2);ax=axs[0];ax.scatter(d.Fuel,d.Output,c=BLUE);ax.set(xlabel='燃料（10^10 Btu）',ylabel='发电量（10^6 kWh）');axs[1].bar(range(1,len(d)+1),out['发电燃料比'],color=ORANGE);axs[1].set(xlabel='电厂原始序号',ylabel='发电量／燃料');cap='美国电厂的燃料与发电量：为何需要多投入比较'
    elif ch==2:
        tab=pd.DataFrame({'变量':cols,'均值':d[cols].mean().values,'最小值':d[cols].min().values,'最大值':d[cols].max().values});fig,axs=plt.subplots(1,2);ax=axs[0];corr=d[cols].corr();im=axs[1].imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r');axs[1].set(xticks=range(4),yticks=range(4),xticklabels=['劳动','燃料','装机','发电'],yticklabels=['劳动','燃料','装机','发电'],ylabel='变量相关系数');
        for i in range(4):
            for j in range(4):axs[1].text(j,i,f'{corr.iloc[i,j]:.2f}',ha='center',va='center',color='white' if abs(corr.iloc[i,j])>.6 else 'black',fontsize=9)
        ax.boxplot([d[c]/d[c].median() for c in cols],tick_labels=['劳动','燃料','装机','发电'],showfliers=True);ax.set(ylabel='变量／各列中位数（仅用于绘图）');cap='读取数据后比较四个变量的相对离散程度'
    else:
        X,Y,_,_=arrays(d,['Labor','Fuel','Capital'],['Output']);out['效率']=[radial(X,Y,x,y,'vrs','output') for x,y in zip(X,Y)];tab=out[['Plant','效率']].head(5).rename(columns={'Plant':'原文件前五家电厂'});fig,axs=plt.subplots(1,2);ax=axs[0];ax.plot(np.arange(1,len(out)+1),out['效率'],'o-',c=BLUE,lw=.8);ax.set(xlabel='原文件行号',ylabel='VRS产出效率');axs[1].bar(np.arange(1,len(out)+1),(1/out['效率']-1)*100,color=ORANGE);axs[1].set(xlabel='原文件行号',ylabel='径向产出扩张空间（%）');cap='从单个电厂求解扩展到19家电厂循环'
    fig.tight_layout();save(ch,out,tab,fig,cap)

FUN={1:lambda:introductory(1),2:lambda:introductory(2),3:ch3,4:ch4,5:ch5,6:ch6,7:ch7,8:ch8,9:ch9,10:ch10,11:lambda:introductory(11),12:ch12,13:ch13}
def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('chapters', nargs='*', type=int)
    parser.add_argument('--output-dir', type=Path, help='Write results/ and images/ here')
    args = parser.parse_args(argv)
    chapters = args.chapters or list(FUN)
    if any(ch not in FUN for ch in chapters):
        parser.error('chapter numbers must be between 1 and 13')
    if args.output_dir is not None:
        configure_output(args.output_dir)
    for ch in chapters:
        FUN[ch]()


if __name__ == '__main__':
    main()
