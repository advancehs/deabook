from run import *
from deabook import DEA,DEAweak
from deabook.constant import RTS_CRS,RTS_VRS1,OPT_LOCAL
checks=[]
d=data('plants');X,Y,_,_=arrays(d,['Labor','Fuel','Capital'],['Output']);z=pd.DataFrame(np.c_[X,Y],columns=['L','F','K','Y'])
for rt in [RTS_CRS,RTS_VRS1]:
 r=DEA.ccr(z,'L F K=Y').optimize(OPT_LOCAL,'highs') if rt==RTS_CRS else DEA.bcc(z,'L F K=Y').optimize(OPT_LOCAL,'highs')
 ref=pd.read_csv(OUT/'ch03_full.csv'); key='crs_output' if rt==RTS_CRS else 'vrs_output'; err=float(np.max(np.abs(r.te.to_numpy()-ref[key])));checks.append({'model':'radial '+rt,'max_error':err});assert err<1e-6
d=data('plants');X,Y,_,_=arrays(d,['Labor','Fuel','Capital'],['Output']);z=pd.DataFrame(np.c_[X,Y],columns=['L','F','K','Y'])
r=DEA.sbm(z,'L F K=Y',gx=[1,1,1],gy=[1],rts=RTS_VRS1).optimize(OPT_LOCAL,'highs')
ref=pd.read_csv(OUT/'ch08_full.csv');err=float(np.max(np.abs(r.rho.to_numpy()-ref['SBM-VRS'])));checks.append({'model':'SBM VRS','max_error':err});assert err<1e-6
d=data('oecd').query('year==2019').reset_index(drop=True);X,Y,B,_=arrays(d,['K','L','E'],['Y'],['CO2']);z=pd.DataFrame(np.c_[X,Y,B],columns=['K','L','E','Y','CO2'])
for rt,label in [(RTS_CRS,'弱处置CRS'),(RTS_VRS1,'弱处置VRS共同因子')]:
 r=DEAweak.ddf(z,'K L E=Y:CO2',gx=[0,0,0],gy=[0],gb=[1],rts=rt).optimize(OPT_LOCAL,'highs');ref=pd.read_csv(OUT/'ch06_full.csv');err=float(np.max(np.abs(1-r.rho.to_numpy()-ref[label])));checks.append({'model':'DDF '+rt,'max_error':err});assert err<1e-6
(OUT/'interface_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
