"""One empirical dataset per chapter; legacy aliases do not change observations."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent

def plants():
    d=pd.read_csv(BASE/'data/us_electric_plants.csv')
    d['dmu']=range(1,len(d)+1)
    return d

def china(scaled=False):
    d=pd.read_excel(ROOT/'deabook/data/china data.xlsx').query('2016<=year<=2019').copy()
    d['id']=pd.factorize(d['省份'],sort=True)[0]+1
    d['t']=d.year-2015
    if scaled:
        cols=['K','L','E','Y','CO2'];d[cols]=d[cols]/d[cols].mean()
    return d.sort_values(['t','id']).reset_index(drop=True)

def china_lower(scaled=True):
    return china(scaled).rename(columns={'K':'capital','L':'labor','E':'energy','Y':'gdp','CO2':'co2'})
