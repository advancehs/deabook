"""One empirical dataset per chapter; legacy aliases do not change observations."""
import pandas as pd
from deabook import data_path

DATA_FILES = {
    'plants': 'us_electric_plants.csv',
    'china': 'china data.xlsx',
    'oecd': 'oecd data.xlsx',
    'power': 'power plant data.xlsx',
    'firms': 'electricityFirms.csv',
}


def read_raw(name):
    """Read a fresh copy of one official dataset through the package locator."""
    filename = DATA_FILES[name]
    path = data_path(filename)
    return pd.read_csv(path) if filename.endswith('.csv') else pd.read_excel(path)


def plants():
    d=read_raw('plants')
    d['dmu']=range(1,len(d)+1)
    return d

def china(scaled=False):
    d=read_raw('china').query('2016<=year<=2019').copy()
    d['id']=pd.factorize(d['省份'],sort=True)[0]+1
    d['t']=d.year-2015
    if scaled:
        cols=['K','L','E','Y','CO2'];d[cols]=d[cols]/d[cols].mean()
    return d.sort_values(['t','id']).reset_index(drop=True)

def china_lower(scaled=True):
    return china(scaled).rename(columns={'K':'capital','L':'labor','E':'energy','Y':'gdp','CO2':'co2'})
