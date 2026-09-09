z = china(scaled=True).query('t==4').copy()
assert cross_validate_with_deabook(z, ['K','L'], ['Y'])
assert validate_efficiency_results(sbm('Y=K L',z))['checks_passed']
assert validate_efficiency_results(ddf('Y:CO2=K L',z),'ddf')['checks_passed']
assert validate_efficiency_results(nddf('Y:CO2=K L',z),'nddf')['checks_passed']
assert SBMModel('Y=K L',z).fit().result['obj'].notna().all()
assert sensitivity_analysis_drop_dmu('Y=K L',z.head(6)).shape[0]==6
np.random.seed(42)
assert sensitivity_analysis_perturbation('Y=K L',z.head(6),n_replicates=2).shape[0]==12
for technology in ('com','seq','global','window 1'):
 for orientation in ('io','oo'):
  r=mpi('Y=K L',china(scaled=True),'id','t',orient=orientation,tech=technology)
  assert r.loc[r.t>1,'mpi'].notna().all()
  if 'effch' in r: assert np.allclose(r.mpi.dropna(),(r.effch*r.techch).dropna())
for technology in ('com','seq','global','window 1'):
 r=mlpi('Y:CO2=K L',china(scaled=True),'id','t',tech=technology)
 assert r.loc[r.t==1,'mlpi'].isna().all()
 assert r.mlpi.notna().any()
 if 'effch' in r: assert np.allclose(r.mlpi.dropna(),(r.effch*r.techch).dropna())
print('Chapter 12 additional checks passed')
