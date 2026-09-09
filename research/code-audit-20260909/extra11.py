z=plants()
for rts in ('crs','vrs'):
 for orientation in ('io','oo'):
  a=dea8('Output=Labor Fuel Capital',z,rts,orientation,None,None)
  b=DEAModel('Output=Labor Fuel Capital',z,orient=orientation,rts=rts).solve()
  assert np.allclose(a.te,b.te)
# Exercise the checked solver with an actual plant model.
instance=DEAModel('Output=Labor Fuel Capital',z)
d,ref=instance._prepare_data()
xcols=['Labor','Fuel','Capital'];ycols=['Output']
lp=instance._build_model(d.loc[d.index[0],xcols],d.loc[d.index[0],ycols],ref[xcols],ref[ycols])
solve_with_check(lp)
print('Chapter 11 additional checks passed')
