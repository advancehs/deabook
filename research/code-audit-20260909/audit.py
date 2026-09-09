import sys,re,json,traceback,contextlib,io
from pathlib import Path
# Optional local dependency cache; ordinary pip environments need no override.
cache=Path('C:/Users/hs/.cache/deabook-book-runtime')
if cache.exists(): sys.path.insert(0,str(cache))
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
sys.stdout.reconfigure(encoding='utf-8')
result=[]
for ch in map(int,sys.argv[1:] or ['11','12','13']):
 p=next((ROOT/'docs/deabook文档').glob(f'{ch:02}_*.md'));s=p.read_text(encoding='utf-8');s=re.sub(r'^> ?', '', s, flags=re.M);ns={'__name__':'book_example'}
 for i,m in enumerate(re.finditer(r'^```python\s*\n(.*?)^```',s,re.M|re.S),1):
  code=m[1];line=s[:m.start()].count('\n')+1
  entry={'chapter':ch,'block':i,'line':line}
  try:
   with contextlib.redirect_stdout(io.StringIO()) as log,contextlib.redirect_stderr(log):exec(compile(code,f'{p.name}:{line}','exec'),ns)
   (ROOT/f'research/code-audit-20260909/ch{ch}_block{i:02}.log').write_text(log.getvalue(),encoding='utf-8')
   entry['status']='pass'
  except Exception as e:entry.update(status='fail',error=str(e),trace=traceback.format_exc(limit=2))
  result.append(entry);print(ch,i,entry['status'],entry.get('error','')[:150],flush=True)
  (ROOT/'research/code-audit-20260909/results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')

 # Calls for definitions not exercised by the chapter sequence.
 extra=ROOT/f'research/code-audit-20260909/extra{ch}.py'
 if extra.exists():
  try:
   exec(compile(extra.read_text(encoding='utf-8'),str(extra),'exec'),ns)
   result.append({'chapter':ch,'block':'additional checks','status':'pass'})
  except Exception as e:
   result.append({'chapter':ch,'block':'additional checks','status':'fail','trace':traceback.format_exc()})
   print(traceback.format_exc(),flush=True)
 (ROOT/'research/code-audit-20260909/results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')

if any(item['status']=='fail' for item in result):
 sys.exit(1)
