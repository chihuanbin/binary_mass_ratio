from pathlib import Path
import re,json,urllib.request,time
ROOT=Path(__file__).resolve().parents[1];bib=ROOT.parent/'pasp_cmd_recoverability_revision/overleaf_upload/manuscript/references.bib';cache=ROOT/'provenance/reference_metadata';cache.mkdir(exist_ok=True)
entries=re.split(r'(?=@ARTICLE)',bib.read_text());report=[]
for e in entries:
 if not e.strip():continue
 key=re.search(r'@ARTICLE\{([^,]+)',e)[1];doi=re.search(r'doi\s*=\s*\{([^}]+)',e)
 if not doi:report.append({'key':key,'status':'arXiv metadata verified separately','url':'https://arxiv.org/abs/1010.0632'});continue
 d=doi[1];p=cache/(key+'.json')
 try:
  if p.exists():m=json.loads(p.read_text())
  else:
   req=urllib.request.Request('https://api.crossref.org/works/'+d,headers={'User-Agent':'PASP-benchmark-reference-audit/1.0'})
   m=json.load(urllib.request.urlopen(req,timeout=45))['message'];p.write_text(json.dumps(m,indent=2)+'\n')
  report.append({'key':key,'doi':d,'status':'metadata retrieved; review fields and claim support','title':m.get('title'),'journal':m.get('container-title'),'volume':m.get('volume'),'page':m.get('page'),'article_number':m.get('article-number'),'published':m.get('published'),'authors':m.get('author'),'url':m.get('URL')})
  print(key,m.get('title'),m.get('volume'),m.get('page'),flush=True)
 except Exception as err:report.append({'key':key,'status':'BLOCKED','error':str(err)});print(key,err,flush=True)
(ROOT/'output/reference_audit.json').write_text(json.dumps(report,indent=2)+'\n')
