from pathlib import Path
import hashlib,json,zipfile,datetime
ROOT=Path(__file__).resolve().parents[1];WORK=ROOT.parent;DEST=ROOT/'release';DEST.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def allowed(p):return p.is_file() and p.name!='.DS_Store' and '__pycache__' not in p.parts and '.pytest_cache' not in p.parts and not p.name.startswith('Lucy1974_original')

def main():
 final=WORK/'manuscript/pasp_submission_final';required=[final/x for x in ['main.tex','main.pdf','references.bib','aastex701.cls','aasjournalv7.bst']]+list((final/'figures').glob('*.pdf'))+list((final/'tables').glob('*.tex'))
 with zipfile.ZipFile(DEST/'pasp_submission_final_overleaf.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in sorted(required):z.write(p,p.relative_to(final))
 oldfiles=[p for p in (WORK/'pasp_cmd_recoverability_revision').rglob('*') if allowed(p)]
 newfiles=[p for p in ROOT.rglob('*') if allowed(p) and DEST not in p.parents and 'page_renders' not in p.parts]
 archfiles=[p for p in (WORK/'arch_code_results').rglob('*') if allowed(p) and ('data' not in p.parts or p.name in ['params.csv','photometric_errors.npz'])]
 sourcefiles=sorted(set(oldfiles+newfiles+archfiles+required));rels={str(p.relative_to(WORK)) for p in sourcefiles}
 original=json.loads((ROOT/'provenance/original_asset_hashes.json').read_text());excluded=[r['path'] for r in original if r['path'] not in rels]
 release_meta={'date_prepared':datetime.date.today().isoformat(),'status':'LOCAL_PREPARED_NOT_PUBLISHED','release_tag':'RELEASE_TAG_TO_BE_FILLED','release_commit':'RELEASE_COMMIT_TO_BE_FILLED','release_date':'RELEASE_DATE_TO_BE_FILLED','doi':'ZENODO_DOI_TO_BE_FILLED','existing_public_repository':'https://github.com/chihuanbin/binary_mass_ratio','excluded_original_raw_assets':excluded,'included_source_assets':len(sourcefiles),'includes_all_new_scientific_outputs':True}
 if (ROOT/'provenance/public_release.json').exists():
  public=json.loads((ROOT/'provenance/public_release.json').read_text())
  release_meta.update(status=public['status'],release_tag=public['tag'],release_commit=public['commit'],release_date=public['release_date'],doi=public['doi'],github_release_url=public['github_release_url'])
 manifest=[{'path':str(p.relative_to(WORK)),'size':p.stat().st_size,'sha256':sha(p)} for p in sourcefiles]
 (DEST/'release_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(DEST/'deposit_metadata.json').write_text(json.dumps(release_meta,indent=2)+'\n')
 with zipfile.ZipFile(DEST/'pasp_final_reproducibility.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in sourcefiles:z.write(p,p.relative_to(WORK))
  z.writestr('revision_pasp_final/provenance/portable_excluded_assets.json',json.dumps(excluded,indent=2)+'\n')
  z.writestr('revision_pasp_final/provenance/portable_release_manifest.json',json.dumps(manifest,indent=2)+'\n')
  z.writestr('revision_pasp_final/provenance/portable_deposit_metadata.json',json.dumps(release_meta,indent=2)+'\n')
 for p in DEST.glob('*.zip'):
  with zipfile.ZipFile(p) as z:assert z.testzip() is None
 digest={p.name:{'sha256':sha(p),'size':p.stat().st_size} for p in DEST.glob('*.zip')};(DEST/'archive_hashes.json').write_text(json.dumps(digest,indent=2)+'\n');print(json.dumps(digest,indent=2))
if __name__=='__main__':main()
