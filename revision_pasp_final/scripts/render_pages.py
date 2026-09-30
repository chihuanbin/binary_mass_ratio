from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[1];pdf=ROOT.parent/'manuscript/pasp_submission_final/main.pdf';out=ROOT/'reports/page_renders';out.mkdir(exist_ok=True)
subprocess.run(['pdftoppm','-png','-r','110',str(pdf),str(out/'final')],check=True)
