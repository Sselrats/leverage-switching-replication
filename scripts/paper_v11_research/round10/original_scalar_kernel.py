"""Load the scalar general-F interval kernel of independent_certificates.py without its command-line driver.
Project layout: this file lives in scripts/paper_v11_research/round10/, the kernel module one level up."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]            # scripts/paper_v11_research
source=ROOT/'independent_certificates.py'
text=source.read_text();head,tail=text.split('\nstart=time.monotonic()\n',1)
old=sys.argv[:];sys.argv=[str(source),str(ROOT/'records')]
ns={'__name__':'supplied_scalar_kernel','__file__':str(source)}
try:exec(compile(head,str(source),'exec'),ns)
finally:sys.argv=old
for name in ('CIR','V','I','lo','hi','enc','Q','mp','iv','M','N','sg'):
 globals()[name]=ns[name]
