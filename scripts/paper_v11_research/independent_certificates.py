#!/usr/bin/env python3
"""Independent interval certificates for local endpoint branches and a base-case
max-pair bridge. mpmath==1.3.0; all accepted inequalities use iv arithmetic.
mp arithmetic ONLY proposes root brackets and preconditions Krawczyk maps.
Usage: python independent_certificates.py OUTDIR endpoints|base
The scalar F and eta are exact rationals. This is not a rectangle engine.
"""
import sys,json,time,math,hashlib
from pathlib import Path
from fractions import Fraction as Q
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from decimal_bounds import decimal_bounds

assert mp.__version__=='1.3.0'
mp.mp.dps=90;mp.iv.dps=80;iv=mp.iv;N=180
OUT=Path(sys.argv[1]);OUT.mkdir(parents=True,exist_ok=True)
def V(x):
 x=Q(x);return iv.mpf(x.numerator)/x.denominator
def exact_endpoint(t):
 s,man,ex,bc=t
 return Q((-1 if s else 1)*int(man))*(Q(2)**int(ex))
def lo(x):return exact_endpoint(x._mpi_[0])
def hi(x):return exact_endpoint(x._mpi_[1])
def I(l,h):return iv.mpf([V(l).a,V(h).b])
def M(q):q=Q(q);return mp.mpf(q.numerator)/q.denominator
def sg(x):return 1 if lo(x)>0 else -1 if hi(x)<0 else 0
def enc(x):
 if hasattr(x,'_mpi_'):
  lower,upper=lo(x),hi(x)
  left,right=decimal_bounds(lower,upper,28)
  return {'lo':str(lower),'hi':str(upper),'decimal_lo':left,'decimal_hi':right}
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v) for v in x]
 return x
def write(n,x):(OUT/n).write_text(json.dumps(enc(x),indent=2))
def vmax(a,b):return I(max(lo(a),lo(b)),max(hi(a),hi(b)))
class CIR:
 def __init__(self,F,eta=Q(3,2)):
  self.F0=Q(F);self.e0=Q(eta);self.F=V(F);self.eta=V(eta)
  self.q0=self.F**(self.F-1)/iv.gamma(self.F)
  self.fmid=M(F);self.emid=M(eta);self.qmid=self.fmid**(self.fmid-1)/mp.gamma(self.fmid)
  self.cache={};self.root_checks=[];self.hcache={};self.seed={}
 def q(self,x):return self.q0*x**self.F*iv.exp(-self.F*x)
 def E(self,x):
  """Integral 1 to x of 1/q, via integrated exponential series.
  Uniform remainder <= length * sup(t^-F) * exponential tail / q0.
  Handles integer F by the logarithmic term exactly.
  """
  t=V(1);s=V(0)
  for n in range(N+1):
   p=Q(n+1)-self.F0
   s+=t*(iv.log(x) if p==0 else (x**V(p)-1)/V(p))
   t*=self.F/(n+1)
  xmax=max(Q(1),hi(x));xmin=min(Q(1),lo(x));length=max(abs(lo(x)-1),abs(hi(x)-1))
  assert xmin>0
  ratio=self.F*V(xmax)/(N+2);assert hi(ratio)<1
  tail=t*V(xmax)**(N+1)/(1-ratio)*V(xmin)**(-self.F)*V(length)/self.q0
  rem=I(0,hi(tail)) if lo(x)>=1 else I(-hi(tail),0) if hi(x)<=1 else I(-hi(tail),hi(tail))
  return s/self.q0+rem
 def theta_A(self,x):
  t=V(1);theta=t;A=x
  for n in range(1,N+1):
   t*=self.F*x/(self.F+n);theta+=t;A+=t*x/(n+1)
  nxt=t*self.F*x/(self.F+N+1);ratio=self.F*x/(self.F+N+2);assert hi(ratio)<1
  rt=nxt/(1-ratio);ra=nxt*x/(N+2)/(1-ratio)
  return theta+I(0,hi(rt)),A+I(0,hi(ra))
 def A(self,x):return self.theta_A(x)[1]
 def G(self,p,x):
  q=self.q(x);theta,_=self.theta_A(x);P=q*theta
  if p==23:return (self.eta-V(Q(5,2)))*P+V(Q(5,2))*q
  i,j=(1,2) if p==12 else (1,3);B=i*self.eta-V(Q(i*i,2));coef=V(Q(j*j-i*i,2))
  return B+((j-i)*self.eta-coef)*P+coef*q
 def K(self,x):return (self.eta-V(Q(5,2)))*self.A(x)+V(Q(5,2))*x
 def L(self,p,x,g):
  i,j=(1,2) if p==12 else (1,3);B=i*self.eta-V(Q(i*i,2));coef=V(Q(j*j-i*i,2))
  return (B-g)*self.E(x)+((j-i)*self.eta-coef)*self.A(x)+coef*x
 def Gn(self,p,x):
  f=self.fmid;e=self.emid;q=self.qmid*x**f*mp.exp(-f*x);P=mp.gammainc(f,0,f*x)/mp.gamma(f)
  if p==23:return (e-mp.mpf('2.5'))*P+mp.mpf('2.5')*q
  i,j=(1,2) if p==12 else (1,3);coef=mp.mpf(j*j-i*i)/2
  return i*e-mp.mpf(i*i)/2+((j-i)*e-coef)*P+coef*q
 def proposal(self,p,level,side):
  xc=self.emid/(mp.mpf('1.5') if p==12 else 2 if p==13 else mp.mpf('2.5'))
  left,right=(mp.mpf('1e-35'),xc) if side=='a' else (xc,mp.mpf(12));fn=lambda x:self.Gn(p,x)-M(level)
  seed=self.seed.get((p,side),xc*(mp.mpf('.5') if side=='a' else mp.mpf('1.7')))
  try:
   r=mp.findroot(fn,(seed,seed*(1+mp.mpf('1e-6'))),tol=mp.mpf('1e-75'),maxsteps=50)
   assert left<r<right
  except Exception:
   fl,fr=fn(left),fn(right);assert fl*fr<0
   for _ in range(150):
    z=(left+right)/2
    if (fn(z)<0)==(side=='a'):left=z
    else:right=z
   r=(left+right)/2
  self.seed[p,side]=r;return r
 def root(self,p,level,side):
  key=(p,lo(level),hi(level),side)
  if key in self.cache:return self.cache[key]
  xs=[self.proposal(p,t,side) for t in (lo(level),hi(level))];scale=10**42
  l0=Q(int(mp.floor(min(xs)*scale)),scale);h0=Q(int(mp.ceil(max(xs)*scale)),scale);rad=Q(1,10**38)
  xc=self.e0/(Q(3,2) if p==12 else 2 if p==13 else Q(5,2));left,right=(Q(1,10**35),xc) if side=='a' else (xc,Q(12))
  desired=(-1,1) if side=='a' else (1,-1)
  for _ in range(110):
   l=max(left,l0-rad);h=min(right,h0+rad);fl=self.G(p,V(l))-level;fh=self.G(p,V(h))-level
   if (sg(fl),sg(fh))==desired and left<l<h<right:
    r=I(l,h);self.cache[key]=r;self.root_checks.append({'pair':p,'side':side,'level':level,'root':r,'residuals':(fl,fh),'branch':(left,right)});return r
   rad*=2
  raise AssertionError('interval bracket verification failed')
 def H(self,p,g):
  if not hasattr(g,'_mpi_'):g=V(g)
  key=(p,lo(g),hi(g))
  if key not in self.hcache:
   a,b=self.root(p,g,'a'),self.root(p,g,'b');self.hcache[key]=self.L(p,b,g)-self.L(p,a,g)
  return self.hcache[key]
 def Hbox(self,p,gl,gh):return I(lo(self.H(p,gh)),hi(self.H(p,gl)))
 def cost(self,g):return vmax(self.H(13,g),self.H(12,g))/2
 def boundary(self,gl,gh):
  g=I(gl,gh);a,b=self.root(13,g,'a'),self.root(13,g,'b');u,be=self.root(12,g,'a'),self.root(12,g,'b');ze=self.root(23,V(0),'b')
  H13=self.Hbox(13,gl,gh);H12=self.Hbox(12,gl,gh);lam=I(lo(self.cost(gh)),hi(self.cost(gl)))
  va=self.L(12,be,g)-self.L(12,a,g);v0=va-self.K(a);t0=v0-H13+lam;tz=va-H13+lam+self.K(ze)-self.K(a)
  tests={'a>0':a,'u>a':u-a,'zeta>u':ze-u,'b>zeta':b-ze,'beta>b':be-b,'a<2eta5':self.eta/V(Q(5,2))-a,'beta>2eta3':be-self.eta/V(Q(3,2)),'T0>0':t0,'V0>0':v0,'2lambda>Tzeta':2*lam-tz}
  return {'gain':(gl,gh),'lambda':lam,'roots':dict(a=a,b=b,u=u,beta=be,zeta=ze),'H13':H13,'H12':H12,'tests':tests,'passed':all(sg(v)==1 for v in tests.values())}
 def Droot(self,l,h,iterations=65):
  fn=lambda x:self.H(13,x)-self.H(12,x);fl,fh=fn(l),fn(h);assert sg(fl)==-1 and sg(fh)==1
  for _ in range(iterations):
   z=(l+h)/2;s=sg(fn(z))
   if s==-1:l=z
   elif s==1:h=z
   else:break
  g=I(l,h);a,b=self.root(13,g,'a'),self.root(13,g,'b');u,be=self.root(12,g,'a'),self.root(12,g,'b');J13=self.E(b)-self.E(a);J12=self.E(be)-self.E(u)
  assert sg(J12-J13)==1
  return {'gain':g,'lambda':I(lo(self.H(13,h)/2),hi(self.H(13,l)/2)),'C2':2/J13-2/J12,'J13':J13,'J12':J12,'difference_derivative':J12-J13,'signs':(fn(l),fn(h))}
 def lower(self,g,ga):
  a,b=self.root(13,g,'a'),self.root(13,g,'b');k=self.G(23,ga);u,be=self.root(12,g-k,'a'),self.root(12,g-k,'b');ze=self.root(23,k,'b')
  ea,eb,eg,ez,ee=[self.E(x) for x in [a,b,ga,ze,be]]
  va=self.L(12,be,g-k)-self.L(12,a,g-k);lam=self.H(13,g)/2
  Qc=va+self.K(ga)-self.K(a)-k*(eg-ea);Z=va+self.K(ze)-self.K(a)-k*(ez-ea)
  kp=(self.eta-V(Q(5,2))*ga)*self.F*self.q(ga)/ga
  A=ea-ee;B=kp*(ee-eg);C=kp*(ee-ez);J13=eb-ea
  J=[[A+J13/2,B],[A+3*J13/2,C]]
  p=(ez-eg)/(ee-eg);s4=(2+p)/(p*(ee-ea));Vu=self.L(12,be,g-k)-self.L(12,u,g-k)
  checks={'a>gamma':a-ga,'u>a':u-a,'zeta>u':ze-u,'b>zeta':b-ze,'beta>b':be-b,'a<2eta5':self.eta/V(Q(5,2))-a,'beta>2eta3':be-self.eta/V(Q(3,2)),'gamma<2eta5':self.eta/V(Q(5,2))-ga,'kprime':kp,'independent_lambda_jacobian':A*(C-B),'Tzeta>Vu':Z-lam-Vu,'C1':s4-2/J13}
  return {'residual':[Qc-lam,Z-3*lam],'Jacobian':J,'checks':checks,'C1':s4-2/J13,'lambda':lam,'s4':s4,'s13':2/J13,'p':p,'roots':dict(a=a,b=b,gamma=ga,u=u,zeta=ze,beta=be)}

def endpoints():
 C=CIR(2);mid=[Q('1.4453304056722230242250485'),Q('.2810797139292023310764377')];rad=Q(1,10**20);box=[I(t-rad,t+rad) for t in mid]
 z=C.lower(*box);zm=C.lower(*map(V,mid));J=z['Jacobian'];Jm=mp.matrix([[M((lo(v)+hi(v))/2) for v in row] for row in J]);R=Jm**-1
 Rv=[[V(Q(mp.nstr(R[i,j],75))) for j in range(2)] for i in range(2)]
 K=[]
 for i in range(2):
  k=V(mid[i])-sum((Rv[i][j]*zm['residual'][j] for j in range(2)),V(0))
  for j in range(2):k+=((1 if i==j else 0)-sum((Rv[i][t]*J[t][j] for t in range(2)),V(0)))*(box[j]-V(mid[j]))
  K.append(k)
 assert all(lo(box[i])<lo(K[i])<=hi(K[i])<hi(box[i]) for i in range(2))
 # A contraction bound makes uniqueness on the box explicit as well.
 norm=max(sum(max(abs(lo((1 if i==j else 0)-sum((Rv[i][t]*J[t][j] for t in range(2)),V(0)))),abs(hi((1 if i==j else 0)-sum((Rv[i][t]*J[t][j] for t in range(2)),V(0))))) for j in range(2)) for i in range(2));assert norm<1
 zk=C.lower(*K);assert all(sg(v)==1 for v in zk['checks'].values())
 root=C.Droot(Q('1.13172670023213'),Q('1.13172670023215'),70);boundary=C.boundary(Q('1.13172670023213'),Q('1.13172670023215'));assert boundary['passed']
 rec={'status':'COMPUTER-ASSISTED PROOF','judgement':'verified','F':'2','eta':'3/2','mpmath':mp.__version__,'iv_dps':iv.dps,'N':N,'lower_root':{'initial_box':box,'Krawczyk_image':K,'contraction_bound':norm,'evaluation':zk},'C1':zk['C1'],'C2':root['C2'],'upper_root':root,'upper_local_bridge':boundary,'root_checks':C.root_checks,'proof':'proofs/R2_R7_local_branches.tex'}
 write('independent_endpoints.json',rec);print('ENDPOINTS PASSED',enc(rec['C1']),enc(rec['C2']),flush=True)
def base():
 C=CIR(Q(124,45));gl,gh=Q('1.145'),Q('1.17');passed=[];failed=[]
 root=C.Droot(Q('1.15'),Q('1.17'))
 def visit(l,h,d=0):
  r=C.boundary(l,h)
  if r['passed']:passed.append(r);return
  z=(l+h)/2;point=C.boundary(z,z)
  if not point['passed']:
   failed.append({'gain':(l,h),'reason':'The midpoint test of this particular entrance construction did not pass; no impossibility claim.','point':point});return
  if d>=7:
   failed.append({'gain':(l,h),'reason':'Interval dependency at maximum depth 7.','box':r});return
  visit(l,z,d+1);visit(z,h,d+1)
 for j in range(100):
  visit(gl+(gh-gl)*j/100,gl+(gh-gl)*(j+1)/100);print('BASE BOX',j,len(passed),len(failed),flush=True)
 intervals=sorted((r['gain'][0],r['gain'][1]) for r in passed);components=[]
 for l,h in intervals:
  if components and l==components[-1][1]:components[-1]=(components[-1][0],h)
  else:components.append((l,h))
 selected=[(l,h) for l,h in components if l<=lo(root['gain']) and hi(root['gain'])<=h]
 rec={'status':'COMPUTER-ASSISTED PROOF OF LOCAL S=0' if selected else 'PARTIAL','judgement':'verified' if selected else 'plausible but unverified','F':'124/45','eta':'3/2','mpmath':mp.__version__,'iv_dps':iv.dps,'N':N,'tested_gain_cover':(gl,gh),'certified_gain_components':components,'D_root':root,'boxes':passed,'failed':failed,'root_checks':C.root_checks,'proof':'proofs/R2_R7_local_branches.tex'}
 if selected:
  l,h=selected[0];rec['selected_gain_component']=(l,h);rec['S_zero_cost_interval']=(hi(C.cost(h)),lo(C.cost(l)))
 if passed:rec['minimum_margins']={k:min(lo(r['tests'][k]) for r in passed) for k in passed[0]['tests']}
 write('independent_base_boundary.json',rec);print('BASE',rec['status'],enc(rec.get('S_zero_cost_interval')),enc(root['lambda']),flush=True)
start=time.monotonic()
try:
 {'endpoints':endpoints,'base':base}[sys.argv[2]]()
except Exception as e:
 import traceback
 write('independent_'+sys.argv[2]+'_FAILED.json',{'status':'FAILED','judgement':'verified','error':repr(e),'traceback':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start});raise
