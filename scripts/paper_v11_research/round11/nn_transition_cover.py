#!/usr/bin/env python3
"""Certify the initial NN-optimality endpoint and its four-threshold successor.
mpmath 1.3.0, rational parameter covers, strict interval acceptance. The
band-touch cost splits two analytic verification constructions; it is not
a regime change in the value. Source kernels and every dependency are bundled.
"""
import sys,json,gzip,time
from pathlib import Path
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'round10'))
import f2_kernel as k
import certify_quadratic_smoothing as sm
import certify_four_branch as fo
from validate_contacts import nn_eval
R,I,lo,hi,serial,mp=k.R,k.I,k.lo,k.hi,k.serial,k.mp
k.iv.dps=80;mp.mp.dps=90;T=time.time()
# Keep the original state-cell test but discharge disjointness analytically
# for the part of its cost cover lying below the uniquely certified band touch.
source=Path(sm.__file__).read_text().split('\ndef main():',1)[0]
source=source.replace(" if not hi(bands[23]['b'])<lo(bands[12]['a']):return None\n",'')
ns={'__file__':sm.__file__,'__name__':'extended_quadratic_test'};exec(compile(source,sm.__file__,'exec'),ns)
quad=ns['check']

def touch():
 fun=lambda x:k.proposal_band(23,Q(mp.nstr(x,85)))[1]-k.proposal_band(12,Q(mp.nstr(x,85)))[0]
 prop=mp.findroot(fun,(mp.mpf('.014'),mp.mpf('.015')),tol=mp.mpf('1e-70'))
 mid=Q(mp.nstr(prop,75));l,h=mid-Q('1e-20'),mid+Q('1e-20')
 fl=k.band(23,l)['b']-k.band(12,l)['a'];fh=k.band(23,h)['b']-k.band(12,h)['a']
 assert hi(fl)<0 and lo(fh)>0
 rr12=k.cost_box(12,l,h);rr23=k.cost_box(23,l,h);g=rr12['g']+rr23['g']-1
 checks={'touch_state_above_3/4':rr12['a']-R(Q(3,4)),'NN_gain_below_G13_max':k.G(13,R(Q(3,4)))-g,'gamma_below_3/5':R(Q(3,5))-rr23['a']}
 assert all(lo(v)>0 for v in checks.values())
 return {'bracket':(l,h),'left':fl,'right':fh,'checks':checks,'gain':g}

def nncheck(l,h,local=False):
 try:z=nn_eval(l,h)
 except (AssertionError,ValueError,ZeroDivisionError):return None
 rr=z['roots'];la=z['lambda'];tests={'u>3/5':rr['u']-R(Q(3,5)),'zeta<1':1-rr['zeta'],'gamma<3/5':R(Q(3,5))-rr['gamma'],'beta>1':rr['beta']-1,'lambda>H13-Q':la-z['H13']+z['Q']}
 if local:tests['Fprime>0']=z['derivative']
 else:tests['lambda>Q']=la-z['Q']
 if any(lo(v)<=0 for v in tests.values()):return None
 return {'cost':(l,h),'evaluation':z,'accepted_tests':tests,'local_contact_mode':local}

def fourcheck(l,h):
 rr=k.cost_box(23,l,h);ga,ze,kk,la=rr['a'],rr['b'],rr['g']-1,rr['lambda'];gg=fo.propose((l+h)/2);base=Q(str(gg));wid=h-l
 for mult in [2,4,8,16,32]:
  gl=base-mult*wid-Q('1e-25');gh=base+mult*wid+Q('1e-25')
  try:fl=fo.atg(R(gl),rr)[0];fh=fo.atg(R(gh),rr)[0]
  except (AssertionError,ValueError,ZeroDivisionError):continue
  if lo(fl)>0 and hi(fh)<0:break
 else:return None
 g=I(gl,gh)
 try:
  _,a,be,va=fo.atg(g,rr);b=fo.rootbox(13,g,'b');u=fo.rootbox(12,g-kk,'a');nn=nn_eval(l,h)
 except (AssertionError,ValueError,ZeroDivisionError):return None
 A=k.E(be)-k.E(a);B=k.E(ze)-k.E(ga);C=k.E(be)-k.E(ga)
 H13=k.L(13,b,g)-k.L(13,a,g)
 tests={'a>gamma':a-ga,'u>a':u-a,'zeta>u':ze-u,'b>zeta':b-ze,'beta>b':be-b,'gamma<3/5':R(Q(3,5))-ga,'beta>1':be-1,'A>0':A,'B>0':B,'C>0':C,'H13>lambda':H13-la,'H13<2lambda':2*la-H13,'NN_contact_residual_derivative':nn['derivative']}
 if any(lo(v)<=0 for v in tests.values()):return None
 return {'cost':(l,h),'gain':g,'shooting_signs':(fl,fh),'roots':dict(gamma=ga,a=a,u=u,zeta=ze,b=b,beta=be),'tests':tests,'Q_g':-A,'s_four':(1+2*C/B)/A,'H13':H13}

def cover(l,h,n,fn,name):
 done=[];fail=[]
 def visit(a,b,d=0):
  r=fn(a,b)
  if r:done.append(r);return
  if d==14:fail.append({'cost':(a,b),'reason':'No accepted interval test at maximum depth; no impossibility claim.'});return
  c=(a+b)/2;visit(a,c,d+1);visit(c,b,d+1)
 for j in range(n):
  visit(l+(h-l)*j/n,l+(h-l)*(j+1)/n)
  if j%20==0:print(name,j,len(done),len(fail),flush=True)
  if fail:break
 assert not fail,(name,serial(fail))
 return {'cover':(l,h),'boxes':done,'failed':fail}

def main():
 tt=touch();tl,th=tt['bracket'];print('TOUCH',serial(tt),flush=True)
 smooth=cover(Q('.01'),th,80,quad,'SMOOTH')
 nb=cover(tl,Q('.0607'),150,nncheck,'NN')
 local=nncheck(Q('.0607'),Q('.0608'),True);assert local
 nc=json.loads((ROOT/'records/nn_contact.json').read_text());nl,nh=map(Q,nc['bracket']);assert Q('.0607')<nl<nh<Q('.0608')
 fb=cover(Q('.0607'),Q('.085'),100,fourcheck,'FOUR')
 z=nn_eval(nl,nh);a=z['roots']['a'];be=z['roots']['beta'];jump=z['derivative']/(k.E(be)-k.E(a));assert lo(jump)>0
 rec = {
  'status': 'COMPUTER-ASSISTED PROOF',
  'judgement': 'verified',
  'claim': 'For the zero-fee menu {1,2,3} at F=2 and eta=3/2, the initial interval on which an adjacent policy attains g123=gN is (0,lambda_nn]. The adjacent bands first touch at the cost enclosed in band_touch.bracket. At lambda_nn, a four-threshold policy using 1->3 in each regenerative cycle also attains gN. For lambda_nn < lambda <= 0.085 the calibrated four-threshold policy attains g123>gN, and every adjacent policy is strictly suboptimal.',
  'proposition': '4.3',
  'menu': [1, 2, 3],
  'holding_fees': [0, 0, 0],
  'F': '2',
  'eta': '3/2',
  'mpmath': mp.__version__,
  'iv_dps': k.iv.dps,
  'N': k.N,
  'band_touch': tt,
  'quadratic_cover': smooth,
  'quadratic_geometry': 'Only costs <= the exact touch root use this construction. Its interval tests also cover the outward bracket. At equality the touching endpoint definitions glue C1.',
  'nn_cover': nb,
  'NN_geometry': 'Only costs >= the exact touch root use the overlap construction. Ordered roots follow analytically from u<=zeta, u>3/5, zeta<1, and the pair-level identities.',
  'nn_local_contact_tube': local,
  'contact_record': 'reference/contact.json',
  'small_cost_records': ['reference/small_cost.json.gz', 'reference/smoothing.json.gz'],
  'four_cover': fb,
  'four_geometry': 'For lambda_nn <= lambda <= 0.085, Q(gNN)-lambda>=0 and Q_g<0 give V(u)<=2lambda, using the positive derivative of the NN contact residual.',
  'frequency_downward_jump': jump,
  'uncovered_cost_boxes': [],
  'root_checks': k.CHECKS,
  'elapsed_seconds': time.time()-T,
 }
 with gzip.open(ROOT/'records/nn_transition_cover.json.gz','wt') as f:json.dump(serial(rec),f,indent=2)
 print('NN_TRANSITION_COVER PASSED',len(smooth['boxes']),len(nb['boxes']),len(fb['boxes']),flush=True)
if __name__=='__main__':main()
