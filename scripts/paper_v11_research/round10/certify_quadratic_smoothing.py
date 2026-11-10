#!/usr/bin/env python3
"""Quadratic smoothing with epsilon=lambda for Proposition 4.3.
Check the inactive generator inequality on the normalized state cells
x=a+t(b-a), 0<=t<=1, in each adjacent-pair band.
"""
from pathlib import Path
import json,time
import f2_kernel as k
from f2_kernel import Q,R,I,lo,hi,iv,serial
ROOT=Path(__file__).resolve().parents[1];t0=time.time();NS=80

def Lc(p,x,g):
 c=(lo(x)+hi(x))/2
 return k.L(p,R(c),g)+(x-R(c))*(k.G(p,x)-g)/k.q(x)
def Gc(p,x):
 c=(lo(x)+hi(x))/2
 return k.G(p,R(c))+(x-R(c))*k.Gp(p,x)
def check(l,h):
 ans=[];bands={p:k.cost_box(p,l,h) for p in (12,23)}
 if not hi(bands[23]['b'])<lo(bands[12]['a']):return None
 for p,other in [(12,3),(23,1)]:
  a,b,g=[bands[p][s] for s in ('a','b','g')];lam=I(l,h);La=k.L(p,a,g);mins=[]
  i,j=k.levels(p)
  for jcell in range(NS):
   tt=I(Q(jcell,NS),Q(jcell+1,NS));x=a+tt*(b-a)
   d=lam-Lc(p,x,g)+La
   th=(d+lam)/(2*lam);th=I(max(Q(0),lo(th)),min(Q(1),hi(th)))
   def reward(ll):return R(Q(3*ll,2))-R(Q(ll*ll,2))*x
   # Bilinear interpolation reaches its extrema at rectangle corners.
   def gap_at(xx,tt):
    ri=R(Q(3*(i-other),2))-R(Q(i*i-other*other,2))*xx
    rj=R(Q(3*(j-other),2))-R(Q(j*j-other*other,2))*xx
    return ri+tt*(rj-ri)
   gap=k.hull(*[gap_at(R(xx),R(tt)) for xx in (lo(x),hi(x)) for tt in (lo(th),hi(th))])
   slope=(g-Gc(p,x))/k.q(x)
   margin=gap-x*slope**2/(4*lam)
   if lo(margin)<=0:return None
   mins.append(lo(margin))
  ans.append(dict(pair=p,inactive=other,minimum_generator_margin=min(mins),cell_margins=mins,bands=bands[p]))
 return dict(l=l,h=h,tests=ans,disjoint_margin=lo(bands[12]['a'])-hi(bands[23]['b']))

def main():
 start=Q('0.0018625115132');end=Q('.01');done=[];failed=[]
 def visit(l,h,depth=0):
  r=check(l,h)
  if r:done.append(r);return
  if depth>=10:failed.append(dict(l=l,h=h,reason='Quadratic smoothing cell test did not pass at maximum depth; no nonexistence conclusion.'));return
  z=(l+h)/2;visit(l,z,depth+1);visit(z,h,depth+1)
 n=100
 for j in range(n):
  visit(start+(end-start)*j/n,start+(end-start)*(j+1)/n)
  if j%10==0:print('SMOOTH',j,len(done),len(failed),flush=True)
  if failed:break
 rec=dict(status='COMPUTER-ASSISTED PROOF' if not failed else 'PARTIAL',judgement='verified' if not failed else 'plausible but unverified',F='2',eta='3/2',mpmath=k.mp.__version__,iv_dps=iv.dps,N=k.N,epsilon='lambda',cost_cover=(start,end),spatial_normalised_cells=NS,boxes=done,failed=failed,all_passed=not failed,root_checks=k.CHECKS,elapsed_seconds=time.time()-t0)
 (ROOT/'records/quadratic_smoothing.json').write_text(json.dumps(serial(rec),indent=2))
 print('SMOOTH_DONE',not failed,len(done),flush=True)
if __name__=='__main__':main()
