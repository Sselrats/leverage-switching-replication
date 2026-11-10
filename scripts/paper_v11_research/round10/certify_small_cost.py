#!/usr/bin/env python3
"""Initial small-cost cover for the adjacent policy in Proposition 4.3. Exact interval arithmetic throughout acceptance.
All costs in the finite cover are rational. The omitted punctured interval at zero
is covered analytically by a curvature/area estimate recorded below.
"""
import json,time
from pathlib import Path
import f2_kernel as k
from f2_kernel import Q,R,I,lo,hi,iv,M,mp,serial
ROOT=Path(__file__).resolve().parents[1]

def numeric_margin(pair,lam):
 a,b,g=k.proposal_band(pair,Q(str(lam)));c=k.xc(pair)
 qq=min(k.qn(a),k.qn(b));e=k.Gn(pair,M(c))-g
 delta=4*a-3 if pair==12 else 3-4*b
 return delta-b/(2*lam)*(e/qq)**2

def margin_box(pair,l,h):
 r=k.cost_box(pair,l,h);a,b,g=r['a'],r['b'],r['g']
 e=k.G(pair,R(k.xc(pair)))-g
 qa,qb=k.q(a),k.q(b)
 # q increases to 1 and decreases thereafter: its minimum on a band is at an endpoint.
 qm=I(min(lo(qa),lo(qb)),min(hi(qa),hi(qb)))
 delta=4*a-3 if pair==12 else 3-4*b
 err=b/(2*I(l,h))*(e/qm)**2
 return dict(pair=pair,l=l,h=h,a=a,b=b,g=g,delta=delta,q_min=qm,error=err,margin=delta-err)

def main():
 t=time.time();anchor=Q('0.0000001');eps=Q('0.01');zero=[]
 for pair in (12,23):
  x=k.xc(pair);X=I(x-eps,x+eps);qc=k.q(X);cv=-k.Gpp(pair,X)
  assert lo(qc)>0 and lo(cv)>0
  qm=R(lo(qc));qM=R(hi(qc));cm=R(lo(cv));cM=R(hi(cv))
  e0=cm*R(eps)**2/2
  CH=4*iv.sqrt(2)/(3*qM*iv.sqrt(cM))
  area=CH*e0**R(Q(3,2));assert lo(area)>2*anchor
  delta=4*R(x-eps)-3 if pair==12 else 3-4*R(x+eps)
  error=R(x+eps)/2/qm**2*(2/CH)**R(Q(4,3))*R(anchor)**R(Q(1,3))
  assert lo(delta-error)>0
  zero.append(dict(pair=pair,window=X,curvature=cv,q=qc,e0=e0,H_lower_at_e0=area,area_coefficient=CH,delta=delta,error_at_anchor=error,margin=delta-error))
 # Floating arithmetic proposes the bracket; signs below are certified.
 rt=mp.findroot(lambda ll:numeric_margin(23,ll),(M('.0018'),M('.0019')),tol=M('1e-50'))
 L,LF=k.near(rt,n=13,pad=2)
 boxes=[]
 def visit(l,h,depth=0):
  r12=margin_box(12,l,h);r23=margin_box(23,l,h)
  disjoint=lo(r12['a'])-hi(r23['b'])
  if lo(r12['margin'])>0 and lo(r23['margin'])>0 and disjoint>0:
   boxes.append(dict(l=l,h=h,test12=r12,test23=r23,disjoint_margin=disjoint));return
  assert depth<55,(l,h,depth)
  z=(l+h)/2;visit(l,z,depth+1);visit(z,h,depth+1)
 # Geometric stepping avoids a low-cost overestimation by the 1/lambda term.
 l=anchor
 while l<L:
  h=min(l*Q(11,10),L);visit(l,h);l=h
  if len(boxes)%20==0:print('BOXES',len(boxes),k.dec(l),flush=True)
 fail=margin_box(23,LF,LF);assert hi(fail['margin'])<0
 record=dict(status='COMPUTER-ASSISTED PROOF',judgement='verified',mpmath=mp.__version__,iv_dps=iv.dps,series_N=k.N,F='2',eta='3/2',anchor=anchor,punctured_zero=zero,certified_L=L,first_certified_failure=LF,maximal_initial_test_endpoint_bracket=(L,LF),finite_boxes=boxes,finite_box_count=len(boxes),failure=fail,root_checks=k.CHECKS,elapsed_seconds=time.time()-t)
 out=ROOT/'records'/'small_cost_test.json';out.write_text(json.dumps(serial(record),indent=2))
 print('RESULT',k.dec(L),k.dec(LF),len(boxes),flush=True)
if __name__=='__main__':main()
