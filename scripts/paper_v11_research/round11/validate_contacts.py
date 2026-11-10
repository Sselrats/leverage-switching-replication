#!/usr/bin/env python3
"""Validated contact systems, mpmath 1.3.0. Exact rational boxes;
mp used only for preconditioners/proposals, strict iv acceptance throughout.
Run from any directory. Reports local contact results, NOT global zero-set topology.
"""
import sys,json,time,hashlib
from pathlib import Path
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'round10'))
import original_scalar_kernel as sk
import f2_kernel as fk
mp=sk.mp;iv=sk.iv;mp.mp.dps=90;iv.dps=80
V,I,lo,hi,sg,enc,M=sk.V,sk.I,sk.lo,sk.hi,sk.sg,sk.enc,sk.M
OUT=ROOT/'records';OUT.mkdir(exist_ok=True)
def write(name,r):
 (OUT/name).write_text(json.dumps(enc(r),indent=2));print(name,r['status'],flush=True)
def krawczyk(fun,mid,radii):
 box=[I(t-r,t+r) for t,r in zip(mid,radii)];z=fun(*box);z0=fun(*map(V,mid));J=z['Jacobian'];n=len(mid)
 Rm=mp.matrix([[M((lo(v)+hi(v))/2) for v in row] for row in J])**-1
 R=[[V(Q(mp.nstr(Rm[i,j],75))) for j in range(n)] for i in range(n)]
 D=[[(1 if i==j else 0)-sum((R[i][t]*J[t][j] for t in range(n)),V(0)) for j in range(n)] for i in range(n)]
 image=[V(mid[i])-sum((R[i][j]*z0['residual'][j] for j in range(n)),V(0))+sum((D[i][j]*(box[j]-V(mid[j])) for j in range(n)),V(0)) for i in range(n)]
 norm=max(sum(max(abs(lo(v)),abs(hi(v))) for v in row) for row in D)
 assert all(lo(box[i])<lo(image[i])<=hi(image[i])<hi(box[i]) for i in range(n)),enc({'box':box,'image':image})
 assert norm<1
 zi=fun(*image)
 return {'initial_box':box,'midpoint':mid,'rational_preconditioner':R,'residual_at_midpoint':z0['residual'],'Jacobian_on_box':J,'I_minus_RJ':D,'contraction_bound':norm,'Krawczyk_image':image,'evaluation':zi}

def base_upper_eval(c,g,ga):
 a,b=c.root(13,g,'a'),c.root(13,g,'b');k=c.G(23,ga);h=g-k
 u,be=c.root(12,h,'a'),c.root(12,h,'b');ze=c.root(23,k,'b')
 ea,eb,eg,eu,ee,ez=[c.E(x) for x in (a,b,ga,u,be,ze)]
 A=ee-ea;C=ee-eg;J=eb-ea;J12=ee-eu;kp=(c.eta-V(Q(5,2))*ga)*c.F*c.q(ga)/ga
 lam=c.H(13,g)/2;va=c.L(12,be,h)-c.L(12,a,h)
 Qc=va+c.K(ga)-c.K(a)-k*(eg-ea);Z=va+c.K(ze)-c.K(a)-k*(ez-ea)
 H12=c.L(12,be,h)-c.L(12,u,h)
 p=(ee-eu)/(ee-eg);Tc=p*(eb-eg);sbr=(2+p)/Tc
 Jmat=[[-A+J/2,kp*C],[-J12+J,kp*J12]]
 # Along the entrance-contact branch Q=lambda(g), k_lambda=(1-2A/J)/C.
 kl=(1-2*A/J)/C
 M12prime=2-2*J12/J-J12*kl
 checks={'gamma>0':ga,'a>gamma':a-ga,'u>a':u-a,'zeta>u':ze-u,'b>zeta':b-ze,'beta>b':be-b,'a<xc23':c.eta/V(Q(5,2))-a,'b>xc12':b-c.eta/V(Q(3,2)),'k>0':k,'kprime>0':kp,'M23':3*lam-Z,'s13_minus_sbranch':2/J-sbr,'minus_M12prime':-M12prime}
 return {'residual':[Qc-lam,H12-2*lam],'Jacobian':Jmat,'lambda':lam,'g':g,'k':k,'roots':dict(a=a,b=b,u=u,beta=be,gamma=ga,zeta=ze),'checks':checks,'p':p,'Tc':Tc,'Nc':2+p,'frequency':sbr,'s13':2/J,'contact_derivatives':{'k_lambda':kl,'M12_lambda':M12prime},'Q':Qc,'Z':Z}

def base_lower():
 c=sk.CIR(Q(124,45));mid=list(map(Q,['1.27817577160732696897867395974719610496575','.212786830745532072006121174569243151450474']))
 r=krawczyk(c.lower,mid,[Q('1e-22'),Q('1e-22')]);assert all(sg(v)==1 for v in r['evaluation']['checks'].values())
 # Parameterized entrance branch derivatives on a rational neighbourhood.
 g,ga=r['initial_box'];z=c.lower(g,ga);rr=z['roots'];ea,eb,eg,ez,ee=[c.E(rr[t]) for t in ['a','b','gamma','zeta','beta']]
 kl=(1-2*(ee-ea)/(eb-ea))/(ee-eg);M23prime=2+(ez-eg)*kl
 assert sg(M23prime)==1
 write('base_lower_contact.json',{'status':'COMPUTER-ASSISTED PROOF','judgement':'verified','claim':'Unique local lower contact; local one-sided zero/positive transition. Global joining is not asserted here.','F':'124/45','eta':'3/2','mpmath':mp.__version__,'iv_dps':iv.dps,'N':sk.N,'contact':r,'M23_lambda':M23prime,'root_checks':c.root_checks})

def base_upper():
 from explore_contacts import Model
 mp.mp.dps=70;c0=Model(mp.mpf(124)/45)
 seed=['.295760536','1.240661025','.414387018','1.97684312','.000273','1.03','1.155015864','.466738999']
 d=c0.contact_fast('upper',seed);print('UPPER_PROPOSAL',json.dumps({k:mp.nstr(v,50) for k,v in d.items()}),flush=True)
 mp.mp.dps=90
 c=sk.CIR(Q(124,45));mid=[Q(mp.nstr(d[t],65)) for t in ['g','gamma']]
 r=krawczyk(lambda g,ga:base_upper_eval(c,g,ga),mid,[Q('1e-22'),Q('1e-20')]);assert all(sg(v)==1 for v in r['evaluation']['checks'].values()),enc(r['evaluation']['checks'])
 write('base_upper_contact.json',{'status':'COMPUTER-ASSISTED PROOF','judgement':'verified','claim':'Unique local upper contact with positive gamma and V(u)=2lambda, not the k=0 entrance contact. Global joining is not asserted here.','F':'124/45','eta':'3/2','mpmath':mp.__version__,'iv_dps':iv.dps,'N':sk.N,'contact':r,'root_checks':c.root_checks})

def nn_eval(l,h):
 b12=fk.cost_box(12,l,h);b23=fk.cost_box(23,l,h);lam=I(l,h);g=b12['g']+b23['g']-1;k=b23['g']-1
 a=fk.hull(fk.root(13,lo(g),'a'),fk.root(13,hi(g),'a'));b=fk.hull(fk.root(13,lo(g),'b'),fk.root(13,hi(g),'b'))
 u,be,ga,ze=b12['a'],b12['b'],b23['a'],b23['b'];ea,eb,eu,ee,eg,ez=[fk.E(x) for x in [a,b,u,be,ga,ze]]
 Qc=fk.L(12,be,b12['g'])-fk.L(12,a,b12['g'])+fk.L(23,ga,b23['g'])-fk.L(23,a,b23['g'])
 H=fk.L(13,b,g)-fk.L(13,a,g)
 D12=ee-eu;B=ez-eg;A=ee-ea;C=ee-eg
 derivative=2*A/D12-2*(ea-eg)/B-1
 checks={'gamma>0':ga,'a>gamma':a-ga,'u>a':u-a,'zeta>u':ze-u,'b>zeta':b-ze,'beta>b':be-b,'Q>0':Qc,'H13>Q':H-Qc,'lambda>H13-Q':lam-H+Qc,'Fprime>0':derivative}
 return {'lambda':lam,'g':g,'g12':b12['g'],'k':k,'residual':Qc-lam,'derivative':derivative,'H13':H,'Q':Qc,'roots':dict(a=a,b=b,u=u,beta=be,gamma=ga,zeta=ze),'checks':checks}

def nn():
 l=Q('.0607725779352778027476');h=Q('.0607725779352778027477')
 left=nn_eval(l,l);right=nn_eval(h,h);assert sg(left['residual'])==-1 and sg(right['residual'])==1
 z=nn_eval(l,h);assert all(sg(v)==1 for v in z['checks'].values())
 write('nn_contact.json',{'status':'COMPUTER-ASSISTED PROOF','judgement':'verified','claim':'Unique NN/direct-1-to-3 contact in the stated bracket, with a local NN-optimal left neighbourhood and four-threshold-optimal right neighbourhood; first global transition requires a connected verification cover.','F':'2','eta':'3/2','mpmath':mp.__version__,'iv_dps':iv.dps,'N':fk.N,'bracket':[l,h],'left':left,'right':right,'evaluation':z,'root_checks':fk.CHECKS})

if __name__=='__main__':
 for name in sys.argv[1:] or ['base_lower','base_upper','nn']:
  try: globals()[name]()
  except Exception as e:
   import traceback
   write(name+'_FAILED.json',{'status':'FAILED','judgement':'verified','error':repr(e),'traceback':traceback.format_exc()});raise
