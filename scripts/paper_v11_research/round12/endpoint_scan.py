#!/usr/bin/env python3
"""Tenth round: EXPLORATORY endpoint scan, not an interval certificate (numerical helper imported by general_family_check.py).
General-F, no-fee CIR menu {1,2,3}. All decimal outputs are numerical proposals.
Usage: python audit_work/endpoint_scan.py [--out records/endpoint_scan.json]
Run audit_work/general_family_check.py for strict interval acceptance.
"""
from __future__ import annotations
import argparse, json, sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import mpmath as mp
mp.mp.dps=40

def m(x):
    if isinstance(x,Fraction):return mp.mpf(x.numerator)/x.denominator
    return mp.mpf(x)

class NumericCIR:
    def __init__(self,F,eta='1.5'):
        self.F=m(F);self.eta=m(eta);self.q0=self.F**(self.F-1)/mp.gamma(self.F)
        self.seeds={}
    def B(self,i):return i*self.eta-m(i*i)/2
    def q(self,x):return self.q0*x**self.F*mp.exp(-self.F*x)
    def P(self,x):return mp.gammainc(self.F,0,self.F*x)/mp.gamma(self.F)
    def G(self,p,x):
        if p==23:return (self.eta-m('2.5'))*self.P(x)+m('2.5')*self.q(x) # K_23, NOT G_23
        i,j=divmod(p,10);c=m(j*j-i*i)/2
        return self.B(i)+((j-i)*self.eta-c)*self.P(x)+c*self.q(x)
    @lru_cache(maxsize=30000)
    def E(self,x):
        if self.F==2:return mp.ei(2*x)-mp.exp(2*x)/(2*x) # additive constant immaterial
        t=m(1);s=m(0)
        for n in range(181):
            p=n+1-self.F
            s+=t*(mp.log(x) if p==0 else mp.expm1(p*mp.log(x))/p)
            t*=self.F/(n+1)
        return s/self.q0
    @lru_cache(maxsize=30000)
    def A(self,x):
        if not x:return m(0)
        t=m(1);s=x
        for n in range(1,181):
            t*=self.F*x/(self.F+n);s+=t*x/(n+1)
        return s
    def K(self,x):return (self.eta-m('2.5'))*self.A(x)+m('2.5')*x
    def L(self,p,x,g):
        i,j=divmod(p,10);c=m(j*j-i*i)/2
        return (self.B(i)-g)*self.E(x)+((j-i)*self.eta-c)*self.A(x)+c*x
    def root(self,p,h,side):
        xc=self.eta/(m('2.5') if p==23 else m(sum(divmod(p,10)))/2)
        l,r=(m('1e-25'),xc) if side=='a' else (xc,max(m(12),4*xc))
        fun=lambda x:self.G(p,x)-h
        if fun(l)*fun(r)>=0:raise ValueError(f'No branch root {p,side,mp.nstr(h,12)}')
        # Safeguarded Newton: do not let special functions see out-of-branch iterates.
        x=self.seeds.get((p,side),(l+r)/2)
        if not l<x<r:x=(l+r)/2
        for _ in range(180):
            fx=fun(x)
            if (fx<0)==(side=='a'):l=x
            else:r=x
            if abs(fx)<m('1e-34') or r-l<m('1e-32'):break
            f=self.F*self.q(x)/x
            slope=((self.eta-m('2.5')*x) if p==23 else ((self.eta-m('1.5')*x) if p==12 else 2*(self.eta-2*x)))*f
            xn=x-fx/slope if slope else (l+r)/2
            x=xn if l<xn<r else (l+r)/2
        self.seeds[p,side]=x;return x
    def pair(self,p,g):
        a,b=self.root(p,g,'a'),self.root(p,g,'b');J=self.E(b)-self.E(a)
        return dict(a=a,b=b,g=g,H=self.L(p,b,g)-self.L(p,a,g),J=J,s=2/J)
    def least(self,g):
        rr=self.pair(13,g);a,b,lam=rr['a'],rr['b'],rr['H']/2
        peak=self.G(12,self.eta/m('1.5'));t0=None
        def at_ga(ga):
            kk=self.G(23,ga) if ga else m(0);h=g-kk
            be=self.root(12,h,'b')
            Q=self.L(12,be,h)+(self.K(ga)-kk*self.E(ga) if ga else m(0))-self.L(13,a,g)
            return Q-lam,be,kk
        if self.B(2)<g<peak:
            t0,be,kk=at_ga(m(0))
        if t0 is not None and t0>=0:
            ga=m(0)
        else:
            xc=self.eta/m('2.5')
            lo=m('1e-18') if g<peak else self.root(23,g-peak,'a')*(1+m('1e-15'))
            # h must remain above B2 on the decreasing G12 branch.
            hi=xc*(1-m('1e-6'))
            if self.G(23,hi)>=g-self.B(2):hi=self.root(23,g-self.B(2),'a')*(1-m('1e-6'))
            fl=at_ga(lo)[0];fh=at_ga(hi)[0]
            if not fl<0<fh:raise ValueError('No positive-gamma chord bracket')
            ga=(lo+hi)/2
            for _ in range(110):
                v,be,kk=at_ga(ga)
                if abs(v)<m('1e-29'):break
                if v<0:lo=ga
                else:hi=ga
                kp=(self.eta-m('2.5')*ga)*self.F*self.q(ga)/ga
                dv=kp*(self.E(be)-self.E(ga))
                gn=ga-v/dv
                ga=gn if lo<gn<hi else (lo+hi)/2
            _,be,kk=at_ga(ga)
        h=g-kk
        def d(x):
            if x<=a:return lam
            if x>=b:return -lam
            return lam-self.L(13,x,g)+self.L(13,a,g)
        def V(x):
            if x>=be:return m(0)
            if ga and x<=ga:return d(x)
            va=self.L(12,be,h)-self.L(12,a,h)
            if x>=a:return self.L(12,be,h)-self.L(12,x,h)
            return va+self.K(x)-self.K(a)-(kk*(self.E(x)-self.E(a)) if kk else 0)
        roots={'gamma':ga,'a':a,'b':b,'beta':be};margins={}
        try:
            u=self.root(12,h,'a');roots['u']=u
            if max(a,ga)<u<be:margins['M_u']=2*lam-V(u)
        except ValueError:pass
        try:
            ze=self.root(23,kk,'b');roots['zeta']=ze
            if ga<ze<min(b,be):margins['M_zeta']=2*lam-(V(ze)-d(ze))
        except ValueError:pass
        for name,x in [('a',a),('b',b)]:
            if ga<x<be:
                margins['M_V_'+name]=2*lam-V(x);margins['M_T_'+name]=2*lam-(V(x)-d(x))
        if not ga:
            margins['M_V_0']=2*lam-V(m(0));margins['M_T_0']=2*lam-(V(m(0))-lam)
        derivatives={}
        if ga>0 and ga<a and be>b:
            J=rr['J'];AA=self.E(be)-self.E(a);CC=self.E(be)-self.E(ga)
            kd=(1-2*AA/J)/CC;hd=-2/J-kd
            derivatives['dk_dlambda']=kd;derivatives['dh_dlambda']=hd
            if 'M_u' in margins:
                derivatives['dM_u_dlambda']=2+(self.E(be)-self.E(roots['u']))*hd
            if 'M_zeta' in margins:
                derivatives['dM_zeta_dlambda']=2+(self.E(roots['zeta'])-self.E(ga))*kd
        return dict(g=g,lambda_=lam,roots=roots,k=kk,T0_entrance_candidate=t0,margins=margins,derivatives=derivatives,canonical_feasible_numerically=min(margins.values(),default=0)>=0)

def serial(z):
    if isinstance(z,dict):return {k:serial(v) for k,v in z.items()}
    if isinstance(z,(list,tuple)):return [serial(v) for v in z]
    if isinstance(z,mp.mpf):return mp.nstr(z,32)
    return z

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'records/endpoint_scan.json');args=p.parse_args()
    rows=[]
    for F,gains in [(Fraction(2),['1.44534','1.44532','1.44','1.4','1.3','1.2','1.14','1.132','1.1318']), (Fraction(124,45),['1.28','1.275','1.26','1.24','1.2','1.17','1.16','1.157','1.156','1.1555','1.1551','1.15502'])]:
        C=NumericCIR(F)
        for gs in gains:
            try:
                z=C.least(m(gs));z['F']=str(F);z['eta']='3/2';rows.append(z)
                print('POINT',str(F),gs,mp.nstr(z['lambda_'],14),'gamma',mp.nstr(z['roots']['gamma'],16),{k:mp.nstr(v,10) for k,v in z['margins'].items()},flush=True)
            except Exception as e:
                rows.append({'F':str(F),'g':gs,'error':repr(e)});print('FAILED POINT',F,gs,repr(e),flush=True)
    record={'status':'EXPLORATORY ONLY','judgement':'plausible but unverified','mpmath':mp.__version__,'mp_dps':mp.mp.dps,'scope':'Numerical proposals, not interval acceptance or a phase diagram. p=23 denotes K_23, not G_23.','rows':rows}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(serial(record),indent=2)+'\n')
if __name__=='__main__':main()
