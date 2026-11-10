"""Validated F=2 pair values. mpmath.mp is used ONLY for root proposals.
All accepted brackets, signs, primitives and reported enclosures use
mpmath 1.3.0 iv arithmetic and exact rational inputs. No float is evidence.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as Fr
import sys, json
# The F=2 elementary interval kernel is the paper's certificate program (scripts/certificate_zero_interval_F2.py).
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
_saved=sys.argv[:];sys.argv=[sys.argv[0],'3/2']
import certificate_zero_interval_F2 as c
sys.argv=_saved
m=c.mm; iv=c.iv;m.mp.dps=85
assert m.__version__=='1.3.0'
R=c.R; lo=c.lo; hi=c.hi; ends=c.ends
_cache_num={};_cache_iv={}

def M(x):
    if isinstance(x,Fr):return m.mpf(x.numerator)/x.denominator
    return m.mpf(x)

def qn(x):return 2*x*x*m.exp(-2*x)
def pn(x):return 1-m.exp(-2*x)*(1+2*x)
def En(x):return m.ei(2*x)-m.exp(2*x)/(2*x) # additive constant irrelevant

def Gn(pair,x):
    if pair==12:return 1+3*x*x*m.exp(-2*x)
    return m.exp(-2*x)*(1+2*x+8*x*x)

def Ln(pair,x,g):
    if pair==12:return (1-g)*En(x)+M('1.5')*x
    return -g*En(x)-1/(2*x)+m.log(x)+4*x

def Gp(pair,x):
    f=4*x*iv.exp(-2*x)
    return ((R(Fr(3,2))-R(Fr(3,2))*x) if pair==12 else (3-4*x))*f

def Gi(pair,x):return c.G12(x) if pair==12 else c.G13(x)
def Li(pair,x,g):return c.L12(x,g) if pair==12 else c.Lg(x,g)
def xc(pair):return Fr(1) if pair==12 else Fr(3,4)

def numerical_band(pair,lam):
    key=(pair,Fr(lam))
    if key in _cache_num:return _cache_num[key]
    l=M(Fr(lam));x=M(xc(pair));dmc=M('1.5') if pair==12 else M(4)
    curvature=dmc*(4*x*m.exp(-2*x));width=(3*l*qn(x)/curvature)**(M(1)/3)
    a=max(x/M(20),x-width);b=x+width;g=Gn(pair,a)
    prev=[(abs(m.log(M(k[1])/l)),v) for k,v in _cache_num.items() if k[0]==pair]
    if prev and min(prev,key=lambda p:p[0])[0]<M('0.8'):
        _,(aa,bb,gg)=min(prev,key=lambda p:p[0]);a,b,g=aa,bb,gg
    def fun(aa,bb,gg):return Gn(pair,aa)-gg,Gn(pair,bb)-gg,Ln(pair,bb,gg)-Ln(pair,aa,gg)-2*l
    try:
        sol=m.findroot(fun,(a,b,g),tol=M('1e-65'),maxsteps=70)
        a,b,g=sol
        if not (0<a<x<b<12 and g>1):raise ValueError('wrong branch')
    except Exception:
        # A bounded scalar fallback supplies proposals, not proof.
        def point(g):
            al,ah=M('1e-40'),x
            bl,bh=x,M(12)
            for _ in range(100):
                z=(al+ah)/2
                if Gn(pair,z)<g:al=z
                else:ah=z
                z=(bl+bh)/2
                if Gn(pair,z)>g:bl=z
                else:bh=z
            aa=(al+ah)/2;bb=(bl+bh)/2
            return aa,bb,Ln(pair,bb,g)-Ln(pair,aa,g)
        gl,gh=M(1)+M('1e-50'),Gn(pair,x)-M('1e-50')
        for _ in range(80):
            gg=(gl+gh)/2;aa,bb,H=point(gg)
            if H>2*l:gl=gg
            else:gh=gg
        sol=m.findroot(fun,(aa,bb,gg),tol=M('1e-65'),maxsteps=30);a,b,g=sol
    _cache_num[key]=(a,b,g);return a,b,g

def rational_near(x,places=35):
    den=10**places
    n=int(m.floor(x*den));return Fr(n,den),Fr(n+1,den)

def root_bracket(pair,g,side,seed=None):
    gm=M(g);x=M(xc(pair))
    if seed is None and pair==13:
        aa,bb=(M('1e-40'),x) if side=='a' else (x,M(12))
        for _ in range(90):
            zz=(aa+bb)/2
            if (Gn(13,zz)<gm)==(side=='a'):aa=zz
            else:bb=zz
        seed=(aa+bb)/2
    if seed is None: seed=x*(M('0.6') if side=='a' else M('1.7'))
    # mp proposal on the specified monotone branch
    if pair==12:
        y=-m.sqrt((gm-1)/3)
        root=-m.lambertw(y,0 if side=='a' else -1).real
    else:
        root=m.findroot(lambda z:Gn(pair,z)-gm,(seed,seed*(1+M('1e-6'))),tol=M('1e-65'),maxsteps=50)
    a,b=rational_near(root,35);a-=Fr(1,10**35);b+=Fr(1,10**35)
    assert (0<a<b<xc(pair)) if side=='a' else (xc(pair)<a<b)
    sa,sb=c.sign(Gi(pair,R(a))-R(g)),c.sign(Gi(pair,R(b))-R(g))
    assert (sa,sb)==((-1,1) if side=='a' else (1,-1)),(pair,g,side,'root signs',sa,sb)
    return iv.mpf([R(a).a,R(b).b])

def at_gain(pair,g,seed=None):
    a=root_bracket(pair,g,'a',None if seed is None else seed[0])
    b=root_bracket(pair,g,'b',None if seed is None else seed[1])
    H=Li(pair,b,R(g))-Li(pair,a,R(g));J=c.E(b)-c.E(a)
    Jg=1/(c.q(b)*Gp(pair,b))-1/(c.q(a)*Gp(pair,a))
    return {'a':a,'b':b,'H':H,'J':J,'Jg':Jg}

def validated_band(pair,lam):
    lam=Fr(lam);key=(pair,lam)
    if key in _cache_iv:return _cache_iv[key]
    a,b,g=numerical_band(pair,lam)
    gl,gh=rational_near(g,24);gl-=Fr(2,10**24);gh+=Fr(2,10**24)
    low=at_gain(pair,gl,(a,b));high=at_gain(pair,gh,(a,b))
    assert lo(low['H']-2*R(lam))>0 and hi(high['H']-2*R(lam))<0,('gain signs',pair,lam)
    aa=c.hull(low['a'],high['a']);bb=c.hull(low['b'],high['b'])
    J=c.E(bb)-c.E(aa)
    Jg=1/(c.q(bb)*Gp(pair,bb))-1/(c.q(aa)*Gp(pair,aa))
    obj={'lambda':str(lam),'pair':pair,'g':iv.mpf([R(gl).a,R(gh).b]),'a':aa,'b':bb,'J':J,'Jg':Jg,
         'H_glo_minus_2lambda':low['H']-2*R(lam),'H_ghi_minus_2lambda':high['H']-2*R(lam)}
    _cache_iv[key]=obj;return obj

def serial(row):return {k:(ends(v) if hasattr(v,'_mpi_') else v) for k,v in row.items()}

def cost_box(l,h):
    low={p:validated_band(p,l) for p in (12,13)}
    high={p:validated_band(p,h) for p in (12,13)}
    out={}
    for p in (12,13):
        out[p]={
          'g':iv.mpf([high[p]['g'].a,low[p]['g'].b]),
          'a':c.hull(low[p]['a'],high[p]['a']),
          'b':c.hull(low[p]['b'],high[p]['b']),
          'J':iv.mpf([low[p]['J'].a,high[p]['J'].b])}
        a,b=out[p]['a'],out[p]['b']
        out[p]['Jg']=1/(c.q(b)*Gp(p,b))-1/(c.q(a)*Gp(p,a))
    return out

if __name__=='__main__':
    for lam in ['0.0000001','0.01','0.81','1.4734','1.5','1.56353865713912']:
        for p in (12,13):
            if p==13 and Fr(lam)>Fr('1.56'):continue
            r=validated_band(p,Fr(lam));print(lam,p,ends(r['g'])['dec'],ends(2/r['J'])['dec'],flush=True)
