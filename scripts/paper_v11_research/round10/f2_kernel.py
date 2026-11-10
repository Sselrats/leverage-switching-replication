#!/usr/bin/env python3
"""Independent elementary F=2 kernel for pairs 12,13,23 at eta=3/2.
mp supplies proposals only. Accepted root brackets use exact rational endpoints,
strict iv signs and analytically specified monotone branches. Every E evaluation
uses a positive series and a geometric tail, not mp.ei or quadrature.
"""
from fractions import Fraction as Q
from functools import lru_cache
import mpmath as mp
assert mp.__version__ == '1.3.0'
iv=mp.iv; iv.dps=70; mp.mp.dps=90
N=140
CHECKS=[]
def M(x):
    if isinstance(x,Q): return mp.mpf(x.numerator)/x.denominator
    return mp.mpf(x)
def R(x):
    x=Q(x); return iv.mpf(x.numerator)/iv.mpf(x.denominator)
def frac(t):
    s,m,e,b=t
    return Q((-1 if s else 1)*m)*Q(2)**e
def lo(x): return frac(x._mpi_[0])
def hi(x): return frac(x._mpi_[1])
def I(a,b):return iv.mpf([R(a).a,R(b).b])
def hull(*xs):return I(min(lo(x) for x in xs),max(hi(x) for x in xs))
def sign(x):return 1 if lo(x)>0 else -1 if hi(x)<0 else 0
def near(x,n=38,pad=2):
    d=10**n; k=int(mp.floor(x*d)); return Q(k-pad,d),Q(k+1+pad,d)
def dec(x,up=False,digits=26):
    z=Q(x)*10**digits;n=z.numerator//z.denominator
    if up and z.denominator!=1:n+=1
    neg=n<0;ss=str(abs(n)).zfill(digits+1)
    return ('-' if neg else '')+ss[:-digits]+'.'+ss[-digits:]
def ends(x):return {'lo':str(lo(x)),'hi':str(hi(x)),'dec':[dec(lo(x)),dec(hi(x),True)]}
def serial(x):
    if hasattr(x,'_mpi_'):return ends(x)
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def q(x):return 2*x*x*iv.exp(-2*x)
def f(x):return 4*x*iv.exp(-2*x)
def P(x):return 1-iv.exp(-2*x)*(1+2*x)
def qp(x):return (1-x)*f(x)
def Upsilon(x):return 1/(2*x)-iv.log(x)
@lru_cache(None)
def Epoint(t):
    x=R(t);z=2*x;term=iv.mpf(1);s=iv.mpf(0)
    for n in range(1,N+1):
        term*=z/n;s+=term/n
    nxt=term*z/(N+1)/(N+1);r=z/(N+2)
    assert hi(r)<1
    tail=nxt/(1-r)
    return iv.log(x)+s-iv.exp(z)/(2*x)+I(0,hi(tail))
def E(x):
    c=(lo(x)+hi(x))/2
    return Epoint(c)+(x-R(c))/q(x)
def levels(pair):return pair//10,pair%10
def xc(pair):
    i,j=levels(pair);return Q(3,i+j)
def B(i):return Q(3*i-i*i,2)
def G(pair,x):
    i,j=levels(pair);d=j-i;m=Q(i+j,2)
    return R(B(i))+R(d*(Q(3,2)-m))*P(x)+R(d*m)*q(x)
def Gp(pair,x):
    i,j=levels(pair);return R(j-i)*(R(Q(3,2))-R(Q(i+j,2))*x)*f(x)
def Gpp(pair,x):
    i,j=levels(pair);d=j-i;m=Q(i+j,2)
    return R(d)*(-R(m)*f(x)+(R(Q(3,2))-R(m)*x)*4*(1-2*x)*iv.exp(-2*x))
def L(pair,x,g):
    i,j=levels(pair);d=j-i;m=Q(i+j,2)
    return (R(B(j))-g)*E(x)+R(d*(Q(3,2)-m))*Upsilon(x)+R(d*m)*x

def qn(x):return 2*x*x*mp.exp(-2*x)
def pn(x):return 1-mp.exp(-2*x)*(1+2*x)
def En(x):return mp.ei(2*x)-mp.exp(2*x)/(2*x)
def Gn(pair,x):
    i,j=levels(pair);d=j-i;m=M(Q(i+j,2))
    return M(B(i))+d*(M('1.5')-m)*pn(x)+d*m*qn(x)
def Ln(pair,x,g):
    i,j=levels(pair);d=j-i;m=M(Q(i+j,2))
    return (M(B(j))-g)*En(x)+d*(M('1.5')-m)*(1/(2*x)-mp.log(x))+d*m*x
@lru_cache(None)
def root(pair,g,side):
    xx=M(xc(pair));gm=M(g)
    l,h=(M('1e-50'),xx) if side=='a' else (xx,M(20))
    assert (Gn(pair,l)-gm)*(Gn(pair,h)-gm)<0
    for _ in range(95):
        z=(l+h)/2
        if (Gn(pair,z)<gm)==(side=='a'):l=z
        else:h=z
    prop=mp.findroot(lambda z:Gn(pair,z)-gm,((l+h)/2,(l+h)/2*(1+M('1e-8'))),tol=M('1e-78'))
    l,h=near(prop)
    fl,fh=G(pair,R(l))-R(g),G(pair,R(h))-R(g)
    assert (0<l<h<xc(pair)) if side=='a' else (xc(pair)<l<h)
    assert (sign(fl),sign(fh))==((-1,1) if side=='a' else (1,-1))
    ans=I(l,h);CHECKS.append({'kind':'state_root','pair':pair,'g':g,'side':side,'bracket':ans,'left_residual':fl,'right_residual':fh})
    return ans
@lru_cache(None)
def at_gain(pair,g):
    a=root(pair,g,'a');b=root(pair,g,'b');J=E(b)-E(a)
    return {'a':a,'b':b,'g':R(g),'H':L(pair,b,R(g))-L(pair,a,R(g)),'J':J}
@lru_cache(None)
def proposal_band(pair,lam):
    x=M(xc(pair));l=M(lam);i,j=levels(pair);curv=(j-i)*M(Q(i+j,2))*4*x*mp.exp(-2*x)
    width=(3*l*qn(x)/curv)**(M(1)/3)
    a=max(x/20,x-width);b=x+width;g=Gn(pair,a)
    fun=lambda aa,bb,gg:(Gn(pair,aa)-gg,Gn(pair,bb)-gg,Ln(pair,bb,gg)-Ln(pair,aa,gg)-2*l)
    sol=mp.findroot(fun,(a,b,g),tol=M('1e-75'),maxsteps=100)
    aa,bb,gg=sol
    assert 0<aa<x<bb<20 and gg>max(M(B(i)),M(B(j)))
    return aa,bb,gg
@lru_cache(None)
def band(pair,lam):
    lam=Q(lam);aa,bb,gg=proposal_band(pair,lam)
    gl,gh=near(gg,27,3)
    rlo=at_gain(pair,gl);rhi=at_gain(pair,gh)
    fl=rlo['H']-2*R(lam);fh=rhi['H']-2*R(lam)
    assert lo(fl)>0 and hi(fh)<0,(pair,lam,serial(fl),serial(fh))
    a=hull(rlo['a'],rhi['a']);b=hull(rlo['b'],rhi['b']);J=E(b)-E(a)
    ans={'lambda':lam,'pair':pair,'a':a,'b':b,'g':I(gl,gh),'J':J,'s':2/J,'gain_left_residual':fl,'gain_right_residual':fh}
    CHECKS.append({'kind':'gain_root',**ans})
    return ans

def cost_box(pair,l,h):
    l=Q(l);h=Q(h);r0=band(pair,l);r1=band(pair,h)
    return {'lambda':I(l,h),'a':hull(r0['a'],r1['a']),'b':hull(r0['b'],r1['b']),'g':I(lo(r1['g']),hi(r0['g'])),'J':I(lo(r0['J']),hi(r1['J']))}
