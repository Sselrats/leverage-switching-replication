"""CIR certificates beyond F=2; mpmath 1.3.0 interval arithmetic.
P/q and its integral: positive series with geometric tails.
E=integral_1^x 1/q: exponential series with an integral remainder bound.
The only special-function kernel added to the original is iv.gamma(F)
for q0=F**(F-1)/Gamma(F); no incomplete gamma or quadrature is trusted.
"""
from validated_pairs_F2 import R,lo,hi,ends,Fr,iv,m,c,BASE,M,rational_near
import json,sys,time

def I(l,h):return iv.mpf([R(Fr(l)).a,R(Fr(h)).b])
def hull(a,b):return c.hull(a,b)
class CIR:
 def __init__(self,Flo,Fhi,elo,ehi,N=100):
  self.F=I(Flo,Fhi);self.eta=I(elo,ehi);self.N=N
  self.q0=self.F**(self.F-1)/iv.gamma(self.F)
  self.B1=self.eta-R(Fr(1,2));self.theta_coeff=[iv.mpf(1)];self.exp_coeff=[iv.mpf(1)]
  for n in range(1,N+2):
   self.theta_coeff.append(self.theta_coeff[-1]*self.F/(self.F+n))
   self.exp_coeff.append(self.exp_coeff[-1]*self.F/n)
  self.Fmid=M((lo(self.F)+hi(self.F))/2);self.emid=M((lo(self.eta)+hi(self.eta))/2)
  self._ecache={};self._acache={}
 def q(self,x):return self.q0*x**self.F*iv.exp(-self.F*x)
 def f(self,x):return self.F*self.q(x)/x
 def theta(self,x):
  z=self.theta_coeff[self.N]
  for n in range(self.N-1,-1,-1):z=z*x+self.theta_coeff[n]
  ratio=self.F*x/(self.F+self.N+2);assert hi(ratio)<1
  tail=self.theta_coeff[self.N+1]*x**(self.N+1)/(1-ratio)
  return z+iv.mpf([0,tail.b])
 def P(self,x):return self.q(x)*self.theta(x)
 def Apoint(self,x):
  z=self.theta_coeff[self.N]/(self.N+1)
  for n in range(self.N-1,-1,-1):z=z*x+self.theta_coeff[n]/(n+1)
  z=z*x;ratio=self.F*x/(self.F+self.N+2);assert hi(ratio)<1
  tail=self.theta_coeff[self.N+1]*x**(self.N+2)/(self.N+2)/(1-ratio)
  return z+iv.mpf([0,tail.b])
 def A(self,x):
  mid=(lo(x)+hi(x))/2
  if mid not in self._acache:self._acache[mid]=self.Apoint(R(mid))
  return self._acache[mid]+(x-R(mid))*self.theta(x)
 def phi(self,r,x):
  L=iv.log(x)
  if lo(r)<=0<=hi(r) or max(abs(lo(r)),abs(hi(r)))<Fr(1,100):
   # (exp(r L)-1)/r = L sum (r L)^n/(n+1)!.
   z=r*L;term=iv.mpf(1);s=term
   for n in range(1,25):term=term*z/(n+1);s+=term
   a=max(abs(lo(z)),abs(hi(z)));assert a<Fr(26)
   nxt=abs(term)*R(a)/26;tail=nxt/(1-R(a)/27)
   return L*(s+iv.mpf([-tail.b,tail.b]))
  return (iv.exp(r*L)-1)/r
 def Epoint(self,x):
  z=iv.mpf(0)
  for n in range(self.N+1):z+=self.exp_coeff[n]*self.phi(n+1-self.F,x)
  tx=I(min(Fr(1),lo(x)),max(Fr(1),hi(x)))
  ratio=self.F*tx/(self.N+2);assert hi(ratio)<1
  bound=abs(x-1)*tx**(-self.F)*self.exp_coeff[self.N+1]*tx**(self.N+1)/(1-ratio)/self.q0
  return z/self.q0+iv.mpf([-bound.b,bound.b])
 def E(self,x):
  mid=(lo(x)+hi(x))/2
  if mid not in self._ecache:self._ecache[mid]=self.Epoint(R(mid))
  return self._ecache[mid]+(x-R(mid))/self.q(x)
 def G(self,p,x):
  if p==13:return self.B1+2*(self.eta-2)*self.P(x)+4*self.q(x)
  if p==12:return self.B1+(self.eta-R(Fr(3,2)))*self.P(x)+R(Fr(3,2))*self.q(x)
  return (self.eta-R(Fr(5,2)))*self.P(x)+R(Fr(5,2))*self.q(x)
 def L(self,p,x,g):
  if p==13:return (self.B1-g)*self.E(x)+2*(self.eta-2)*self.A(x)+4*x
  return (self.B1-g)*self.E(x)+(self.eta-R(Fr(3,2)))*self.A(x)+R(Fr(3,2))*x
 def K(self,x):return (self.eta-R(Fr(5,2)))*self.A(x)+R(Fr(5,2))*x
 def Gn(self,p,x):
  F=self.Fmid;e=self.emid;P=m.gammainc(F,0,F*x)/m.gamma(F);q=F**(F-1)/m.gamma(F)*x**F*m.exp(-F*x)
  if p==13:return e-M('.5')+2*(e-2)*P+4*q
  if p==12:return e-M('.5')+(e-M('1.5'))*P+M('1.5')*q
  return (e-M('2.5'))*P+M('2.5')*q
 def root(self,p,level,side):
  xc=self.eta/(2 if p==13 else R(Fr(3,2)) if p==12 else R(Fr(5,2)))
  left,right=(Fr(1,1000),lo(xc)) if side=='a' else (hi(xc),Fr(12))
  # The initial interval must work uniformly for F, eta and level.
  fun=lambda x:self.G(p,x)-level
  sl,sr=c.sign(fun(R(left))),c.sign(fun(R(right)))
  assert (sl,sr)==((-1,1) if side=='a' else (1,-1)),('no initial bracket',p,side,sl,sr)
  # Two one-sided sign searches: interval parameter uncertainty must not
  # terminate bisection at its first indeterminate midpoint.  The returned
  # endpoints retain opposite strict signs uniformly over every parameter.
  orient=-sl
  L,H=left,right
  for _ in range(58):
   mid=(L+H)/2
   if hi(orient*fun(R(mid)))<0:L=mid
   else:H=mid
  lower=L
  L,H=left,right
  for _ in range(58):
   mid=(L+H)/2
   if lo(orient*fun(R(mid)))>0:H=mid
   else:L=mid
  upper=H
  assert hi(orient*fun(R(lower)))<0 and lo(orient*fun(R(upper)))>0
  return I(lower,upper)
 def cost(self,g):
  a=self.root(13,R(g),'a');b=self.root(13,R(g),'b')
  return (self.L(13,b,R(g))-self.L(13,a,R(g)))/2
 def critical13(self):
  b=self.root(13,self.B1,'b')
  # E is anchored at 1; its coefficient vanishes at g=B1; A(0)=0.
  return self.L(13,b,self.B1)/2,b
 def box(self,gl,gh,gamma):
  g=I(gl,gh);ga=R(gamma);k=self.G(23,ga)
  a=self.root(13,g,'a');b=self.root(13,g,'b');u=self.root(12,g-k,'a');be=self.root(12,g-k,'b');ze=self.root(23,k,'b')
  clo,chi=self.cost(gh),self.cost(gl);lam=iv.mpf([clo.a,chi.b])
  va=self.L(12,be,g-k)-self.L(12,a,g-k);vu=self.L(12,be,g-k)-self.L(12,u,g-k)
  vg=va+self.K(ga)-self.K(a)-k*(self.E(ga)-self.E(a))
  tz=va-lam+self.K(ze)-self.K(a)-k*(self.E(ze)-self.E(a))
  lowmin=min(lo(vg-lam),Fr(0));highmax=max(hi(vu),hi(tz),hi(lam),Fr(0));margin=2*lo(lam)-(highmax-lowmin)
  order=[a-ga,u-a,ze-u,b-ze,be-b,self.eta/R(Fr(5,2))-a,be-self.eta/R(Fr(3,2))]
  passed=all(lo(v)>0 for v in order) and margin>0
  return {'g_lo':str(gl),'g_hi':str(gh),'gamma':str(gamma),'root_order_margins':[ends(v) for v in order],'roots':{k:ends(v) for k,v in {'a':a,'b':b,'u':u,'beta':be,'zeta':ze}.items()},'q0':ends(self.q0),'lambda':ends(lam),'Vgamma':ends(vg),'Vu':ends(vu),'Tzeta':ends(tz),'margin_lower':str(margin),'margin_lower_dec':c.dec(margin),'passed':passed}

def single(F='124/45',g='1.3',gamma='0.1'):
 C=CIR(Fr(F),Fr(F),Fr('1.5'),Fr('1.5'))
 r=C.box(Fr(g),Fr(g),Fr(gamma));print(json.dumps(r,indent=1))
if __name__=='__main__':single(*sys.argv[1:])
