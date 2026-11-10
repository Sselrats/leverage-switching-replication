"""Numerical proposals only: full-menu four-threshold branch, fees zero."""
from pathlib import Path
import sys, math, json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import root

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/cir_premium_bregman_20260927'))
from check_derivations import Pair

COST=-math.log1p(-.0018)
LAM=3.1*COST/.04

class Full:
    def __init__(self,eta):
        self.eta=eta
        self.pairs={p:Pair(.04*eta,*ij,COST) for p,ij in [(12,(1,2)),(13,(1,3)),(23,(2,3))]}
        self.sol={p:v.solve() for p,v in self.pairs.items()}
        self.q=self.pairs[13].q; self.P=self.pairs[13].P
        self.B2=2*eta-2
    def G(self,p,x):return self.pairs[p].G(x)
    def K(self,x):return self.G(23,x)-self.B2
    def intL(self,p,lo,hi,g):
        return quad(lambda x:(self.G(p,x)-g)/self.q(x),lo,hi,epsabs=1e-10)[0]
    def intK(self,lo,hi,k):
        return quad(lambda x:(self.K(x)-k)/self.q(x),lo,hi,epsabs=1e-10)[0]
    def E(self,x):return quad(lambda z:1/self.q(z),1,x,epsabs=1e-10)[0]
    def A(self,x):return quad(lambda z:self.P(z)/self.q(z),1,x,epsabs=1e-10)[0]
    def solve(self):
        gamma,zeta=(self.sol[23][k] for k in ['u_dim','w_dim'])
        k=self.sol[23]['g_dim']-self.B2
        def equations(v):
            a,be=np.exp(v[:2]); g=v[2]
            return [self.G(13,a)-g,self.G(12,be)-(g-k),
                    self.intL(12,a,be,g-k)+self.intK(a,gamma,k)-LAM]
        init=[math.log(self.sol[13]['u_dim']),math.log(self.sol[12]['w_dim']),self.sol[13]['g_dim']]
        sol=root(equations,init,tol=1e-11)
        assert sol.success and max(abs(x) for x in equations(sol.x))<1e-8,sol
        a,be=np.exp(sol.x[:2]); g=sol.x[2]
        aa,b=self.pairs[13].band(g); u,bebe=self.pairs[12].band(g-k)
        assert abs(a-aa)+abs(be-bebe)<1e-8
        H13=self.intL(13,a,b,g); Vu=self.intL(12,u,be,g-k)
        orders=[a-gamma,u-a,zeta-u,b-zeta,be-b,2*self.eta/5-gamma,be-2*self.eta/3]
        margins={'H13_above_lambda':H13-LAM,'2lambda_above_H13':2*LAM-H13,'2lambda_above_Vu':2*LAM-Vu}
        p=(self.E(zeta)-self.E(gamma))/(self.E(be)-self.E(gamma))
        Te=(1-p)*self.A(gamma)+p*self.A(be)-self.A(zeta)
        Tl=self.A(zeta)-self.A(gamma)
        Tret=self.A(a)-self.E(a)-self.A(be)+self.E(be)
        Tup=self.A(zeta)-self.A(a)
        T=Te+(1-p)*Tl+p*(Tret+Tup)
        ell=(2*Te+3*(1-p)*Tl+p*(Tret+3*Tup))/T
        return {'eta':self.eta,'gain_dim':g,'gain':.04*g,'ell':ell,
                'roots':dict(gamma=gamma,a=a,u=u,zeta=zeta,b=b,beta=be),
                'order_margins':orders,'margins':margins,
                'passes_numeric_sufficiency':min(orders+list(margins.values()))>0,
                'pair13_gain':self.sol[13]['g'],'pair12_gain':self.sol[12]['g'],
                'deletion_loss':.04*g-self.sol[13]['g']}

    def solve_upper(self):
        # The 12 contact band is calibrated; 3 returns directly to 1 at b.
        u,beta=(self.sol[12][key] for key in ['u_dim','w_dim'])
        g12=self.sol[12]['g_dim']
        def equations(v):
            a,b,gamma=np.exp(v[:3]); g=v[3]; k=g-g12
            Q=self.intL(12,a,beta,g12)+self.intK(a,gamma,k)
            H13=self.intL(13,a,b,g)
            return [self.G(13,a)-g,self.G(13,b)-g,self.K(gamma)-k,H13-Q-LAM]
        init=[math.log(self.sol[13]['u_dim']),math.log(self.sol[13]['w_dim']),
              math.log(self.sol[23]['u_dim']),self.sol[13]['g_dim']+.03]
        sol=root(equations,init,tol=1e-11)
        assert sol.success and max(abs(x) for x in equations(sol.x))<1e-8,sol
        a,b,gamma=np.exp(sol.x[:3]); g=sol.x[3]; k=g-g12
        aa,bb=self.pairs[13].band(g)
        gg,zeta=self.pairs[23].band(k+self.B2)
        assert abs(a-aa)+abs(b-bb)+abs(gamma-gg)<1e-8
        H13=self.intL(13,a,b,g)
        area23=self.intK(gamma,zeta,k)
        margins={'H13_above_lambda':H13-LAM,'2lambda_above_H13':2*LAM-H13,
                 '2lambda_above_area23':2*LAM-area23}
        orders=[a-gamma,u-a,zeta-u,b-zeta,beta-b,2*self.eta/5-gamma,beta-2*self.eta/3]
        p=(self.E(beta)-self.E(u))/(self.E(beta)-self.E(gamma))
        Te=(1-p)*self.A(beta)+p*self.A(gamma)-self.A(u)
        Tupper=self.A(u)-self.E(u)-self.A(beta)+self.E(beta)
        Tlow3=self.A(b)-self.A(gamma)
        Tlow1=self.A(u)-self.E(u)-self.A(b)+self.E(b)
        T=Te+(1-p)*Tupper+p*(Tlow3+Tlow1)
        ell=(2*Te+(1-p)*Tupper+p*(3*Tlow3+Tlow1))/T
        return {'eta':self.eta,'branch':'12-calibrated','gain_dim':g,'gain':.04*g,'ell':ell,
                'roots':dict(gamma=gamma,a=a,u=u,zeta=zeta,b=b,beta=beta),
                'order_margins':orders,'margins':margins,
                'passes_numeric_sufficiency':min(orders+list(margins.values()))>0,
                'pair13_gain':self.sol[13]['g'],'deletion_loss':.04*g-self.sol[13]['g']}

    def reduced_test(self,gamma):
        g=self.sol[13]['g_dim']; a,b=(self.sol[13][key] for key in ['u_dim','w_dim'])
        k=self.K(gamma)
        u,beta=self.pairs[12].band(g-k)
        _,zeta=self.pairs[23].band(k+self.B2)
        va=self.intL(12,a,beta,g-k); vu=self.intL(12,u,beta,g-k)
        vg=va+self.intK(a,gamma,k)
        tz=va-LAM+self.intK(a,zeta,k)
        m=min(vg-LAM,0); mm=max(vu,tz,LAM)
        return {'gamma':gamma,'margin':2*LAM-(mm-m),
                'orders':[gamma,a-gamma,u-a,zeta-u,b-zeta,beta-b],
                'g':g,'vu':vu,'vg':vg,'tz':tz}

    def solve_ladder(self):
        gamma,zeta=(self.sol[23][key] for key in ['u_dim','w_dim'])
        u,beta=(self.sol[12][key] for key in ['u_dim','w_dim'])
        k=self.sol[23]['g_dim']-self.B2
        g=self.sol[12]['g_dim']+k
        a,b=self.pairs[13].band(g)
        Q=self.intL(12,a,beta,g-k)+self.intK(a,gamma,k)
        H13=self.intL(13,a,b,g)
        margins={'Q':Q,'lambda_minus_Q':LAM-Q,'H13_minus_Q':H13-Q,
                 'lambda_minus_H13_plus_Q':LAM-H13+Q}
        orders=[a-gamma,u-a,zeta-u,b-zeta,beta-b,2*self.eta/5-gamma,beta-2*self.eta/3]
        ell=self.sol[12]['ell']+self.sol[23]['ell']-2
        return {'eta':self.eta,'branch':'adjacent ladder','gain_dim':g,'gain':.04*g,'ell':ell,
                'roots':dict(gamma=gamma,a=a,u=u,zeta=zeta,b=b,beta=beta),
                'order_margins':orders,'margins':margins,
                'passes_numeric_sufficiency':min(orders+list(margins.values()))>0,
                'pair13_gain':self.sol[13]['g'],'deletion_loss':.04*g-self.sol[13]['g']}

def main():
    records=[]
    for eta in [.95,1.4999,1.5,1.5001,2.05]:
        try:r=Full(eta).solve()
        except Exception as e:r={'eta':eta,'error':str(e)}
        print(json.dumps(r),flush=True); records.append(r)
    r=Full(2.05).solve_upper(); print(json.dumps(r),flush=True); records.append(r)
    r=Full(2.05).solve_ladder(); print(json.dumps(r),flush=True); records.append(r)
    low=Full(.95)
    tests=[]
    for gamma in [.02,.03,.04,.05,.06,.08,.1,.12,.14]:
        try:tests.append(low.reduced_test(gamma))
        except Exception:pass
    print(json.dumps(tests),flush=True)
    (ROOT/'artifacts/cir_fullmenu_endpoints_20260927/exploration.json').write_text(
        json.dumps({'status':'diagnostic only','proportional_cost':.0018,'log_cost':COST,
                    'lambda':LAM,'records':records,'lower_endpoint_tests':tests},indent=2)+'\n')

if __name__=='__main__':main()
