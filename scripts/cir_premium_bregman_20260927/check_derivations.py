"""Fixed algebra, convex counterexamples and deterministic CIR pair diagnostics.

No market observations, Monte Carlo, or certified floating-point quadrature.
Run: uv run --no-project --with sympy==1.14.0 --with scipy==1.17.1 python ...
"""
from pathlib import Path
import hashlib
import json
import math
import warnings

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gammainc, gammaln
from scipy.stats import norm
import sympy as s

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/cir_premium_bregman_20260927"


class Pair:
    def __init__(self, b, i=0, j=3, d=0.001):
        self.theta, self.kappa, self.xi2 = 0.04, 3.1, 0.09
        self.F = 2*self.kappa*self.theta/self.xi2
        self.eta, self.i, self.j = b/self.theta, i, j
        self.mid = (i+j)/2
        self.lam = self.kappa*d/self.theta
        self.Bi = i*self.eta-i*i/2
        self.Bj = j*self.eta-j*j/2

    def q(self, x):
        return math.exp((self.F-1)*math.log(self.F)-gammaln(self.F)
                        +self.F*math.log(x)-self.F*x)

    def P(self, x):
        return gammainc(self.F, self.F*x)

    def G(self, x):
        return self.Bi+(self.j-self.i)*(
            (self.eta-self.mid)*self.P(x)+self.mid*self.q(x))

    def A(self, x):
        return self.i+(self.j-self.i)*self.P(x)

    def Gx(self, x):
        return self.F*self.q(x)/x*(self.j-self.i)*(self.eta-self.mid*x)

    def band(self, g):
        xc = self.eta/self.mid
        u = brentq(lambda x:self.G(x)-g, 1e-12, xc, xtol=1e-14)
        right = max(2*xc, 2.)
        while self.G(right)>g:
            right *= 2
            assert right < 100, "Numerical tail bracket not obtained"
        w = brentq(lambda x:self.G(x)-g, xc, right, xtol=1e-14)
        return u, w

    def solve(self):
        if self.eta<=0:
            return {"active":False,"g":self.theta*self.Bi,"ell":self.i,
                    "curvature":0.}
        base = max(self.Bi, self.Bj)
        top = self.G(self.eta/self.mid)
        # An inactive pair is possible only when the lower label has larger mean.
        if self.Bi>self.Bj:
            xc = self.eta/self.mid
            wr = max(2., 2*xc)
            while self.G(wr)>base:
                wr *= 2
            w0 = brentq(lambda x:self.G(x)-base,xc,wr)
            critical = quad(lambda x:(self.G(x)-base)/self.q(x),
                            1e-12,w0,epsabs=1e-10,epsrel=1e-10)[0]/2
            if self.lam>=critical:
                return {"active":False,"g":self.theta*base,"ell":self.i,
                        "curvature":0.}

        def cost_error(g):
            u,w = self.band(g)
            return quad(lambda x:(self.G(x)-g)/self.q(x),u,w,
                        epsabs=1e-11,epsrel=1e-11)[0]-2*self.lam
        # Keep endpoint roots numerically separated from the maximum.
        g = brentq(cost_error, base+1e-9, top-1e-9, xtol=1e-13)
        u,w = self.band(g)
        J = quad(lambda x:1/self.q(x),u,w,epsabs=1e-11,epsrel=1e-11)[0]
        ell = quad(lambda x:self.A(x)/self.q(x),u,w,
                   epsabs=1e-11,epsrel=1e-11)[0]/J
        curv = ((ell-self.A(u))**2/(self.q(u)*self.Gx(u))
                -(self.A(w)-ell)**2/(self.q(w)*self.Gx(w)))/J/self.theta
        # Independently integrate the invariant holding density.
        def holding(x):
            return quad(lambda y:1/self.q(y),x,w,epsabs=1e-10)[0]/J
        occupation = self.P(u)+quad(
            lambda x:holding(x)*self.F*self.q(x)/x,u,w,epsabs=1e-10)[0]
        ell_occ = self.i+(self.j-self.i)*occupation
        assert abs(ell_occ-ell)<1e-9
        return {"active":True,"g":self.theta*g,"ell":ell,"curvature":curv,
                "u_dim":u,"w_dim":w,"u":u*self.theta,"w":w*self.theta,
                "g_dim":g,"J_dim":J,"ell_occupation":ell_occ}


def main():
    identities = {}
    def check(name, expr):
        value = s.simplify(expr)
        assert value==0, (name,value)
        identities[name] = str(value)
    g0,g1,b0,b1,ell,c0 = s.symbols("g0 g1 b0 b1 ell c0")
    check("affine_policy_bregman",g1-(ell*b1+c0)
          -(g1-g0-ell*(b1-b0)+(g0-(ell*b0+c0))))
    Au,Aw,Gu,Gw,qu,qw,J = s.symbols("Au Aw Gu Gw qu qw J",nonzero=True)
    up,wp = (ell-Au)/Gu,(ell-Aw)/Gw
    check("pair_curvature",((Aw-ell)*wp/qw-(Au-ell)*up/qu)/J
          -((ell-Au)**2/(qu*Gu)-(Aw-ell)**2/(qw*Gw))/J)
    # A kink with an interior supporting slope has two different coefficients.
    l,r,chosen = 1.,3.,1.5
    for delta in [-.1,.1]:
        value = max(l*delta,r*delta)-chosen*delta
        expected = ((chosen-l) if delta<0 else (r-chosen))*abs(delta)
        assert abs(value-expected)<1e-15
    # Continuous derivative but unbounded curvature need not give quadratic regret.
    ratios = [h**1.5/h**2 for h in [1e-2,1e-4,1e-6]]
    assert ratios[0]<ratios[1]<ratios[2]
    nominal = Pair(.06).solve()
    derivative_checks=[]
    for b in [.035,.05,.06,.07,.085]:
        sol=Pair(b).solve(); h=1e-5
        left,right=Pair(b-h).solve(),Pair(b+h).solve()
        num_slope=(right['g']-left['g'])/(2*h)
        num_curvature=(right['ell']-left['ell'])/(2*h)
        assert abs(num_slope-sol['ell'])<1e-6
        assert abs(num_curvature-sol['curvature'])<1e-4
        derivative_checks.append({"b":b,"analytic_slope":sol['ell'],
            "finite_difference_slope":num_slope,"analytic_curvature":sol['curvature'],
            "finite_difference_curvature":num_curvature})
    qinv=2*3.1/(2*3.1*.04-.09); z=norm.ppf(.975)
    rows=[]
    for T in [30,64,100,192,1000]:
        rad=z/math.sqrt(T*qinv)
        ends=[Pair(.06-rad).solve(),Pair(.06+rad).solve()]
        regrets=[x['g']-nominal['g']-nominal['ell']*delta
                 for x,delta in zip(ends,[-rad,rad])]
        endpoint=max(regrets)
        # Independent grid comparisons challenge endpoint maximization.
        grid=[Pair(b).solve()['g']-nominal['g']-nominal['ell']*(b-.06)
              for b in np.linspace(.06-rad,.06+rad,13)]
        assert min(grid)>-1e-10 and max(grid)<=endpoint+1e-10
        rows.append({"T_planning_only":T,"I_proxy":T*qinv,"radius":rad,
            "linear_W3":3*rad,"pair_endpoint_regret":endpoint,
            "endpoint_regrets":regrets,"endpoint_active":[x['active'] for x in ends],
            "unrestricted_frictionless_95_benchmark":z*z/(2*T),
            "unrestricted_expected_loss_asymptotic":1/(2*T)})
    result={"status":"passed","scope":"Symbolic identities and floating-point diagnostics; "
            "not interval certificates or observed-data results.",
        "versions":{"sympy":s.__version__,"scipy":scipy.__version__},
        "identities":identities,"kink_coefficients":[chosen-l,r-chosen],
        "continuous_derivative_nonquadratic_ratios":ratios,
        "parameters":{"kappa":3.1,"theta":.04,"xi2":.09,"bhat":.06,
                      "menu":[0,3],"fees":[0,0],"physical_log_cost":.001},
        "nominal_pair":nominal,"derivative_checks":derivative_checks,
        "q_inv":qinv,"z95":z,"planning_rows":rows,
        "linear_theta_z_years_eps01":9*z*z*.04/.01**2,
        "linear_theta_chebyshev_years_eps01":9*.04/(.05*.01**2),
        "linear_qinv_proxy_z_years_eps01":9*z*z/(qinv*.01**2),
        "unrestricted_quadratic_planning_years":{str(e):z*z/(2*e) for e in [.03,.01]}}
    frozen=json.loads((OUT/'input_manifest.json').read_text())['protected_files']
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
               for p,h in frozen.items())
    result['protected_inputs_unchanged']=len(frozen)
    (OUT/'derivation_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
