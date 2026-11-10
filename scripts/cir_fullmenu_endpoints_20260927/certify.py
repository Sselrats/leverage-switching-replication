"""Three exact continuum gains and one deployed-band interval loss certificate.

Acceptance uses outward interval arithmetic and strict rational signs.
SciPy supplies proposals only. c=9/5000 is proportional, so log cost is
the exact transcendental -log(4991/5000), enclosed using iv.log.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/paper_v11_research'))
from general_CIR_certificate import CIR,I,R,lo,hi,ends,iv
from explore import Full,COST

OUT=ROOT/'artifacts/cir_fullmenu_endpoints_20260927'
THETA=Q(1,25)
LOG_COST=-iv.log(R(Q(4991,5000)))
LAM=R(Q(155,2))*LOG_COST

def enc(x):return ends(x)
def hull(x,y):return I(min(lo(x),lo(y)),max(hi(x),hi(y)))
def lower_min(xs):return min(lo(x) for x in xs)
def positive(values):
    assert all(lo(x)>0 for x in values.values()),{k:enc(v) for k,v in values.items()}

class Model(CIR):
    def __init__(self,eta,eta_hi=None):
        eta_hi=Q(eta) if eta_hi is None else Q(eta_hi)
        eta=Q(eta); self.eta_exact=(eta+eta_hi)/2; self.eta_width=eta_hi-eta
        super().__init__(Q(124,45),Q(124,45),eta,eta_hi,N=140)
        self.numeric=Full(float(self.eta_exact))
        self.B2=2*self.eta-2

    def rootbox(self,p,g,side):
        return self.root(p,g,side)

    def primitive(self,p,x,g):
        return self.K(x)-g*self.E(x) if p==23 else self.L(p,x,g)

    def area(self,p,u,w,g):
        return self.primitive(p,w,g)-self.primitive(p,u,g)

    def pair(self,p):
        proposal=self.numeric.sol[p]
        center=Q(str(proposal['g_dim']))-(2*self.eta_exact-2 if p==23 else 0)
        radius=Q(1,10**11)+8*self.eta_width
        gl,gh=center-radius,center+radius
        roots=[]; areas=[]
        for g in [gl,gh]:
            u=self.rootbox(p,R(g),'a'); w=self.rootbox(p,R(g),'b')
            roots.append((u,w));areas.append(self.area(p,u,w,R(g)))
        signs={'low_gain_area_minus_cost':areas[0]-2*LAM,
               'cost_minus_high_gain_area':2*LAM-areas[1]}
        positive(signs)
        g=I(gl,gh); u=hull(roots[0][0],roots[1][0]); w=hull(roots[0][1],roots[1][1])
        J=self.E(w)-self.E(u); ap=self.A(w)-self.A(u)
        i,j={12:(1,2),13:(1,3),23:(2,3)}[p]
        ell=i+(j-i)*ap/J
        total=g+self.B2 if p==23 else g
        return {'g':g,'total':total,'u':u,'w':w,'ell':ell,'signs':signs}

    def values(self,g,k,gamma):
        a=self.rootbox(13,g,'a'); b=self.rootbox(13,g,'b')
        u=self.rootbox(12,g-k,'a'); beta=self.rootbox(12,g-k,'b')
        zeta=self.rootbox(23,k,'b')
        Va=self.area(12,a,beta,g-k)
        shoot=Va+self.area(23,a,gamma,k)
        H=self.area(13,a,b,g)
        Vmax=self.area(12,u,beta,g-k)
        H23=self.area(23,gamma,zeta,k)
        roots=dict(gamma=gamma,a=a,u=u,zeta=zeta,b=b,beta=beta)
        orders=dict(gamma_positive=gamma,a_minus_gamma=a-gamma,u_minus_a=u-a,
                    zeta_minus_u=zeta-u,b_minus_zeta=b-zeta,beta_minus_b=beta-b,
                    x23_minus_gamma=2*self.eta/5-gamma,beta_minus_x12=beta-2*self.eta/3)
        positive(orders)
        return roots,orders,{'Q':shoot,'H13':H,'Vmax':Vmax,'H23':H23}

    def low(self):
        pair=self.pair(13)
        gamma=R(Q(3,50)); k=self.G(23,gamma)
        roots,orders,vals=self.values(pair['g'],k,gamma)
        # Existing compatible-middle certificate, permitting a constant shift.
        Qv=vals['Q']; Tv=Qv-LAM+vals['H23']
        minimum=I(min(Q(0),lo(Qv-LAM)),min(Q(0),hi(Qv-LAM)))
        maximum=I(max(lo(vals['Vmax']),lo(Tv),lo(LAM)),
                  max(hi(vals['Vmax']),hi(Tv),hi(LAM)))
        margin=2*LAM-(maximum-minimum)
        positive({'switching_span_margin':margin})
        return pair['g'],{'branch':'optimal outer pair with compatible middle certificate',
            'gain_dim':enc(pair['g']),'gain_physical':enc(R(THETA)*pair['g']),
            'mean_leverage':enc(pair['ell']),'pair_root_signs':{k:enc(v) for k,v in pair['signs'].items()},
            'roots':{k:enc(v) for k,v in roots.items()},'orders':{k:enc(v) for k,v in orders.items()},
            'middle_shift_switching_margin':enc(margin)}

    def nominal(self):
        pair23=self.pair(23); k=pair23['g']; gamma=pair23['u']
        proposal=self.numeric.solve()
        center=Q(str(proposal['gain_dim'])); gl,gh=center-Q(1,10**8),center+Q(1,10**8)
        signs={}
        for key,g in [('low',gl),('high',gh)]:
            _,_,vals=self.values(R(g),k,gamma)
            signs[key]=vals['Q']-LAM
        positive({'low_shoot':signs['low'],'minus_high_shoot':-signs['high']})
        g=I(gl,gh); roots,orders,vals=self.values(g,k,gamma)
        margins={'H13_above_lambda':vals['H13']-LAM,
                 '2lambda_above_H13':2*LAM-vals['H13'],
                 '2lambda_above_Vmax':2*LAM-vals['Vmax']}
        positive(margins)
        ga,a,ze,be=(roots[x] for x in ['gamma','a','zeta','beta'])
        Ea,Ega,Eze,Ebe=(self.E(x) for x in [a,ga,ze,be])
        Aa,Aga,Aze,Abe=(self.A(x) for x in [a,ga,ze,be])
        J=Ebe-Ea; prob=(Eze-Ega)/(Ebe-Ega)
        Te=(1-prob)*Aga+prob*Abe-Aze
        Tlower=Aze-Aga; Tret=Aa-Ea-Abe+Ebe; Tup=Aze-Aa
        Tcycle=prob*J
        ell_cycle=(2*Te+3*(1-prob)*Tlower+prob*(Tret+3*Tup))/Tcycle
        # Independent differentiation of the shooting equation in eta.
        kprime=pair23['ell']-2
        ell_shoot=1+(Abe+Aga-2*Aa+kprime*(Ebe-Ega))/J
        assert max(lo(ell_cycle),lo(ell_shoot))<=min(hi(ell_cycle),hi(ell_shoot))
        positive({'minus_Q_g':J,'gamma_simple':self.f(ga)*(self.eta-R(Q(5,2))*ga),
                  'zeta_simple':-self.f(ze)*(self.eta-R(Q(5,2))*ze)})
        return g,{'branch':'23-calibrated four-threshold optimum',
            'gain_dim':enc(g),'gain_physical':enc(R(THETA)*g),
            'shooting_signs':{k:enc(v) for k,v in signs.items()},
            'pair23_signs':{k:enc(v) for k,v in pair23['signs'].items()},
            'roots':{k:enc(v) for k,v in roots.items()},'orders':{k:enc(v) for k,v in orders.items()},
            'strict_margins':{k:enc(v) for k,v in margins.items()},
            'Q_g':enc(-J),'mean_leverage_cycle':enc(ell_cycle),'mean_leverage_shooting':enc(ell_shoot),
            'local_no_kink':'Strict margins and simple roots imply an analytic full-optimal branch on an open premium neighborhood; neighborhood width not certified.'}

    def high(self):
        p12=self.pair(12); p23=self.pair(23); k=p23['g']; g=p12['g']+k
        roots,orders,vals=self.values(g,k,p23['u'])
        margins={'Q_positive':vals['Q'],'lambda_minus_Q':LAM-vals['Q'],
            'H13_minus_Q':vals['H13']-vals['Q'],'lambda_minus_H13_plus_Q':LAM-vals['H13']+vals['Q']}
        positive(margins)
        ell=p12['ell']+p23['ell']-2
        return g,{'branch':'adjacent ladder optimum, gain g12+g23-B2',
            'gain_dim':enc(g),'gain_physical':enc(R(THETA)*g),'mean_leverage':enc(ell),
            'pair12_signs':{k:enc(v) for k,v in p12['signs'].items()},
            'pair23_signs':{k:enc(v) for k,v in p23['signs'].items()},
            'roots':{k:enc(v) for k,v in roots.items()},'orders':{k:enc(v) for k,v in orders.items()},
            'strict_margins':{k:enc(v) for k,v in margins.items()},
            'active_equalities':'Vmax=H12=2lambda and H23=2lambda hold at the exact pair roots; interval endpoint comparisons are not used to assert equality.'}

    def contact_quantities(self):
        p13=self.pair(13); p23=self.pair(23)
        roots,orders,vals=self.values(p13['g'],p23['g'],p23['u'])
        a,ga,be=(roots[key] for key in ['a','gamma','beta'])
        J=self.E(be)-self.E(a)
        ell=1+(self.A(be)+self.A(ga)-2*self.A(a)
            +(p23['ell']-2)*(self.E(be)-self.E(ga)))/J
        jump=ell-p13['ell']
        return p13,roots,orders,vals,J,ell,jump

def certify_contact():
    elo,ehi=Q('1.39631804'),Q('1.39631806')
    signs={}
    for name,eta in [('low',elo),('high',ehi)]:
        m=Model(eta); vals=m.contact_quantities()[3]
        signs[name]=vals['Q']-LAM
    positive({'minus_phi_low':-signs['low'],'phi_high':signs['high']})
    m=Model(elo,ehi)
    pair,roots,orders,vals,J,ell,jump=m.contact_quantities()
    margins={'2lambda_above_Vmax':2*LAM-vals['Vmax'],'mean_leverage_jump':jump,
             'phi_eta_derivative':J*jump,'minus_Q_g':J}
    positive(margins)
    return {'status':'certified local kink at a premium deletion boundary',
        'eta_bracket':[str(elo),str(ehi)],'physical_premium_bracket':[str(elo*THETA),str(ehi*THETA)],
        'phi_definition':'Q(g13(eta), k23(eta); eta) - lambda, k23=g23-B2',
        'phi_endpoint_signs':{k:enc(v) for k,v in signs.items()},
        'gain_physical':enc(R(THETA)*pair['g']),
        'mean_leverage_pair13':enc(pair['ell']),'mean_leverage_four':enc(ell),
        'leverage_jump':enc(jump),'strict_margins':{k:enc(v) for k,v in margins.items()},
        'roots':{k:enc(v) for k,v in roots.items()},'orders':{k:enc(v) for k,v in orders.items()},
        'equalities':'At the unique Phi root in this bracket, Q=lambda; exact pair equations give H13=H23=2lambda. Equality is analytic, not inferred from overlapping intervals.',
        'local_branch_identification':'Full optimum is pair13 on a left neighborhood and the 23-calibrated four-threshold branch on a right neighborhood; only local neighborhoods, not a global premium phase diagram.',
        'loss_expansion':'For the exact pair13 policy calibrated at b_star: regret(b_star+h) = jump*h + O(h^2), regret(b_star-h)=O(h^2), as h down to zero.'}

def main():
    cases={}; gains={}; models={}
    for name,eta,method in [('lower',Q(19,20),'low'),('nominal',Q(3,2),'nominal'),('upper',Q(41,20),'high')]:
        model=Model(eta); models[name]=model
        gains[name],rec=getattr(model,method)(); cases[name]=rec
        print(json.dumps({'case':name,'gain':rec['gain_physical']['dec']}),flush=True)
    nominal=models['nominal']; pair=nominal.pair(13)
    proposal=nominal.numeric.sol[13]
    a=Q(round(proposal['u_dim']*10**13),10**13)
    b=Q(round(proposal['w_dim']*10**13),10**13)
    J=nominal.E(R(b))-nominal.E(R(a))
    gain0=R(THETA)*(nominal.area(13,R(a),R(b),R(Q(0)))-2*LAM)/J
    ell=1+2*(nominal.A(R(b))-nominal.A(R(a)))/J
    radius=R(Q(11,500))
    regret_low=R(THETA)*gains['lower']-gain0+ell*radius
    regret_high=R(THETA)*gains['upper']-gain0-ell*radius
    worst=I(max(lo(regret_low),lo(regret_high)),max(hi(regret_low),hi(regret_high)))
    nominal_gap=R(THETA)*gains['nominal']-gain0
    exact_deletion=R(THETA)*(gains['nominal']-pair['g'])
    positive({'nominal_deletion':exact_deletion,'worst_positive':worst})
    assert hi(worst)<Q(7,1000)
    out={'status':'certified','scope':'Known variance law; prescribed premium interval; three continuum optimal gains via global verification and attained policies; no observed data or variance-error guarantee.',
        'parameters':{'kappa':'31/10','theta':'1/25','xi2':'9/100','F':'124/45','fees':[0,0,0],
            'full_menu':[1,2,3],'candidate_menu':[1,3],'proportional_switching_cost':'9/5000',
            'log_cost_definition':'-log(4991/5000)','log_cost':enc(LOG_COST),'lambda':enc(LAM),
            'nominal_premium':'3/50','premium_interval':['19/500','41/500']},
        'cases':cases,'candidate':{'a_dim':str(a),'b_dim':str(b),'a_physical':str(a*THETA),
            'b_physical':str(b*THETA),'nominal_gain':enc(gain0),'mean_leverage':enc(ell),
            'nominal_optimal_pair_gain':enc(R(THETA)*pair['g'])},
        'nominal_full_menu_gap':enc(nominal_gap),'nominal_exact_deletion_value':enc(exact_deletion),
        'endpoint_regrets':{'lower':enc(regret_low),'upper':enc(regret_high)},
        'worst_interval_regret':enc(worst),'W_radius':'11/250',
        'C1_comparison_note':'W*r=0.044 alone is not a valid C1 certificate here: nominal deletion is positive. It is the premium term only.',
        'local_loss_expansion_note':'Positive nominal gap plus a signed linear term from the leverage difference; this is not kink regret at a nominal full optimizer.',
        'trusted_arithmetic':'mpmath 1.3.0 iv, 50 decimal digits, positive series through index 140 with bounded tails; strict rational sign acceptance; floating proposals only'}
    out['premium_deletion_boundary']=certify_contact()
    paths=[Path(__file__),ROOT/'scripts/cir_fullmenu_endpoints_20260927/explore.py',
        ROOT/'scripts/cir_premium_bregman_20260927/check_derivations.py',
        ROOT/'scripts/paper_v11_research/general_CIR_certificate.py',
        ROOT/'scripts/paper_v11_research/validated_pairs_F2.py',ROOT/'scripts/certificate_zero_interval_F2.py']
    out['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (OUT/'full_menu_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'certified','regrets':{k:v['dec'] for k,v in out['endpoint_regrets'].items()},
        'nominal_gap':out['nominal_full_menu_gap']['dec'],'worst':out['worst_interval_regret']['dec']},indent=2))

if __name__=='__main__':main()
