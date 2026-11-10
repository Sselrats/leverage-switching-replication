#!/usr/bin/env python3
"""Shooting functions for the four-threshold policy in Proposition 4.3."""
from functools import lru_cache

import f2_kernel as k
from f2_kernel import R, lo, hi, mp

@lru_cache(None)
def propose(l):
 gamma,zeta,g23=k.proposal_band(23,l);kk=g23-1;la=k.M(l)
 def eq(a,be,g):
  va=k.Ln(12,be,g-kk)-k.Ln(12,a,g-kk)
  qq=va+k.Ln(23,gamma,1)-k.Ln(23,a,1)-kk*(k.En(gamma)-k.En(a))
  return k.Gn(13,a)-g,k.Gn(12,be)-(g-kk),qq-la
 return mp.findroot(eq,(mp.mpf('.48'),mp.mpf('1.35'),mp.mpf('1.46')),tol=mp.mpf('1e-75'),maxsteps=100)[2]
def rootbox(pair,g,side):return k.hull(k.root(pair,lo(g),side),k.root(pair,hi(g),side))
def atg(g,rr):
 ga,ze,kk,lam=rr['a'],rr['b'],rr['g']-1,rr['lambda']
 a=rootbox(13,g,'a');be=rootbox(12,g-kk,'b')
 va=k.L(12,be,g-kk)-k.L(12,a,g-kk)
 qq=va+k.L(23,ga,R(1))-k.L(23,a,R(1))-kk*(k.E(ga)-k.E(a))
 return qq-lam,a,be,va

