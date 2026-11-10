#!/usr/bin/env python3
"""Rational enclosures for pair costs and finite-menu path thresholds.

The graph weights in Theorem 3.1 are enclosed using fractions.Fraction,
positive series with explicit remainder bounds, and rational bisection.
The program checks each increasing path and rounds displayed decimals outward.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import platform
from fractions import Fraction as Q
from pathlib import Path
from typing import Dict, List, Tuple

Interval = Tuple[Q, Q]
TAIL_TOL = Q(1, 10**44)
ROOT_STEPS = 64

def need(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def series_bounds(F: Q, x: Q, integral: bool = False) -> Interval:
    """Enclose D(x)=Theta(x)-1 or A(x)=int_0^x D(t)dt.

    t_n=(Fx)^n/(F+1)_n, n>=1. For all n>=N the term ratios
    are <=r=Fx/(F+N+1). The same r bounds integral-term ratios.
    The first omitted term divided by 1-r encloses the whole tail.
    """
    need(F > 1 and x >= 0, 'Invalid F or x')
    if x == 0:
        return Q(0), Q(0)
    term = Q(1)
    total = Q(0)
    for n in range(1, 10001):
        term *= F*x/(F+n)
        total += term*x/(n+1) if integral else term
        r = F*x/(F+n+1)
        if r < 1:
            first_omitted = term*r
            if integral:
                first_omitted *= x/(n+2)
            tail = first_omitted/(1-r)
            if tail <= TAIL_TOL:
                return total, total+tail
    raise ArithmeticError('Series failed to meet the tail tolerance')


def root_bracket(F: Q, eta: Q, m: Q) -> Interval:
    """Enclose the unique root D(z)=eta/(m-eta)."""
    need(m > eta > 0, 'Root requires m > eta > 0')
    target = eta/(m-eta)
    lo, hi = Q(0), Q(1)
    while series_bounds(F, hi)[0] <= target:
        hi *= 2
        need(hi < 1024, 'Root outside supported computational range')
    for _ in range(ROOT_STEPS):
        mid = (lo+hi)/2
        dl, du = series_bounds(F, mid)
        if du < target:
            lo = mid
        elif dl > target:
            hi = mid
        else:
            raise ArithmeticError('Root sign undecided: tighten TAIL_TOL')
    need(series_bounds(F, lo)[1] < target, 'Lower root sign is not strict')
    need(series_bounds(F, hi)[0] > target, 'Upper root sign is not strict')
    return lo, hi


def edge_bounds(F: Q, eta: Q, a: Q, b: Q) -> Tuple[Interval, Interval]:
    """Enclose w_ab = (b-a)*sigma((a+b)/2).

    For psi(z)=eta*z-(m-eta)A(z), concavity implies
    psi(lo)<=sigma<=psi(lo)+psi'(lo)*(hi-lo).
    All errors in A, D, and the root are enclosed rationally.
    """
    m, d = (a+b)/2, (a+b)/2-eta
    need(0 <= a < b and d > 0, 'Invalid pair')
    lo, hi = root_bracket(F, eta, m)
    al, au = series_bounds(F, lo, integral=True)
    dl, _ = series_bounds(F, lo)
    lower = eta*lo-d*au
    upper = eta*lo-d*al+(eta-d*dl)*(hi-lo)
    need(0 < lower <= upper, 'Invalid sigma interval')
    return ((b-a)*lower, (b-a)*upper), (lo, hi)


def paths(indices: Tuple[int, ...]):
    """Increasing paths beginning at the minimum label of a submenu."""
    for k in range(1, len(indices)):
        for chosen in itertools.combinations(indices[1:], k):
            yield (indices[0],) + chosen


def path_bound(path: Tuple[int, ...], edges: Dict[Tuple[int, int], Interval]) -> Interval:
    k = len(path)-1
    intervals = [edges[a,b] for a,b in zip(path, path[1:])]
    return (sum((p[0] for p in intervals), Q(0))/(k+1),
            sum((p[1] for p in intervals), Q(0))/(k+1))


def max_bounds(intervals: List[Interval]) -> Interval:
    need(bool(intervals), 'Cannot maximize an empty list')
    return max(v[0] for v in intervals), max(v[1] for v in intervals)


def out_decimal(x: Q, places: int, upper: bool) -> str:
    scale = 10**places
    n = -((-x.numerator*scale)//x.denominator) if upper else (x.numerator*scale)//x.denominator
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n//scale}.{n%scale:0{places}d}'


def emit(v: Interval, places: int = 14) -> dict:
    return {'lower_fraction':str(v[0]), 'upper_fraction':str(v[1]),
            'outward_decimal':[out_decimal(v[0],places,False),out_decimal(v[1],places,True)]}


def menu_record(name: str, F: Q, eta: Q, levels: List[Q]) -> dict:
    need(len(levels)>=2 and all(levels[i]<levels[i+1] for i in range(len(levels)-1)), 'Unordered menu')
    need(eta < (levels[0]+levels[1])/2, 'Theorem 3.1 assumptions fail')
    edges = {}
    pairs = []
    for i,j in itertools.combinations(range(len(levels)),2):
        value, root = edge_bounds(F,eta,levels[i],levels[j])
        edges[i,j] = value
        pairs.append({'indices':[i,j], 'levels':[str(levels[i]),str(levels[j])],
                      'weight_2c':emit(value), 'positive_root':emit(root,22)})
    all_paths = list(paths(tuple(range(len(levels)))))
    path_values = [path_bound(p,edges) for p in all_paths]
    full = tuple(range(len(levels)))
    full_value = path_bound(full,edges)
    shorter = max_bounds([v for p,v in zip(all_paths,path_values) if p != full])
    gap = full_value[0]-shorter[1], full_value[1]-shorter[0]
    threshold = max_bounds(path_values)
    cbar = max_bounds([(edges[0,j][0]/2,edges[0,j][1]/2) for j in range(1,len(levels))])
    proper = []
    for r in range(2,len(levels)):
        for idx in itertools.combinations(range(len(levels)),r):
            bounds=max_bounds([path_bound(p,edges) for p in paths(idx)])
            proper.append({'indices':list(idx), 'threshold':emit(bounds)})
    proper_max=max_bounds([(Q(v['threshold']['lower_fraction']),Q(v['threshold']['upper_fraction'])) for v in proper])
    data={'name':name, 'F':str(F), 'eta':str(eta), 'levels':[str(x) for x in levels],
          'W':str((1+1/F)*eta**2), 'pairs':pairs, 'cbar':emit(cbar),
          'anchored_paths':[{'indices':list(p),'value':emit(v)} for p,v in zip(all_paths,path_values)],
          'full_path':emit(full_value), 'threshold':emit(threshold),
          'max_shorter_anchored_path':emit(shorter), 'full_minus_shorter':emit(gap),
          'proper_submenus_with_at_least_two_levels':proper, 'max_proper_threshold':emit(proper_max)}
    if name == 'four_level_counterexample':
        need(gap[0] > Q('0.0631052289'), 'Counterexample gap is insufficient')
        data['checks']={'full_minus_all_at_most_two_edge_anchored_paths_gt_0.0631052289':True}
    if name == 'geometric_R100':
        cost = Q(21,20)
        need(proper_max[1] < cost < full_value[0], 'Cost does not separate full and proper menus')
        data['lambda']=str(cost)
        data['checks']={'all_ten_nontrivial_proper_submenus_inactive':True,
                        'four_singletons_inactive_by_definition':True,
                        'full_menu_active_by_theorem_3_1':True,
                        'full_path_unique_threshold_maximizer':bool(gap[0]>0)}
    return data


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('exact_records.json'))
    args=parser.parse_args()
    cases=[('three_level_example',Q(2),Q(61,50),[Q(1),Q(3,2),Q(3)]),
           ('four_level_counterexample',Q(2),Q(44973,40000),[Q(1),Q(5,4),Q(11),Q(111)]),
           ('geometric_R100',Q(2),Q(1),[Q(1),Q(101),Q(10001),Q(1000001)])]
    records=[]
    for case in cases:
        record=menu_record(*case)
        records.append(record)
        print(record['name'])
        for key in ['cbar','threshold','max_shorter_anchored_path','full_minus_shorter','max_proper_threshold']:
            print(' ',key,record[key]['outward_decimal'])
    result={'description':__doc__, 'python_version':platform.python_version(),
            'arithmetic':'fractions.Fraction (all acceptance tests exact)',
            'root_bisections':ROOT_STEPS,'series_tail_tolerance':str(TAIL_TOL),
            'program_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'records':records}
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Saved',args.output)

if __name__=='__main__':
    main()
