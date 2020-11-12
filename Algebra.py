#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Filename:   Algebra.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2018
Disclaimer: This code is presented "as is" and it has been written to
            generate random models of exams for the subject of Linear
            at ESEIAAT
License:    This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.

 	        See https://www.gnu.org/licenses/
"""

import random
import math
import copy
import collections
import itertools
from fractions import gcd
from functools import reduce
from sympy import *
from sympy.solvers.solveset import linsolve
from sympy.vector import Vector
from fractions import gcd
from functools import reduce
from sympy import Basic, Function, Symbol
from sympy.printing.printer import Printer
from sympy.printing.latex import print_latex
from sympy.core.basic import Basic
from itertools import permutations

var('x y z t')
ddict = collections.OrderedDict([(x**2,1),(y**2,2),(z**2,3),(t**2,4),
                                ((x,y),5),((x,z),6),((x,t),7),
                                ((y,z),8),((y,t),9),((z,t),10),
                                (x,11),(y,12), (z,13), (t,14)])

class Impresora(Printer):
    """
       La funció latex() del sympy té la mania d'escriure les variables x, y, z i t
       en l'ordre t, x, y i z. L'única manera que, de moment, he trobat per resoldre
       aquest inconvenient és definir la classe Impresorai la funció mylatex().
       Ho he trobat a StackOverflow.
    """
    printmethod = 'impresora'
    #
    #
    #
    def _print_Add(self,expr):
        expr_args=expr.args

        def new_place(el):
            if el in ddict:
                return ddict[el]
            else:
                return len(ddict)+1

        def get_place(el):
            if el.is_integer:
                return new_place(el)
            elif el.is_symbol:
                return new_place(el)
            if len(el.args) == 2:
                if el.args[0].is_symbol and el.args[1].is_symbol:
                    return new_place(el.args)
                q = el.args[len(el.args)-1]
                if q.is_symbol:
                    return new_place(q)
                elif q.args[0].is_symbol:
                    return new_place(q)
                else:
                    return 0
            elif len(el.args) == 3:
                return new_place(el.args[1:])
            else:
                return 0

        def write_coeff(el):
            if el.is_integer:
                if el > 0:
                    return "+%s" % el
                else:
                    return "%s" % el
            elif el.is_symbol:
                return "+%s" % el
            elif len(el.args) == 2 and el.args[0].is_symbol and el.args[1].is_symbol:
                return "+%s" % latex(el)
            elif len(el.args) > 0:
                if el.args[len(el.args)-1].is_symbol:
                    if el.args[0].is_rational:
                        if el.args[0] > 0:
                            return "+%s" % latex(el)
                        else:
                            return "%s" % latex(el)
                    else:
                        return "%s" % latex(el)
                else:
                    return "%s" % latex(el)
            else:
                return "%s" % el
        list_place = [get_place(a) for a in expr.args]
        expr_args = list(zip(*sorted(zip(list_place,expr_args))))[1]
        to_print = [write_coeff(a) for a in expr_args]
        to_print[0] = str(latex(expr_args[0]))
        return "".join(a for a in to_print)

def mylatex(e):
    return (Impresora().doprint(e))

def mcd_llista(list):
    x = reduce(gcd, list)
    return x

def mcm_llista(list):
    x = reduce(gcd, list)
    p = 1
    for k in list:
        p *= k
    return p // x

def norma(m):
    f, c = m.shape
    n = 0
    for i in range(f):
        for j in range(c):
            if abs(m[i,j]) > n:
                n = abs(m[i,j])
    return n

def nzeros(m):
    f, c = m.shape
    z = 0
    for i in range(f):
        for j in range(c):
            if m[i,j] == 0:
                z += 1
    return z

def simplifica(v):
    m = 1
    for i in range(len(v)):
        if isinstance(v[i],Rational):
            m *= v [i].q
    v = [m * x for x in v]
    mcd = mcd_llista(v)
    v = [x // mcd for x in v]
    if v[0] < 0:
        v = [-x for x in v]
    return v

def matriu_aleatoria(f=3,c=3,maxim=5,nuls=True):
    m = Matrix(f,c,lambda i, j : random.randint(-maxim,maxim))
    if not nuls:
        values = [i for i in range(1,maxim + 1)] + [-i for i in range(1,maxim + 1)]
        for i in range(f):
            for j in range(c):
                if m[i,j] == 0:
                    m[i,j] = values[random.randint(0,2 * maxim - 1)]
    return m

def llista_aleatoria(l=3,maxim=5,nuls=True,positius=False):
    if positius:
        c = [random.randint(1,maxim) for i in range(l)]
        values = [i for i in range(1,maxim + 1)]
        m = maxim
    else:
        c = [random.randint(-maxim,maxim) for i in range(l)]
        values = [i for i in range(1,maxim + 1)] + [-i for i in range(1,maxim + 1)]
        m = 2 * maxim - 1
    if not nuls:
        for i in range(l):
            if c[i] == 0:
                c[i] = values[random.randint(0,m)]
    return c

def matriu_amb_rang(f=3,c=3,r=3,maxim=5,nuls=True):
    trobat = False
    while not trobat:
        m = matriu_aleatoria(f,c,maxim,nuls)
        if m.rank() != r:
            continue
        trobat = True
    return m

def mti(i,j):
    values = (-1,1)
    if i > j:
        return 0
    elif i == j:
        return values[random.randint(0,1)]
    else:
        return random.randint(-1,2)

def mts(i,j,values):
    if i < j:
        return 0
    elif i == j:
        return  values[i]
    else:
        return random.randint(-2,1)

def matriu_invertible(ordre=3,maxim=5,mzeros=-1,unitaria=False):
    opcions = []
    for i in range(1,13):
        opcions += [i for j in range(0,2**(13-i))] + [-i for j in range(0,2**(13-i))]
    unitats = (-1,1)
    els = len(opcions)
    random.shuffle(opcions)
    trobat = False
    while not trobat:
        if unitaria:
            values = [unitats[random.randint(0,1)] for i in range(ordre)]
        else:
            values = [opcions[random.randint(0,els-1)] for i in range(ordre)]
        ti = Matrix(ordre,ordre,lambda i, j : mti(i,j))
        ts = Matrix(ordre,ordre,lambda i, j : mts(i,j,values))
        m = ti * ts
        if norma(m) > maxim:
            continue
        if mzeros >= 0 and nzeros(m) > mzeros:
            continue
        trobat = True
    return m

def matriu_diagonalitzable(ordre=3,maxim=5,mzeros=-1,mvaps=3,vapsnuls=False):
    trobat = False
    while not trobat:
        random.seed()
        c = matriu_invertible(ordre,maxim=3,mzeros=0,unitaria=True)
        v = [-i for i in range(1,mvaps+1)] + [i for i in range(1,mvaps+1)]
        if vapsnuls:
            vaps = [random.randint(-mvaps,mvaps) for i in range(ordre)]
        else:
            vaps = [v[random.randint(0,2*mvaps-1)] for i in range(ordre)]
        if not vapsnuls and 0 in vaps:
            continue
        if len(set(vaps)) == 1:
            continue
        d = diag(*vaps)
        a = c * d * c**(-1)
        if norma(a) > maxim:
            continue
        if mzeros >= 0 and nzeros(a) > mzeros:
            continue
        trobat = True
    return a, vaps

def matriu_gram(ordre=3,maxim=5,mzeros=-1):
    trobat = False
    while not trobat:
        c = matriu_invertible(ordre,maxim=5)
        g = c.T * c
        if mzeros >= 0 and nzeros(g) > mzeros:
            continue
        if norma(g) > maxim:
            continue
        trobat = True
    return g

def matriu_latex(m,format=None,det=False):
    f, c = m.shape
    tipus = "pmatrix"
    if det:
        tipus = "vmatrix"
    if format is None:
        text = "\\begin{TIPUS}{*{%d}r} LINES\end{TIPUS}" % c
    else:
        text = "\\begin{TIPUS}{%s} LINES\end{TIPUS}" % format
    text = text.replace('TIPUS',tipus)
    lines = []
    for i in range(f):
        line = []
        for j in range(c):
            line.append(latex(m[i,j]))
        lines.append(" & ".join(map(str,line)))
    return (text.replace('LINES',"\\\\ ".join(lines)))

def vector_latex(v,unitari=False):
    s = list(map(latex,v))
    s = ",".join(s)
    r = ""
    if unitari:
        r = f"\\deufrac{{1}}{{{longitud_latex(v)}}}"
    return r + f"({s})"

def matriu_a_llista(m):
    l = []
    f, c = m.shape
    for i in range(f):
        for j in range(c):
            l.append(m[i,j])
    return l

def matriu_inversa_latex(m,format=None):
    f,c = m.shape
    d = m.det()
    n = d * m**(-1)
    l = matriu_a_llista(n)
    l.append(d)
    mcd = abs(mcd_llista(l))
    s = 1
    if d < 0:
        d *= -1
        s = -1
    d /= mcd
    for i in range(f):
        for j in range(c):
            n[i,j] = s * n[i,j] // mcd
    text = matriu_latex(n,format)
    if d > 1:
        text = ("\\deufrac{1}{%d}" % d) + text
    return text

def combinacio_lineal(v,l):
    s = set(list(map(len,v)))
    if len(s) != 1:
        return None
    s = s.pop()
    r =  []
    for i in range(s):
        t = sum([l[k]*v[k][i] for k in range(len(l))])
        r.append(t)
    return r

def text_vectors_matriu(m,simplificar=True):
    f, c = m.shape
    vecs = []
    for i in range(c):
        v = matriu_a_llista(m.col(i))
        if simplificar:
            v = simplifica(v)
        v = ",".join(map(str,v))
        vecs.append(f"({v})")
    return ",".join(vecs)

def vectors_matriu(m,simplificar=True):
    f, c = m.shape
    vecs = []
    for i in range(c):
        v = []
        for j in range(f):
            v.append(m[j,i])
        if simplificar:
            v = simplifica(v)
        vecs.append(v)
    return vecs

def vectors_fila_matriu(m,simplificar=True):
    f, c = m.shape
    vecs = []
    for i in range(f):
        v = []
        for j in range(c):
            v.append(m[i,j])
        if simplificar:
            v = simplifica(v)
        vecs.append(v)
    return vecs

def producte_vectorial(v1,v2,simplificar=True):
    v = [v1[1] * v2[2] - v1[2] * v2[1],- v1[0] * v2[2] + v1[2] * v2[0],v1[0] * v2[1] - v1[1] * v2[0]]
    if simplificar:
        v = simplifica(v)
    return v

def producte_escalar(v1,v2):
    p = 0
    for i in range(len(v1)):
        p += v1[i] * v2[i]
    return p

def sistema_equacions(m,b=None,simplificar=True):
    f, c = m.shape
    if b is None:
        b = zeros(f,1)
    if c <= 4:
        x, y, z, t = symbols('x y z t')
        incg = [x,y,z,t]
    else:
        x1, x2, x3, x4, x5, x6, x7 = symbols('x1 x2 x3 x4 x5 x6 x7')
        incg = [x1,x2,x3,x4,x5,x6,x7]
    incg = Matrix(c,1,incg[0:c])
    eqs = []
    if simplificar:
        n = m * incg - b
        for i in range(f):
            eqs.append(simplificar_equacio(n[i]))
    else:
        a = m * incg
        for i in range(f):
            eqs.append(f"{mylatex(a[i])} &= {latex(b[i])}")
    eqs = " \\\\ ".join(eqs)
    return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"

def solucio_general_sistema(u,s):
    f = len(u)
    eqs = " \\\\ ".join([f"{latex(u[i])} &= {mylatex(s[i,0])}" for i in range(f)])
    return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"

def equacio_del_pla_vectorial(n,m=None,r=None):
    x, y, z = symbols('x y z')
    if m is None:
        v1, v2 = vectors_matriu(n)
    else:
        v1, v2 = vectors_matriu(m**(-1) * n)
    v = producte_vectorial(v1,v2)
    eq = mylatex (v[0] * x + v[1] * y + v[2] * z)
    if r is not None:
        for k, v in r.items():
            eq = eq.replace(k,v)
    eq += ' = 0'
    return eq

def equacio_del_pla_afi(w,p):
    t = producte_escalar(w,p)
    x, y, z = symbols('x y z')
    return f"{mylatex(w[0] * x + w[1] * y + w[2] * z)} = {t}"

def equacio_continua(p,v):
    v = simplifica(v)
    x, y, z = symbols('x y z')
    incg = [x,y,z]
    eq = []
    for i in range(3):
        if abs(v[i]) == 1:
            eq.append(latex(incg[i]-p[i]))
        else:
            eq.append(f"\\frac{{{latex(incg[i]-p[i])}}}{{{v[i]}}}")
    return " = ".join(eq)

def simplificar_equacio(eq,r=None,amp=True,positiu=None):
    d = eq.as_coefficients_dict()
    l = d.values()
    mcd = mcd_llista(l)
    eq = 0
    for k in d.keys():
        d[k] = d[k] // mcd
        eq += d[k] * k
    if positiu is not None and d[positiu] < 0:
        eq *= -1
        d[1] *= -1
    eq -= d[1]
    if amp:
        str = f"{mylatex(eq)} &= {-d[1]}"
    else:
        str = f"{mylatex(eq)} = {-d[1]}"
    if r is None:
        return str
    for k,v in r.items():
        str = str.replace(k,v)
    return str

def eliminacio_dos_parametres(n,m):
    f, c = n.shape
    if f <= 4:
        x, y, z, t = symbols('x y z t')
        incg = [x,y,z,t]
    else:
        x1, x2, x3, x4, x5, x6, x7 = symbols('x1 x2 x3 x4 x5 x6 x7')
        incg = [x1,x2,x3,x4,x5,x6,x7]
    incg = incg[0:f]
    a = n.col_insert(2,Matrix(f,1,incg) - m)
    all = list(range(f))
    options = permutations(all,2)
    for o in options:
        if n[o,:].det() != 0:
            break
    for e in o:
        all.remove(e)
    eqs = []
    for e in all:
        f = [el for el in o]
        f.append(e)
        f.sort()
        d = a[f,:].det()
        eq = simplificar_equacio(d)
        eqs.append(eq)
    if len(eqs) == 1:
        return eqs[0].replace('&','')
    eqs = " \\\\ ".join(eqs)
    return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"

def equacions_parametriques(a,b,r=None):
    u,v =  symbols('u v')
    f, c = a.shape
    if f <= 4:
        x, y, z, t = symbols('x y z t')
        incg = [x,y,z,t]
    else:
        x1, x2, x3, x4, x5, x6, x7 = symbols('x1 x2 x3 x4 x5 x6 x7')
        incg = [x1,x2,x3,x4,x5,x6,x7]
    m = a * Matrix(2,1,[u,v]) + b
    eqs = []
    for i in range(f):
        str = f"{latex(incg[i])} &= {latex(m[i,0])}"
        for k,v in r.items():
            str = str.replace(k,v)
        eqs.append(str)
    eqs = " \\\\ ".join(eqs)
    return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"

def multiplicar_matriu(a):
    l = matriu_a_llista(a)
    d = []
    for k in l:
        if not k.is_integer:
            d.append(k.q)
    if len(d) == 0:
        return a
    if len(d) == 1:
        return d[0] * a
    return mcm_llista(d) * a

def solucio_equacio_matricial(a,x):
    f,c = a.shape
    z1, z2, z3, z4, z5, z6, z7 = symbols('z1 z2 z3 z4 z5 z6 z7')
    incg = [z1,z2,z3,z4,z5,z6,z7]
    n = a.nullspace(simplify=True)
    if len(n) == 0:
        return matriu_latex(x)
    elif len(n) == 1:
        n = multiplicar_matriu(n[0])
        s = x.col(0) + z1 * n
        for i in range(1,c):
            s = s.col_insert(i,x.col(i) + incg[i] * n)
        return matriu_latex(s,format='*{%d}c' % c)
    else:
        print('No implementat')
        sys.exit(0)

def nucli(m):
    n = m.nullspace()
    v = []
    for i in range(len(n)):
        v.append(vectors_matriu(n[i])[0])
    return v

def base_ortogonal_nucli(m):
    v1, v2 = nucli(m)
    t1 = producte_escalar(v1,v1)
    t2 = producte_escalar(v1,v2)
    u2 = simplifica([t1*v2[i]-t2*v1[i] for i in range(len(v1))])
    return v1,u2

def longitud_latex(v):
    s = 0
    for i in range(len(v)):
        s += v[i]**2
    return latex(sqrt(s))
