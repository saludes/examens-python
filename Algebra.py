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

var('x y z t')
ddict = collections.OrderedDict([(x,1),(y,2), (z,3), (t,4)])

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
            elif len(el.args)>0:
                if el.args[len(el.args)-1].is_symbol:
                    return new_place(el.args[len(el.args)-1])
                else:
                    return 0
            else:
                return 0

        def write_coeff(el):
            if el.is_integer:
                if el > 0:
                    return "+%s" % el
                else:
                    return "%s" % el
            elif el.is_symbol:
                return "+%s" %el
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
    s = 1
    mcd = mcd_llista(v)
    v = [x / mcd for x in v]
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

def matriu_latex(m,format=None):
    f, c = m.shape
    if format is None:
        text = "\\begin{pmatrix}{*{%d}r} LINES\end{pmatrix}" % c
    else:
        text = "\\begin{pmatrix}{%s} LINES\end{pmatrix}" % format
    lines = []
    for i in range(f):
        line = []
        for j in range(c):
            line.append(latex(m[i,j]))
        lines.append(" & ".join(map(str,line)))
    return (text.replace('LINES',"\\\\ ".join(lines)))

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

def mylatex(e):
    return (Impresora().doprint(e))

def sistema_equacions(m,b):
    f, c = m.shape
    if c <= 4:
        x, y, z, t = symbols('x y z t')
        incg = [x,y,z,t]
    else:
        x1, x2, x3, x4, x5, x6, x7 = symbols('x1 x2 x3 x4 x5 x6 x7')
        incg = [x1,x2,x3,x4,x5,x6,x7]
    incg = Matrix(c,1,incg[0:c])
    n = m * incg - b
    eqs = []
    for i in range(f):
        eqs.append(simplificar_equacio(n[i]))
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

def simplificar_equacio(eq,r=None):
    d = eq.as_coefficients_dict()
    l = d.values()
    mcd = mcd_llista(l)
    eq = 0
    for k in d.keys():
        d[k] = d[k] // mcd
        eq += d[k] * k
    eq -= d[1]
    str = f"{mylatex(eq)} &= {-d[1]}"
    if r is None:
        return str
    for k,v in r.items():
        str = str.replace(k,v)
    return str

def eliminacio_dos_parametres(n,m):
    a = n.col_insert(2,Matrix(4,1,[x,y,z,t]) - m)
    options = ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
    for o in options:
        if n[o,:].det() != 0:
            break
    all = [0,1,2,3]
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
    n = multiplicar_matriu(n[0])
    s = x.col(0) + z1 * n
    for i in range(1,c):
        s = s.col_insert(i,x.col(i) + incg[i] * n)
    return matriu_latex(s,format='*{%d}c' % c)
