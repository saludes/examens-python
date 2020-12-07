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
    return reduce(gcd, list)

def mcm_llista(list):
    mcm = list[0]
    for i in list[1:]:
        mcm = mcm * i // gcd(mcm, i)
    return mcm

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

def norma_maxim(m):
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

def primer_no_nul(list):
    if list is None or len(list) == 0:
        return None
    for i in range(len(list)):
        if list[i] != 0:
            return i
    return None

class Vector:
    #
    #
    #
    def __init__(self,c):
        self.dimensio = len(c)
        self.components = list(c)
    #
    #
    #
    @classmethod
    def nul(cls,dim):
        l = [0 for i in range(dim)]
        return cls(l)
    #
    #
    #
    @classmethod
    def aleatori(cls,l=3,maxim=5,nuls=True,positius=False):
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
        return cls(c)
    #
    #
    #
    def __repr__(self):
        s = []
        r = False
        for x in self.components:
            if not (isinstance(x,int) or isinstance(x,Integer)):
                r = True
                s.append(latex(x))
            else:
                s.append(f"{x}")
        s = ",".join(s)
        if r:
            return f"\\left({s}\\right)"
        return f"({s})"
    #
    #
    #
    def __add__(self,other):
        if not isinstance(other,Vector):
            return None
        if self.dimensio != other.dimensio:
            return None
        c1 = self.components
        c2 = other.components
        c = []
        for i in range(self.dimensio):
            c.append(c1[i] + c2[i])
        return Vector(c)
    #
    #
    #
    def __sub__(self,other):
        if not isinstance(other,Vector):
            return None
        if self.dimensio != other.dimensio:
            return None
        c1 = self.components
        c2 = other.components
        c = []
        for i in range(self.dimensio):
            c.append(c1[i] - c2[i])
        return Vector(c)
    #
    #
    #
    def __mul__(self,other):
        types = [Rational,float,int,Float,Pow,Add,Mul]
        for t in types:
            if isinstance(other,t):
                return Vector([other * x for x in self.components])
        if isinstance(other,Matriu):
            if self.dimensio != other.files:
                return None
            u = Matrix(1,self.dimensio,self.components)
            v = u * other.matriu
            f, c = v.shape
            r = []
            for i in range(c):
                r.append(v[0,i])
            return Vector(r)
        return None
    #
    #
    #
    __rmul__ = __mul__
    #
    #
    #
    def __getitem__(self,i):
        try:
            return self.components[i]
        except:
            return None
    #
    #
    #
    def dot(self,other):
        if not isinstance(other,Vector):
            return None
        if self.dimensio != other.dimensio:
            return None
        s = 0
        for i in range(self.dimensio):
            s += self.components[i] * other.components[i]
        return s
    #
    #
    #
    def length(self):
        s = 0
        for i in range(self.dimensio):
            s += self.components[i]**2
        return sqrt(s)
    #
    #
    #
    def simplificar(self):
        m = 1
        for i in range(self.dimensio):
            if isinstance(self.components[i],Rational):
                m *= self.components[i].q
            elif isinstance(self.components[i],int):
                pass
            else:
                return
        v = [m * x for x in self.components]
        mcd = mcd_llista(v)
        v = [x // mcd for x in v]
        if v[0] < 0:
            v = [-x for x in v]
        self.components = v
    #
    #
    #
    def cross(self,other,simplificar=False):
        if not isinstance(other,Vector):
            return None
        if self.dimensio != other.dimensio or self.dimensio != 3:
            return None
        v1 = self.components
        v2 = other.components
        v = [v1[1] * v2[2] - v1[2] * v2[1],- v1[0] * v2[2] + v1[2] * v2[0],v1[0] * v2[1] - v1[1] * v2[0]]
        vec = Vector(v)
        if simplificar:
            vec.simplificar()
        return vec
    #
    #
    #
    def latex(self,unitari=False):
        s = str(self)
        r = ""
        if unitari:
            l = latex(sqrt(self.length()))
            r = f"\\deufrac{{1}}{{{l}}}"
        return r + f"{s}"

class Matriu:
    #
    #
    #
    def __init__(self,matrix=eye(3)):
        f, c = matrix.shape
        self.files = f
        self.columnes = c
        self.matriu = matrix
        self.vaps = None
    #
    #
    #
    def set_vaps(self,vaps):
        self.vaps = vaps
    #
    #
    #
    @classmethod
    def aleatoria(cls,f=3,c=3,maxim=5,nuls=True):
        m = Matrix(f,c,lambda i, j : random.randint(-maxim,maxim))
        if not nuls:
            values = [i for i in range(1,maxim + 1)] + [-i for i in range(1,maxim + 1)]
            for i in range(f):
                for j in range(c):
                    if m[i,j] == 0:
                        m[i,j] = values[random.randint(0,2 * maxim - 1)]
        return cls(m)
    #
    #
    #
    def norma_maxim(self):
        return norma_maxim(self.matriu)
    #
    #
    #
    def nzeros(self):
        return nzeros(self.matriu)
    #
    #
    #
    def rank(self):
        return self.matriu.rank()
    #
    #
    #
    def determinant(self):
        return self.matriu.det()
    #
    #
    #
    def __add__(self,other):
        if not isinstance(other,Matriu):
            return None
        if self.files != other.files:
            return None
        if self.columnes != other.columnes:
            return None
        m = self.matriu + other.matriu
        return Matriu(m)
    #
    #
    #
    def __sub__(self,other):
        if not isinstance(other,Matriu):
            return None
        if self.files != other.files:
            return None
        if self.columnes != other.columnes:
            return None
        m = self.matriu - other.matriu
        return Matriu(m)
    #
    #
    #
    def __mul__(self,other):
        if isinstance(other,Matriu):
            if self.columnes != other.files:
                return None
            m = self.matriu * other.matriu
            return Matriu(m)
        if isinstance(other,Vector):
            if self.columnes != other.dimensio:
                return None
            u = Matrix(other.dimensio,1,other.components)
            v = self.matriu * u
            f, c = v.shape
            c = []
            for i in range(f):
                c.append(v[i,0])
            return Vector(c)
        return None
    #
    #
    #
    def __repr__(self):
        m = 1
        l = []
        for i in range(self.files):
            for j in range(self.columnes):
                if isinstance(self.matriu[i,j],Rational):
                    m * self.matriu[i,j].q
                elif isinstance(self.matriu[i,j],int):
                    pass
                else:
                    return matriu_latex(self.matriu)
        for i in range(self.files):
            for j in range(self.columnes):
                l.append(m * self.matriu[i,j])
        mcd = mcd_llista(l)
        f = Rational(mcd,m)
        p, q = f.p, f.q
        s = ""
        if q > 1:
            s = f"\\deufrac{{1}}{{{q}}}"
        l = [p * x for x in l]
        m = Matrix(self.files,self.columnes,l)
        return s + matriu_latex(m)
    #
    #
    #
    @classmethod
    def amb_rang(cls,f=3,c=3,r=3,maxim=5,nuls=True):
        trobat = False
        while not trobat:
            m = cls.aleatoria(f,c,maxim,nuls)
            if m.rank() != r:
                continue
            trobat = True
        return cls(m.matriu)
    #
    #
    #
    def inversa(self):
        return Matriu(self.matriu**(-1))
    #
    #
    #
    @classmethod
    def invertible(cls,ordre=3,maxim=5,mzeros=-1,unitaria=False):
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
        return cls(m)
    #
    #
    #
    @classmethod
    def diagonalitzable(cls,ordre=3,maxim=5,mzeros=-1,mvaps=3,vapsnuls=False):
        trobat = False
        while not trobat:
            random.seed()
            c = Matriu.invertible(ordre=ordre,maxim=3,mzeros=0,unitaria=True)
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
            a = c.matriu * d * c.matriu**(-1)
            if norma_maxim(a) > maxim:
                continue
            if mzeros >= 0 and nzeros(a) > mzeros:
                continue
            trobat = True
        m = cls(a)
        m.set_vaps(vaps)
        return m
    #
    #
    #
    @classmethod
    def gram(cls,ordre=3,maxim=5,mzeros=-1):
        trobat = False
        while not trobat:
            c = matriu_invertible(ordre,maxim=5)
            g = c.T * c
            if mzeros >= 0 and nzeros(g) > mzeros:
                continue
            if norma_maxim(g) > maxim:
                continue
            trobat = True
        return cls(g)
    #
    #
    #
    @classmethod
    def matriu_fila(cls,v):
        if not isinstance(v,Vector):
            return None
        m = Matrix(1,v.dimensio,v.components)
        return cls(m)
    #
    #
    #
    @classmethod
    def matriu_columna(cls,v):
        if not isinstance(v,Vector):
            return None
        m = Matrix(v.dimensio,1,v.components)
        return cls(m)
    #
    #
    #
    @classmethod
    def from_vectors_fila(cls,vecs):
        if len(vecs) == 0:
            return None
        if not isinstance(vecs[0],Vector):
            return None
        c = vecs[0].dimensio
        f = len(vecs)
        l = []
        for v in vecs:
            if not isinstance(v,Vector):
                return None
            if v.dimensio != c:
                return None
            l += v.components
        m = Matrix(f,c,l)
        return cls(m)
    #
    #
    #
    @classmethod
    def from_vectors_columna(cls,vecs):
        if len(vecs) == 0:
            return None
        if not isinstance(vecs[0],Vector):
            return None
        c = vecs[0].dimensio
        f = len(vecs)
        l = []
        for v in vecs:
            if not isinstance(v,Vector):
                return None
            if v.dimensio != c:
                return None
            l += v.components
        m = Matrix(f,c,l)
        return cls(m.T)
    #
    #
    #
    def vectors_columna(self,simplificar=False):
        vecs = []
        m = self.matriu
        for i in range(self.columnes):
            v = []
            for j in range(self.files):
                v.append(m[j,i])
            m = Vector(v)
            if simplificar:
                m.simplificar()
            vecs.append(m)
        return vecs
    #
    #
    #
    def vectors_fila(self,simplificar=False):
        vecs = []
        m = self.matriu
        for i in range(self.files):
            v = []
            for j in range(self.columnes):
                v.append(m[i,j])
            u = Vector(v)
            if simplificar:
                u.simplificar()
            vecs.append(u)
        return vecs
    #
    #
    #
    def nucli(self):
        n = self.matriu.nullspace()
        vecs = []
        for i in range(len(n)):
            m = Matriu(n[i])
            vecs += m.vectors_columna()
        return vecs


class EquacioLineal:
    #
    #
    #
    def __init__(self,eq,amp=True):
        self.equacio = eq
        self.amp = amp
        d = self.equacio.as_coefficients_dict()
        self.unknowns = d.keys()
    #
    #
    #
    @classmethod
    def coeficients(cls,a,b,amp=True):
        if not isinstance(a,Vector):
            return None
        if a.dimensio <= 4:
            x, y, z, t = symbols('x y z t')
            unknowns = [x,y,z,t]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
        eq = 0
        for k in range(a.dimensio):
            eq += a[k] * unknowns[k]
        return cls(eq - b,amp)
    #
    #
    #
    def set_coeficient_positiu(self,k):
        d = self.equacio.as_coefficients_dict()
        if d[k] < 0:
            self.equacio *= -1
    #
    #
    #
    def __repr__(self):
        d = self.equacio.as_coefficients_dict()
        l = list(d.values())
        m = 1
        for i in range(len(l)):
            if isinstance(l[i],Rational):
                m *= l[i].q
            elif isinstance(l[i],int):
                pass
            else:
                return ""
        v = [m * x for x in l]
        mcd = mcd_llista(v)
        factor = m // mcd
        eq = 0
        for k in d.keys():
            d[k] = factor * d[k]
            eq += d[k] * k
        eq -= d[1]
        if t in self.unknowns:
            str = mylatex(eq)
        else:
            str = latex(eq)
        if self.amp:
            str = f"{str} &= {-d[1]}"
        else:
            str = f"{str} = {-d[1]}"
        return str
    #
    #
    #
    def __add__(self,other):
        if not isinstance(other,EquacioLineal):
            return None
        eq = self.equacio + other.equacio
        return EquacioLineal(eq,self.amp or other.amp)
    #
    #
    #
    def __sub__(self,other):
        if not isinstance(other,EquacioLineal):
            return None
        eq = self.equacio - other.equacio
        return EquacioLineal(eq,self.amp or other.amp)
    #
    #
    #
    def __mul__(self,other):
        types = [Rational,float,int,Float,Pow,Add,Mul]
        for t in types:
            if isinstance(other,t):
                return EquacioLineal(other * self.equacio,self.amp)
        return None
    #
    #
    #
    __rmul__ = __mul__


class SistemaEquacions:
    #
    #
    #
    def __init__(self,a,b):
        if not isinstance(a,Matriu):
            return None
        if not isinstance(b,Vector):
            return None
        if a.files != b.dimensio:
            return None
        self.A = a
        self.B = b
        self.solucio = None
        eq = []
        files = a.vectors_fila()
        for k in range(self.A.files):
            eq.append(EquacioLineal.coeficients(files[k],b[k]))
        self.equacions = eq
        self.nombre = len(eq)
        if self.A.columnes <= 4:
            x, y, z, t = symbols('x y z t')
            unknowns = [x,y,z,t]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
        self.unknowns = unknowns[0:self.A.columnes]
    #
    #
    #
    @classmethod
    def from_equacions(cls,eqs):
        if len(eqs) <= 4:
            x, y, z, t = symbols('x y z t')
            unknowns = [x,y,z,t]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
        t = []
        vecs = []
        for e in eqs:
            d = e.as_coefficients_dict()
            c = []
            for i in unknowns:
                c.append(d[i])
            t.append(-d[1])
            vecs.append(Vector(c))
        a = Matriu.from_vectors_fila(vecs)
        b = Vector(t)
        return cls(a,b)
    #
    #
    #
    def __repr__(self):
        eqs = list(map(str,self.equacions))
        if len(eqs) == 0:
            return ""
        if len(eqs) == 1:
            return eqs[0].replace('&','')
        eqs = " \\\\ ".join(eqs)
        return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"
    #
    #
    #
    def resol(self):
        system = (self.A.matriu,Matrix(self.B.dimensio,1,self.B.components))
        self.solucio = list(linsolve(system,*self.unknowns))[0]
    #
    #
    #
    def solucio_latex(self):
        if len(self.solucio) == 0:
            return ""
        eqs = []
        for i in range(self.A.columnes):
            if self.solucio[i] == self.unknowns[i]:
                continue
            d = self.solucio[i].as_coefficients_dict()
            l = list(d.values())
            m = []
            for k in range(len(l)):
                if isinstance(l[k],Rational):
                    m.append(l[k].q)
                elif isinstance(l[k],int):
                    pass
                else:
                    return ""
            if len(m) == 0:
                factor = 1
            else:
                print(m)
                factor = mcm_llista(m)
                print(factor)
            eq = 0
            for k in d.keys():
                d[k] = factor * d[k]
                eq += d[k] * k
            if factor == 1:
                eqs.append(f"{self.unknowns[i]} &= {mylatex(eq)}")
            else:
                eqs.append(f"{self.unknowns[i]} &= \\frac{{{mylatex(eq)}}}{{{factor}}}")
        if len(eqs) == 1:
            return eqs[0].replace('','')
        eqs = " \\\\ ".join(eqs)
        return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"

class EquacioParametrica:
    #
    #
    #
    def __init__(self,eq,amp=True):
        x, y, z, t = symbols('x y z t')
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        unknowns = [x1,x2,x3,x4,x5,x6,x7,x8] + [x,y,z,t]
        self.equacio = eq
        self.amp = amp
        d = eq.as_coefficients_dict()
        self.b = 0
        self.coefs = {}
        for k in d.keys():
            if k == 1:
                self.b = d[k]
            elif k in unknowns:
                self.unknown = k
            else:
                self.coefs[k] = d[k]
    #
    #
    #
    @classmethod
    def coeficients(cls,a,b,p=0,total=1,amp=True):
        if not isinstance(a,Vector):
            return None
        if p >= total:
            return None
        if total <= 4:
            x, y, z, t = symbols('x y z t')
            unknowns = [x,y,z,t]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
        t1, t2, t3, t4, t5, t6, t7, t8 = symbols('t1 t2 t3 t4 t5 t6 t7 t8')
        parameters = [t1,t2,t3,t4,t5,t6,t7,t8]
        eq = - unknowns[p]
        cls.unknown = unknowns[p]
        cls.params = parameters[0:a.dimensio]
        for k in range(a.dimensio):
            eq += a[k] * parameters[k]
        eq += b
        return cls(eq,amp)
    #
    #
    #
    def __repr__(self):
        t1, t2, t3, t4, t5, t6, t7, t8 = symbols('t1 t2 t3 t4 t5 t6 t7 t8')
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        x, y, z, t = symbols('x y z t')
        if self.amp:
            s = f"{self.unknown} &= "
        else:
            s = f"{self.unknown} = "
        d = []
        if isinstance(self.b,Rational):
            d.append(self.b.q)
        for k,v in self.coefs.items():
            if isinstance(v,Rational):
                d.append(v.q)
        if len(d) == 0:
            factor = 1
        else:
            factor = mcm_llista(d)
        if factor != 1:
            s += f"\\frac{{NUMERADOR}}{{{factor}}}"
        else:
            s += "NUMERADOR"
        e = 0
        first = None
        for k in self.params:
            try:
                self.coefs[k]
            except:
                continue
            if first is None and self.coefs[k] != 0:
                first = self.coefs[k]
            e += factor * self.coefs[k] * k
        f = ""
        if first is not None and self.b != 0:
            if first > 0:
                f = f"{factor*self.b} + "
            else:
                f = f"{factor*self.b} "
        s = s.replace('NUMERADOR',f + latex(e))
        return s

class EquacionsParametriques:
    #
    #
    #
    def __init__(self,a,b,amp=True):
        if not isinstance(a,Matriu):
            return None
        if not isinstance(b,Vector):
            return None
        if a.files != b.dimensio:
            return None
        self.A = a
        self.B = b
        eq = []
        files = a.vectors_fila()
        for k in range(self.A.files):
            eq.append(EquacioParametrica.coeficients(files[k],b[k],k,self.A.files,amp))
        self.equacions = eq
        self.nombre = len(eq)
    #
    #
    #
    def __repr__(self):
        l = list(map(str,self.equacions))
        eqs = " \\\\ ".join(l)
        return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"
    #
    #
    #
    def eliminar_parametres(self):
        L, U, _ = self.A.matriu.LUdecomposition()
        r = U.rank()
        v = Vector([e.unknown for e in self.equacions])
        t = Matrix(self.nombre,1,(v - self.B).components)
        t = (L**(-1) * t)[r:]
        return SistemaEquacions.from_equacions(t)


class PlaVectorial:
    #
    #
    #
    def __init__(self,u1,u2):
        if not isinstance(u1,Vector):
            return None
        if not isinstance(u2,Vector):
            return None
        self.u1 = u1
        self.u2 = u2
    #
    #
    #
    def __repr__():
        return f"{self.u1}, {self.u2}"
    #
    #
    #
    def equacio_implicita(self):
        w = self.u1.cross(self.u2)
        return EquacioLineal.coeficients(w,0,False)
    #
    #
    #
    def base_ortogonal(self):
        v1 = self.u1
        u2 = self.u2
        v2 = v1.dot(v1) * u2 - v1.dot(u2) * v1
        v2.simplificar()
        return [v1,v2]

class PlaAfi:
    #
    #
    #
    def __init__(self,u1,u2,p):
        if not isinstance(u1,Vector):
            return None
        if not isinstance(u2,Vector):
            return None
        if not isinstance(p,Vector):
            return None
        self.u1 = u1
        self.u2 = u2
        self.p = p
    #
    #
    #
    def equacio_implicita(self):
        w = self.u1.cross(self.u2)
        return EquacioLineal.coeficients(w,w.dot(self.p),False)
    #
    #
    #
    def base_ortogonal(self):
        v1 = self.u1
        u2 = self.u2
        v2 = v1.dot(v1) * u2 - v1.dot(u2) * v1
        v2.simplificar()
        return [v1,v2]

class RectaVectorial:
    #
    #
    #
    def __init__(self,u):
        if not isinstance(u,Vector):
            return None
        self.u = u
        self.set_equacions_implicites()
    #
    #
    #
    def equacio_continua(self):
        self.u.simplificar()
        x, y, z = symbols('x y z')
        incg = [x,y,z]
        eq = []
        for i in range(3):
            if self.u[i] == 1:
                eq.append(latex(incg[i]))
            else:
                eq.append(f"\\frac{{{latex(incg[i])}}}{{{v[i]}}}")
        return " = ".join(eq)
    #
    #
    #
    def set_equacions_implicites(self):
        a = Matriu.matriu_fila(self.u)
        l = a.nucli()
        a = Matriu.from_vectors_fila(l)
        b = Vector.nul(len(l))
        self.implicites = SistemaEquacions(a,b)
    #
    #
    #
    def equacions_implicites(self):
        m = Matriu.amb_rang(f=2,c=2,r=2,maxim=3,nuls=False)
        a = m * self.implicites.A
        b = Vector.nul(2)
        return SistemaEquacions(a,b)

class RectaAfi:
    #
    #
    #
    def __init__(self,u,p):
        if not isinstance(u,Vector):
            return None
        if not isinstance(p,Vector):
            return None
        self.u = u
        self.p = p
        self.set_equacions_implicites()
    #
    #
    #
    def equacio_continua(self):
        self.u.simplificar()
        x, y, z = symbols('x y z')
        incg = [x,y,z]
        eq = []
        for i in range(3):
            if self.u[i] == 1:
                eq.append(latex(incg[i]-p[i]))
            else:
                eq.append(f"\\frac{{{latex(incg[i]-p[i])}}}{{{v[i]}}}")
        return " = ".join(eq)
    #
    #
    #
    def set_equacions_implicites(self):
        a = Matriu.matriu_fila(self.u)
        l = a.nucli()
        a = Matriu.from_vectors_fila(l)
        b = a * self.p
        self.implicites = SistemaEquacions(a,b)
    #
    #
    #
    def equacions_implicites(self):
        m = Matriu.amb_rang(f=2,c=2,r=2,maxim=3,nuls=False)
        a = m * self.implicites.A
        b = m * self.implicites.B
        return SistemaEquacions(a,b)

class SubespaiVectorial:
    #
    #
    #
    def __init__(self,vecs):
        self.generadors = [v.simplificar() for v in vecs]
        a = Matriu.from_vectors_columna(vecs)
        if a is None:
            return None
        L, U, _ = a.matriu.LUdecomposition()
        files = Matriu(U).vectors_fila()
        self.base = []
        for f in files:
            k = primer_no_nul(f.components)
            if k is not None:
                self.base.append(vecs[k])
        self.dimensio = len(self.base)
    #
    #
    #
    def equacions_implicites(self):
        if len(self.base) == 0:
            return None
        a = Matriu.from_vectors_columna(self.base)
        b = Vector.nul(self.base[0].dimensio)
        p = EquacionsParametriques(a,b)
        return (p.eliminar_parametres())
    #
    #
    #
    def suplementari(self):
        a = Matriu.from_vectors_fila(self.base)
        return SubespaiVectorial(a.nucli())
    #
    #
    #
    def base_ortogonal(self):
        if self.dimensio == 0:
            return None
        if len(self.base) == 1:
            return self.base
        L = []
        for b in self.base:
            a = Matriu.from_vectors_columna([b])
            L.append(a.matriu)
        Q = GramSchmidt(L)
        base = []
        for m in Q:
            m = Matriu(m)
            v = m.vectors_columna()[0]
            v.simplificar()
            base.append(v)
        return base
