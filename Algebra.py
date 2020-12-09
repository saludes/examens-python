#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Filename:   Algebra.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2020
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
    """
    Retorna el màxim comú divisor d'una llista d'enters
    """
    return reduce(gcd, list)

def mcm_llista(list):
    """
    Retorna el mínim comú múltiple d'una llista d'enters
    """
    if len(list) == 0:
        return 1
    mcm = list[0]
    for i in list[1:]:
        mcm = mcm * i // gcd(mcm, i)
    return mcm

def mti(i,j):
    """
    Funció auxiliar per crear una matriu triangular superior i uns
    o menys uns a la diagonal
    Retorna zero si el coeficients està per sota de la diagonal principal
    """
    values = (-1,1)
    if i > j:
        return 0
    elif i == j:
        return values[random.randint(0,1)]
    else:
        return random.randint(-1,2)

def mts(i,j,values):
    """
    Funció auxiliar per crear una matriu triangular inferior i uns
    o menys uns a la diagonal
    Retorna zero si el coeficients està per sobre de la diagonal principal
    """
    if i < j:
        return 0
    elif i == j:
        return  values[i]
    else:
        return random.randint(-2,1)

def norma_maxim(m):
    """
    Retorna el màxim en valor absolut d'entre els coeficients d'una matriu
    del tipus Matrix del sympy
    """
    f, c = m.shape
    n = 0
    for i in range(f):
        for j in range(c):
            if abs(m[i,j]) > n:
                n = abs(m[i,j])
    return n

def nzeros(m):
    """
    Retorna el nombre de zeros d'una matriu del tipus Matrix del sympy
    """
    f, c = m.shape
    z = 0
    for i in range(f):
        for j in range(c):
            if m[i,j] == 0:
                z += 1
    return z

def matriu_latex(m,format=None):
    """
    Retorna l'expressió en latex d'una matriu del tipus Matrix del sympy
    """
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
    """
    Retorna l'índel del primer coeficient no nul d'una llista
    """
    if list is None or len(list) == 0:
        return None
    for i in range(len(list)):
        if list[i] != 0:
            return i
    return None

class Vector:
    """
    Classe que ens permetrà representar vectors i punts
    Atributs:
        dimensio: el nombre de components o longitud del vector
        components: llista amb les components del vector
    """
    #
    #
    #
    def __init__(self,c):
        """
        Constructor.
        Paràmetres:
           c: Llista de nombres
        """
        self.dimensio = len(c)
        self.components = list(c)
    #
    #
    #
    @classmethod
    def nul(cls,dim):
        """
        Retorna el vector nul de longitud dim
        """
        l = [0 for i in range(dim)]
        return cls(l)
    #
    #
    #
    @classmethod
    def aleatori(cls,l=3,maxim=5,nuls=True,positius=False):
        """
        Retorna un vector aleatori de longitud l
        Paràmetres:
           l: longitud del vector
           maxim: Nombre màxin que pot contenir en valor absolut
           nuls: Si pot contenir el valor 0 o no
           positius: Si els coeficients han de set tots positius
        """
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
    def factor_comu(self):
        """
        En cas que totes les components del vector sigui enteres o racionals,
        torna el racional que es pot treure factor comú i el vector simplificat.
        Si hi ha components que no són ni enteres ni racionals, torna 1 i el mateix
        vector
        """
        d = []
        other = False
        for k in self.components:
            if isinstance(k,Rational):
                d.append(k.q)
            elif (isinstance(k,int) or isinstance(k,Integer)):
                pass
            else:
                other = True
        if other:
            return (1,Vector(self.components))
        mcm = mcm_llista(d)
        v = [mcm * k for k in self.components]
        mcd = mcd_llista(v)
        v = [k // mcd for k in v]
        f = Rational(mcd,mcm)
        if f < 0:
            f *= -1
            v = [-k for k in self.components]
        return (f,Vector(v))
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex del vector
        """
        d = []
        other = False
        for k in self.components:
            if isinstance(k,Rational):
                d.append(k.q)
            elif (isinstance(k,int) or isinstance(k,Integer)):
                pass
            else:
                other = True
        if other:
            s = ",".join([latex(k) for k in self.components])
            return f"\\left({s}\\right)"
        mcm = mcm_llista(d)
        v = [mcm * k for k in self.components]
        mcd = mcd_llista(v)
        v = [k // mcd for k in v]
        f = Rational(mcd,mcm)
        if f == 1:
            s = ",".join([latex(k) for k in v])
            return f"({s})"
        elif f == -1:
            s = ",".join([latex(-k) for k in v])
            return f"({s})"
        else:
            s = ",".join([latex(k) for k in v])
            return f"\\deufrac{{{f.p}}}{{{f.q}}}({s})"
    #
    #
    #
    def __add__(self,other):
        """
        Defineix la suma de vectors. Per exemple:
        u1 = Vector([3,2,1,3])
        u2 = Vector[-2,4,-3,1])
        v = u1 + u2
        """
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
        """
        Defineix la resta de vectors.Per exemple:
        u1 = Vector([3,2,1,3])
        u2 = Vector[-2,4,-3,1])
        v = u1 - u2
        """
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
        """
        Defineix el producte d'un escalar per un vector i el
        producte d'una vector per una matriu. Per exemple,
        u1 = Vector([3,2,1,3])
        u2 = Vector[-2,4,-3,1])
        v = 5*u1 - 4*u2
        """
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
        """
        Permet indexar els elements d'un vector. Per exemple:
        v = Vector([3,2,1,3,5])
        x = v[2] + v[3]
        """
        try:
            return self.components[i]
        except:
            return None
    #
    #
    #
    def __setitem__(self,i,value):
        """
        Permet assignar valors mitjançant índexs. Per exemple:
        v = Vector.aleatori(l=4)
        v[0] = 3
        """
        self.components[i] = value
    #
    #
    #
    def dot(self,other):
        """
        Retorna el producte escalar de dos vectors: Per exemple:
        u = Vector.aleatori(l=3)
        v = Vector.aleatori(l=3)
        p = u.dot(v)
        """
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
        """
        Retorna la longitud o mòdul del vector
        """
        s = 0
        for i in range(self.dimensio):
            s += self.components[i]**2
        return sqrt(s)
    #
    #
    #
    def simplificar(self):
        """
        Simplifica el vector, és a dir, converteix les seves components en una
        llista d'enters amb mcd igual a 1.
        Només té sentit si totes les components del vector són nombres enters
        o racionals
        """
        d = []
        for i in range(self.dimensio):
            if isinstance(self.components[i],Rational):
                d.append(self.components[i].q)
            elif isinstance(self.components[i],int):
                pass
            elif isinstance(self.components[i],Integer):
                pass
            else:
                return
        mcm = mcm_llista(d)
        v = [mcm * x for x in self.components]
        mcd = mcd_llista(v)
        v = [x // mcd for x in v]
        if v[0] < 0:
            v = [-x for x in v]
        self.components = v
    #
    #
    #
    def cross(self,other,simplificar=False):
        """
        Retorna un nou vector que és el producte vectorial de dos vectors:
        Si simplificar és True, simplifica el vector resultant
        Per exemple:
        u = Vector.aleatori(l=3)
        v = Vector.aleatori(l=3)
        p = u.cross(v)
        """
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
        """
        Retorna l'expressió en latex del vector. Si unitari és True,
        divideix per la seva longitud
        """
        if unitari:
            f,v = self.factor_comu()
            l = v.length()
            if l == 1:
                return str(v)
            return f"\\deufrac{{1}}{{{latex(l)}}}{str(v)}"
        return str(self)
    #
    #
    #
    def components_en_base(self,base=None):
        """
        Retorna un nou vector amb les components del vector (donem per fet
        que estan en la base canònica) en la base "base"
        """
        if base is None:
            return Vector(self.components)
        if not isinstance(base,Base):
            return None
        if self.dimensio != base.dimensio:
            return None
        c = base.matriu()
        return c.inversa() * self
    #
    #
    #
    def reordena_aleatoriament(self):
        """
        Retorna un nou vector amb les components reordenades aleatoriament
        """
        p = list(range(self.dimensio))
        random.shuffle(p)
        r = [self.components[i] for i in p]
        return Vector(r)
    #
    #
    #
    def nzeros(self):
        """
        Retorna el nombre de zeros del vector
        """
        n = 0
        for i in self.components:
            if i == 0:
                n += 1
        return n

class Base:
    """
    Classe que ens permetrà representar bases de R^n
    Atributs:
        vecs: una llista de n vectors de R^n
        dimensio: el valor de n
        unitaria: En funcio de si volem imprimir els seus
                  vectors unitaris o no
    """
    #
    #
    #
    def __init__(self,vecs,unitaria=False):
        """
        Contructor.
        Paràmetres:
           vecs: llista de vectors
           unitaria: valor de unitaria.
                     Serveix per si volem imprimir o fer servir els seus
                     vectors com a unitaris
        """
        self.unitaria = unitaria
        if len(vecs) == 0:
            return None
        if not isinstance(vecs[0],Vector):
            return None
        d = vecs[0].dimensio
        for v in vecs:
            if not isinstance(v,Vector):
                return None
            if v.dimensio != d:
                return None
        self.vecs = vecs
        self.dimensio = d
        m = Matriu.from_vectors_columna(self.vecs)
        if m.columnes != d:
            return None
        if m.rank() != d:
            return None
    #
    #
    #
    def matriu(self):
        """
        Retorna la matriu de la classe Matriu que té per coulumnes
        els vectors de la base
        """
        if not self.unitaria:
            return Matriu.from_vectors_columna(self.vecs)
        unitaris = [(1 / v.length()) * v for v in self.vecs]
        return Matriu.from_vectors_columna(unitaris)
    #
    #
    #
    @classmethod
    def from_matriu(cls,m):
        """
        Crea una nova base a partir d'una matriu de la classe Matriu
        Si la matriu no és quadrada o no té rang màxim, retorna None
        """
        if m.files != m.columnes:
            return None
        if m.rank() != m.columnes:
            return None
        return cls(m.vectors_columna())
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex de la base
        """
        if not self.unitaria:
            base = ",".join(map(str,self.vecs))
            return f"\\{{{base}\\}}"
        base = ",".join([v.latex(True) for v in self.vecs])
        return f"\\left\\{{{base}\\right\\}}"

class Matriu:
    """
    Classe que ens permetrà representar matrius. El problema de la classe
    Matrix del sympy és que només es poden multiplicar per elements del tipus
    Matrix.
    Ens insteressa poder multiplicar Matrius per Vectors
    Atributs:
        dimensio: nombre de files de la matriu
        columnes: nombre de columnes de la matriu
        matriu: matriu de la classe Matrix del sympy

        Només s'utilitzen quen generem una matriu diagonalitzble
          vaps: llista de vectors propis de la matriu
          veps: llista de vectors propis de la matriu
    """
    #
    #
    #
    def __init__(self,matrix=eye(3)):
        """
        Constructor.
        Paràmetres:
          matrix: matriu del tipus Matrix de sympy
                  per defecte és la matriu unitat d'ordre 3
        """
        f, c = matrix.shape
        self.files = f
        self.columnes = c
        self.matriu = matrix
        self.vaps = None
        self.veps = None
    #
    #
    #
    def set_vaps(self,vaps):
        """
        Assigna un llista de valors propis a la variable self.vaps
        """
        self.vaps = vaps
    #
    #
    #
    def set_veps(self,veps):
        """
        Assigna un llista de vectors propis simplificats a la variable self.veps
        """
        for v in veps:
            v.simplificar()
        self.veps = veps
    #
    #
    #
    @classmethod
    def aleatoria(cls,f=3,c=3,maxim=5,nuls=True):
        """
        Genera una matriu aleatoria.
        Paràmetres:
          f: nombre de files de la matriu
          c: nombre de columnes de la matriu
          maxim: tots els elements tindran valor absolut menor que "maxim"
          nuls: la matriu pot contenir coeficients nuls o no
        """
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
        """
        Retorna la norma del màxim de la matriu
        """
        return norma_maxim(self.matriu)
    #
    #
    #
    def nzeros(self):
        """
        Retorna el nombre de zeros de la matriu
        """
        return nzeros(self.matriu)
    #
    #
    #
    def rank(self):
        """
        Retorna el rang de la matriu
        """
        return self.matriu.rank()
    #
    #
    #
    def determinant(self):
        """
        Retorna el determiant de la matriu
        """
        return self.matriu.det()
    #
    #
    #
    def polinomi_caracteristic(self):
        """
        Retorna el polinomi característic de la matriu
        """
        lamda = symbols('lamda')
        p = self.matriu.charpoly(lamda)
        return p
    #
    #
    #
    def __add__(self,other):
        """
        Defineix la suma de matrius. Per exemple:
        a1 = Matriu.aleatoria()
        a2 = Matriu.aleatoria())
        b = a1 + a2
        """
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
        """
        Defineix la resta de matrius. Per exemple:
        a1 = Matriu.aleatoria()
        a2 = Matriu.aleatoria())
        b = a1 - a2
        """
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
    def __rmul__(self,other):
        """
        Defineix el producte d'un escalar per un vector i el
        producte d'una vector per una matriu. Per exemple:
        a = Matriu.aleatoria(f=4,c=3,maxim=7,nuls=False)
        b = 4 * a
        """
        types = [Rational,float,int,Float,Pow,Add,Mul]
        for t in types:
            if isinstance(other,t):
                return Matriu(other * self.matriu)
        return None
    #
    #
    #
    def __mul__(self,other):
        """
        Defineix el producte de matrius. Per exemple:
        a1 = Matriu.aleatoria()
        a2 = Matriu.aleatoria())
        b = a1 * a2
        """
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
        """
        Retorna l'expressió en latex de la matriu
        """
        m = 1
        l = []
        for i in range(self.files):
            for j in range(self.columnes):
                if isinstance(self.matriu[i,j],Rational):
                    l.append(self.matriu[i,j].q)
                elif isinstance(self.matriu[i,j],int):
                    pass
                else:
                    return matriu_latex(self.matriu)
        m = mcm_llista(l)
        l = []
        for i in range(self.files):
            for j in range(self.columnes):
                l.append(m * self.matriu[i,j])
        mcd = mcd_llista(l)
        f = Rational(mcd,m)
        p, q = f.p, f.q
        s = ""
        if q > 1:
            s = f"\\deufrac{{1}}{{{q}}}"
        m = Matrix(self.files,self.columnes,l)
        return s + matriu_latex(m)
    #
    #
    #
    @classmethod
    def diagonal(cls,vals):
        """
        Genera una matriu diagonal amb valors "vals" a la diagonal
        """
        d = len(vals)
        if d == 0:
            return None
        m = diag(*vals)
        return cls(m)
    #
    #
    #
    @classmethod
    def amb_rang(cls,f=3,c=3,r=3,maxim=5,nuls=True):
        """
        Retorna una matriu aleatoria amb rang r.
        Paràmetres:
          f: nombre de files de la matriu
          c: nombre de columnes de la matriu
          r: rang de la matriu
          maxim: tots els elements tindran valor absolut menor que "maxim"
          nuls: la matriu pot contenir coeficients nuls o no
        """
        trobat = False
        while not trobat:
            m = cls.aleatoria(f,c,maxim,nuls)
            trobat = m.rank() == r
        return cls(m.matriu)
    #
    #
    #
    def inversa(self):
        """
        Retorna una nova matriu que és la inversa de l'actual
        """
        return Matriu(self.matriu**(-1))
    #
    #
    #
    @classmethod
    def invertible(cls,ordre=3,maxim=5,mzeros=-1,unitaria=False):
        """
        Retorna una matriu quadrada aleatoria invertible.
        Paràmetres:
          ordre: nombre de files i columnes de la matriu
          maxim: tots els elements tindran valor absolut menor que "maxim"
          mzeros: si nzeros >= 0, nombre màxim de zeros que tindrà la matrius
                  si nzeros < 0, el nombre de zeros no està limitat
          unitaria: si volem que el determinant sigui 1 o -1
        """
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
            if norma_maxim(m) > maxim:
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
        """
        Retorna una matriu quadrada aleatoria diagonalitzable.
        Paràmetres:
          ordre: nombre de files i columnes de la matriu
          maxim: tots els elements tindran valor absolut menor o igual que "maxim"
          mzeros: si nzeros >= 0, nombre màxim de zeros que tindrà la matrius
                  si nzeros < 0, el nombre de zeros no està limitat
          mvaps: tots els valors propis tindran valor absolut menor o igual que "mvaps"
          vapsnuls: si hi pot aparèixer el valor propi nul
        """
        trobat = False
        while not trobat:
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
        m.set_veps(c.vectors_columna())
        return m
    #
    #
    #
    @classmethod
    def gram(cls,ordre=3,maxim=5,mzeros=-1):
        """
        Retorna una matriu quadrada aleatoria que serà d'un producte escalar.
        Paràmetres:
          ordre: nombre de files i columnes de la matriu
          maxim: tots els elements tindran valor absolut menor o igual que "maxim"
          mzeros: si nzeros >= 0, nombre màxim de zeros que tindrà la matrius
                  si nzeros < 0, el nombre de zeros no està limitat
        """
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
        """
        Retorna una nova matriu fila a partir de les components del vector v
        """
        if not isinstance(v,Vector):
            return None
        m = Matrix(1,v.dimensio,v.components)
        return cls(m)
    #
    #
    #
    @classmethod
    def matriu_columna(cls,v):
        """
        Retorna una nova matriu columna a partir de les components del vector v
        """
        if not isinstance(v,Vector):
            return None
        m = Matrix(v.dimensio,1,v.components)
        return cls(m)
    #
    #
    #
    @classmethod
    def from_vectors_fila(cls,vecs):
        """
        Retorna una nova matriu a partir d'una llista de vectors.
        Les components dels vectors seran les files de la nova matriu
        """
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
        """
        Retorna una nova matriu a partir d'una llista de vectors.
        Les components dels vectors seran les columnes de la nova matriu
        """
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
        """
        Retorna una llista amb els vectors columna de la matriu
        """
        vecs = []
        m = self.matriu
        for i in range(self.columnes):
            v = [m[j,i] for j in range(self.files)]
            u = Vector(v)
            if simplificar:
                u.simplificar()
            vecs.append(u)
        return vecs
    #
    #
    #
    def vectors_fila(self,simplificar=False):
        """
        Retorna una llista amb els vectors fila de la matriu
        """
        vecs = []
        m = self.matriu
        for i in range(self.files):
            v = [m[i,j] for j in range(self.columnes)]
            u = Vector(v)
            if simplificar:
                u.simplificar()
            vecs.append(u)
        return vecs
    #
    #
    #
    def nucli(self):
        """
        Retorna una llista de vectors que formen una base del nucli de la matriu
        """
        n = self.matriu.nullspace()
        vecs = []
        for i in range(len(n)):
            m = Matriu(n[i])
            vecs += m.vectors_columna()
        return vecs
    #
    #
    #
    def reordena_aleatoriament_columnes(self):
        """
        Retorna una nova matriu amb les columnes reordenades aleatòriament
        """
        c = self.vectors_columna()
        p = list(range(self.columnes))
        random.shuffle(p)
        d = [c[i] for i in p]
        return Matriu.from_vectors_columna(d)
    #
    #
    #
    def reordena_aleatoriament_files(self):
        """
        Retorna una nova matriu amb les files reordenades aleatòriament
        """
        c = self.vectors_fila()
        p = list(range(self.files))
        random.shuffle(p)
        d = [c[i] for i in p]
        return Matriu.from_vectors_fila(d)
    #
    #
    #
    def inserta_fila(self,pos,fila):
        """
        Retorna una nova matriu amb la fila "fila" insertada a la posició "pos"
        """
        if not isinstance(fila,Vector):
            return None
        if fila.dimensio != self.columnes:
            return None
        if pos > self.files:
            pos = self.files
        files = self.vectors_fila()
        files.insert(pos,fila)
        return Matriu.from_vectors_fila(files)
    #
    #
    #
    def inserta_columna(self,pos,columna):
        """
        Retorna una nova matriu amb la columna "columna" insertada a la posició "pos"
        """
        if not isinstance(fila,Vector):
            return None
        if columna.dimensio != self.files:
            return None
        if pos > self.columnes:
            return None
        columnes = self.vectors_fila()
        columnes.insert(pos,columna)
        return Matriu.from_vectors_columna(columnes)

class EquacioLineal:
    """
    Classe per treballar amb equacions lineals
    """
    #
    #
    #
    def __init__(self,eq,amp=True,prime=0):
        """
        Constructor.
        Paràmentres:
           eq: expressió lineal que ha de contenir tots els termes, aleshores
               l'equació serà "eq = 0". Només guardem la part "eq"
           amp: quan escrivim l'equació en latex ha d'aparèixer &= o només =
           prime: nombre de primes que s'han de posar a les incògnites

         Per exemple:
         x, y, z, t = symbols('x y z t')
         eq = 2*x-3*y+4*z-3*t-4
         e = EquacioLineal(eq)

         Atributs:
           equacio: terme de l'esquerra en la equacio "eq = 0"
           unknowns: incògnites que apareixen a l'equació
           amp: True o False
           prime: nombre de primes que escriurem a l'equaqció
        """
        self.equacio = eq
        self.amp = amp
        d = self.equacio.as_coefficients_dict()
        self.unknowns = d.keys()
        self.prime = prime
    #
    #
    #
    @classmethod
    def coeficients(cls,a,b,amp=True,prime=0):
        """
        Retorna una nova equació amb coeficients de les incògnites el vector "a" i
        terme independent b
        """
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
        return cls(eq - b,amp,prime)
    #
    #
    #
    def set_coeficient_positiu(self,k):
        """
        Si el coeficient de "k" és negatiu, canvia de signe tota l'equació, de
        manera que el coeficients de "k" passa a ser positiu
        """
        d = self.equacio.as_coefficients_dict()
        if d[k] < 0:
            self.equacio *= -1
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex de l'equació.
        Si els coeficients són enters o racionals, treu el denominador comú
        """
        other = False
        d = self.equacio.as_coefficients_dict()
        l = list(d.values())
        m = []
        for k in l:
            if isinstance(k,Rational):
                m.append(k.q)
            elif isinstance(k,int) or isinstance(k,Integer):
                pass
            else:
                other = True
        if not other:
            mcm = mcm_llista(m)
            v = [mcm * x for x in l]
            mcd = mcd_llista(v)
            factor = Rational(mcm,mcd)
            eq = 0
            for k in d.keys():
                d[k] = factor * d[k]
                eq += d[k] * k
        else:
            eq = self.equacio
        eq -= d[1]
        if t in self.unknowns:
            str = mylatex(eq)
        else:
            str = latex(eq)
        if self.amp:
            str = f"{str} &= {-d[1]}"
        else:
            str = f"{str} = {-d[1]}"
        if self.prime > 0:
            s = self.prime * "'"
            for i in self.unknowns:
                str = str.replace(latex(i),latex(i) + s)
        return str
    #
    #
    #
    def __add__(self,other):
        """
        Serveix per sumar equacions
        """
        if not isinstance(other,EquacioLineal):
            return None
        eq = self.equacio + other.equacio
        return EquacioLineal(eq,self.amp or other.amp)
    #
    #
    #
    def __sub__(self,other):
        """
        Serveix per restar equacions
        """
        if not isinstance(other,EquacioLineal):
            return None
        eq = self.equacio - other.equacio
        return EquacioLineal(eq,self.amp or other.amp)
    #
    #
    #
    def __mul__(self,other):
        """
        Serveix per multiplicar un escalar per una equació
        """
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
    """
    Classe per treballar amb sistemes d'equacions lineals
    """
    #
    #
    #
    def __init__(self,a,b):
        """
        Constructor.
        Paràmetres:
          a: matriu dels coeficients de les incògnites
          b: Termes independents

        Atributs:
          A: matriu dels coeficients de les incógnites
          B: vector de termes independents
          equacions: llista de EquacioLineal
          nombre: nombre d'equacions
          solucio: solucio del sistema d'equacions
          unknowns: llista d'incògnites
        """
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
        """
        Retorna un sistema d'equacions amb equacions "eqs"
        Paràmetre:
          eqs: llista de EquacioLineal
        """
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
        """
        Retorna l'expressió en latex del sistema d'equacions
        """
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
        """
        Resol el sistema d'equacions utilitzant la funció linsolve del sympy.
        El resultat és una llista d'expressions on hi poden aparèixer les
        incògnites del sistema com a paràmetres
        """
        system = (self.A.matriu,Matrix(self.B.dimensio,1,self.B.components))
        self.solucio = list(linsolve(system,*self.unknowns))[0]
    #
    #
    #
    def solucio_latex(self):
        """
        Retorna l'expressió en latexs de la solució del sistema d'equacions
        """
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
                factor = mcm_llista(m)
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
    """
    Classe per treballar amb equacions paramètriques
    """
    #
    #
    #
    def __init__(self,eq,amp=True):
        """
        Contructor.
        Paràmetres:
          eq: equació paramètrica. Ha de ser del tipus
              -x + 2*t1 - 3*t2 + t3 - 4
             amb el signe menys a la incògnita
          amp: True o False en funció si hem d'escriure &= o només = en la representació
               en latex de l'equació

        Atributs:
          equacio: l'equació paramètrica
          b: terme independent de l'equaqció
          coefs: coeficients dels paràmetres
          unknown: incògina de l'equació paramètrica
        """
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
        """
        Genera una equació amb coeficients dels paràmtres el vector "a", terme
        independent b i incògnita número p d'un total de "total". Per exemple:

           e = EquacioParametrica(Vector([3,-2,1],5,1,4) genera l'equació
               - y + 3*t1 - 2*t2 + t3 + 5

           e = EquacioParametrica(Vector([3,-2,1,3,7],-5,3,7) genera l'equació
               - x6 + 3*t1 - 2*t2 + t3 + 3*t4 + 7*t5 - 5
        """
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
        """
        Retorna l'expessió en latex de l'equació paramètrica
        """
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
    """
    Classe per treballar amb sistemes d'equacions paramètriques
    """
    #
    #
    #
    def __init__(self,a,b,amp=True):
        """
        Contructor.
        Genera les equacions X = b + a.T on X són les inognites i T els paràmetres

        Paràmetres:
          a: matriu dels coeficients dels paràmetres
          b: vector de termes independents
          amp: True o False en funció si hem d'escriure &= o només = en la representació
               en latex del sistema

        Atributs:
          A: matriu dels coeficients dels paràmetres
          B: vector dels termes independents
          equacions: llista de EquacioParametrica
          noimbre: nombre d'equacions
        """
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
        """
        Retorna la representació en latex del sistema d'equacions paramètriques
        """
        l = list(map(str,self.equacions))
        eqs = " \\\\ ".join(l)
        return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"
    #
    #
    #
    def eliminar_parametres(self):
        """
        Retorna el SistemaEquacions que s'obté en eliminar els paràmetres dels
        sistema
        """
        L, U, _ = self.A.matriu.LUdecomposition()
        r = U.rank()
        v = Vector([e.unknown for e in self.equacions])
        t = Matrix(self.nombre,1,(v - self.B).components)
        t = (L**(-1) * t)[r:]
        return SistemaEquacions.from_equacions(t)

class PlaVectorial:
    """
    Classe per treballar amb plans vectorials
    """
    #
    #
    #
    def __init__(self,u1,u2):
        """
        Constructor.
        Paràmetres:
          u1, u2: generadors del pla
        """
        if not isinstance(u1,Vector):
            return None
        if not isinstance(u2,Vector):
            return None
        if u1.dimensio != 3 or u2.dimensio != 3:
            return None
        m = Matriu.from_vectors_columna([u1,u2])
        if m.rank() != 2:
            return None
        self.u1 = u1
        self.u2 = u2
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex dels dos vectors generadors
        """
        return f"{self.u1}, {self.u2}"
    #
    #
    #
    def equacio_implicita(self,base=None,prime=0):
        """
        Retorna l'equació implícita del pla en la base "base" en format latex.
        Normalment, si la base no és la canònica posarem prime > 0 perquè el resultat
        sigui de l'estil 2x'-3y'+4z'=0
        """
        if base is None:
            w = self.u1.cross(self.u2)
            return EquacioLineal.coeficients(w,0,False)
        if not isinstance(base,Base):
            return None
        v1 = self.u1.components_en_base(base)
        v2 = self.u2.components_en_base(base)
        w = v1.cross(v2)
        return EquacioLineal.coeficients(w,0,False,prime)
    #
    #
    #
    @classmethod
    def from_matriu(cls,m):
        """
        Genera el pla vectorial generat per les columnes de la matriu "m"
        """
        if not isinstance(m,Matriu):
            return None
        if m.files != 3 or m.columnes != 2 or m.rank() != 2:
            return None
        v = m.vectors_columna()
        return cls(v[0],v[1])
    #
    #
    #
    @classmethod
    def amb_associat(cls,w):
        """
        Genera el pla vectorial que té vector perpendicular "w"
        """
        a = Matriu.from_vectors_fila([w])
        l = a.nucli()
        return cls(l[0],l[1])
    #
    #
    #
    def base_ortogonal(self):
        """
        Retorna una base orogonal del pla vectorial
        """
        v1 = self.u1
        u2 = self.u2
        v2 = v1.dot(v1) * u2 - v1.dot(u2) * v1
        v2.simplificar()
        return Base([v1,v2])

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
    @classmethod
    def amb_associat(cls,w,p):
        a = Matriu.from_vectors_fila([w])
        l = a.nucli()
        return cls(l[0],l[1],p)
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
                eq.append(latex(incg[i]-self.p[i]))
            else:
                eq.append(f"\\frac{{{latex(incg[i]-self.p[i])}}}{{{self.u[i]}}}")
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

class TransformacioLineal:
    #
    #
    #
    def __init__(self,matriu,base=None):
        if matriu.files != matriu.columnes:
            return None
        self.dimensio = matriu.files
        self.base = base
        self.matriu = matriu
        if self.base is None:
            self.canonica = Matriu(self.matriu.matriu)
        else:
            if not isinstance(base,Base):
                return None
            c = base.matriu()
            if c is None:
                return None
            if c.files != self.dimensio:
                return None
            if c.coulumnes != self.dimensio:
                return None
            m = c.matriu * self.matriu.matriu * c.matriu**(-1)
            self.canonica = m
    #
    #
    #
    def latex(self,base=None):
        if self.dimensio <= 3:
            x, y, z = symbols('x y z')
            u, v, w = symbols('u v w')
            o = [x,y,z]
            d = [u,v,w]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            u1, u2, u3, u4, u5, u6, u7, u8 = symbols('u1 u2 u3 u4 u5 u6 u7 u8')
            o = [x1, x2, x3, x4, x5, x6, x7, x8]
            d = [u1, u2, u3, u4, u5, u6, u7, u8]
        o = o[0:self.dimensio]
        d = d[0:self.dimensio]
        if base is None:
             o = " \\\\ ".join([latex(x) for x in o])
             d = " \\\\ ".join([latex(x) for x in d])
             m = self.canonica
        else:
            if not isinstance(base,Base):
                return None
            c = base.matriu()
            o = " \\\\ ".join([latex(x) + "'" for x in o])
            d = " \\\\ ".join([latex(x) + "'" for x in d])
            m = c.matriu**(-1) * self.canonica.matriu * c.matriu
            m = Matriu(m)
        s = "\\begin{pmatrix}{c} " + d + "\\end{pmatrix} = \n"
        s += f"{m}\n"
        s += "\\begin{pmatrix}{c} " + o + "\\end{pmatrix}"
        return s
    #
    #
    #
    def __repr__(self):
        return self.latex()
    #
    #
    #
    def transforma(self,vec,base=None):
        if not isinstance(vec,Vector):
            return None
        if vec.dimensio != self.dimensio:
            return None
        if base is None:
            return self.canonica * vec
        else:
            if not isinstance(base,Base):
                return None
            c = base.matriu()
            m = c.matriu**(-1) * self.canonica.matriu * c.matriu
            m = Matriu(m)
            return m * vec

class TransformacioAfi:
    #
    #
    #
    def __init__(self,terme,matriu):
        if not isinstance(terme,Vector):
            return None
        self.terme = terme
        if not isinstance(matriu,Matriu):
            return None
        if matriu.files != matriu.columnes:
            return None
        if terme.dimensio != matriu.columnes:
            return None
        self.matriu = matriu
        self.dimensio = terme.dimensio
    #
    #
    #
    def __repr__(self):
        if self.dimensio <= 3:
            x, y, z = symbols('x y z')
            u, v, w = symbols('u v w')
            o = [x,y,z]
            d = [u,v,w]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            u1, u2, u3, u4, u5, u6, u7, u8 = symbols('u1 u2 u3 u4 u5 u6 u7 u8')
            o = [x1, x2, x3, x4, x5, x6, x7, x8]
            d = [u1, u2, u3, u4, u5, u6, u7, u8]
        o = o[0:self.terme.dimensio]
        d = d[0:self.terme.dimensio]
        o = " \\\\ ".join([latex(x) for x in o])
        d = " \\\\ ".join([latex(x) for x in d])
        m = Matriu.matriu_columna(self.terme)
        s = "\\begin{pmatrix}{c} " + d + "\\end{pmatrix} = \n"
        s += f"{m} + \n"
        s += f"{self.matriu}\n"
        s += "\\begin{pmatrix}{c} " + o + "\\end{pmatrix}"
        return s
    #
    #
    #
    def transforma(self,p):
        if not isinstance(p,Vector):
            return None
        if self.dimensio != p.dimensio:
            return None
        return self.terme + self.matriu * p
