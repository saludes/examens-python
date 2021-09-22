#!/usr/bin/python3
# -*- coding: utf-8
"""
Filename:   Algebra.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2020--2021
Disclaimer: This program is provided "as is", without warranty of any kind,
            either expressed or implied, including, but not linmited to, the
            implied warranties of merchantability and fitness for a particular
            purpose.
            It has been written to generate random models of exams for the
            subject of Linear Algebra at ESEIAAT, Technic University of Catalonia
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
import re
from sympy import *
from sympy.solvers.solveset import linsolve
from sympy.vector import Vector
from functools import reduce
from math import gcd
from sympy import Basic, Function, Symbol
from sympy.printing.printer import Printer
from sympy.printing.latex import print_latex
from sympy.core.basic import Basic
from itertools import permutations

var('p q u v')
ddict = collections.OrderedDict([(p**2,1),(q**2,2),(u**2,3),(v**2,4),
                                ((p,q),5),((p,u),6),((p,v),7),
                                ((q,u),8),((q,v),9),((u,v),10),
                                ((q,p),5),((u,p),6),((v,p),7),
                                ((u,q),8),((v,q),9),((v,u),10),
                                (p,11),(q,12), (u,13), (v,14)])

class Impresora(Printer):
    """
    La funció latex() del sympy té la mania d'escriure les variables x, y, z i t
    en l'ordre t, x, y i z. L'única manera que, de moment, he trobat per resoldre
    aquest inconvenient és definir la classe Impresora i la funció mylatex().
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
            p, q, r, s = symbols('p q u v')
            if isinstance(el,Pow):
                return new_place(el)
            if el.is_integer:
                return new_place(el)
            elif el.is_symbol:
                return new_place(el)
            if len(el.args) == 2:
                if el.args[0].is_symbol and el.args[1].is_symbol:
                    k = -1
                    for i in range(0,2):
                        if el.args[i] not in [p,q,r,s]:
                            k = i
                    if k >= 0:
                        k = (k+1) % 2
                        return new_place(el.args[k])
                    return new_place(el.args)
                q = el.args[len(el.args)-1]
                if q.is_symbol or q.is_integer or q.is_rational:
                    return new_place(q)
                elif q.args[0].is_symbol:
                    return new_place(q)
                else:
                    return 0
            elif len(el.args) == 3:
                k = -1
                for i in range(0,2):
                    if el.args[i+1] not in [p,q,r,s]:
                        k = i
                if k >= 0:
                    k = (k+1) % 2
                    return new_place(el.args[k+1])
                return new_place(el.args[1:])
            else:
                return 0

        def write_coeff(el):
            if isinstance(el,Pow):
                return " + %s" % latex(el)
            if el.is_integer:
                if el > 0:
                    return " + %s" % el
                else:
                    return " %s" % el
            elif el.is_symbol:
                return " + %s" % el
            elif len(el.args) == 2 and el.args[0].is_symbol and el.args[1].is_symbol:
                return " + %s" % latex(el)
            elif len(el.args) > 0:
                if el.args[len(el.args)-1].is_symbol:
                    if el.args[0].is_rational or el.args[0].is_integer:
                        if el.args[0] > 0:
                            return " + %s" % latex(el)
                        else:
                            return " %s" % latex(el)
                    else:
                        return " %s" % latex(el)
                else:
                    if el.args[0].is_rational or el.args[0].is_integer:
                        if el.args[0] > 0:
                            return " + %s" % latex(el)
                    return " %s" % latex(el)
            else:
                return " %s" % el
        list_place = [get_place(a) for a in expr.args]
        expr_args = list(zip(*sorted(zip(list_place,expr_args))))[1]
        to_print = [write_coeff(a) for a in expr_args]
        to_print[0] = str(latex(expr_args[0]))
        return "".join(a for a in to_print)

def mylatex(e):
    p, q, r, s = symbols('p q u v')
    e = str(e)
    e = e.replace('x','p')
    e = e.replace('y','q')
    e = e.replace('z','u')
    e = e.replace('t','v')
    e = sympify(e)
    if not isinstance(e,Add):
        e = latex(e)
    else:
        e = Impresora().doprint(e)
    e = e.replace('p','x')
    e = e.replace('q','y')
    e = e.replace('u','z')
    e = e.replace('v','t')
    return e

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

def matriu_latex(m,format=None,ampliada=False,tipus="p"):
    """
    Retorna l'expressió en latex d'una matriu del tipus Matrix del sympy
    Parametres:
        format:   format de les columnes de la matriu. Per defecte "r"
        ampliada: si es vol separar amb una línia vertical l'última
                  columna de la matriu
    """
    f, c = m.shape
    vert = ""
    if ampliada:
        cols = c - 1
        vert = "|r"
    else:
        cols = c
    if format is None:
        text = "\\begin{TIPUSmatrix}{*{%d}r%s} LINES\\end{TIPUSmatrix}" % (cols,vert)
    else:
        text = "\\begin{TIPUSmatrix}{%s} LINES\\end{TIPUSmatrix}" % format
    text = text.replace('TIPUS',tipus)
    lines = []
    for i in range(f):
        line = []
        for j in range(c):
            line.append(latex(m[i,j]))
        lines.append(" & ".join(map(str,line)))
    return (text.replace('LINES',"\\\\ ".join(lines)))

def matriu_mathematica(m):
    """
    Retorna l'expressió en Mathematica d'una matriu del tipus Matrix del sympy
    """
    r = re.compile(r'sqrt\((\d+)\)')
    f, c = m.shape
    lines = []
    slines = "{LINES}"
    for i in range(f):
        line = []
        sline = "{LINE}"
        for j in range(c):
            line.append(m[i,j])
        sline = sline.replace('LINE',",".join(map(str,line)))
        sline = r.sub(r'Sqrt[\g<1>]',sline)
        lines.append(sline)
    return (slines.replace('LINES',",".join(lines)))

def primer_no_nul(list):
    """
    Retorna l'índex del primer coeficient no nul d'una llista
    """
    if list is None or len(list) == 0:
        return None
    for i in range(len(list)):
        if list[i] != 0:
            return i
    return None

def vectors_latex(l):
    """
    Retorna la llista de vectors l escrita en latex
    Paràmetres:
        l: llista de vectors o punts
    """
    return ",".join([str(k) for k in l])

def vaps_veps_amb_signe(result,signe=1):
    """
    Donada una matriu A del sympy, result ha de ser el resultat de la funció
    r = A.eigenvects(). Aleshores aquesta funció retorna la llista valors propis
    positius (signe > 0) o negatius (signe < 0) i els seus vectors propis.
    Paràmetres:
        result: resultat de la funció eigenvects() del sympy
        signe: positiu o negatiu, en funció de quins valors i vectors propis es volen
    """
    vaps = []
    veps = []
    for l, m, us in result:
        if l * signe <= 0:
            continue
        for k in range(m):
            vaps.append(l)
            n1 = Matriu(us[k])
            vep = n1.vectors_columna()[0]
            vep.simplificar()
            veps.append(vep)
    return (vaps,veps)

def vaps_veps(result):
    """
    Retorna la llista valors propis i els seus vectors propis.
    Paràmetres:
        result: resultat de la funció eigenvects() del sympy
    """
    vaps = []
    veps = []
    for l, m, us in result:
        for k in range(m):
            vaps.append(l)
            n1 = Matriu(us[k])
            vep = n1.vectors_columna()[0]
            vep.simplificar()
            veps.append(vep)
    return (vaps,veps)

class Radicals(object):
    """
    Classe per treure factor comú en expressions on hi apareixen arrels quadrades
    """
    def __init__(self):
        self.quadrats = []
        self.enters = []
        self.fraccions = []
    #
    #
    #
    def busca_quadrats(self,el):
        """
        Afegeix els termes que apareixen dins d'arrels quadrades a la llista
        self.quadrats
        Paràmetres:
            el: expressió del sympy
        """
        if isinstance(el,Pow) and el.args[1] == Rational(1,2):
            if el.args[0] not in self.quadrats:
                self.quadrats.append(el.args[0])
        elif isinstance(el,int) or isinstance(el,Integer):
            if el not in self.enters:
                self.enters.append(el)
        else:
            for k in el.args:
                self.busca_quadrats(k)
    #
    #
    #
    def busca_fraccions(self,el):
        """
        Afegeix els termes que apareixen als denominadors a la llista
        self.fraccions
        Paràmetres:
            el: expressió del sympy
        """
        if isinstance(el,Pow):
            return
        if isinstance(el,int) or isinstance(el,Integer):
            return
        if isinstance(el,Rational):
            if el.q not in self.fraccions:
                self.fraccions.append(el.q)
        for k in el.args:
            self.busca_fraccions(k)
    #
    #
    #
    def mcd(self):
        """
        Retorna el màxim comú divisor dels elements de la llista self.quadrats
        """
        if len(self.quadrats) == 0:
            a = 1
        else:
            a = mcd_llista(self.quadrats)
        if len(self.enters) == 0:
            b = 1
        else:
            b = mcd_llista(self.enters)
        return (a,b)
    #
    #
    #
    def mcm(self):
        """
        Retorna el mínim comú múltiple dels elements de la llista self.faccions
        """
        if len(self.fraccions) == 0:
            return 1
        return mcm_llista(self.fraccions)

class Vector(object):
    """
    Classe que ens permetrà representar vectors i punts
    Atributs:
        dimensio: el nombre de components o longitud del vector
        components: llista amb les components del vector
    """
    #
    #
    #
    def __new__(cls,*args):
        if len(args) == 0:
            return None
        if len(args) == 1:
            if isinstance(args[0],list) or isinstance(args[0],tuple):
                return super(Vector,cls).__new__(cls)
        else:
            if isinstance(args[0],list) or isinstance(args[0],tuple):
                return None
            return super(Vector,cls).__new__(cls)
        return None
    #
    #
    #
    def __init__(self,*args):
        """
        Constructor.
        Paràmetres:
           c: Una única llista de nombres o una llista de paràmetres
              que han de ser nombres
        Exemples:
           u = Vector([2,3,1,2])
           v = Vector(3,1,-2)
        """
        if len(args) == 1:
            c = args[0]
        else:
            c = list(args)
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
            v = [-k for k in v]
        return (f,Vector(v))
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex del vector
        """
        l = []
        s = []
        square = False
        other = False
        for k in self.components:
            if isinstance(k,Rational):
                l.append(k.q)
                s.append(False)
            elif isinstance(k,int) or isinstance(k,Integer):
                l.append(1)
                s.append(False)
            elif isinstance(k**2,Rational):
                square = True
                k2 = k**2
                l.append(k2.q)
                s.append(True)
            else:
                other = True
        if other:
            s = ",".join([latex(k) for k in self.components])
            return f"\\left({s}\\right)"
        if square:
            for k in range(len(l)):
                if not s[k]:
                    l[k] = l[k]**2
        m = mcm_llista(l)
        if square:
            m = sqrt(m)
        l = [m*k for k in self.components]
        s = ""
        if m != 1:
            s = f"\\deufrac{{1}}{{{latex(m)}}}"
        r = ",".join([latex(k) for k in l])
        return f"{s}({r})"
    #
    #
    #
    def __add__(self,other):
        """
        Defineix la suma de vectors.
        Paràmetres:
            other: un altre vector
        Exemple:
            u1 = Vector(3,2,1,3)
            u2 = Vector(-2,4,-3,1)
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
        Defineix la resta de vectors.
        Paràmetres:
            other: un altre vector
        Exemple:
            u1 = Vector(3,2,1,3)
            u2 = Vector(-2,4,-3,1)
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
    def __neg__(self):
        """
        Retorna l'oposat d'un vector.
        Exemple:
            u1 = Vector(3,2,1,3)
            u2 = - u1
        """
        p = [-k for k in self.components]
        return Vector(p)
    #
    #
    #
    def __mul__(self,other):
        """
        Defineix el producte d'un escalar per un vector i el
        producte d'una vector per una matriu.
        Paràmetres:
            other: un escalar o una matriu (classe Matriu)
        Exemples:
            u1 = Vector(3,2,1,3)
            u2 = Vector(-2,4,-3,1)
            v = 5 * u1 - 4 * u2
            a = Matriu.aleatoria(f=4,c=2)
            u = v * a
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
        Permet indexar els elements d'un vector.
        Paràmetres:
            i : índex
        Exemple:
            v = Vector(3,2,1,3,5)
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
        Permet assignar valors mitjançant índexs.
        Paràmetres:
            i     : índex
            valor : valor vàlid del sympy
        Exemple:
            v = Vector.aleatori(l=4)
            v[0] = 3
        """
        try:
            self.components[i] = value
        except:
            pass
    #
    #
    #
    def __eq__(self,other):
        """
        Determina si dos vectors són iguals
        Paràmetres:
            other: un altre vector
        Exemple:
            u = Vector(1,2,3)
            v = Vector(1,2,3)
            if u == v:
                pass
        """
        if self.dimensio != other.dimensio:
            return False
        for k in range(self.dimensio):
            if self[k] != other[k]:
                return False
        return True
    #
    #
    #
    def dot(self,other):
        """
        Retorna el producte escalar de dos vectors.
        Paràmetres:
            other: un altre vector
        Exemple:
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
    def normalitzar(self):
        """
        Converteix el vector en unitari
        Exemple:
            u = Vector(1,2,3)
            u.normalitzar()
        """
        l = self.length()
        self.components = [k / l for k in self.components]
    #
    #
    #
    def normalitzat(self):
        """
        Retorna el vector unitari en la direcció i sentit del vector
        Exemple:
            u = Vector(1,2,3)
            v = u.normalitzat()
        """
        l = self.length()
        return Vector([k / l for k in self.components])
    #
    #
    #
    def maxim(self):
        """
        Retorna el màxim dels valors absoluts de les seves components
        Exemple:
            u = Vector.aleatori(l=5,maxim=4,nuls=False)
            m = u.maxim()
        """
        m = 0
        for k in self.components:
            v = abs(k)
            if v > m:
                m = v
        return m
    #
    #
    #
    def simplificar(self):
        """
        Simplifica el vector, és a dir, converteix les seves components en una
        llista d'enters amb mcd igual a 1.
        Només té sentit si totes les components del vector són nombres enters
        o racionals
        Exemple:
            u = Vector(-1,2,2,Rational(3,2))
            u.simplificar()
        """
        if self.length() == 0:
            return
        d = []
        for i in range(self.dimensio):
            k = self.components[i]
            if isinstance(k,Rational):
                d.append(k.q)
            elif isinstance(k,int):
                pass
            elif isinstance(k,Integer):
                pass
            else:
                return
        mcm = mcm_llista(d)
        v = [mcm * x for x in self.components]
        mcd = mcd_llista(v)
        v = [simplify(k // mcd) for k in v]
        if v[0] < 0:
            v = [-k for k in v]
        self.components = list(v)
    #
    #
    #
    def radsimplificar(self):
        """
        Simplifica el vector quan alguna de les seves components té radicals
        Exemple:
            u = Vector(sqrt(2),-3*sqrt(2),-sqrt(2))
            u.radsimplificar()
        """
        r = Radicals()
        for k in self.components:
            r.busca_quadrats(k)
        a, b = r.mcd()
        if abs(b) != 1:
            for k in range(self.dimensio):
                self.components[k] /= b
        if len(r.quadrats) <=  1:
            for k in range(self.dimensio):
                n = sqrt(a) * self.components[k]
                self.components[k] = simplify(n.expand())
            self.simplificar()
            return
        for k in range(self.dimensio):
            n = sqrt(a) * self.components[k] / a
            self.components[k] = simplify(n.expand())
        r = Radicals()
        for k in self.components:
            r.busca_fraccions(k)
        m = r.mcm()
        if m != 1:
            for k in range(self.dimensio):
                self.components[k] *= m
        self.simplificar()
    #
    #
    #
    def cross(self,other,simplificar=False):
        """
        Retorna un nou vector que és el producte vectorial de dos vectors:
        Si simplificar és True, simplifica el vector resultant
        Paràmetres:
            other: un altre vector
            simplificar: si es vol simplificar el resultat o no
        Exemple:
            u = Vector.aleatori(l=3)
            v = Vector.aleatori(l=3)
            w = u.cross(v,simplificar=True)
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
        Retorna l'expressió en latex del vector.
        Paràmetres:
            unitari: si unitari és True, el retorna dividit per la seva longitud
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
        Paràmetres:
            base: una base (classe Base)
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
    #
    #
    #
    def punt(self):
        return Punt(self.components)

class Punt(Vector):
    """
    Classe per treballar amb punts.
    Internament un punt és el mateix que un vector
    """
    def __new__(cls,*args):
        """
        Constructor.
        Paràmetres:
           c: Llista de nombres
        """
        if len(args) == 0:
            return None
        if len(args) == 1:
            if isinstance(args[0],list) or isinstance(args[0],tuple):
                return super(Vector,cls).__new__(cls)
        else:
            for k in args:
                if isinstance(k,list) or isinstance(k,tuple):
                    return None
            return super(Vector,cls).__new__(cls)
        return None
    #
    #
    #
    def __init__(self,*args):
        Vector.__init__(self,*args)
    #
    #
    #
    def coordenades_en_referencia(self,ref):
        """
        Retorna les coordenades del punt en la referència "ref"
        Paràmetres:
            ref: referència de la classe ReferenciaAfi
        """
        if not isinstance(ref,ReferenciaAfi):
            return None
        if self.dimensio != ref.dimensio:
            return None
        op = self - ref.origen
        return Punt(op.components_en_base(ref.base).components)

class Base(object):
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
    def __new__(cls,vecs,unitaria=False):
        """
        Contructor.
        Paràmetres:
           vecs: llista de vectors
           unitaria: True o False.
                     Serveix per si volem imprimir o fer servir els seus
                     vectors com a unitaris
        Exemple:
            base = Base([Vector(1,1),Vector(-1,1)],unitaria=True)
            print(base)
            print(base.vecs)
        """
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
        m = Matriu.from_vectors_columna(vecs)
        if m.files != d:
            return None
        if m.rank() != len(vecs):
            return None
        return super(Base,cls).__new__(cls)
    #
    #
    #
    def __init__(self,vecs,unitaria=False):
        self.unitaria = unitaria
        self.vecs = vecs
        self.dimensio = vecs[0].dimensio
    #
    #
    #
    def es_unitaria(self):
        """
        Retorna si la base és unitària. Notem que els vectors no es guarden
        com a unitaris.
        """
        return self.unitaria
    #
    #
    #
    def es_ortogonal(self):
        """
        Retorna si la base és ortogonal
        """
        for i in range(self.dimensio):
            for j in range(i+1,self.dimensio):
                if self.vecs[i].dot(self.vecs[j]) != 0:
                    return False
        return True
    #
    #
    #
    def set_unitaria(self):
        """
        Si la matriu és ortogonal, la passa a unitària
        """
        if not self.es_ortogonal():
            return
        self.unitaria = True
    #
    #
    #
    def te_orientacio_positiva(self):
        """
        Retorna si té orientació positiva
        """
        m = Matriu.from_vectors_columna(self.vecs)
        return m.determinant() > 0
    #
    #
    #
    def orientacio_positiva(self):
        """
        Fa que tingui orientació positiva canviant, si cal, de signe l'últim
        vector
        """
        if not self.te_orientacio_positiva():
            self.vecs[-1] = - self.vecs[-1]
    #
    #
    #
    def matriu(self):
        """
        Retorna la matriu de la classe Matriu que té per columnes
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
    #
    #
    #
    def vectors(self,unitaris=False):
        """
        Retorna els vectors la base
        Paràmetres:
            unitaris: si és True, els divideix per la seva longitud
        """
        if not unitaris:
            return self.vecs
        return [v/v.length() for v in self.vecs]
    #
    #
    #
    def vectors_latex(self):
        """
        Retorna l'expressió en latex dels vector de la base, sense les claus
        inicial i final
        """
        if not self.unitaria:
            base = ",".join(map(str,self.vecs))
        else:
            base = ",".join([v.latex(unitari=True) for v in self.vecs])
        return f"{base}"
    #
    #
    #
    def vector_de_components(self,vec):
        """
        Retorna un nou vector expressat en la base canònica del vector
        que en aquesta base té components "vec"
        Paràmetres:
            vec: vector
        """
        if not isinstance(vec,Vector):
            return None
        if vec.dimensio != self.dimensio:
            return None
        if not self.unitaria:
            c = self.matriu()
        else:
            unitaris = [(1 / v.length()) * v for v in self.vecs]
            c = Matriu.from_vectors_columna(unitaris)
        return c * vec
    #
    #
    #
    @classmethod
    def ortogonal(cls,ordre=3,maxim=5,unitaria=False):
        """
        Retorna una base ortogonal "aleatòria"
        Paràmetres:
            ordre: dimensió
            maxim: màxim per a les components dels vectors de la base
            unitaria: si és True, la base serà ortonormal
        """
        trobat = False
        while not trobat:
            m = Matriu.invertible(ordre=ordre,maxim=4,mzeros=0,unitaria=True)
            L = []
            for v in m.vectors_columna():
                a = Matriu.from_vectors_columna([v])
                L.append(a.matriu)
            Q = GramSchmidt(L)
            base = []
            for m in Q:
                m = Matriu(m)
                v = m.vectors_columna()[0]
                v.simplificar()
                base.append(v)
            m = 0
            for v in base:
                a = v.maxim()
                if a > m:
                    m = a
            trobat = m <= maxim
        return cls(base,unitaria)
    #
    #
    #
    def quadrats_longituds(self):
        """
        Retorna els quadrats de les longituds dels vectors de la base sense tenir
        en compte si la base és unitària
        """
        return [v.length()**2 for v in self.vecs]

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

        Només s'utilitzen quan generem una matriu diagonalitzble
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
        Paràmetres:
            vaps: llista de nombres
        """
        self.vaps = vaps
    #
    #
    #
    def set_veps(self,veps):
        """
        Assigna un llista de vectors propis simplificats a la variable self.veps
        Paràmetres:
            vaps: llista de vectors
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
    @classmethod
    def transformacio_elemental(cls,ordre,i,j,s,t):
        a = eye(ordre)
        a[i,:] = s * a[i,:] + t * a[j,:]
        return cls(a)
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
    def max_diagonal(self):
        """
        Retorna el màxim en valor absolut dels coeficients de la diagonal
        Si la matriu no és quadrada retorna None
        """
        if self.files == 0 or self.columnes == 0:
            return 0
        max = 0
        m = min(self.files,self.columnes)
        for i in range(m):
            k = abs(self.matriu[i,i])
            if k > max:
                max = k
        return max
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
    def rang(self):
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
    def det(self):
        """
        Retorna el determiant de la matriu
        """
        return self.matriu.det()
    #
    #
    #
    def clona(self):
        return Matriu(self.matriu[:,:])
    #
    #
    #
    def latex(self,format=None,tipus='p'):
        return matriu_latex(self.matriu,format,tipus)
    #
    #
    #
    def polinomi_caracteristic(self):
        """
        Retorna el polinomi característic de la matriu
        """
        lamda = symbols('lamda')
        p = (-1)**(self.files) * self.matriu.charpoly(lamda)
        return latex(p.as_expr())
    #
    #
    #
    def __getitem__(self,tup):
        """
        Permet indexar els elements d'una matriu.
        Paràmetres:
            tup: tupla d'índexs
        Exemple:
        m = Matriu.aleatoria()
        k = m[2,1]
        """
        i, j = tup
        try:
            return self.matriu[i,j]
        except:
            return None
    #
    #
    #
    def __setitem__(self,tup,value):
        """
        Permet assignar valors mitjançant índexs.
        Paràmetres:
            tup: tupla d'índexs
            valor: valor assignat
        Exemple:
            m = Matriu.aleatoria()
            m[2,1] = 5
        """
        i, j = tup
        self.matriu[i,j] = value
    #
    #
    #
    def __add__(self,other):
        """
        Defineix la suma de matrius.
        Paràmetres:
            other: una altra matriu
        Exemple:
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
        Defineix la resta de matrius.
        Paràmetres:
            other: una altra matriu
        Exemple:
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
    def __neg__(self):
        """
        Retorna l'oposada de la matriu
        """
        return Matriu(-self.matriu)
    #
    #
    #
    def __rmul__(self,other):
        """
        Defineix el producte d'un escalar per una matriu
        Paràmetrers:
            other: un escalar
        Exemple:
            a = Matriu.aleatoria(f=4,c=3,maxim=7,nuls=False)
            c = 4 * a
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
        Defineix el producte de matrius i el producte d'una matriu per un vector
        Paràmetrers:
            other: un vector (classe Vector) o un matriu (classe Matriu)
        Exemple:
            a1 = Matriu.aleatoria()
            a2 = Matriu.aleatoria()
            b = a1 * a2
            v = Vector(2,1,-1)
            w = a1 * v
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
    def transposada(self):
        """
        Retorna la transposada de la matriu
        """
        return Matriu(self.matriu.T)
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex de la matriu
        """
        l = []
        s = []
        square = False
        for i in range(self.files):
            for j in range(self.columnes):
                k = self.matriu[i,j]
                if isinstance(k,Rational):
                    l.append(k.q)
                    s.append(False)
                elif isinstance(k,int) or isinstance(k,Integer):
                    l.append(1)
                    s.append(False)
                elif isinstance(k**2,Rational):
                    square = True
                    k2 = k**2
                    l.append(k2.q)
                    s.append(True)
                else:
                    return matriu_latex(self.matriu)
        if square:
            for k in range(len(l)):
                if not s[k]:
                    l[k] = l[k]**2
        m = mcm_llista(l)
        if square:
            m = sqrt(m)
        l = []
        for i in range(self.files):
            for j in range(self.columnes):
                l.append(m * self.matriu[i,j])
        s = ""
        if m != 1:
            s = f"\\deufrac{{1}}{{{latex(m)}}}"
        m = Matrix(self.files,self.columnes,l)
        return s + matriu_latex(m)
    #
    #
    #
    def __eq__(self,other):
        """
        Compara si dues matriu són iguals
        Paràmetres:
            other: una altra matriu
        Exemple:
            a1 = Matriu.aleatoria()
            a2 = Matriu.aleatoria()
            if a1 == a2:
                pass
        """
        if self.files != other.files:
            return False
        if self.columnes != other.columnes:
            return False
        for i in range(self.files):
            for j in range(self.columnes):
                if self.matriu[i,j] != other.matriu[i,j]:
                    return False
        return True
    #
    #
    #
    @classmethod
    def diagonal(cls,vals):
        """
        Retorna una matriu diagonal amb valors "vals" a la diagonal
        Paràmetres:
            vals: llista d'escalars o vector (class Vector o Punt)
        """
        if isinstance(vals,Vector):
            vals = vals.components
        if not isinstance(vals,list) and not isinstance(vals,tuple):
            return None
        d = len(vals)
        if d == 0:
            return None
        m = diag(*vals)
        return cls(m)
    #
    #
    #
    @classmethod
    def amb_rang(cls,f=3,c=3,r=3,maxim=5,mzeros=-1):
        """
        Retorna una matriu aleatoria amb rang r.
        Paràmetres:
          f: nombre de files de la matriu
          c: nombre de columnes de la matriu
          r: rang de la matriu
          maxim: tots els elements tindran valor absolut menor que "maxim"
          nuls: la matriu pot contenir coeficients nuls o no
        """
        values = [-1,1,2]
        trobat = False
        while not trobat:
            mc = cls.invertible(f,maxim=2,mzeros=2,unitaria=True)
            md = cls.invertible(c,maxim=2,mzeros=2,unitaria=True)
            m = Matrix.zeros(f,c)
            for k in range(r):
                m[k,k] = values[random.randint(0,2)]
            n = mc.matriu * m * md.matriu
            if mzeros >= 0 and nzeros(n) > mzeros:
                continue
            trobat = norma_maxim(n) <= maxim
        return cls(n)
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
    def adjunta(self):
        """
        Retorna una nova matriu que és l'adjunta de l'actual
        """
        return Matriu(self.matriu.adjugate().T)
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
    def diagonalitzable(cls,ordre=3,maxim=5,mzeros=-1,mvaps=3,vapsnuls=False,vapsrepetits=True):
        """
        Retorna una matriu quadrada aleatoria diagonalitzable.
        Paràmetres:
            ordre: nombre de files i columnes de la matriu
            maxim: tots els elements tindran valor absolut menor o igual que "maxim"
            mzeros: si nzeros >= 0, nombre màxim de zeros que tindrà la matrius
                    si nzeros < 0, el nombre de zeros no està limitat
            mvaps: tots els valors propis tindran valor absolut menor o igual que "mvaps"
            vapsnuls: si hi pot aparèixer el valor propi nul
            vapsrepetits: si hi pot aparèixer valors propis repetits
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
            if not vapsrepetits and len(set(vaps)) != ordre:
                continue
            vaps.sort()
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
        Retorna una matriu quadrada aleatoria que serà d'un producte escalar,
        és a dir, una matriu de Gram
        Paràmetres:
            ordre: nombre de files i columnes de la matriu
            maxim: tots els elements tindran valor absolut menor o igual que "maxim"
            mzeros: si nzeros >= 0, nombre màxim de zeros que tindrà la matrius
                    si nzeros < 0, el nombre de zeros no està limitat
        """
        trobat = False
        while not trobat:
            c = Matriu.invertible(ordre,maxim=5)
            g = c.transposada() * c
            if mzeros >= 0 and g.nzeros() > mzeros:
                continue
            if g.norma_maxim() > maxim:
                continue
            trobat = True
        return cls(g.matriu)
    #
    #
    #
    @classmethod
    def matriu_fila(cls,v):
        """
        Retorna una nova matriu fila a partir de les components del vector v
        Paràmetres:
            v: vector o punt
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
        Paràmetres:
            v: vector o punt
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
        Paràmetres:
            v: llista de vectors o punts
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
        Paràmetres:
            v: llista de vectors o punts
        """
        if not (isinstance(vecs,list) or isinstance(vecs,tuple)):
            return None
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
    @classmethod
    def identitat(cls,ordre):
        return cls(eye(ordre))
    #
    #
    #
    def vectors_columna(self,simplificar=False):
        """
        Retorna una llista amb els vectors columna de la matriu
        Paràmetres:
            simplificar: si és True retornarà els vectors simplificats
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
        Paràmetres:
            simplificar: si és True retornarà els vectors simplificats
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
        for v in vecs:
            v.simplificar()
        return vecs
    #
    #
    #
    def es_simetrica(self):
        """
        Retorna True si és simètrica
        """
        return self == self.transposada()
    #
    #
    #
    def es_diagonal(self):
        """
        Retorna True si és una matriu diagonal
        """
        for i in range(self.files):
            for j in range(self.columnes):
                if i != j and self.matriu[i,j] != 0:
                    return False
        return True
    #
    #
    #
    def intercanvia_columnes(self,i,j):
        """
        Retorna una matriu amb les columnes i i j permutades
        Paràmetres:
            i, j: índexs de les columnes
        """
        if i >= self.columnes or j >= self.columnes:
            return None
        if i == j:
            return Matriu(self.matriu)
        c = self.vectors_columna()
        k = c[i]
        c[i] = c[j]
        c[j] = k
        return Matriu.from_vectors_columna(c)
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
        Paràmetrers:
            fila: nova fila de la matriu
            pos: posició que ha d'ocupar la nova fila
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
        Paràmetrers:
            columna: nova columna de la matriu
            pos: posició que ha d'ocupar la nova columna
        """
        if not isinstance(columna,Vector):
            return None
        if columna.dimensio != self.files:
            return None
        if pos > self.columnes:
            pos = self.columnes
        columnes = self.vectors_columna()
        columnes.insert(pos,columna)
        return Matriu.from_vectors_columna(columnes)
    #
    #
    #
    def factor_comu(self):
        """
        Retorna quin factor comú podem treure de la matriu
        """
        d = []
        for i in range(self.files):
            for j in range(self.columnes):
                k = self[i,j]
            if (isinstance(k,int) or isinstance(k,Integer)):
                d.append(k)
            else:
                return 1
        return (mcd_llista(d))
    #
    #
    #
    def simplificar(self):
        """
        Simplifica la matriu, és a dir, converteix les seves entrades en una
        llista d'enters amb mcd igual a 1.
        Només té sentit si totes les components del vector són nombres enters
        o racionals
        """
        d = []
        for i in range(self.columnes):
            for j in range(self.files):
                if isinstance(self.matriu[i,j],Rational):
                    d.append(self.matriu[i,j].q)
                elif isinstance(self.matriu[i,j],int):
                    pass
                elif isinstance(self.matriu[i,j],Integer):
                    pass
                else:
                    return
        mcm = mcm_llista(d)
        m = mcm * self.matriu
        d = []
        for i in range(self.columnes):
            for j in range(self.files):
                d.append(m[i,j])
        mcd = mcd_llista(d)
        for i in range(self.columnes):
            for j in range(self.files):
                self.matriu[i,j] = m[i,j] // mcd
        if self.matriu[0,0] < 0:
            self.matriu = - self.matriu
    #
    #
    #
    def submatriu(self,files,columnes):
        """
        Retorna la submatriu determinada per les files "files" i les
        columnes "columnes".
        Paràmetrers:
            files: llista de files
            columnes: llista de columnes
        """
        if not (isinstance(files,list) or isinstance(files,tuple)):
            return None
        if not (isinstance(columnes,list) or isinstance(columnes,tuple)):
            return None
        if max(columnes) >= self.columnes or min(columnes) < 0:
            return None
        if max(files) >= self.files or min(files) < 0:
            return None
        m = self.matriu[files,columnes]
        return Matriu(m)
    #
    #
    #
    def sistema_propi(self):
        """
        Retorna el sistema d'equacions en format latex corresponent al
        càlcul dels valors propis de la matriu
        """
        if self.columnes <= 4:
            x, y, z, t = symbols('x y z t')
            unknowns = [x,y,z,t]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
        unknowns = unknowns[0:self.columnes]
        A = self * Vector(unknowns)
        eqs = [f"{latex(A[i])} &= \\lambda {unknowns[i]}" for i in range(A.dimensio)]
        eqs = " \\\\ ".join(eqs)
        return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"

class EquacioLineal:
    """
    Classe per treballar amb equacions lineals.
    Atributs:
        equacio: terme de l'esquerra en la equacio "eq = 0"
        unknowns: incògnites que apareixen a l'equació
        amp: True o False
        prime: nombre de primes que escriurem a l'equació
    """
    #
    #
    #
    def __init__(self,eq,amp=False,prime=0):
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
        """
        x, y, z, t = symbols('x y z t')
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        incognites = [x,y,z,t,x1,x2,x3,x4,x5,x6,x7,x8]
        self.equacio = eq
        self.amp = amp
        d = self.equacio.as_coefficients_dict()
        self.unknowns = []
        for k in incognites:
            if diff(eq,k) != 0:
                self.unknowns.append(k)
        self.prime = prime
    #
    #
    #
    @classmethod
    def coeficients(cls,a,b,amp=False,prime=0):
        """
        Retorna una nova equació amb coeficients de les incògnites el vector "a" i
        terme independent b
        Paràmetres:
            a: Vector amb els coeficients de les incògnites
            b: terme independents
            amp: si és True l'equació s'escriurà amb el &= per al LaTeX
            prime: nombre de primes que s'escriuran a les incògnites
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
        eq = simplify((eq - b).expand())
        return cls(eq,amp,prime)
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
        t = symbols('t')
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        incgs = {x1, x2, x3, x4, x5, x6, x7, x8}
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
            terme = 0
            for k in d.keys():
                if k.is_symbol and k not in self.unknowns:
                    terme += factor * d[k] * k
                else:
                    d[k] = factor * d[k]
                    eq += d[k] * k
        else:
            eq = self.equacio
        eq -= d[1]
        if t in self.unknowns:
            str = mylatex(eq)
        else:
            if len(incgs & set(self.unknowns)) > 0:
                str = latex(simplify(eq.expand()))
            else:
                str = mylatex(simplify(eq.expand()))
        if self.amp:
            str = f"{str} &= {latex(-d[1] - terme)}"
        else:
            str = f"{str} = {latex(-d[1] - terme)}"
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
        Suma d'equacions
        Paràmetres:
            other: una altra equació
        """
        if not isinstance(other,EquacioLineal):
            return None
        if self.prime != other.prime:
            return None
        eq = self.equacio + other.equacio
        return EquacioLineal(eq,self.amp or other.amp)
    #
    #
    #
    def __sub__(self,other):
        """
        Resta d'equacions
        Paràmetres:
            other: una altra equació
        """
        if not isinstance(other,EquacioLineal):
            return None
        if self.prime != other.prime:
            return None
        eq = self.equacio - other.equacio
        return EquacioLineal(eq,self.amp or other.amp)
    #
    #
    #
    def __mul__(self,other):
        """
        Producte d'un escalar per una equació
        Paràmetres:
            other: un escalar
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
    #
    #
    #
    def to_sistema_equacions(self):
        x, y, z, t = symbols('x y z t')
        unknowns = [x,y,z,t]
        n = -1
        for u in self.unknowns:
            m = unknowns.index(u)
            if m > n:
                n = m
        if n >= 0:
            return SistemaEquacions.from_equacions([self.equacio],n+1,self.prime)
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
        for u in self.unknowns:
            m = unknowns.index(u)
            if m > n:
                n = m
        return SistemaEquacions.from_equacions([self.equacio],n+1,self.prime)


class SistemaEquacions:
    """
    Classe per treballar amb sistemes d'equacions lineals
    Atributs:
      A: matriu dels coeficients de les incógnites
      B: vector de termes independents
      equacions: llista de EquacioLineal
      nombre: nombre d'equacions
      solucio: solucio del sistema d'equacions
      unknowns: llista d'incògnites
      prime: nombre de primes que escriurem a l'equació
    """
    #
    #
    #
    def __new__(cls,a,b,unknowns=None,prime=0):
        """
        Constructor.
        Paràmetres:
          a: matriu dels coeficients de les incògnites
          b: Termes independents
        """
        if not isinstance(a,Matriu):
            return None
        if not isinstance(b,Vector):
            return None
        if a.files != b.dimensio:
            return None
        if unknowns is not None:
            if not (insinatance(unknowns,list) or insinatance(unknowns,tuple)):
                return None
            if len(unknowns) != a.columnes:
                return None
        return super(SistemaEquacions,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a,b,unknowns=None,prime=0):
        self.A = a
        self.B = b
        self.solucio = None
        eq = []
        files = a.vectors_fila()
        for k in range(self.A.files):
            eq.append(EquacioLineal.coeficients(files[k],b[k],amp=True,prime=0))
        self.equacions = eq
        self.nombre = len(eq)
        if unknowns is not None:
            self.unknowns = unknowns
        else:
            if self.A.columnes <= 4:
                x, y, z, t = symbols('x y z t')
                unknowns = [x,y,z,t]
            else:
                x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
                unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
            self.unknowns = unknowns[0:self.A.columnes]
        self.prime = prime
    #
    #
    #
    @classmethod
    def from_equacions(cls,eqs,nombre,prime=0):
        """
        Retorna un sistema d'equacions amb equacions "eqs"
        Paràmetres:
            eqs: llista de EquacioLineal
            nombre: nombre d'incògnites
            prime: nombre de primes que s'escriuran a les incògnites
        """
        if nombre <= 4:
            x, y, z, t = symbols('x y z t')
            unknowns = [x,y,z,t]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            unknowns = [x1,x2,x3,x4,x5,x6,x7,x8]
        unknowns = unknowns[0:nombre]
        t = []
        vecs = []
        for e in eqs:
            d = e.as_coefficients_dict()
            c = []
            for k in unknowns:
                c.append(d[k])
            t.append(-d[1])
            vecs.append(Vector(c))
        a = Matriu.from_vectors_fila(vecs)
        b = Vector(t)
        return cls(a,b,prime)
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex del sistema d'equacions
        """
        p = ""
        eqs = list(map(str,self.equacions))
        if self.prime > 0:
            p = self.prime * "'"
        if len(eqs) == 0:
            return ""
        if len(eqs) == 1:
            for k in self.unknowns:
                eqs[0] = eqs[0].replace(latex(k),latex(k) + p)
            return eqs[0].replace('&','')
        eqs = " \\\\ ".join(eqs)
        for k in self.unknowns:
            eqs = eqs.replace(latex(k),latex(k) + p)
        return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"
    #
    #
    #
    def matriu_incognites(self):
        """
        Retorna la matriu dels coeficients de les incògnites expressada
        en LaTeX
        """
        return f"{self.A}"
    #
    #
    #
    def matriu_ampliada(self):
        """
        Retorna la matriu ampliada del sistema d'equacions expressada
        en LaTeX
        """
        c = self.A.inserta_columna(self.A.columnes,self.B)
        return matriu_latex(c.matriu,ampliada=True)
    #
    #
    #
    def resol(self,unknowns=None):
        """
        Resol el sistema d'equacions utilitzant la funció linsolve del sympy.
        El resultat és una llista d'expressions on hi poden aparèixer les
        incògnites del sistema com a paràmetres
        """
        system = (self.A.matriu,Matrix(self.B.dimensio,1,self.B.components))
        X = Vector(self.unknowns)
        system = self.A * X - self.B
        if unknowns is not None:
            s = solve(system.components,*unknowns)
        else:
            s = solve(system.components,*self.unknowns)
        self.solucio = []
        for k in self.unknowns:
            try:
                self.solucio.append(s[k])
            except:
                self.solucio.append(k)
        return self.solucio
    #
    #
    #
    def solucio_latex(self,linia=False):
        """
        Retorna l'expressió en LaTex de la solució del sistema d'equacions
        Paràmetres:
            linia: si és True escriu la solució en una línia, en cas contrari
                   ho fa com un sistema d'equacions
        """
        if self.solucio is None:
            self.resol()
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
        if not linia:
            eqs = " \\\\ ".join(eqs)
            return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"
        eqs = " $, $".join(eqs)
        eqs = eqs.replace('&','')
        return f"$ {eqs} $"

class EquacioParametrica:
    """
    Classe per treballar amb equacions paramètriques
    Atributs:
        equacio: l'equació paramètrica
        b: terme independent de l'equaqció
        coefs: coeficients dels paràmetres
        unknown: incògina de l'equació paramètrica
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
               en LaTeX de l'equació
        """
        x, y, z, t = symbols('x y z t')
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        unknowns = [x1,x2,x3,x4,x5,x6,x7,x8] + [x,y,z,t]
        self.equacio = eq
        self.amp = amp
        d = eq.as_coefficients_dict()
        self.b = 0
        self.coefs = {}
        self.params = []
        for k in d.keys():
            if k == 1:
                self.b = d[k]
            elif k in unknowns:
                self.unknown = k
            else:
                if k not in self.params:
                    self.params.append(k)
                self.coefs[k] = d[k]
        self.params.sort(key=latex)
    #
    #
    #
    @classmethod
    def coeficients(cls,a,b,p=0,total=1,amp=True):
        """
        Genera una equació amb coeficients dels paràmtres el vector "a", terme
        independent b i incògnita número p d'un total de "total".
        Paràmetres:
            a: vector amb els coeficients dels paràmetres
            b: terme independent
            p: índex que representa la incògnita
            total: nombre total d'incònites
        Exemple:
           e = EquacioParametrica(Vector(3,-2,1),5,1,4) genera l'equació
                        - y + 3*t1 - 2*t2 + t3 + 5
           e = EquacioParametrica(Vector(3,-2,1,3],-5,3,7) genera l'equació
               - x3 + 3*t1 - 2*t2 + t3 + 3*t4 + 7*t5 - 5
        Observació:
            si hi ha un màxim de quatre incògnites, son (x,y,z,t)
            si n'hi ha més, són (x1,x2,x3,x4,...)
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
        for k in range(a.dimensio):
            eq += a[k] * parameters[k]
        eq += b
        return cls(eq,amp)
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expessió en LaTeX de l'equació paramètrica
        """
        t1, t2, t3, t4, t5, t6, t7, t8 = symbols('t1 t2 t3 t4 t5 t6 t7 t8')
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        x, y, z, t = symbols('x y z t')
        if self.amp:
            s = f"{latex(self.unknown)} &= "
        else:
            s = f"{latex(self.unknown)} = "
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
                first = k
            e += factor * self.coefs[k] * k
        f = ""
        if first is not None and self.b != 0:
            if self.coefs[first] > 0:
                f = f"{factor*self.b} + "
            else:
                f = f"{factor*self.b} "
        s = s.replace('NUMERADOR',f + latex(e))
        return s

class EquacionsParametriques(object):
    """
    Classe per treballar amb sistemes d'equacions paramètriques
    Atributs:
        A: matriu dels coeficients dels paràmetres
        B: vector dels termes independents
        equacions: llista de EquacioParametrica
        nombre: nombre d'equacions
    """
    #
    #
    #
    def __new__(cls,a,b,amp=True):
        """
        Contructor.
        Genera les equacions X = b + a.T on X són les incògnites i T els paràmetres

        Paràmetres:
          a: matriu dels coeficients dels paràmetres
          b: vector de termes independents
          amp: True o False en funció si s'ha d'escriure &= o només = en la representació
               en LaTeX del sistema
        """
        if not isinstance(a,Matriu):
            return None
        if not isinstance(b,Vector):
            return None
        if a.files != b.dimensio:
            return None
        return super(EquacionsParametriques,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a,b,amp=True):
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
        Retorna la representació en LaTeX del sistema d'equacions paramètriques
        """
        l = list(map(str,self.equacions))
        eqs = " \\\\ ".join(l)
        return f"\\left.\\aligned {eqs} \\endaligned\\;\\right\\}}"
    #
    #
    #
    def eliminar_parametres(self,prime=0):
        """
        Retorna el SistemaEquacions que s'obté en eliminar els paràmetres dels
        sistema
        """
        L, U, _ = self.A.matriu.LUdecomposition()
        r = U.rank()
        v = Vector([e.unknown for e in self.equacions])
        t = Matrix(self.nombre,1,(v - self.B).components)
        t = (L**(-1) * t)[r:]
        return SistemaEquacions.from_equacions(t,self.A.files,prime)

class PlaVectorial(object):
    """
    Classe per treballar amb plans vectorials
    """
    #
    #
    #
    def __new__(cls,u1,u2):
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
        return super(PlaVectorial,cls).__new__(cls)
    #
    #
    #
    def __init__(self,u1,u2):
        self.u1 = u1
        self.u2 = u2
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en LaTeX dels dos vectors generadors
        """
        return f"{self.u1}, {self.u2}"
    #
    #
    #
    def equacio_implicita(self,base=None,prime=0):
        """
        Retorna l'equació implícita del pla en la base "base" en format LaTeX.
        Normalment, si la base no és la canònica es posa prime > 0 perquè el resultat
        sigui de l'estil 2x'-3y'+4z'= 0
        Paràmetres:
            base: una base (classe Base). Si és None, serà la canònica
            prime: nombre de primes que s'escriuran a les incògnites
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
        Crea el pla vectorial generat per les columnes de la matriu "m"
        Paràmetres:
            m: matriu. Ha de tenir 3 files, dues columnes i rang 2
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
        Paràmetres:
            v: vector no nul de dimensió 3
        """
        if not isinstance(w,Vector):
            return None
        if w.dimensio != 3:
            return None
        if w.length() == 0:
            return None
        a = Matriu.from_vectors_fila([w])
        l = a.nucli()
        return cls(l[0],l[1])
    #
    #
    #
    def base_ortogonal(self):
        """
        Retorna una base ortogonal del pla vectorial
        """
        v1 = self.u1
        u2 = self.u2
        v2 = v1.dot(v1) * u2 - v1.dot(u2) * v1
        v2.simplificar()
        return Base([v1,v2])

class RectaVectorial(object):
    """
    Classe per treballar amb rectes vectorials, dimensió 2 o 3
    Atributs:
        u: generador de la recta vectorial
    """
    #
    #
    #
    def __new__(cls,u):
        """
        Constructor.
        Paràmetres:
            u: vector amb dimensió 2 o 3
        """
        if not isinstance(u,Vector):
            return None
        if u.dimensio not in [2,3]:
            return None
        return super(RectaVectorial,cls).__new__(cls)
    #
    #
    #
    def __init__(self,u):
        self.u = u
    #
    #
    #
    def equacions_implicites(self,base=None,prime=0,aleatori=False):
        """
        Retorna l'equació implícita (dimensió 2) o el sistema d'equacions implícites
        (dimensió 3) en la base "base".
        Paràmetres:
            base: Base en la que calculem les equacions implícites
            prime: Quantes primes volem posar a les equacions
            aleatori: només s'aplica a dimensió 3 i genera unes equacions implícites no trivials
        """
        if base is None:
            a = Matriu.matriu_fila(self.u)
        else:
            a = Matriu.matriu_fila(self.u.components_en_base(base))
        l = a.nucli()
        if len(l) == 1:
            return EquacioLineal.coeficients(l[0],0,False,prime)
        a = Matriu.from_vectors_fila(l)
        b = Vector.nul(len(l))
        if aleatori:
            m = Matriu.amb_rang(f=2,c=2,r=2,maxim=3,mzeros=0)
            a = m * a
        return SistemaEquacions(a,b,prime=prime)
    #
    #
    #
    def equacio_continua(self,base=None,prime=0):
        """
        Retorna l'expressió en LaTeX de l'equació contínua de la recta vectorial
        en la base "base".
        Paràmetres:
            base: base del pla o de l'espai vectorial (classe Base)
            prime: nombre de primes amb el que s'escriuran les incògnites
        """
        if base is None:
            v = Vector(self.u.components)
        else:
            v = self.u.components_en_base(base)
        v.simplificar()
        x, y, z = symbols('x y z')
        incg = [x,y,z]
        eq = []
        p = ""
        if prime > 0:
            p = prime * "'"
        for i in range(v.dimensio):
            if v[i] == 1:
                eq.append(latex(incg[i]) + p)
            else:
                eq.append(f"\\frac{{{latex(incg[i]) + p}}}{{{v[i]}}}")
        return " = ".join(eq)

class ReferenciaAfi(object):
    """
    Classe per treballar amb referències afins de P^n
    Atributs:
        origen: origen de la referència (classe Punt)
        base: base de la referència (classe Base)
        dimensio: dimensió de l'espai corresponent
    """
    #
    #
    #
    def __new__(cls,origen,base):
        """
        Constructor.
        Paràmetres:
            origen: origen de la referència
            base: base de la referència
        """
        if not isinstance(origen,Punt):
            return None
        if not isinstance(base,Base):
            return None
        if origen.dimensio != base.dimensio:
            return None
        return super(ReferenciaAfi,cls).__new__(cls)
    #
    #
    #
    def __init__ (self,origen,base):
        self.origen = origen
        self.base = base
        self.dimensio = origen.dimensio
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en LaTeX de la referència
        """
        return f"\\left\\{{{self.origen};{self.base.vectors_latex()}\\right\\}}"
    #
    #
    #
    @classmethod
    def aleatoria(cls,dimensio=3,unitaria=False):
        """
        Retorna una referència aleatòria
        Paràmetres:
            dimensio: dimensió de l'espai corresponent
            unitaria: si és True la matriu del canvi de base tindrà determinant 1 o -1
        """
        origen = Punt.aleatori(l=dimensio,maxim=3,nuls=False)
        m = Matriu.invertible(ordre=3,maxim=3,mzeros=0,unitaria=unitaria)
        base = Base.from_matriu(m)
        return cls(origen,base)
    #
    #
    #
    def punt_de_coordenades(self,punt):
        """
        Retorna un nou punt expressat en la referència canònica del
        punt que en aquesta referencia té coordenades "punt"
        Paràmetres:
            punt: coordenades d'un punt en la referència actual
        """
        if not isinstance(punt,Punt):
            return None
        if punt.dimensio != self.dimensio:
            return None
        if not self.base.unitaria:
            c = self.base.matriu()
        else:
            unitaris = [(1 / v.length()) * v for v in self.base.vecs]
            c = Matriu.from_vectors_columna(unitaris)
        p = self.origen + c * punt
        return Punt(p.components)
    #
    #
    #
    def canvi_coordenades(self,prime1=0,prime2=1):
        """
        Restorna en format latex l'expressió del canvi de coordenades de la referència
        actual a la referència canònica
        Paràmetres:
            prime1: primes que s'escriuran a les coordenades en la referència canònica
            prime2: primes que s'escriuran a les coordenades en la referència actual
        """
        if self.dimensio <= 3:
            x, y, z = symbols('x y z')
            coords = [x,y,z]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            coords = [x1, x2, x3, x4, x5, x6, x7, x8]
        coords = coords[0:self.dimensio]
        p1 = ""
        p2 = ""
        if prime1 > 0:
            p1 = prime1 * "'"
        if prime2 > 0:
            p2 = prime2 * "'"
        o = " \\\\ ".join([latex(k) + p1 for k in coords])
        d = " \\\\ ".join([latex(k) + p2 for k in coords])
        m = Matriu.matriu_columna(self.origen)
        s = "\\begin{pmatrix}{c} " + o + "\\end{pmatrix} = \n"
        s += f"{m} + \n"
        s += f"{self.base.matriu()}\n"
        s += "\\begin{pmatrix}{c} " + d + "\\end{pmatrix}"
        return s
    #
    #
    #
    def referencia_inversa(self):
        """
        Retorna una referència que es correspon amb el canvi de coordenades de
        de la referència canònica a l'actual
        """
        p = Punt.nul(self.dimensio)
        q = p.coordenades_en_referencia(self)
        m = self.base.matriu().inversa()
        b = Base.from_matriu(m)
        return ReferenciaAfi(q,b)
    #
    #
    #
    def vectors(self,unitaris=False):
        """
        Retorna els vectors de la base de la referència
        Paràmetres:
            unitaris: si és True els retorna dividits per la seva longitud
        """
        return self.base.vectors(unitaris)

class PlaAfi(object):
    """
    Classe per treballar amb plans afins.
    Atributs:
        u1, u2: vectors directors del pla
        p: punt de pas
    """
    #
    #
    #
    def __new__(cls,p,u1,u2,ref=None):
        """
        Constructor.
        Paràmetres:
            p: punt de pas
            u1, u2: generadors del pla
            ref: referència en que estan expressats u1, u2 i p
        """
        if not isinstance(u1,Vector):
            return None
        if not isinstance(u2,Vector):
            return None
        if not isinstance(p,Punt):
            return None
        if u1.dimensio != 3 or u2.dimensio != 3:
            return None
        if p.dimensio != 3:
            return None
        if ref is not None:
            if not isinstance(ref,ReferenciaAfi):
                return None
            if ref.dimensio != 3:
                return None
        m = Matriu.from_vectors_columna([u1,u2])
        if m.rank() != 2:
            return None
        return super(PlaAfi,cls).__new__(cls)
    #
    #
    #
    def __init__(self,p,u1,u2,ref=None):
        if ref is None:
            self.u1 = u1
            self.u2 = u2
            self.p = p
        else:
            self.u1 = ref.base.vector_de_components(u1)
            self.u2 = ref.base.vector_de_components(u2)
            self.p = ref.punt_de_coordenades(p)
    #
    #
    #
    @classmethod
    def amb_associat(cls,w,p):
        """
        Genera el pla afí que té vector perpendicular "w" i passa pel punt p
        Paràmetres:
            w: vector associat
            p: punt de pas
        """
        if not isinstance(w,Vector):
            return None
        if not isinstance(p,Punt):
            return None
        a = Matriu.from_vectors_fila([w])
        l = a.nucli()
        return cls(p,l[0],l[1])
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'equació vectorial del pla en LaTeX
        """
        return f"(x,y,z)={self.p}+t_1{self.u1}+t_2{self.u2}"
    #
    #
    #
    def equacio_implicita(self,ref=None,prime=0):
        """
        Retorna l'expressió de l'equació implícita del pla en la referència "ref"
        Paràmetres:
            ref: referència afí
            prime: nombre de primes que s'escriuran a les incògnites
        """
        if ref is None:
            w = self.u1.cross(self.u2)
            return EquacioLineal.coeficients(w,w.dot(self.p),False)
        if not isinstance(ref,ReferenciaAfi):
            return None
        c = self.p.coordenades_en_referencia(ref)
        v1 = self.u1.components_en_base(ref.base)
        v2 = self.u2.components_en_base(ref.base)
        w = v1.cross(v2)
        return EquacioLineal.coeficients(w,w.dot(c),False,prime)
    #
    #
    #
    def base_ortogonal(self):
        """
        Retorna una base orogonal (vectors directors perpendiculars) del pla
        """
        v1 = self.u1
        u2 = self.u2
        v2 = v1.dot(v1) * u2 - v1.dot(u2) * v1
        v2.simplificar()
        return [v1,v2]
    #
    #
    #
    def associat(self,base=None):
        """
        Retorna un vector perpendicular al pla expressat en la base "base"
        Paràmetres:
            base: base de l'espai vectorial. Si és None, serà la canònica
        """
        w = self.u1.cross(self.u2,simplificar=True)
        if base is None:
            return w
        return w.components_en_base(base)
    #
    #
    #
    def distancia(self,other):
        """
        Retorna la distància entre el pla actual i un punt, una recta o un altre pla
        Paràmetres:
            other: un punt (classe Punt), una recta (classe RectaAfi) o
            un pla (class PlaAfi)
        """
        if isinstance(other,Punt):
            v = self.p - other
            w = self.u1.cross(self.u2,simplificar=True)
            return abs(v.dot(w))/w.length()
        if isinstance(other,PlaAfi):
            w1 = self.associat()
            w2 = other.associat()
            w = w1.cross(w2)
            if w.length() == 0:
                return self.distancia(other.p)
            return 0
        if isinstance(other,RectaAfi):
            return other.distancia(self)
        return None

class RectaAfi(object):
    """
    Classe per treballar amb rectes afins, dimensió 2 o 3
    Atributs:
        u: generador de la recta vectorial
        p: punt de pas
    """
    #
    #
    #
    def __new__(cls,p,u,ref=None):
        """
        Constructor.
        Paràmetres:
            p: punt de pas
            u: vector director de la recta
            ref: referència en que estan expressats u i p
        """
        if not isinstance(u,Vector):
            return None
        if not isinstance(p,Punt):
            return None
        if p.dimensio != u.dimensio:
            return None
        if p.dimensio not in [2,3]:
            return None
        if ref is not None:
            if not isinstance(ref,ReferenciaAfi):
                return None
            if ref.dimensio != p.dimensio:
                return None
        return super(RectaAfi,cls).__new__(cls)
    #
    #
    #
    def __init__(self,p,u,ref=None):
        if ref is None:
            self.u = u
            self.p = p
        else:
            self.u = ref.base.vector_de_components(u)
            self.p = ref.punt_de_coordenades(p)
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'equació vectorial de la recta en LaTeX
        """
        if self.u.dimensio == 2:
            return f"(x,y)={self.p}+t{self.u}"
        else:
            return f"(x,y,z)={self.p}+t{self.u}"
    #
    #
    #
    def equacio_continua(self,ref=None,prime=0):
        """
        Retorna l'expressió en latex de l'equqció contínua de la recta afí
        en la referència "ref".
        Paràmetres:
            ref: referència afí. Si és None, serà la canònica
            prime: nombre de primes que s'escriuran a les incògnites
        """
        if ref is None:
            v = Vector(self.u.components)
            q = Punt(self.p.components)
        else:
            v = self.u.components_en_base(ref.base)
            q = self.p.coordenades_en_referencia(ref)
        v.simplificar()
        x, y, z = symbols('x y z')
        incg = [x,y,z]
        eq = []
        p = ""
        if prime > 0:
            p = prime * "'"
        for i in range(v.dimensio):
            if v[i] == 1:
                eq.append(latex(incg[i] - q[i]))
            else:
                eq.append(f"\\frac{{{latex(incg[i] - q[i])}}}{{{v[i]}}}")
        eq = " = ".join(eq)
        for i in incg:
            eq = eq.replace(latex(i),latex(i) + p)
        return eq
    #
    #
    #
    def equacions_implicites(self,ref=None,prime=0,aleatori=False):
        """
        Retorna l'equació implícita (dimensió 2) o el sistema d'equacions implícites
        (dimensió 3) de la recta afí en la referència "ref".
        Paràmetres:
           ref: referència en la que calculem les equacions implícites
           prime: nombre de primes que s'escriran a les incògnites
           aleatori: només s'aplica a dimensió 3 i genera unes equacions implícites
           amb tots els coeficients de les incògnites no nuls
        """
        if ref is None:
            v = Vector(self.u.components)
            a = Matriu.matriu_fila(self.u)
            q = Punt(self.p.components)
        else:
            v = self.u.components_en_base(ref.base)
            a = Matriu.matriu_fila(v)
            q = self.p.coordenades_en_referencia(ref)
        l = a.nucli()
        if len(l) == 1:
            return EquacioLineal.coeficients(l[0],l[0].dot(q),False,prime)
        a = Matriu.from_vectors_fila(l)
        b = a * q
        if aleatori:
            trobat = False
            while not trobat:
                m = Matriu.amb_rang(f=2,c=2,r=2,maxim=3,mzeros=0)
                aux = m * a
                trobat = aux.nzeros() == 0
            a = m * a
            b = m * b
        return SistemaEquacions(a,b,prime=prime)
    #
    #
    #
    def distancia(self,other):
        """
        Retorna la distància entre la recta actual i un punt, una recta o un pla
        Paràmetres:
            other: un punt (classe Punt), una recta (classe RectaAfi) o un plan
            (class PlaAfi)
        """
        if isinstance(other,Punt):
            v = self.p - other
            w = self.u.cross(v)
            return w.length()/self.u.length()
        if isinstance(other,RectaAfi):
            w = self.u.cross(other.u,simplificar=True)
            if w.length() > 0:
                u = self.p - other.p
                return abs(u.dot(w) / w.length())
            return self.distancia(other.p)
        if isinstance(other,PlaAfi):
            w = other.associat()
            if w.dot(self.u) == 0:
                return other.distancia(self.p)
            return 0
        return None
    #
    #
    #
    def punt(self,t):
        """
        Retorna el punt de la recta amb paràmetre t
        Paràmetres:
            t: escalar
        """
        q = self.p + t * self.u
        return Punt(q.components)

class SubespaiVectorial(object):
    """
    Classe per treballar amb subespais vectorials
    Atributs:
       generadors: generadors del subespai
       base: base del subespai
       dimensio: dimensio del subespai
       espai: n si és un suespai de R^n
    """
    #
    #
    #
    def __new__(cls,vecs,basern=None):
        """
        Constructor.
        Retorna el subespai vectorial generat per una llista de vectors
        Paràmetres:
           vecs: llista de vectors
           basern: base en la que estan expressats els vectors
        """
        a = Matriu.from_vectors_columna(vecs)
        if a is None:
            return None
        if basern is not None:
            if not isinstance(basern,Base):
                return None
            if basern.dimensio != a.files:
                return None
        return super(SubespaiVectorial,cls).__new__(cls)
    #
    #
    #
    def __init__(self,vecs,basern=None):
        if basern is None:
            for k in vecs:
                k.simplificar()
            self.generadors = [k for k in vecs]
        else:
            v = [basern.vector_de_components(k) for k in vecs]
            for k in v:
                k.simplificar()
            self.generadors = [k for k in v]
        a = Matriu.from_vectors_columna(self.generadors)
        L, U, _ = a.matriu.LUdecomposition()
        files = Matriu(U).vectors_fila()
        self.base = []
        for f in files:
            k = primer_no_nul(f.components)
            if k is not None:
                self.base.append(self.generadors[k])
        self.dimensio = len(self.base)
        self.espai = self.generadors[0].dimensio
    #
    #
    #
    @classmethod
    def from_equacions_implicites(cls,eqs,basern=None):
        """
        Retorna el subespai vectorial que té equacions implícites "eqs"
        Paràmetres:
            eqs: equacions implícites (classe SistemaEquacions)
            basern: Base de R^n
        """
        if basern is None:
            return cls(eqs.A.nucli())
        c = [basern.vector_de_components(u) for u in eqs.A.nucli()]
        s = cls(c)
        eqs = s.equacions_implicites()
        return cls.from_equacions_implicites(eqs)
    #
    #
    #
    @classmethod
    def from_equacio_implicita(cls,eq,basern=None):
        s = eq.to_sistema_equacions()
        return cls.from_equacions_implicites(s,basern)
    #
    #
    #
    def es_zero(self):
        """
        Retorna True si és el subespai {0}
        """
        return len(self.base) == 0
    #
    #
    #
    def es_total(self):
        """
        Retorna True si és el subespai R^n
        """
        return len(self.base) == self.espai
    #
    #
    #
    def equacions_implicites(self,basern=None,prime=0):
        """
        Retorna unes equacions implícites del subespai
        Paràmetres:
            basern: base en la que s'escriuran les equacions implícites
            prime: nombre de primes a les incògnites
        """
        if self.es_total():
            return None
        b = Vector.nul(self.espai)
        if basern is None:
            a = Matriu.from_vectors_columna(self.base)
        else:
            c = [v.components_en_base(basern) for v in self.base]
            a = Matriu.from_vectors_columna(c)
        p = EquacionsParametriques(a,b)
        return (p.eliminar_parametres(prime))
    #
    #
    #
    def suplementari_ortogonal(self):
        """
        Retorna el suplementari ortogonal
        """
        a = Matriu.from_vectors_fila(self.base)
        return SubespaiVectorial(a.nucli())
    #
    #
    #
    def base_ortogonal(self):
        """
        Retorna una base ortogonal del subespai
        """
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
    #
    #
    #
    def amplia_base(self,unitaria=False):
        """
        Retorna una base ortogonal amb orientació positiva de R^n
        que comença amb una base del subespai
        Paràmetres:
            unitaria: si és True, retorna una base ortonormal
        """
        if self.es_total():
            b = Base(self.base_ortogonal(),unitaria)
        else:
            h = self.suplementari_ortogonal()
            b = Base(self.base_ortogonal() + h.base_ortogonal(),unitaria)
        b.orientacio_positiva()
        return b
    #
    #
    #
    def amplia_base_suplementari(self,unitaria=False):
        """
        Retorna una base ortogonal amb orientació positiva de R^n
        que comença amb una base del suplementari ortogonal del subespai
        Paràmetres:
            unitaria: si és True, retorna una base ortonormal
        """
        h = self.suplementari_ortogonal()
        b = Base(h.base_ortogonal() + self.base_ortogonal(),unitaria)
        b.orientacio_positiva()
        return b

class VarietatLineal(object):
    """
    Classe per treballar amb varietats lineal
    Atributs:
       punt: punt de pas
       subespai: SubespaiVectorial
    """
    #
    #
    #
    def __new__(cls,p,s,ref=None):
        """
        Constructor.
        Retorna la varietat lineal que passa per p i té subespai director s
        Paràmetres:
           p: punt de pas
           s: subespai director
           ref: referència afí en la que estan expressats el punt p
        """
        if not isinstance(p,Punt):
            return None
        if not isinstance(s,SubespaiVectorial):
            return None
        if p.dimensio != s.espai:
            return None
        if ref is not None:
            if not isinstance(ref,ReferenciaAfi):
                return None
            if ref.dimensio != p.dimensio:
                return None
        return super(VarietatLineal,cls).__new__(cls)
    #
    #
    #
    def __init__(self,p,s,ref=None):
        self.subespai = s
        if ref is None:
            self.punt = p
        else:
            self.punt = ref.punt_de_coordenades(p)
    #
    #
    #
    def es_un_punt(self):
        """
        Retorna True si la varietat lineal és un punt
        """
        return self.subespai.es_zero()
    #
    #
    #
    def es_total(self):
        """
        Retorna True si la varietat lineal és P_n
        """
        return self.subespai.es_total()
    #
    #
    #
    def equacions_implicites(self,ref=None,prime=0):
        """
        Retorna unes equacions implícites de la varietat lineal
        Paràmetres:
            ref: referència en la que s'escriuran les equacions implícites
            prime: nombre de primes de les incògnites
        """
        if self.es_total():
            return None
        if ref is None:
            q = self.punt
            base = self.subespai.base
        else:
            q = self.punt.coordenades_en_referencia(ref)
            base = [v.components_en_base(ref.base) for v in self.subespai.base]
        a = Matriu.from_vectors_columna(base)
        p = EquacionsParametriques(a,q)
        return (p.eliminar_parametres(prime))
    #
    #
    #
    def base_ortogonal(self):
        """
        Retorna una base ortogonal del subespai director de la varietat lineal
        """
        return self.subespai.base_ortogonal()
    #
    #
    #
    def varietat_ortogonal(self,p):
        """
        Retorna la varietat ortogonal a l'actual que passa pel punt p
        """
        if not isinstance(p,Punt):
            return None
        if p.dimensio != self.subespai.espai:
            return None
        s = self.subespai.suplementari_ortogonal()
        return VarietatLineal(p,s)

class TransformacioLineal(object):
    """
    Classe per treballar amb transformacions lineals T:R^n ----> R^n
    Atributs:
       dimensio: n
       canonica: matriu de la transformació en la base canònica
    """
    #
    #
    #
    def __new__(cls,matriu,base=None):
        """
        Constructor.
        Retorna una nova transformació lineal
        Paràmetres:
            matriu: matriu de la transformació lineal en la base "base"
            base: base de R^n. Si és None serà la canònica
        """
        if not isinstance(matriu,Matriu):
            return None
        if matriu.files != matriu.columnes:
            return None
        if base is not None:
            if not isinstance(base,Base):
                return None
            if base.dimensio != matriu.files:
                return None
        return super(TransformacioLineal,cls).__new__(cls)
    #
    #
    #
    def __init__(self,matriu,base=None):
        self.dimensio = matriu.files
        if base is None:
            self.canonica = Matriu(matriu.matriu)
        else:
            c = base.matriu()
            m = c.matriu * matriu.matriu * c.matriu**(-1)
            self.canonica = Matriu(m)
    #
    #
    #
    def determinant(self):
        """
        Retorna el deterinant de la transformació lineal
        """
        return self.canonica.determinant()
    #
    #
    #
    def es_rotacio(self):
        """
        Ens diu si és una rotació tridimensional o no
        """
        if self.dimensio != 3:
            return False
        d = self.determinant()
        if abs(N(d) - 1.0) > 10**(-8):
            return False
        m = self.canonica
        u = Matriu()
        k = m.transposada() * m
        for i in range(k.files):
            for j in range(k.columnes):
                p = N(k[i,j])
                if abs(p - u[i,j]) > 10**(-8):
                    return False
        return True
    #
    #
    #
    def eix_angle_rotacio(self,radiants=False):
        """
        Retorna l'eix i l'angle de rotació
        Paràmetres:
            radiants: si és True retorna l'angle en radiants, en cas contrari,
            ho fa en graus
        """
        if not self.es_rotacio():
            return None
        m = self.canonica - Matriu()
        e = m.nucli()[0]
        e = Vector([simplify(k.expand()) for k in e.components])
        q = Quaternion.from_rotation_matrix(self.canonica.matriu)
        v, alpha = q.to_axis_angle()
        alpha = alpha.expand()
        v = [simplify(k.expand()) for k in v]
        u = Vector(v)
        l = e.length()
        u = l * u
        if e[0] != u[0]:
            alpha = 2 * pi - alpha
        if not radiants:
            alpha = alpha * 180/pi
        return (e,alpha)
    #
    #
    #
    def angles_euler(self,radiants=False):
        """
        Retorna els angles d'Euler de la rotació
        Paràmetres:
            radiants: si és True retorna els angles en radiants, en cas contrari,
            ho fa en graus
        """
        if not self.es_rotacio():
            return None
        m = self.canonica
        if abs(m[2,2]) != 1:
            theta = acos(m[2,2])
            psi = atan2(m[2,0],m[2,1])
            phi = atan2(m[0,2],-m[1,2])
        else:
            phi = 0
            theta = acos(m[2,2])
            psi = atan2(m[2,2] * m[1,0],m[0,0])
        if not radiants:
            theta *= 180 / pi
            psi *= 180 / pi
            phi *= 180 / pi
        return (psi,theta,phi)
    #
    #
    #
    @classmethod
    def gir(cls,angle,radiants=False):
        """
        Retorna el gir d'angle "angle" en dimensió 2
        Paràmetres:
            angle: angle de rotació
            radiants: si l'angle està en radiants, ha de ser True
        """
        if not radiants:
            angle *= pi / 180
        m = Matriu.from_vectors_columna([Vector(cos(angle),sin(angle)),Vector(-sin(angle),cos(angle))])
        return cls(m)
    #
    #
    #
    @classmethod
    def rotacio(cls,eix,angle,radiants=False):
        """
        Retorna la rotació d'angle "angle" al voltant del vector "eix"
        Paràmetres:
            eix: vector al voltant del qual fem la rotació
            angle: angle de rotació
            radiants: si l'angle està en radiants, ha de ser True
        """
        if not radiants:
            angle *= pi / 180
        v = Vector(eix.components)
        v.normalitzar()
        q = Quaternion.from_axis_angle(v.components,angle)
        m = Matriu(q.to_rotation_matrix())
        for i in range(m.files):
            for j in range(m.columnes):
                m[(i,j)] = simplify(m[(i,j)].expand())
        return cls(m)
    #
    #
    #
    @classmethod
    def projeccio_ortogonal(cls,s):
        """
        Retorna la projecció ortogonal sobre el subespai "s"
        Paràmetres:
            s: subespai vectorial (classe SubespaiVectorial)
        """
        if not isinstance(s,SubespaiVectorial):
            return None
        n = s.espai
        m = s.dimensio
        uns = [1 for k in range(m)]
        zeros = [0 for k in range(n-m)]
        d = uns + zeros
        b = Base(s.base + s.suplementari_ortogonal().base)
        m = Matriu.diagonal(d)
        return cls(m,b)
    #
    #
    #
    @classmethod
    def simetria(cls,s):
        """
        Retorna la simetria respecte al subespai "s"
        Paràmetres:
            s: subespai vectorial (classe SubespaiVectorial)
        """
        if not isinstance(s,SubespaiVectorial):
            return None
        n = s.espai
        m = s.dimensio
        uns = [1 for k in range(m)]
        menys = [-1 for k in range(n-m)]
        d = uns + menys
        b = Base(s.base + s.suplementari_ortogonal().base)
        m = Matriu.diagonal(d)
        return cls(m,b)
    #
    #
    #
    def matriu_en_base(self,base):
        """
        Retorna la matriu de la transformacio lineal en la base "base"
        Paràmetres:
            base: base de R^n (classe Base)
        """
        if base is None:
            return self.canonica
        c = base.matriu()
        m = c.matriu**(-1) * self.canonica.matriu * c.matriu
        return Matriu(m)
    #
    #
    #
    def latex(self,base=None,prime=0):
        """
        Retorna l'expressió en latex de la transformació lineal en la base "base"
        Si base és None, serà en la base canònica
        Paràmetres:
            base: base de R^n
            prime: nombre de primes que s'han d'escriure
        """
        if self.dimensio <= 4:
            x, y, z, t = symbols('x y z t')
            u, v, w, r = symbols('u v w r')
            o = [x,y,z,t]
            d = [u,v,w,r]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            u1, u2, u3, u4, u5, u6, u7, u8 = symbols('u1 u2 u3 u4 u5 u6 u7 u8')
            o = [x1, x2, x3, x4, x5, x6, x7, x8]
            d = [u1, u2, u3, u4, u5, u6, u7, u8]
        o = o[0:self.dimensio]
        d = d[0:self.dimensio]
        p = prime * "'"
        o = " \\\\ ".join([latex(x)+p for x in o])
        d = " \\\\ ".join([latex(x)+p for x in d])
        if base is None:
             m = self.canonica
        else:
            if not isinstance(base,Base):
                return None
            c = base.matriu()
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
        """
        Retorna el latex de l'expressió en la base canònica de la transformació lineal
        """
        return self.latex()
    #
    #
    #
    def __eq__(self,other):
        """
        Permet comparar transfonacions lineals amb ==
        Paràmetres:
            other: una altra transformació lineal
        """
        return self.canonica == other.canonica
    #
    #
    #
    def __add__(self,other):
        """
        Suma de transformacions lineals
        Paràmetres:
            other: una altra transformació lineal
        """
        if other.dimensio != self.dimensio:
            return None
        m = self.canonica + other.canonica
        return TransformacioLineal(m)
    #
    #
    #
    def __mul__(self,other):
        """
        Composició de transformacions lineal o càlcul de imatges
        Paràmetres:
            other: transformació lineal o vector
        """
        if isinstance(other,TransformacioLineal):
            if other.dimensio != self.dimensio:
                return None
            m = self.canonica * other.canonica
            return TransformacioLineal(m)
        if isinstance(other,Vector):
            if other.dimensio != self.dimensio:
                return None
            return self.canonica * other
        return None
    #
    #
    #
    def es_simetrica(self):
        """
        Retorna si la transformació lineal és simètrica
        """
        return self.canonica == self.canonica.transposada()
    #
    #
    #
    def transforma(self,vec,base=None):
        """
        Calcula el transformat (image) del vector "vec".
        Paràmetres:
            vec: vector
            base: si no és None, vec seran les compoents del vectors en aquesta
            base. El transformat o imatge també estarà expressat en aquesta base
        """
        if not isinstance(vec,Vector):
            return None
        if vec.dimensio != self.dimensio:
            return None
        if base is None:
            t = self.canonica * vec
        else:
            if not isinstance(base,Base):
                return None
            c = base.matriu()
            m = c.matriu**(-1) * self.canonica.matriu * c.matriu
            m = Matriu(m)
            t = m * vec
        for k in range(t.dimensio):
            t[k] = simplify(t[k].expand())
        return t
    #
    #
    #
    def polinomi_caracteristic(self):
        """
        Retorna el polinomi característic de la transformació lineal
        """
        return self.canonica.polinomi_caracteristic()

class TransformacioAfi:
    """
    Classe per treballar amb transformacions afins T:P^n ----> P^n, on
    T(p) = T + A(p)
    Atributs:
        transformacio: Transformació lineal donada per la matriu A
        translacio: translació de la transformació afí en la referència canònica
        dimensio: n
    """
    #
    #
    #
    def __new__(cls,p,t):
        """
        Contructor:
        Retorna una transformació afí
        Paràmetres:
           p: Punt que representa la translació en la referencia canònica,
              és a dir, el transformat del punt zero
           t: Transformació lineal
        """
        if not isinstance(p,Vector):
            return None
        if not isinstance(t,TransformacioLineal):
            return None
        if p.dimensio != t.dimensio:
            return None
        return super(TransformacioAfi,cls).__new__(cls)
    #
    #
    #
    def __init__(self,p,t):
        self.translacio = p
        self.transformacio = t
        self.dimensio = self.transformacio.dimensio
    #
    #
    #
    @classmethod
    def gir(cls,angle,origen,radiants=False):
        """
        Retorna el gir d'angle "angle" al voltant del punt "origen"
        Paràmetres:
            origen: centre del gir (classe Punt)
            angle: angle de rotació
            radiants: si és True, l'angle ha d'estar expressat en radiants
        """
        if not isinstance(origen,Punt):
            return None
        if origen.dimensio != 2:
            return None
        t = TransformacioLineal.gir(angle,radiants)
        p = - t.transforma(origen) + origen
        return cls(p,t)
    #
    #
    #
    @classmethod
    def simetria(cls,v):
        """
        Retorna la simetria respecte a la variatat lineal "v"
        Paràmetres:
            v: varietat lineal (classe VarietatLineal)
        """
        if not isinstance(v,VarietatLineal):
            return None
        t = TransformacioLineal.simetria(v.subespai)
        p = - t.transforma(v.punt) + v.punt
        return cls(p,t)
    #
    #
    #
    @classmethod
    def projeccio_ortogonal(cls,v):
        """
        Retorna la projecció ortogonal sobre la variatat lineal "v"
        Paràmetres:
            v: varietat lineal (classe VarietatLineal)
        """
        if not isinstance(v,VarietatLineal):
            return None
        t = TransformacioLineal.projeccio_ortogonal(v.subespai)
        p = - t.transforma(v.punt) + v.punt
        return cls(p,t)
    #
    #
    #
    @classmethod
    def moviment_helicoidal(cls,recta,angle,radiants=False,alpha=0):
        """
        Retorna el moviment helicoidal de P_3 que consisteix en la rotació
        d'angle "angle" al voltant de la recta "recta" seguit d'una translació
        de vector alpha * vector director de la recta.
        Paràmetres:
            recta: eix de rotació
            angle: angle de rotació
            radiants: si és True, l'angle ha d'estar expressat en radiants
            alpha: factor de la translació
        """
        if not isinstance(recta,RectaAfi):
            return None
        if not recta.u.dimensio == 3:
            return None
        rotacio = TransformacioLineal.rotacio(recta.u,angle,radiants)
        p = - rotacio.transforma(recta.p) + recta.p + alpha * recta.u
        return cls(p,rotacio)
    #
    #
    #
    @classmethod
    def rotacio(cls,recta,angle,radiants=False):
        """
        Retorna la rotació d'angle "angle" al voltant de la recta "recta"
        Paràmetres:
            recta: eix de rotació
            angle: angle de rotació
            radiants: si és True, l'angle ha d'estar expressat en radiants
        """
        return cls.moviment_helicoidal(recta,angle,radiants)
    #
    #
    #
    def __repr__(self):
        """
        Restorna en format latex l'expressió de la transformació afí en la
        referència canònica
        """
        return self.latex()
    #
    #
    #
    def latex(self,ref=None,prime=0):
        """
        Retorna l'expressió en latex de la transformació afí en la referència
        "ref". Si ref és None, serà en la referència canònica
        Paràmetres:
            base: referència de P^n
            prime: nombre de primes que s'escriuran
        """
        if self.dimensio <= 4:
            x, y, z, t = symbols('x y z t')
            u, v, w, r = symbols('u v w r')
            o = [x,y,z,t]
            d = [u,v,w,r]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            u1, u2, u3, u4, u5, u6, u7, u8 = symbols('u1 u2 u3 u4 u5 u6 u7 u8')
            o = [x1, x2, x3, x4, x5, x6, x7, x8]
            d = [u1, u2, u3, u4, u5, u6, u7, u8]
        o = o[0:self.dimensio]
        d = d[0:self.dimensio]
        p = prime * "'"
        o = " \\\\ ".join([latex(x)+p for x in o])
        d = " \\\\ ".join([latex(x)+p for x in d])
        z = Punt.nul(self.dimensio)
        tz = self.transforma(z,ref)
        if ref is None:
            m = self.transformacio.canonica
        else:
            m = self.transformacio.matriu_en_base(ref.base)
        s = "\\begin{pmatrix}{c} " + d + "\\end{pmatrix} = \n"
        if tz.length() > 0:
            s += f"{Matriu.matriu_columna(tz)} + \n"
        s += f"{m}\n"
        s += "\\begin{pmatrix}{c} " + o + "\\end{pmatrix}"
        return s
    #
    #
    #
    def __eq__(self,other):
        """
        Permet comparar transfonacions afins amb ==
        Paràmetres:
            other: una altra transformació afí
        """
        return self.canonica == other.canonica and self.translacio == other.translacio
    #
    #
    #
    def __add__(self,other):
        """
        Suma de transformacions afins
        Paràmetres:
            other: una altra transformació afi
        """
        if other.dimensio != self.dimensio:
            return None
        m = self.canonica + other.canonica
        t = TansformacioLineal(m)
        tr = self.translacio + other.translacio
        return TransformacioAfi(tr,t)
    #
    #
    #
    def __mul__(self,other):
        """
        Composició de transformacions afins o càlcul de imatges de punts
        Paràmetres:
            other: transformació afi o punt
        """
        if isinstance(other,TransformacioAfi):
            if other.dimensio != self.dimensio:
                return None
            t = self.translacio + self.transformacio * other.translacio
            t = Punt(list(map(simplify,t.components)))
            tr = self.transformacio * other.transformacio
            return TransformacioAfi(t,tr)
        if isinstance(other,Punt):
            if other.dimensio != self.dimensio:
                return None
            t = self.translacio + self.transformacio * other
            return Punt(list(map(simplify,t.components)))
        return None
    #
    #
    #
    def transforma(self,p,ref=None):
        """
        Calcula el transformat o imatge del punt "p". p seran les coordenades
        d'aquest punt en la referència "ref" i el transformat també estarà expressat
        en aquesta referència
        Paràmetres:
            punt: punt (classe Punt)
            ref: referència afí. Si és None, serà la referència canònica
        """
        if not isinstance(p,Vector):
            return None
        if self.dimensio != p.dimensio:
            return None
        if ref is None:
            return self.translacio + self.transformacio.transforma(p)
        q = ref.punt_de_coordenades(p)
        p = Punt(self.transforma(q).components)
        return p.coordenades_en_referencia(ref)

class FormaQuadratica(object):
    """
    Classe per treballar amb formes quadràtiques
    Atributs:
        matriu: matriu de la forma quadrètica en la base canònica
        dimensio: n
        vaps: valors propis
        base: base oronormal en la que diagonalitza
    """
    #
    #
    #
    def __new__(cls,matriu,base=None):
        """
        Contructor:
        Retorna una forma quadràtica
        Paràmetres:
           matriu: matriu simètrica de la forma quadratica en la base "base"
           base: base ortonormal de R^n. Si és None, serà la canònica
        """
        if not isinstance(matriu,Matriu):
            return None
        if matriu.files != matriu.columnes:
            return None
        if matriu != matriu.transposada():
            return None
        if base is not None:
            if not isinstance(base,Base):
                return None
            if not base.es_ortogonal():
                return None
            if not base.es_unitaria():
                return None
        return super(FormaQuadratica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,matriu,base=None):
        self.dimensio = matriu.columnes
        if base is None:
            self.matriu = matriu
        else:
            c = base.matriu()
            self.matriu = c * matriu * c.transposada()
        e = self.matriu.matriu.eigenvects()
        vaps = []
        veps = []
        va, ve = vaps_veps_amb_signe(e,signe=1)
        vaps += va
        veps += ve
        va, ve = vaps_veps_amb_signe(e,signe=-1)
        vaps += va
        veps += ve
        s = SubespaiVectorial(veps)
        self.base = s.amplia_base(unitaria=True)
        self.vaps = Vector(vaps + (self.dimensio - len(vaps)) * [0])
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'expressió en latex de la forma quadràtica com a polinomi
        de segon grau en la base canònica
        """
        return self.latex()
    #
    #
    #
    def latex(self,base=None,prime=0):
        """
        Retorna l'expressió en latex de la forma quadràtica com a polinomi
        de segon grau en la base "base"
        Paràmetres:
            base: base de R^n. No cal que sigui ortormal. Si és None, serà la
            base canònica
            prime: nombre de primes que s'escriuran a les variables del polinomi
            segon grau
        """
        if self.dimensio <= 4:
            x, y, z, t = symbols('x y z t')
            s = [x,y,z,t]
        else:
            x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
            s = [x1, x2, x3, x4, x5, x6, x7, x8]
        s = Vector(s[0:self.dimensio])
        m = Matriu.matriu_columna(s)
        if base is None:
            q = self.matriu
        else:
            c = base.matriu()
            q = c.transposada() * self.matriu * c
        r = m.transposada() * q * m
        if self.dimensio <= 4:
            expr = mylatex(r[0,0].expand())
        else:
            expr = latex(r[0,0].expand())
        if prime == 0:
            return expr
        p = prime * "'"
        for k in s.components:
            expr = expr.replace(latex(k),latex(k)+p)
        return expr
    #
    #
    #
    @classmethod
    def aleatoria(cls,ordre=3,maxim=20,vapsnonuls=2):
        """
        Retorna una forma quadràtica aleatòria
        Paràmetres:
            ordre: ordre de la matriu simètrica de la forma quadràtica
            maxim: nombre màxim dels coeficients de la matriu
            vapsnonuls: nombre mínim de valrs propis no nuls
        """
        trobat = False
        while not trobat:
            vaps = Vector.aleatori(l=ordre,maxim=3)
            if len(set(vaps.components)) == 1:
                continue
            nonuls = [x for x in vaps.components if x != 0]
            if len(nonuls) < vapsnonuls:
                continue
            d = Matriu.diagonal(vaps)
            b = Base.ortogonal(ordre=ordre,maxim=5)
            c = b.matriu()
            m = c * d * c.transposada()
            if m.norma_maxim() > maxim:
                continue
            trobat = not m.es_diagonal()
        return cls(m)
    #
    #
    #
    def signatura(self):
        """
        Retorna la signatura o índexs d'inèrcia de la forma quadràtica
        """
        r = 0
        s = 0
        for k in self.vaps.components:
            if k > 0:
                r += 1
            elif k < 0:
                s += 1
        return (r,s)
    #
    #
    #
    def classificacio(self):
        """
        Retorna la classificació de la forma quadràtica
        """
        r, s = self.signatura()
        if r == self.dimensio:
            return "definida positiva"
        if s == self.dimensio:
            return "definida negativa"
        if s == 0:
            return "semidefinida positiva"
        if r == 0:
            return "semidefinida negativa"
        return "indefinida"
    #
    #
    #
    def rank(self):
        """
        Retorna el rang de la forma quadràtica
        """
        r, s = self.signatura()
        return r + s
    #
    #
    #
    def polinomi_caracteristic(self):
        """
        Retorna el polinomi característic de la forma quadràtica
        """
        return self.matriu.polinomi_caracteristic()
    #
    #
    #
    def expressio_euclidiana(self,prime=1):
        """
        Retorna l'expressió euclidiana reduïda de la forma de polinomi
        expressat en LaTeX
        Paràmetres:
            prime: nombre de primes que s'escriuran a les variables
        """
        return self.latex(self.base,prime=prime)

class Conica(object):
    """
    Classe per treballar amb còniques. L'objectiu no és classificar còniques,
    sinó generar coniques a partir dels elements característics o de manera
    aleatòria.
    Atributs:
        ref: referència afí
        matriu: matriu projectiva de la cònica en la referència "ref"
        canonica: matriu projectiva de la cònica en la referència canònica
    """
    #
    #
    #
    def __new__(cls,matriu,ref=None):
        """
        Retorna una nova cònica
        Paràmetres:
            matriu: matriu projectiva de la cònica en la referència "ref"
            Ha de ser 3x3
            ref: referència afí. Si és None, serà la referència canònica
        """
        if not isinstance(matriu,Matriu):
            return None
        if matriu.files != matriu.columnes:
            return None
        if matriu.files != 3:
            return None
        if not matriu.es_simetrica():
            return None
        if ref is not None:
            if ref.dimensio != 2:
                return None
            if not ref.base.es_unitaria():
                return None
            if not ref.base.es_ortogonal():
                return None
        return super(Conica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,matriu,ref=None):
        self.ref = ref
        self.matriu = matriu
        if ref is None:
            self.canonica = matriu
        else:
            a = ref.base.matriu()
            a = a.inserta_columna(2,ref.origen)
            a = a.inserta_fila(2,Vector([0,0,1]))
            b = a.inversa()
            self.canonica = b.transposada() * matriu * b
            self.canonica.simplificar()
    #
    #
    #
    @classmethod
    def from_equacio(cls,eq):
        """
        Retorna i classifica la cònica a partir de la seva equació.
        Només per a el·lipses, hipèrboles i paràboles
        """
        x, y = symbols('x y')
        unknowns = [x,y]
        a = diff(eq,x,2) / 2
        b = diff(eq,x,y) / 2
        c = diff(eq,y,2) / 2
        d = (diff(eq,x).subs({x:1,y:0}) - 2 * a) / 2
        e = (diff(eq,y).subs({x:0,y:1}) - 2 * c) / 2
        f = eq.subs({x:0,y:0})
        Q = Matrix([[a,b],[b,c]])
        L = Matrix(2,1,[d,e])
        l = Vector(d,e)
        vs = Q.eigenvects()
        (t1,t2), veps = vaps_veps(vs)
        s = list(linsolve((Q,-L),*unknowns))
        if t1 * t2 < 0:
            #
            # Hipèrbola
            #
            s = Punt(list(s[0]))
            fp = eq.subs({x:s[0],y:s[1]})
            a2 = -fp/t1
            b2 = -fp/t2
            if a2 > 0:
                u = veps[0]
            else:
                a2, b2 = b2, a2
                u = veps[1]
            return Hiperbola(a2,-b2,s,u)
        if t1 * t2 > 0:
            #
            # El·lipse
            #
            s = Punt(list(s[0]))
            fp = eq.subs({x:s[0],y:s[1]})
            if fp == 0:
                return None
            a2 = -fp/t1
            b2 = -fp/t2
            if a2 < 0:
                return None
            if a2 >= b2:
                u = veps[0]
            else:
                a2, b2 = b2, a2
                u = veps[1]
            return Ellipse(a2,b2,s,u)
        if len(s) > 0:
            return None
        if t1 == 0 and t2 == 0:
            return None
        #
        # Paràbola
        #
        if t1 == 0:
            t1, t2 = t2, t1
            veps.reverse()
        veps[0].simplificar()
        veps[1].simplificar()
        b = Base(veps,unitaria=True)
        b.orientacio_positiva()
        es = t1 * veps[0].dot(Vector(x,y)) + veps[0].dot(l)
        vertex = Punt(solve([es,eq],x,y)[0])
        ref = ReferenciaAfi(vertex,b)
        ep = veps[1].dot(l)/veps[1].length()
        p = - ep/t1
        focus = ref.punt_de_coordenades(Punt(0,p/2))
        return Parabola(vertex,focus)
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'equació en latex de l'equació de la cònica en la referència
        canònica
        """
        x, y = symbols('x y')
        m = Matriu.matriu_columna(Vector([x,y,1]))
        r = m.transposada() * self.canonica * m
        return mylatex(r[0,0].expand()) + " = 0"
        #
        #
        #
    def equacio(self):
        """
        Retorna l'equació en latex de l'equació de la quàdrica
        """
        x, y = symbols('x y')
        m = Matriu.matriu_columna(Vector([x,y,1]))
        r = m.transposada() * self.canonica * m
        return r[0,0].expand()
    #
    #
    #
    @classmethod
    def ellipse(cls,maxim=30):
        """
        Retorna una el·lipse aleatòria
        Paràmetres:
            maxim: nombre màxim de la matriu projectiva de la cònica
        """
        trobat = False
        while not trobat:
            e = Ellipse.aleatoria()
            trobat = e.canonica.norma_maxim() <= maxim
        return e
    #
    #
    #
    @classmethod
    def hiperbola(cls,maxim=30):
        """
        Retorna una hipèrbola aleatòria
        Paràmetres:
            maxim: nombre màxim de la matriu projectiva de la cònica
        """
        trobat = False
        while not trobat:
            h = Hiperbola.aleatoria()
            trobat = h.canonica.norma_maxim() <= maxim
        return h
    #
    #
    #
    @classmethod
    def parabola(cls,maxim=30):
        """
        Retorna una paràbola aleatòria
        Paràmetres:
            maxim: nombre màxim de la matriu projectiva de la cònica
        """
        trobat = False
        while not trobat:
            p = Parabola.aleatoria()
            trobat = p.canonica.norma_maxim() <= maxim
        return p
    #
    #
    #
    @classmethod
    def aleatoria(cls,maxim=30):
        """
        Retorna una el·lipse, hipèrbola o paràbola aleatòries
        Paràmetres:
            maxim: nombre màxim de la matriu projectiva de la cònica
        """
        r = random.randint(0,2)
        if r == 0:
            return Conica.ellipse(maxim)
        if r == 1:
            return Conica.hiperbola(maxim)
        return Conica.parabola(maxim)
    #
    #
    #
    def tipus(self):
        """
        Retorna el tipus de cònica
        """
        if isinstance(self,Ellipse):
            return "El·lipse"
        if isinstance(self,Hiperbola):
            return "Hipèrbola"
        if isinstance(self,Parabola):
            return "Paràbola"
        return ""
    #
    #
    #
    def referencia_principal(self):
        """
        Retorna la referencia principal
        """
        return self.ref
    #
    #
    #
    def vectors(self,unitaris=False):
        """
        Retorna els vectors de la base de la referència principal
        Paràmetres:
            unitaris: si és True, els retorna unitaris
        """
        return self.ref.vectors(unitaris)

class Ellipse(Conica):
    """
    Classe per treballar amb el·lipses
    """
    #
    #
    #
    def __new__(cls,a2,b2,centre,eix):
        """
        Constructor.
        Retorna una el·lipse
        Paràmetres:
            a2: semieix major al quadrat
            b2: semieix menor al quadrat
            centre: centre de l'el·lipse
            eix: direcció de l'eix principal (de les x')
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix,Vector):
            return None
        if centre.dimensio != 2:
            return None
        if eix.dimensio != 2:
            return None
        if eix.length() == 0:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        if a2 < b2:
            return None
        return super(Conica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,centre,eix):
        s = SubespaiVectorial([eix])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = gcd(a2,b2)
        t = a2 * b2 // g
        a2 = a2 // g
        b2 = b2 // g
        m = Matriu.diagonal(Vector([b2,a2,-t]))
        Conica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna una el·lipse aleatòria
        """
        eix = Vector.aleatori(l=2,maxim=3,nuls=False)
        centre = Punt.aleatori(l=2,maxim=3,nuls=False)
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36,40,45,48,50,60,80,100]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            trobat = a2 != b2
        if b2 > a2:
            a2, b2 = b2, a2
        return cls(a2,b2,centre,eix)
    #
    #
    #
    def centre(self):
        """
        Retorna el centre de la el·lipse
        """
        return (self.ref.origen)
    #
    #
    #
    def semieix_major(self):
        """
        Retorna el semieix major
        """
        l1 = self.matriu[0,0]
        f = - self.matriu[2,2]
        return sqrt(Rational(f,l1))
    #
    #
    #
    def semieix_menor(self):
        """
        Retorna el semieix menor
        """
        l2 = self.matriu[1,1]
        f = - self.matriu[2,2]
        return sqrt(Rational(f,l2))
    #
    #
    #
    def semidistancia_focal(self):
        """
        Retorna la simidistància focal
        """
        a2 = self.semieix_major()**2
        b2 = self.semieix_menor()**2
        return sqrt(a2 - b2)
    #
    #
    #
    def focus(self):
        """
        Retorna els focus de l'el·lipse
        """
        c = self.semidistancia_focal()
        f = [Punt(c,0),Punt(-c,0)]
        return list(map(self.ref.punt_de_coordenades,f))
    #
    #
    #
    def vertexs(self):
        """
        Retorna els vèrtex de l'el·lipse
        """
        a = self.semieix_major()
        b = self.semieix_menor()
        v = [Punt(a,0),Punt(-a,0),Punt(0,b),Punt(0,-b)]
        return list(map(self.ref.punt_de_coordenades,v))
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda de l'el·lipse en format LaTeX
        """
        a2 = self.semieix_major()**2
        b2 = self.semieix_menor()**2
        if a2 == 1:
            return f"x'^2 + \\frac{{y'^2}}{{{b2}}} = 1"
        if b2 == 1:
            return f"\\frac{{x'^2}}{{{a2}}} + y'^2 = 1"
        return f"\\frac{{x'^2}}{{{a2}}} + \\frac{{y'^2}}{{{b2}}} = 1"
    #
    #
    #
    def to_asy(self,x=13,y=10):
        """
        Retorna una expressió per fer servir amb el programa Asymtote
        Paràmetres:
            x, y: nombres enters. El gràfic es representarà en una quadricula
            de límits (-x,x) i (-y,y)
        """
        a2 = self.semieix_major()**2
        b2 = self.semieix_menor()**2
        centre = self.centre()
        vector = self.vectors()[0]
        return f"Ellipse({centre},{vector},{a2},{b2},x={x},y={y})"

class Hiperbola(Conica):
    """
    Classe per treballar amb hipèrboles
    """
    #
    #
    #
    def __new__(cls,a2,b2,centre,eix):
        """
        Constructor.
        Retorna una hipèrbola
        Paràmetres:
            a2: semieix real al quadrat
            b2: semieix imaginari al quadrat
            centre: centre de la hipèrbola
            eix: direcció de l'eix principal (de les x')
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix,Vector):
            return None
        if centre.dimensio != 2:
            return None
        if eix.dimensio != 2:
            return None
        if eix.length() == 0:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        return super(Conica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,centre,eix):
        s = SubespaiVectorial([eix])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = gcd(a2,b2)
        t = a2 * b2 // g
        a2 = a2 // g
        b2 = b2 // g
        m = Matriu.diagonal(Vector([b2,-a2,-t]))
        Conica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna una hipèrbola aleatòria
        """
        eix = Vector.aleatori(l=2,maxim=3,nuls=False)
        centre = Punt.aleatori(l=2,maxim=3,nuls=False)
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36,40,45,48,50,60,80,100]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            trobat = a2 != b2
        return cls(a2,b2,centre,eix)
    #
    #
    #
    def centre(self):
        """
        Retorna el centre de la el·lipse
        """
        return (self.ref.origen)
    #
    #
    #
    def semieix_real(self):
        """
        Retorna el semieix real
        """
        l1 = self.matriu[0,0]
        f = - self.matriu[2,2]
        return sqrt(Rational(f,l1))
    #
    #
    #
    def semieix_imaginari(self):
        """
        Retorna el semieix imaginari
        """
        l2 = - self.matriu[1,1]
        f = - self.matriu[2,2]
        return sqrt(Rational(f,l2))
    #
    #
    #
    def semidistancia_focal(self):
        """
        Retorna la simidistància focal
        """
        a2 = self.semieix_real()**2
        b2 = self.semieix_imaginari()**2
        return sqrt(a2 + b2)
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda de la hipèrbola en format LaTeX
        """
        a2 = self.semieix_real()**2
        b2 = self.semieix_imaginari()**2
        if a2 == 1:
            return f"x'^2 - \\frac{{y'^2}}{{{b2}}} = 1"
        if b2 == 1:
            return f"\\frac{{x'^2}}{{{a2}}} - y'^2 = 1"
        return f"\\frac{{x'^2}}{{{a2}}} - \\frac{{y'^2}}{{{b2}}} = 1"
    #
    #
    #
    def focus(self):
        """
        Retorna els focus de la hipèrbola
        """
        c = self.semidistancia_focal()
        f = [Punt(c,0),Punt(-c,0)]
        return list(map(self.ref.punt_de_coordenades,f))
    #
    #
    #
    def vertexs(self):
        """
        Retorna els vèrtexs de la hipèrbola
        """
        a = self.semieix_real()
        v = [Punt(a,0),Punt(-a,0)]
        return list(map(self.ref.punt_de_coordenades,v))
    #
    #
    #
    def vectors_directors_asimptotes(self):
        """
        Retorna els vectors directors de les asímptotes expressats en la
        base canònica
        """
        l = self.ref.base.vecs[0].length()
        a = l * self.semieix_real()
        b = l * self.semieix_imaginari()
        v1 = self.ref.base.vector_de_components(Vector([a,b]))
        v2 = self.ref.base.vector_de_components(Vector([a,-b]))
        v1.radsimplificar()
        v2.radsimplificar()
        return (v1,v2)
    #
    #
    #
    def equacio_continua_asimptota(self,v):
        """
        Retorna l'equació contínua de l'asímptota amb direccio v en format LaTeX
        Paràmetres:
            v: vector director de l'asímptota
        """
        x, y = symbols('x y')
        centre = self.centre()
        if v[0] == 1:
            str = latex(x - centre[0])
        else:
            str = f"\\frac{{ {latex(x - centre[0])} }}{{ {v[0]} }}"
        str += " = "
        if v[1] == 1:
            str += latex(y - centre[1])
        else:
            str += f"\\frac{{ {latex(y - centre[1])} }}{{ {v[1]} }}"
        return str
    #
    #
    #
    def to_asy(self,x=10,y=10):
        """
        Retorna una expressió per fer servir amb el programa Asymtote
        Paràmetres:
            x, y: nombres enters. El gràfic es representarà en una quadricula
            de límits (-x,x) i (-y,y)
        """
        a2 = self.semieix_real()**2
        b2 = self.semieix_imaginari()**2
        centre = self.centre()
        vector = self.vectors()[0]
        return f"Hiperbola({centre},{vector},{a2},{b2},x={x},y={y})"

class Parabola(Conica):
    """
    Classe per treballar amb paràboles
    """
    #
    #
    #
    def __new__(cls,vertex,focus):
        """
        Constructor.
        Retorna una paràbola
        Paràmetres:
            vertex: vèrtex
            focus: focus de la paràbola
        """
        if not isinstance(vertex,Punt):
            return None
        if not isinstance(focus,Punt):
            return None
        if vertex.dimensio != 2:
            return None
        if focus.dimensio != 2:
            return None
        if vertex == focus:
            return None
        return super(Conica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,vertex,focus):
        eix = focus - vertex
        p = 2 * eix.length()
        m = Matriu(Matrix(3,3,[1,0,0,0,0,-p,0,-p,0]))
        s = SubespaiVectorial([eix])
        base = s.amplia_base_suplementari(unitaria=True)
        r = ReferenciaAfi(vertex,base)
        Conica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna una paràbola aleatòria
        """
        trobat = False
        while not trobat:
            vertex = Punt.aleatori(l=2,maxim=5,nuls=False)
            focus = Punt.aleatori(l=2,maxim=5)
            trobat = vertex[0] != focus[0] and vertex[1] != focus[1]
        return cls(vertex,focus)
    #
    #
    #
    def parametre(self):
        """
        Retorna el paràmetre de la paràbola
        """
        return - self.matriu[2,1]
    #
    #
    #
    def vertex(self):
        """
        Retorna el vèrtex de la paràbola
        """
        return (self.ref.origen)
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda de la paràbola en format LaTeX
        """
        p = self.parametre()
        if 2 * p == 1:
            return f"y' = x'^2"
        return f"y' = \\frac{{x'^2}}{{{ latex(2 * p) }}}"
    #
    #
    #
    def focus(self):
        """
        Retorna el focus de la paràbola
        """
        p2 = self.parametre() / 2
        return self.ref.punt_de_coordenades(Punt([0,p2]))
    #
    #
    #
    def to_asy(self,x=10,y=10):
        """
        Retorna una expressió per fer servir amb el programa Asymtote
        Paràmetres:
            x, y: nombres enters. El gràfic es representarà en una quadricula
            de límits (-x,x) i (-y,y)
        """
        return f"Parabola({self.vertex()},{self.focus()},x={x},y={y})"

class Quadrica(object):
    """
    Classe per treballar amb quàdriques. L'objectiu no és classificar quàdriques,
    sinó generar-les a partir dels elements característics o de manera
    aleatòria.
    Atributs:
        ref: Referència afí
        matriu: matriu projectiva de la quàdrica en la referència "ref"
        canonica: matriu projectiva de la quàdrica en la referència canònica
    """
    #
    #
    #
    def __new__(cls,matriu,ref=None):
        """
        Retorna una nova quàdrica
        Paràmetres:
            matriu: matriu projectiva de la quàdrica en la referència "ref"
            Ha de ser 4x4
            ref: referència afí. Si és None, serà la referència canònica
        """
        if not isinstance(matriu,Matriu):
            return None
        if matriu.files != matriu.columnes:
            return None
        if matriu.files != 4:
            return None
        if not matriu.es_simetrica():
            return None
        if ref is not None:
            if ref.dimensio != 3:
                return None
            if not ref.base.es_unitaria():
                return None
            if not ref.base.es_ortogonal():
                return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,matriu,ref=None):
        self.ref = ref
        self.matriu = matriu
        if ref is None:
            self.canonica = matriu
        else:
            a = ref.base.matriu()
            a = a.inserta_columna(3,ref.origen)
            a = a.inserta_fila(3,Vector([0,0,0,1]))
            b = a.inversa()
            self.canonica = b.transposada() * self.matriu * b
        self.canonica.simplificar()
    #
    #
    #
    def __repr__(self):
        """
        Retorna l'equació en latex de l'equació de la quàdrica en la referència
        canònica en LaTeX
        """
        x, y, z = symbols('x y z')
        m = Matriu.matriu_columna(Vector([x,y,z,1]))
        r = m.transposada() * self.canonica * m
        return mylatex(r[0,0].expand()) + " = 0"
    #
    #
    #
    def equacio(self):
        """
        Retorna l'equació en latex de l'equació de la quàdrica
        """
        x, y, z = symbols('x y z')
        m = Matriu.matriu_columna(Vector([x,y,z,1]))
        r = m.transposada() * self.canonica * m
        return r[0,0].expand()

    #
    #
    #
    @classmethod
    def ellipsoide(cls,maxim=30,diagonal=15):
        """
        Retorna un el·lipsoide de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            e = Ellipsoide.aleatoria()
            trobat = e.canonica.norma_maxim() <= maxim and e.canonica.nzeros() < 3 and e.canonica.max_diagonal() < diagonal
        return e
    #
    #
    #
    @classmethod
    def hiperboloideunafulla(cls,maxim=30,diagonal=15):
        """
        Retorna un hiperboloide d'una fulla de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            h = HiperboloideUnaFulla.aleatoria()
            trobat = h.canonica.norma_maxim() <= maxim and h.canonica.nzeros() < 3 and h.canonica.max_diagonal() < diagonal
        return h
    #
    #
    #
    @classmethod
    def hiperboloideduesfulles(cls,maxim=30,diagonal=15):
        """
        Retorna un hiperboloide de dues fulles de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            h = HiperboloideDuesFulles.aleatoria()
            trobat = h.canonica.norma_maxim() <= maxim and h.canonica.nzeros() < 3 and h.canonica.max_diagonal() < diagonal
        return h
    #
    #
    #
    @classmethod
    def con(cls,maxim=30,diagonal=15):
        """
        Retorna un con de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            c = Con.aleatoria()
            trobat = c.canonica.norma_maxim() <= maxim and c.canonica.nzeros() < 3 and c.canonica.max_diagonal() < diagonal
        return c
    #
    #
    #
    @classmethod
    def cilindreelliptic(cls,maxim=30,diagonal=15):
        """
        Retorna un cilindre el·líptic de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            c = CilindreElliptic.aleatoria()
            trobat = c.canonica.norma_maxim() <= maxim and c.canonica.nzeros() < 3 and c.canonica.max_diagonal() < diagonal
        return c
    #
    #
    #
    @classmethod
    def cilindrehiperbolic(cls,maxim=30,diagonal=15):
        """
        Retorna un cilindre hiperbòlic de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            c = CilindreHiperbolic.aleatoria()
            trobat = c.canonica.norma_maxim() <= maxim and c.canonica.nzeros() < 3 and c.canonica.max_diagonal() < diagonal
        return c
    #
    #
    #
    @classmethod
    def paraboloideelliptic(cls,maxim=30,diagonal=15):
        """
        Retorna un paraboloide el·líptic de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            p = ParaboloideElliptic.aleatoria()
            trobat = p.canonica.norma_maxim() <= maxim and p.canonica.nzeros() < 3 and p.canonica.max_diagonal() < diagonal
        return p
    #
    #
    #
    @classmethod
    def paraboloidehiperbolic(cls,maxim=30,diagonal=15):
        """
        Retorna un paraboloide hiperbòlic de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            p = ParaboloideHiperbolic.aleatoria()
            trobat = p.canonica.norma_maxim() <= maxim and p.canonica.nzeros() < 3 and p.canonica.max_diagonal() < diagonal
        return p
    #
    #
    #
    @classmethod
    def cilindreparabolic(cls,maxim=30,diagonal=15):
        """
        Retorna un cilindre parabòlic de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        trobat = False
        while not trobat:
            c = CilindreParabolic.aleatoria()
            trobat = c.canonica.norma_maxim() <= maxim and c.canonica.nzeros() < 3 and c.canonica.max_diagonal() < diagonal
        return c
    #
    #
    #
    @classmethod
    def aleatoria(cls,maxim=30,diagonal=15):
        """
        Retorna una quàdrica de manera aleatòria
        Paràmetres:
            maxim: valor màxim de la matriu projectiva de la quàdrica
            diagonal: valor màxim de la diagonal de la matriu projectiva de la quàdrica
        """
        r = random.randint(0,9)
        if r == 0:
            return Quadrica.ellipsoide(maxim,diagonal)
        if r == 1:
            return Quadrica.hiperboloideunafulla(maxim,diagonal)
        if r == 2:
            return Quadrica.hiperboloideduesfulles(maxim,diagonal)
        if r == 3:
            return Quadrica.con(maxim,diagonal)
        if r == 4:
            return Quadrica.cilindreelliptic(maxim,diagonal)
        if r == 5:
            return Quadrica.cilindrehiperbolic(maxim,diagonal)
        if r == 6:
            return Quadrica.paraboloideelliptic(maxim,diagonal)
        if r == 7:
            return Quadrica.paraboloidehiperbolic(maxim,diagonal)
        return Quadrica.cilindreparabolic(maxim,diagonal)
    #
    #
    #
    def tipus(self):
        """
        Retorna el tipus de quàdrica
        """
        if isinstance(self,Ellipsoide):
            return "El·lipsoide"
        if isinstance(self,HiperboloideUnaFulla):
            return "Hiperboloide d'una fulla"
        if isinstance(self,HiperboloideDuesFulles):
            return "Hiperboloide de dues fulles"
        if isinstance(self,Con):
            return "Con"
        if isinstance(self,CilindreElliptic):
            return "Cilindre el·líptic"
        if isinstance(self,CilindreHiperbolic):
            return "Cilindre hiperbòlic"
        if isinstance(self,ParaboloideElliptic):
            return "Paraboloide el·líptic"
        if isinstance(self,ParaboloideHiperbolic):
            return "Paraboloide hiperbòlic"
        if isinstance(self,CilindreParabolic):
            return "Cilindre parabòlic"
        return ""
    #
    #
    #
    def referencia_principal(self):
        """
        Retorna la referencia principal
        """
        return self.ref
    #
    #
    #
    def vectors(self,unitaris=False):
        """
        Retorna els vectors de la base de la referència principal
        Paràmetres:
            unitaris: si és True, els retorna unitaris
        """
        return self.ref.vectors(unitaris)
    #
    #
    #
    @classmethod
    def from_equacio(cls,eq):
        """
        Retorna i classifica la quadrica a partir de la seva equació.
        Només per als nou tipus de quàdriques definides en aquesta llibreria
        """
        x, y, z = symbols('x y z')
        unknowns = [x,y,z]
        #
        # Termes quadràtics
        #
        a = diff(eq,x,2) / 2
        b = diff(eq,x,y) / 2
        d = diff(eq,x,z) / 2
        c = diff(eq,y,2) / 2
        e = diff(eq,y,z) / 2
        f = diff(eq,z,2) / 2
        #
        # Termes lineals
        #
        l1 = (diff(eq,x).subs({x:1,y:0,z:0}) - 2 * a) / 2
        l2 = (diff(eq,y).subs({x:0,y:1,z:0}) - 2 * c) / 2
        l3 = (diff(eq,z).subs({x:0,y:0,z:1}) - 2 * f) / 2
        #
        # Terme independent
        #
        ti = eq.subs({x:0,y:0,z:0})
        #
        # Matrius
        #
        Q = Matrix([[a,b,d],[b,c,e],[d,e,f]])
        L = Matrix(3,1,[l1,l2,l3])
        l = Vector(l1,l2,l3)
        #
        # Preparatius
        #
        vs = Q.eigenvects()
        vaps, veps = vaps_veps_amb_signe(vs,1)
        positius = list(zip(vaps,veps))
        vaps, veps = vaps_veps_amb_signe(vs,-1)
        negatius = list(zip(vaps,veps))
        positius.sort(key=lambda item: abs(item[0]))
        negatius.sort(key=lambda item: abs(item[0]))
        s = list(linsolve((Q,-L),*unknowns))
        #
        # Classificació
        #
        if len(positius) == 3 or len(negatius) == 3:
            #
            # Ellipsoide
            #
            s = Punt(list(s[0]))
            tip = eq.subs({x:s[0],y:s[1],z:s[2]})
            if tip == 0:
                return None
            if len(positius) == 3 and tip > 0:
                return None
            if len(negatius) == 3 and tip < 0:
                return None
            v = positius
            if len(v) == 0:
                v = negatius
            a2 = -tip / v[0][0]
            b2 = -tip / v[1][0]
            c2 = -tip / v[2][0]
            return Ellipsoide(a2,b2,c2,s,v[0][1],v[1][1])
        if len(positius) + len(negatius) == 3:
            #
            # Rang 3
            #
            s = Punt(list(s[0]))
            tip = eq.subs({x:s[0],y:s[1],z:s[2]})
            v1 = positius
            v2 = negatius
            signe = 1
            if len(v2) == 2:
                v1, v2 = v2, v1
                signe = -1
            a2 = signe * v1[0][0]
            b2 = signe * v1[1][0]
            c2 = signe * v2[0][0]
            tip *= signe
            if tip == 0:
                #
                # Con
                #
                if a2.is_integer and b2.is_integer and c2.is_integer:
                    g = mcm_llista([a2,b2,c2])
                else:
                    g = a2 * b2 * c2
                a2 = g / a2
                b2 = g / b2
                c2 = g / c2
                return Con(a2,b2,-c2,s,v1[0][1],v1[1][1])
            elif tip > 0:
                #
                # Hiperboloide de dues fulles
                #
                a2 = tip / a2
                b2 = tip / b2
                c2 = tip / c2
                return HiperboloideDuesFulles(a2,b2,-c2,s,v1[0][1],v1[1][1])
            elif tip < 0:
                #
                # Hiperboloide d'una fulla
                #
                a2 = - tip / a2
                b2 = - tip / b2
                c2 = - tip / c2
                return HiperboloideUnaFulla(a2,b2,-c2,s,v1[0][1],v1[1][1])
        if len(positius) + len(negatius) == 2:
            #
            # Rang 2
            #
            if len(s) > 0:
                s = Punt(list(s[0].subs({x:0,y:0,z:0})))
                tip = eq.subs({x:s[0],y:s[1],z:s[2]})
                if len(positius) == 2 or len(negatius) == 2:
                    #
                    # Clindre e·líptic
                    #
                    if tip == 0:
                        return None
                    if len(positius) == 2 and tip > 0:
                        return None
                    if len(negatius) == 2 and tip < 0:
                        return None
                    v = positius
                    if len(v) == 0:
                        v = negatius
                    a2 = -tip / v[0][0]
                    b2 = -tip / v[1][0]
                    return CilindreElliptic(a2,b2,s,v[0][1],v[1][1])
                if len(positius) == 1:
                    #
                    # Clindre hiperbòlic
                    #
                    if tip == 0:
                        return None
                    v1 = positius
                    v2 = negatius
                    a2 = -tip / positius[0][0]
                    b2 = -tip / negatius[0][0]
                    if a2 < 0:
                        a2, b2 = b2, a2
                        v1, v2 = v2, v1
                    return CilindreHiperbolic(a2,-b2,s,v1[0][1],v2[0][1])
            if len(positius) == 2 or len(negatius) == 2:
                #
                # Paraboloide el·líptic
                #
                v = positius
                if len(v) == 0:
                    v = negatius
                v1, v2 = v[0][1],v[1][1]
                v1.simplificar()
                v2.simplificar()
                v3 = v1.cross(v2,simplificar=True)
                ep = v3.dot(l)/v3.length()
                t1, t2 = v[0][0],v[1][0]
                es1 = t1 * v1.dot(Vector(x,y,z)) + v1.dot(l)
                es2 = t2 * v2.dot(Vector(x,y,z)) + v2.dot(l)
                vertex = Punt(solve([es1,es2,eq],x,y,z)[0])
                a2 =  -2 * ep / t1
                b2 =  -2 * ep / t2
                if a2 < 0:
                    a2, b2 = -a2,-b2
                    v2 = -v2
                a2 /= v3.length()
                b2 /= v3.length()
                return ParaboloideElliptic(a2,b2,vertex,v1,v2)
            if len(positius) == 1:
                #
                # Paraboloide hiperbòlic
                #
                v1 = positius
                v2 = negatius
                vec1, vec2 = v1[0][1],v2[0][1]
                vec1.simplificar()
                vec2.simplificar()
                vec3 = vec1.cross(vec2,simplificar=True)
                ep = vec3.dot(l)/vec3.length()
                t1, t2 = v1[0][0],v2[0][0]
                es1 = t1 * vec1.dot(Vector(x,y,z)) + vec1.dot(l)
                es2 = t2 * vec2.dot(Vector(x,y,z)) + vec2.dot(l)
                vertex = Punt(solve([es1,es2,eq],x,y,z)[0])
                a2 =  -2 * ep / t1
                b2 =  -2 * ep / t2
                if a2 < 0:
                    a2, b2 = b2, a2
                    vec1, vec2 = vec2, vec1
                    vec3 *= -1
                a2 /= vec3.length()
                b2 /= vec3.length()
                return ParaboloideHiperbolic(a2,-b2,vertex,vec1,vec2)
        if len(positius) + len(negatius) == 1:
            #
            # Cilindre parabòlic
            #
            if len(s) > 0:
                return None
            v = positius
            if len(v) == 0:
                v = negatius
            t1 = v[0][0]
            v1 = v[0][1]
            v2 = v1.cross(l,simplificar=True)
            v3 = v1.cross(v2,simplificar=True)
            es = t1 * v1.dot(Vector(x,y,z)) + v1.dot(l)
            vertex = solve([es,eq],x,y,z)[0]
            k = 0
            while True:
                v = Punt([item.subs({x:k,y:k,z:k}) for item in vertex])
                f = v.factor_comu()
                if f[0].q == 1:
                    break
                v = Punt([item.subs({x:-k,y:-k,z:-k}) for item in vertex])
                f = v.factor_comu()
                if f[0].q == 1:
                    break
                k += 1
            ep = v3.dot(l)/v3.length()
            a2 =  -2 * ep / t1
            return CilindreParabolic(a2/2,v,v1,v2)

class Ellipsoide(Quadrica):
    """
    Classe per treballar amb el·lipsoides
    """
    #
    #
    #
    def __new__(cls,a2,b2,c2,centre,eix1,eix2):
        """
        Constructor.
        Paràmetres:
           a2, b2, c2: semieixos al quadrat
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if centre.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        if c2 <= 0:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,c2,centre,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = mcm_llista([a2,b2,c2])
        t = g
        a2 = g // a2
        b2 = g // b2
        c2 = g // c2
        m = Matriu.diagonal(Vector([a2,b2,c2,-g]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un el·lipsoide de manera aleatòria aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            centre = Punt.aleatori(l=3,maxim=4)
            trobat = centre.length() > 0
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            a = random.randint(0,len(c) - 1)
            c2 = c[a]
            trobat = a2 != b2 or a2 != c2
        a2, b2, c2 = sorted([a2,b2,c2])[::-1]
        return cls(a2,b2,c2,centre,eix1,eix2)
    #
    #
    #
    def centre(self):
        """
        Retorna el centre de l'el·lipsoide
        """
        return (self.ref.origen)
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos de l'el·lipsoide
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = self.matriu[2,2]
        f = - self.matriu[3,3]
        return (sqrt(Rational(f,a)),sqrt(Rational(f,b)),sqrt(Rational(f,c)))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat de l'el·lipsoide
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = self.matriu[2,2]
        f = - self.matriu[3,3]
        return (Rational(f,a),Rational(f,b),Rational(f,c))
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda de l'el·lipsoide
        """
        a2, b2, c2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "x'^2"
        else:
            str = f"\\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " + y'^2"
        else:
            str += f" + \\frac{{y'^2}}{{ {latex(b2)} }}"
        if c2 == 1:
            str += " + z'^2 = 1"
        else:
            str += f" + \\frac{{z'^2}}{{ {latex(c2)} }} = 1"
        return str

class HiperboloideUnaFulla(Quadrica):
    """
    Classe per treballar amb hiperboloides d'una fulla
    """
    #
    #
    #
    def __new__(cls,a2,b2,c2,centre,eix1,eix2):
        """
        Constructor.
        Paràmetres:
           a2, b2, c2: semieixos al quadrat
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if centre.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        if c2 <= 0:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,c2,centre,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = mcm_llista([a2,b2,c2])
        t = g
        a2 = g // a2
        b2 = g // b2
        c2 = g // c2
        m = Matriu.diagonal(Vector([a2,b2,-c2,-g]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un hiperboloide d'una fulla de manera aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            centre = Punt.aleatori(l=3,maxim=4)
            trobat = centre.length() > 0
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            a = random.randint(0,len(c) - 1)
            c2 = c[a]
            trobat = a2 != b2 or a2 != c2
        a2, b2 = sorted([a2,b2])[::-1]
        return cls(a2,b2,c2,centre,eix1,eix2)
    #
    #
    #
    def centre(self):
        """
        Retorna el centre de l'el·lipsoide
        """
        return (self.ref.origen)
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos de l'el·lipsoide
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = - self.matriu[2,2]
        f = - self.matriu[3,3]
        return (sqrt(Rational(f,a)),sqrt(Rational(f,b)),sqrt(Rational(f,c)))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat de l'el·lipsoide
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = - self.matriu[2,2]
        f = - self.matriu[3,3]
        return (Rational(f,a),Rational(f,b),Rational(f,c))
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda de l'hiperboloide d'una fulla
        """
        a2, b2, c2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "x'^2"
        else:
            str = f"\\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " + y'^2"
        else:
            str += f" + \\frac{{y'^2}}{{ {latex(b2)} }}"
        if c2 == 1:
            str += " - z'^2 = 1"
        else:
            str += f" - \\frac{{z'^2}}{{ {latex(c2)} }} = 1"
        return str

class HiperboloideDuesFulles(Quadrica):
    """
    Classe per treballar amb hiperboloides de dues fulles
    """
    #
    #
    #
    def __new__(cls,a2,b2,c2,centre,eix1,eix2):
        """
        Constructor.
        Paràmetres:
           a2, b2, c2: semieixos al quadrat
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if centre.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        if c2 <= 0:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,c2,centre,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = mcm_llista([a2,b2,c2])
        t = g
        a2 = g // a2
        b2 = g // b2
        c2 = g // c2
        m = Matriu.diagonal(Vector([a2,b2,-c2,g]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un hiperboloide de dues fulles de manera aleatòria aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            centre = Punt.aleatori(l=3,maxim=4)
            trobat = centre.length() > 0
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            a = random.randint(0,len(c) - 1)
            c2 = c[a]
            trobat = a2 != b2 or a2 != c2
        a2, b2 = sorted([a2,b2])[::-1]
        return cls(a2,b2,c2,centre,eix1,eix2)
    #
    #
    #
    def centre(self):
        """
        Retorna el centre de l'hiperboloide de dues fulles
        """
        return (self.ref.origen)
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos de l'hiperboloide de dues fulles
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = - self.matriu[2,2]
        f = self.matriu[3,3]
        return (sqrt(Rational(f,a)),sqrt(Rational(f,b)),sqrt(Rational(f,c)))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat de l'hiperboloide de dues fulles
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = - self.matriu[2,2]
        f = self.matriu[3,3]
        return (Rational(f,a),Rational(f,b),Rational(f,c))
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda de l'hiperboloide de dues fulles
        """

        a2, b2, c2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "x'^2"
        else:
            str = f"\\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " + y'^2"
        else:
            str += f" + \\frac{{y'^2}}{{ {latex(b2)} }}"
        if c2 == 1:
            str += " - z'^2 = -1"
        else:
            str += f" - \\frac{{z'^2}}{{ {latex(c2)} }} = -1"
        return str

class Con(Quadrica):
    """
    Classe per treballar amb cons
    """
    #
    #
    #
    def __new__(cls,a2,b2,c2,centre,eix1,eix2):
        """
        Constructor.
        Retorna un Con
        Paràmetres:
           a2, b2, c2: semieixos al quadrat
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if centre.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        if c2 <= 0:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,c2,centre,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = mcm_llista([a2,b2,c2])
        a2 = g // a2
        b2 = g // b2
        c2 = g // c2
        m = Matriu.diagonal(Vector([a2,b2,-c2,0]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un hiperboloide de dues fulles de manera aleatòria aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            centre = Punt.aleatori(l=3,maxim=4)
            trobat = centre.length() > 0
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            a = random.randint(0,len(c) - 1)
            c2 = c[a]
            trobat = a2 != b2 or a2 != c2
        a2, b2 = sorted([a2,b2])[::-1]
        return cls(a2,b2,c2,centre,eix1,eix2)
    #
    #
    #
    def centre(self):
        """
        Retorna el centre o vèrtex del con
        """
        return (self.ref.origen)
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos del con
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = - self.matriu[2,2]
        m = mcm_llista([a,b,c])
        return (sqrt(Rational(m,a)),sqrt(Rational(m,b)),sqrt(Rational(m,c)))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat del con
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        c = - self.matriu[2,2]
        m = mcm_llista([a,b,c])
        return (Rational(m,a),Rational(m,b),Rational(m,c))
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda del con
        """
        a2, b2, c2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "x'^2"
        else:
            str = f"\\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " + y'^2"
        else:
            str += f" + \\frac{{y'^2}}{{ {latex(b2)} }}"
        if c2 == 1:
            str += " - z'^2 = 0"
        else:
            str += f" - \\frac{{z'^2}}{{ {latex(c2)} }} = 0"
        return str

class CilindreElliptic(Quadrica):
    """
    Classe per treballar amb cilindres el·líptics
    """
    #
    #
    #
    def __new__(cls,a2,b2,centre,eix1,eix2):
        """
        Constructor.
        Paràmetres:
           a2, b2: semieixos al quadrat
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if centre.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,centre,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = mcd_llista([a2,b2])
        t = a2 * b2 // g
        a2 = a2 // g
        b2 = b2 // g
        m = Matriu.diagonal(Vector([b2,a2,0,-t]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un cilindre el·líptic de manera aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            centre = Punt.aleatori(l=3,maxim=4)
            trobat = centre.length() > 0
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            trobat = a2 != b2
        a2, b2 = sorted([a2,b2])[::-1]
        return cls(a2,b2,centre,eix1,eix2)
    #
    #
    #
    def centre(self):
        """
        Retorna un centre del cilindre el·líptic
        """
        return (self.ref.origen)
    #
    #
    #
    def centres(self):
        """
        Retorna la recta de centres
        """
        return RectaAfi(self.ref.origen,self.ref.base[2])
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos del cilindre el·líptic
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        f = - self.matriu[3,3]
        return (sqrt(Rational(f,a)),sqrt(Rational(f,b)))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat del cilindre el·líptic
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        f = - self.matriu[3,3]
        return (Rational(f,a),Rational(f,b))
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda del cilindre el·líptic
        """
        a2, b2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "x'^2"
        else:
            str = f"\\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " + y'^2 = 1"
        else:
            str += f" + \\frac{{y'^2}}{{ {latex(b2)} }} = 1"
        return str

class CilindreHiperbolic(Quadrica):
    """
    Classe per treballar amb cilindres hiperbòlic
    """
    #
    #
    #
    def __new__(cls,a2,b2,centre,eix1,eix2):
        """
        Constructor.
        Paràmetres:
           a2, b2: semieixos al quadrat
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(centre,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if centre.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,centre,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(centre,base)
        g = mcd_llista([a2,b2])
        t = a2 * b2 // g
        a2 = a2 // g
        b2 = b2 // g
        m = Matriu.diagonal(Vector([b2,-a2,0,-t]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un cilindre hiperbòlic de manera aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            centre = Punt.aleatori(l=3,maxim=4)
            trobat = centre.length() > 0
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36,40,45,48,50,60]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            trobat = a2 != b2
        return cls(a2,b2,centre,eix1,eix2)
    #
    #
    #
    def centre(self):
        """
        Retorna un centre del cilindre hiperbòlic
        """
        return (self.ref.origen)
    #
    #
    #
    def centres(self):
        """
        Retorna la recta de centres
        """
        return RectaAfi(self.ref.origen,self.ref.base[2])
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos del cilindre hiperbòlic
        """
        a = self.matriu[0,0]
        b = - self.matriu[1,1]
        f = - self.matriu[3,3]
        return (sqrt(Rational(f,a)),sqrt(Rational(f,b)))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat del cilindre hiperbòlic
        """
        a = self.matriu[0,0]
        b = - self.matriu[1,1]
        f = - self.matriu[3,3]
        return (Rational(f,a),Rational(f,b))
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda del cilindre hiperbòlic
        """
        a2, b2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "x'^2"
        else:
            str = f"\\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " - y'^2 = 1"
        else:
            str += f" - \\frac{{y'^2}}{{ {latex(b2)} }} = 1"
        return str

class ParaboloideElliptic(Quadrica):
    """
    Classe per treballar amb paraboloides el·líptics
    """
    #
    #
    #
    def __new__(cls,a2,b2,vertex,eix1,eix2):
        """
        Constructor.
        Retorna el paraboloide el·líptic amb vèrtex "vertex" i semieixos l*a2
        i l*b2, on l és la longitud de eix3
        Paràmetres:
           a2, b2: nombres enters positius
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(vertex,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if vertex.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        a2 = sympify(a2)
        b2 = sympify(b2)
        if not a2.is_integer or not b2.is_integer:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,vertex,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(vertex,base)
        q = base.quadrats_longituds()
        t = a2 * b2
        a2 = t // a2
        b2 = t // b2
        t *= sqrt(q[2])
        m = Matriu(Matrix([[a2,0,0,0],[0,b2,0,0],[0,0,0,-t/2],[0,0,-t/2,0]]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un paraboloide el·líptic de manera aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            vertex = Punt.aleatori(l=3,maxim=4)
            trobat = vertex.length() > 0
        c = [1,2,3,4,5,8,10,12,16,18,20,25,36,40,45,48,50,60]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            trobat = a2 != b2
        a2, b2 = sorted([a2,b2])[::-1]
        return cls(a2,b2,vertex,eix1,eix2)
    #
    #
    #
    def vertex(self):
        """
        Retorna el vèrtex del paraboloide el·líptic
        """
        return (self.ref.origen)
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos del paraboloide el·líptic
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        t = - 2 * self.matriu[2,3]
        return (sqrt(t/a),sqrt(t/b))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat del paraboloide el·líptic
        """
        a = self.matriu[0,0]
        b = self.matriu[1,1]
        t = - 2 * self.matriu[2,3]
        return (t/a,t/b)
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda del paraboloide el·líptic
        """
        a2, b2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "z' = x'^2"
        else:
            str = f"z' = \\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " + y'^2"
        else:
            str += f" + \\frac{{y'^2}}{{ {latex(b2)} }}"
        return str

class ParaboloideHiperbolic(Quadrica):
    """
    Classe per treballar amb paraboloides hiperbòlics
    """
    #
    #
    #
    def __new__(cls,a2,b2,vertex,eix1,eix2):
        """
        Constructor.
        Retorna el paraboloide hiperbòlic amb vèrtex "vertex" i semieixos l*a2
        i l*b2, on l és la longitud de eix3
        Paràmetres:
           a2, b2: nombres enters positius
           eix1: direcció de l'eix principal (de les x')
           eix2: en pricipi és la direcció de l'eix de les y', però si no és
                 perpendicular a eix1, es calcula el perpendicular que està
                 en el mateix pla vectorial que <eix1,eix2> aplicant el
                 mètode de Gram-Schmidt
        """
        if not isinstance(vertex,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if not isinstance(eix2,Vector):
            return None
        if vertex.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix2.dimensio != 3:
            return None
        if eix1.length() == 0:
            return None
        if eix2.length() == 0:
            return None
        m = Matriu.from_vectors_columna([eix1,eix2])
        if m.rank() != 2:
            return None
        if a2 <= 0:
            return None
        if b2 <= 0:
            return None
        a2 = sympify(a2)
        b2 = sympify(b2)
        if not a2.is_integer or not b2.is_integer:
            return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,a2,b2,vertex,eix1,eix2):
        s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(vertex,base)
        g = mcd_llista([a2,b2])
        t = a2 * b2 // g
        a2 = a2 // g
        b2 = b2 // g
        q = base.quadrats_longituds()
        if random.randint(0,1):
            t /= sqrt(q[2])
        else:
            t *= sqrt(q[2])
        m = Matriu(Matrix([[b2,0,0,0],[0,-a2,0,0],[0,0,0,-t/2],[0,0,-t/2,0]]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un paraboloide hiperbòlic de manera aleatòria
        """
        trobat = False
        while not trobat:
            eix1 = Vector.aleatori(l=3,maxim=3)
            eix2 = Vector.aleatori(l=3,maxim=3)
            if eix1.nzeros() > 1 or eix2.nzeros() > 1:
                continue
            trobat = Matriu.from_vectors_columna([eix1,eix2]).rank() == 2
        trobat = False
        while not trobat:
            vertex = Punt.aleatori(l=3,maxim=4)
            trobat = vertex.length() > 0
        c = [1,2,3,4,5,8,10,12,16]
        trobat = False
        while not trobat:
            a = random.randint(0,len(c) - 1)
            a2 = c[a]
            a = random.randint(0,len(c) - 1)
            b2 = c[a]
            trobat = a2 != b2
        return cls(a2,b2,vertex,eix1,eix2)
    #
    #
    #
    def vertex(self):
        """
        Retorna el vèrtex del paraboloide hiperbòlic
        """
        return (self.ref.origen)
    #
    #
    #
    def semieixos(self):
        """
        Retorna els semieixos del paraboloide hiperbòlic
        """
        a = self.matriu[0,0]
        b = - self.matriu[1,1]
        t = - 2 * self.matriu[2,3]
        return (sqrt(t,a),sqrt(t,b))
    #
    #
    #
    def semieixos_quadrats(self):
        """
        Retorna els semieixos al quadrat del paraboloide hiperbòlic
        """
        a = self.matriu[0,0]
        b = - self.matriu[1,1]
        t = - 2 * self.matriu[2,3]
        return (t/a,t/b)
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda del paraboloide hiperbòlic
        """
        a2, b2 = self.semieixos_quadrats()
        if a2 == 1:
            str = "z' = x'^2"
        else:
            str = f"z' = \\frac{{x'^2}}{{ {latex(a2)} }}"
        if b2 == 1:
            str += " - y'^2"
        else:
            str += f" - \\frac{{y'^2}}{{ {latex(b2)} }}"
        return str

class CilindreParabolic(Quadrica):
    """
    Classe per treballar amb cilindres parabòlics
    """
    #
    #
    #
    def __new__(cls,p,vertex,eix1,eix2=None):
        """
        Constructor.
        Retorna el cilindre parabòlic amb vèrtex "vertex", eix de les x' amb
        direcció "eix", eix de les y' amb direcció "eix2" i paràmetre de
        la paràbola "p"
        Paràmetres:
           p: paràmetre de l'equació reduïda z' = \\frac{x'^2}{2 * p}
           vertex: un vèrtex del cilindre parabòlic
           eix1: direcció de l'eix de les x'
           eix2: direcció de l'eix de les y'. Si és None, el tria el programa
        """
        if not isinstance(vertex,Punt):
            return None
        if not isinstance(eix1,Vector):
            return None
        if vertex.dimensio != 3:
            return None
        if eix1.dimensio != 3:
            return None
        if eix1.nzeros() >= 2:
            return None
        if p == 0:
            return None
        if eix2 is not None:
            if not isinstance(eix2,Vector):
                return None
            if eix2.dimensio != 3:
                return None
            w = eix1.cross(eix2,simplificar=True)
            if w.length() == 0:
                return None
        return super(Quadrica,cls).__new__(cls)
    #
    #
    #
    def __init__(self,p,vertex,eix1,eix2=None):
        if eix2 is None:
            s = SubespaiVectorial([eix1])
        else:
            s = SubespaiVectorial([eix1,eix2])
        base = s.amplia_base(unitaria=True)
        r = ReferenciaAfi(vertex,base)
        m = Matriu(Matrix([[1,0,0,0],[0,0,0,0],[0,0,0,-p],[0,0,-p,0]]))
        Quadrica.__init__(self,m,r)
    #
    #
    #
    @classmethod
    def aleatoria(cls):
        """
        Retorna un cilindre parabòlic de manera aleatòria
        """
        vertex = Punt.aleatori(l=3,maxim=4,nuls=False)
        trobat = False
        while not trobat:
            eix = Vector.aleatori(l=3,maxim=3)
            trobat = eix.nzeros() < 2 and eix.length() < 12
        s = SubespaiVectorial([eix])
        base = s.amplia_base(unitaria=True)
        q = base.quadrats_longituds()
        c = [k for k in range(1,3)] + [-k for k in range(1,3)]
        c += [Rational(1,k) for k in range(1,4)] + [-Rational(1,k) for k in range(1,4)]
        a = random.randint(0,len(c) - 1)
        p = 2 * c[a] * sqrt(q[2])
        return cls(p,vertex,eix)
    #
    #
    #
    def parametre(self):
        """
        Retorna el paràmetre de la paràbola
        """
        return - self.matriu[2,3]
    #
    #
    #
    def vertex(self):
        """
            Retorna l'origen de la referència principal
        """
        return (self.ref.origen)
    #
    #
    #
    def equacio_reduida(self):
        """
        Retorna l'equacio reduïda del cilindre parabòlic
        """
        p = self.parametre()
        if 2 * p == 1:
            return f"z' = x'^2"
        if 2 * p == -1:
            return f"z' = - x'^2"
        if p > 0:
            return f"z' = \\frac{{x'^2}}{{ {latex(2 * p)} }}"
        return f"z' = -\\frac{{x'^2}}{{ {latex(-2 * p)} }}"

class RectaRegressio(object):
    """
    Classe per treballar amb rectes de regressió
    Atributs:
        punts: llista de punts
        A: matriu dels coeficients de les incògnites
        B: vector de termes independents
        solucio: solució del sistema d'equacions A^tAX = A^tB
    """
    #
    #
    #
    def __new__(cls,punts):
        """
        Retorna un objecte RectaRegressio definit a partir d'una llista de punts
        Paràmetres:
            punts: llista de punts de dimensió 2
        """
        if not isinstance(punts,list) and not isinstance(punts,tuple):
            return None
        if len(punts) <= 2:
            return None
        for k in punts:
            if not isinstance(k,Punt):
                return None
            if k.dimensio != 2:
                return None
        return super(RectaRegressio,cls).__new__(cls)
    #
    #
    #
    def __init__(self,punts):
        self.punts = punts
        a = Vector([k[0] for k in self.punts])
        b = Vector(len(self.punts) * [1])
        A = Matriu.from_vectors_columna([a,b])
        B = Vector([k[1] for k in self.punts])
        s = SistemaEquacions(A.transposada() * A,A.transposada() * B)
        s.resol()
        self.A = A
        self.B = B
        self.solucio = Vector(s.solucio)
    #
    #
    #
    def equacio(self):
        """
        Retorna l'equació de la recta de regressió expressada en LaTeX
        """
        x = symbols('x')
        n = Vector.nul(dim=2)
        if self.solucio == n:
            return "y=0"
        a, b = self.solucio.factor_comu()
        eq = a.p * (b[0]*x + b[1])
        if a.q == 1:
            return f"y = {mylatex(eq)}"
        return f"y = \\frac{{ {mylatex(eq)} }}{{ {a.q} }}"
    #
    #
    #
    def error_quadratic(self):
        """
        Retorna d'error quadràtic
        """
        e = self.A * self.solucio - self.B
        return e.length()
    #
    #
    #
    @classmethod
    def aleatoria(cls,l=4,max=4):
        """
        Retorna un problema aleatori amb l punts
        Paràmetres:
            l: nombre de punts
            max: valor màxim de les ys
        """
        if l > 13:
            return None
        p = [-k for k in range(l)] + [k for k in range(l)]
        random.shuffle(p)
        xs = p[0:l]
        xs.sort()
        ys = [random.randint(-max,max) for k in range(l)]
        punts = [Punt(xs[k],ys[k]) for k in range(l)]
        return cls(punts)
    #
    #
    #
    def llista_punts(self):
        """
        Retorna la llista de punts en format LaTeX
        """
        l = list(map(lambda item: f"${item}$",self.punts))
        return ", ".join(l[0:-1]) + f" i {l[-1]}"
    #
    #
    #
    def taula_punts(self):
        """
        Retorna una taula en format LaTeX dels punts
        """
        format = "c|" + len(self.punts) * "r"
        xs = " & ".join([f"${p[0]}$" for p in self.punts])
        ys = " & ".join([f"${p[1]}$" for p in self.punts])
        return f"\\begin{{tabular}}{{{format}}} {xs} \\\\ \\hline {ys} \\end{{tabular}}"
