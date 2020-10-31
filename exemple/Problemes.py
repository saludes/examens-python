#!/usr/bin/python3
# -*- coding: utf-8 -*-

import random
from sympy import *
from Algebra import *

class Problemes:
    def __init__(self):
        pass
    #
    #
    #
    def problema1(self):
        m = matriu_invertible(mzeros=0)
        matriu = matriu_latex(m)
        inversa = matriu_inversa_latex(m)
        return {'MATRIU' : matriu, 'INVERSA' : inversa}
    #
    #
    #
    def problema2(self):
        r = {'x' : "x'", 'y' : "y'", 'z' : "z'"}
        m = matriu_invertible(maxim=2,mzeros=1,unitaria=True)
        trobat = False
        while not trobat:
            n = matriu_aleatoria(f=3,c=2,maxim=2,nuls=False)
            if n.rank() < 2:
                continue
            trobat = True
        vecs = text_vectors_matriu(n)
        base = text_vectors_matriu(m)
        eq = equacio_del_pla_vectorial(n,m,r)
        return {'VECTORS' : vecs, 'BASE' : base, 'EQUACIO' : eq}

    def problema3(self):
        x = matriu_aleatoria(f=3,c=1,maxim=4,nuls=False)
        m = matriu_invertible(maxim=3,mzeros=1)
        p = [0,1,2]
        random.shuffle(p)
        q = random.randint(0,3)
        trobat = False
        while not trobat:
            c = matriu_aleatoria(f=1,c=2,maxim=2,nuls=False)
            row = c[0,0] * m.row(p[0])  +  c[0,1] * m.row(p[1])
            if nzeros(row) > 0:
                continue
            trobat = True
        m = m.row_insert(q,row)
        b = m * x
        sistema = sistema_equacions(m,b)
        solucio = f"$x={latex(x[0])}$, $y={latex(x[1])}$, $z={latex(x[2])}$."
        return {'SISTEMA' : sistema, 'SOLUCIO' : solucio}

    def problema4(self):
        m = matriu_amb_rang(f=3,c=3,r=3,maxim=5,nuls=False)
        p, v, q = vectors_matriu(m)
        eq = equacio_continua(p,v)
        punt = text_vectors_matriu(Matrix(3,1,q))
        sol = equacio_del_pla_afi(v,q)
        return {'RECTA' : eq, 'PUNT' : punt, 'SOLUCIO' : sol}
    #
    #
    #
    def problemes(self):
        return [self.problema1,self.problema2,self.problema3,self.problema4]
