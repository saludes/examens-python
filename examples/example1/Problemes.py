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
        m = Matriu.invertible()
        inversa = m.inversa()
        return {'MATRIU' : f"{m}", 'INVERSA' : f"{inversa}"}
    #
    #
    #
    def problema2(self):
        m = Matriu.invertible(maxim=2,mzeros=1,unitaria=True)
        trobat = False
        while not trobat:
            n = Matriu.aleatoria(f=3,c=2,maxim=2,nuls=False)
            if n.rank() < 2:
                continue
            trobat = True
        base = Base.from_matriu(m)
        p = PlaVectorial.from_matriu(n)
        e = p.equacio_implicita(base,prime=1)
        return {'VECTORS' : f"{p}", 'BASE' : f"{base}", 'EQUACIO' : f"{e}"}

    def problema3(self):
        x = Vector.aleatori(l=3,maxim=4,nuls=False)
        m = Matriu.invertible(maxim=3,mzeros=1)
        r = m.reordena_aleatoriament_files()
        q = random.randint(0,3)
        trobat = False
        while not trobat:
            c = Vector.aleatori(l=3,maxim=2,nuls=False)
            c[2] = 0
            c = c * r
            trobat = c.nzeros() == 0
        m = m.inserta_fila(q,c)
        b = m * x
        s = SistemaEquacions(m,b)
        solucio = f"$x={latex(x[0])}$, $y={latex(x[1])}$, $z={latex(x[2])}$."
        return {'SISTEMA' : f"{s}", 'SOLUCIO' : solucio}

    def problema4(self):
        m = Matriu.amb_rang(f=3,c=3,r=3,maxim=5,mzeros=0)
        p0, v, q0 = m.vectors_columna(m)
        r = RectaAfi(v,p0)
        p = PlaAfi.amb_associat(v,q0)
        return {'RECTA' : f"{r.equacio_continua()}", 'PUNT' : f"{q0}", 'SOLUCIO' : f"{p.equacio_implicita()}"}
    #
    #
    #
    def problemes(self):
        return [self.problema1,self.problema2,self.problema3,self.problema4]
