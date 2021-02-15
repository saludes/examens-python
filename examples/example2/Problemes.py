#!/usr/bin/python3
# -*- coding: utf-8 -*-

import random
from sympy import *
from Algebra import *

class Problemes:
    def __init__(self):
        random.seed()
    #
    #
    #
    def problema1(self):
        a = Matriu.diagonalitzable(ordre=3,maxim=5,mzeros=0,mvaps=3,vapsnuls=False,vapsrepetits=False)
        p = a.polinomi_caracteristic()
        base = Base(a.veps)
        d = Matriu.diagonal(a.vaps)
        c = base.matriu()
        return {'MATRIU' : f"{a}", 'DIAGONAL' : f"{d}",
                'POLINOMI' : f"{p}", 'CANVI' : f"{c}"}
    #
    #
    #
    def problema2(self):
        m = symbols('m')
        trobat = False
        while not trobat:
            a = Matriu.amb_rang(f=3,c=3,r=2,maxim=4,mzeros=0)
            v = a.nucli()[0]
            trobat = v.nzeros() == 0 and v.length() < 5
        aa = a.adjunta()
        x = Vector.aleatori(l=3,maxim=3,nuls=False)
        b = a * x
        trobat = False
        while not trobat:
            i = random.randint(0,2)
            j = random.randint(0,2)
            trobat = aa[i,j] != 0
        valorm = a[i,j]
        a[i,j] = m
        k = random.randint(0,2)
        valorb = b[k]
        b[k] = m + valorb - valorm
        s = SistemaEquacions(a,b)
        sistema = f"{s}"
        a[i,j] = valorm
        b[k] = valorb
        s = SistemaEquacions(a,b)
        solucio = s.solucio_latex()
        return {'EQUACIONS' : sistema, 'VALOR' : f"{valorm}",
                'SOLUCIO' : f"{solucio}"}
    #
    #
    #
    def problema3(self):
        p0 = Punt.aleatori(l=3,maxim=3,nuls=False)
        v = Vector.aleatori(l=3,maxim=1,nuls=False)
        v.simplificar()
        angles = [120,240]
        angle = angles[random.randint(0,1)]
        r = TransformacioLineal.rotacio(v,angle)
        psi, theta, phi = r.angles_euler()
        recta = RectaAfi(p0,v)
        factors = [-3,-2,-1,1,2,3]
        factor = factors[random.randint(0,5)]
        h = TransformacioAfi.moviment_helicoidal(recta,angle,alpha=factor)
        return {'EQUACIO' : recta.equacio_continua(), 'ANGLE' : f"{angle}",
                'TRANSLACIO' : f"{factor * v}", 'PSI' : f"{psi}" , 'VECTOR' : f"{v}",
                'THETA' : f"{theta}", 'PHI' : f"{phi}", 'HELICOIDAL' : f"{h}" }
    #
    #
    def problema4(self):
        u = Vector.aleatori(l=3,maxim=3,nuls=False)
        v = Vector(0,0,0)
        while Matriu.from_vectors_fila([u,v]).rank() < 2:
            v = Vector.aleatori(l=3,maxim=2,nuls=False)
        w = u.cross(v,simplificar=True)
        p0 = Punt.aleatori(l=3,maxim=3,nuls=False)
        q0 = Punt.aleatori(l=3,maxim=3,nuls=False)
        r1 = RectaAfi(p0,u)
        r2 = RectaAfi(q0,v)
        e = EquacioLineal.coeficients(w,w.dot(p0))
        return {'RECTAA' : r1.equacio_continua(),
                'RECTAB' : f"{r2.equacions_implicites(aleatori=True)}",
                'DISTANCIA' : latex(r1.distancia(r2)), 'PLA' : f"{e}"}
    #
    #
    #
    def problema5(self):
        q = Conica.aleatoria(maxim=20)
        o = q.ref.origen
        return {'EQUACIO' : f"{q}", 'PRINCIPAL' : f"{q.ref}",
                'REDUIDA' : q.equacio_reduida(), 'TIPUS' : q.tipus(),
                'CONICA' : q.to_asy(), 'MINIMY' : f"{o[1]-10}",
                'MAXIMY' : f"{o[1]+10}", 'MINIMX' : f"{o[0]-13}", 'MAXIMX' : f"{o[0]+13}"}
    #
    #
    #
    def problemes(self):
        return [self.problema1,self.problema2,self.problema3,
                self.problema4,self.problema5]
