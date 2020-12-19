#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Filename:   test.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2020
Disclaimer: This code is presented "as is" and it has been written to
            generate random models of exams for the subject of Linear
            Algebra at ESEIAAT, Technic University of Catalonia
License:    This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.

 	        See https://www.gnu.org/licenses/
"""

import sys
from sympy import *
from Algebra import *


c = Hiperbola.aleatoria()
v1, v2 = c.vectors_directors_asimptotes()
print(v1)
print(v2)
sys.exit(0)

#
# Tests amb la classe Vector
#
v1 = Vector([2,3,1])
v2 = Vector.aleatori(l=3,maxim=7,nuls=False,positius=True)
print(v1.dot(v2))
print(v1.cross(v2))
v = 3*v1 - 2*v2
print(v)
w = Vector([2,-Rational(1,6),-1,Rational(1,3)])
print(w)
w[0] = Rational(1,2)
print(w.latex(unitari=True))
w.simplificar()
print(w)
v = w.reordena_aleatoriament()
print(v)

v = Vector([sqrt(2) + 2 * sqrt(5),- 2 * sqrt(5) + sqrt(2)])
v.radsimplificar()
print(v)
v = Vector([-6 + 2 * sqrt(2),2 + 6 * sqrt(2)])
v.radsimplificar()
print(v)
v = Vector([-30 * sqrt(20) + 4 * sqrt(10),6 *sqrt(10) + 20 * sqrt(20)])
v.radsimplificar()
print(v)
v = Vector([-30 * sqrt(20),20 * sqrt(20)])
v.radsimplificar()
print(v)
v = Vector([15 - sqrt(15),-15 - sqrt(15)])
v.radsimplificar()
print(v)
print()

#
# Tests amb la classe Punt
#
p = Punt([3,2,1,3])
v = Vector.aleatori(l=4,maxim=4)
print(p + 3*v)
print()

#
# Tests amb la classe Base
#
b = Base([Vector([1,2]),Vector([-1,1])])
print(b)
print(b.matriu())
print(b.vectors_latex())
u = b.vector_de_components(Vector([2,3]))
print(u)
print()

#
# Tests amb la classe Matriu
#
m = Matriu()
print(m)
m = Matriu.aleatoria(f=3,c=5,maxim=5,nuls=True)
print(m)
print(m.nzeros(),m.norma_maxim(),m.rank())
m = Matriu.aleatoria(f=4,c=4,maxim=3,nuls=False)
print(m.polinomi_caracteristic())
m1 = Matriu.aleatoria(f=2,c=3,maxim=6,nuls=False)
m2 = Matriu.aleatoria(f=2,c=3,maxim=6,nuls=False)
print(m1 - m2)
print(3*m1 - 4*m2)
print(m1 * m2.transposada())
print(m2.transposada() * m1)
print(Matriu.diagonal([1,3,-2,1]))
print(Matriu.amb_rang(f=3,c=4,r=2,maxim=3,nuls=False))
m = Matriu.invertible(ordre=3,maxim=5,mzeros=0)
print(m.inversa())
m = Matriu.diagonalitzable(ordre=3,maxim=5,mzeros=0,mvaps=3)
print(m)
print(m.vaps,m.veps)
c = Matriu.from_vectors_columna(m.veps)
b = Base.from_matriu(c)
print(b)
m = Matriu.matriu_fila(Vector([1,-1,2]))
n = Matriu.matriu_columna(Vector([1,2,3]))
print (n * m)
v = (Vector([1,2,-1]),Vector([-1,1,2]),Vector([2,1,2]))
print(Matriu.from_vectors_fila(v))
m = Matriu.from_vectors_columna(v)
print(m)
print(m.vectors_columna())
m = Matriu.aleatoria(f=2,c=4,maxim=3,nuls=False)
print(m)
print(m.nucli())
n = m.reordena_aleatoriament_columnes()
print(n)
n = m.reordena_aleatoriament_files()
print(n)
n = m.inserta_fila(1,Vector([1,2,3,4]))
print(n)
p = n.inserta_columna(10,Vector([1,2,3]))
print(p)
print()

#
# Tests amb la classe EquacioLineal
#
x, y, z, t = symbols('x y z t')
e = EquacioLineal(2*x-y+4*z-6*t-5)
print(e)
e = EquacioLineal.coeficients(Vector([2,5,-2,0,3,1]),4,prime=1)
print(e)
e1 = EquacioLineal.coeficients(Vector([2,5,-2]),-3)
e2 = EquacioLineal.coeficients(Vector([1,-1,2]),4)
print(3*e1 - 5*e2)

#
# Tests amb la classe SistemaEquacions
#
m = Matriu.aleatoria(f=2,c=4,maxim=3,nuls=False)
b = Vector([4,-3])
s = SistemaEquacions(m,b)
print(s)
print(s.matriu_incognites())
print(s.matriu_ampliada())
print(s.solucio_latex())
print()

#
# Tests amb la classe EquacioParametrica
#
t1, t2, t3, t4 = symbols('t1 t2 t3 t4')
eq = -x + 2*t1 - 3*t2 + t3 - 4
e = EquacioParametrica(eq)
print(e)
e = EquacioParametrica.coeficients(Vector([3,2,5,-3]),5,3,7)
print(e)
print()

#
# Tests amb la classe EquacionsParametriques
#
m = Matriu.aleatoria(f=5,c=2,maxim=3,nuls=False)
b = Vector([4,-3,0,-1,0])
e = EquacionsParametriques(m,b)
print(e)
s = e.eliminar_parametres()
print(s)
print()

#
# Tests amb la classe PlaVectorial
#
p = PlaVectorial(Vector([1,2,-1]),Vector([2,-1,3]))
print(p)
print(p.equacio_implicita())
m = Matriu.invertible(ordre=3,maxim=3,mzeros=0)
b = Base.from_matriu(m)
print(p.equacio_implicita(base=b,prime=1))
p = PlaVectorial.amb_associat(Vector([1,-2,4]))
print(p.base_ortogonal())

#
# Tests amb la classe RectaVectorial
#
r = RectaVectorial(Vector([1,2]))
print(r.equacions_implicites())
print(r.equacio_continua())
r = RectaVectorial(Vector([3,2,-2]))
print(r.equacions_implicites(aleatori=True))
print(r.equacio_continua())
print()

#
# Tests amb la classe ReferenciaAfi
#
p = Punt([1,2,-1])
b = Base([Vector([1,1,0]),Vector([-1,1,1]),Vector([1,-1,2])])
r = ReferenciaAfi(p,b)
print(r)
q = r.punt_de_coordenades(Punt([1,5,2]))
print(q)
print(p.coordenades_en_referencia(r))
b = Base([Vector([1,1,0]),Vector([-1,1,1]),Vector([1,-1,2])],unitaria=True)
r = ReferenciaAfi(p,b)
print(r.canvi_coordenades())
p = Punt([3,-2])
b = Base([Vector([1,2]),Vector([-2,1])],unitaria=True)
r = ReferenciaAfi(p,b)
print(r.canvi_coordenades())
r2 = r.referencia_inversa()
print(r2)
print(r2.canvi_coordenades(1,0))
print()

#
# Tests amb la classe PlaAfi
#
p = Punt([2,-3,-1])
u1 = Vector([1,-3,-1])
u2 = Vector([2,1,1])
p = PlaAfi(p,u1,u2)
print(p)
print(p.equacio_implicita())
p0 = Punt([1,2,-1])
b = Base([Vector([1,1,0]),Vector([-1,1,1]),Vector([1,-1,2])])
r = ReferenciaAfi(p0,b)
print(p.equacio_implicita(ref=r,prime=1))
p = Punt([1,1,1])
u1 = Vector([1,0,-1])
u2 = Vector([1,-1,0])
p = PlaAfi(p,u1,u2,r)
print(p.equacio_implicita(ref=r,prime=1))
print()

#
# Tests amb la classe RectaAfi
#
p = Punt([3,2,-1])
u = Vector([4,-3,2])
p0 = Punt([1,2,-1])
b = Base([Vector([1,1,0]),Vector([-1,1,1]),Vector([1,-1,2])])
ref = ReferenciaAfi(p0,b)
r = RectaAfi(p,u,ref)
print(r)
print(r.equacio_continua(ref=ref,prime=1))
print(r.equacions_implicites())
print(r.equacions_implicites(ref=ref,prime=1,aleatori=True))
print()

#
# Tests amb la classe SubespaiVectorial
#
u1 = Vector([1,2,1,1])
u2 = Vector([1,-1,1,-1])
s = SubespaiVectorial([u1,u2])
print(s.base)
print(s.equacions_implicites())
v = s.suplementari_ortogonal()
print(v.base)
print(v.equacions_implicites())
print(v.base_ortogonal())
m = Matriu.invertible(ordre=4,maxim=3,mzeros=2,unitaria=True)
br4 = Base.from_matriu(m)
s = SubespaiVectorial([u1,u2],basern=br4)
print(s.equacions_implicites(basern=br4,prime=1))
m = Matriu.invertible(ordre=4,maxim=1,mzeros=4,unitaria=True)
cr4 = Base.from_matriu(m)
print(s.equacions_implicites(basern=cr4,prime=2))
print()

#
# Tests amb la classe VarietatLineal
#
u1 = Vector([1,1])
u2 = Vector([1,-1])
s = SubespaiVectorial([u1])
b = Base([u1,u2])
p = Punt([3,1])
r = ReferenciaAfi(p,b)
v = VarietatLineal(p,s)
print(v.equacions_implicites())
print (v.equacions_implicites(ref=r,prime=1))
w = v.varietat_ortogonal(Punt([0,0]))
print (w.equacions_implicites())
print()

#
# Tests amb la classe TransformacioLineal
#
a = Matriu.aleatoria(f=2,c=2)
b = Matriu.aleatoria(f=2,c=2)
t1 = TransformacioLineal(a)
t2 = TransformacioLineal(b)
print(t1)
print(t2)
print(t1 * t2)
s = SubespaiVectorial([Vector([2,1,1]),Vector([1,3,-1])])
t = TransformacioLineal.projeccio_ortogonal(s)
print(t)
t = TransformacioLineal.simetria(s)
print(t)
b = t.base
print(t.latex(b,prime=1))
vecs = [Vector.aleatori(l=4,maxim=2,nuls=False) for k in range(3)]
s =  SubespaiVectorial(vecs)
t = TransformacioLineal.simetria(s)
print(t)
v = Vector([1,-1,1])
t = TransformacioLineal.rotacio(v,60)
e, angle = t.eix_angle_rotacio()
print(angle)
print(e)
print(t)
print()

#
# Tests amb la classe TransformacioAfi
#
u1 = Vector([1,-1,0])
u2 = Vector([1,0,-1])
s = SubespaiVectorial([u1,u2])
p = Punt([2,1,3])
v = VarietatLineal(p,s)
t = TransformacioAfi.projeccio_ortogonal(v)
print (t)
q = p + 5*Vector([1,1,1])
print(q,t.transforma(q))

#
# Tests amb la classe FormaQuadratica
#
q = FormaQuadratica.aleatoria()
print(q.polinomi_caracteristic())
