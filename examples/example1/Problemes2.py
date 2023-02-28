#!/usr/bin/python3
# -*- coding: utf-8 -*-

import random
from sympy import *
from Algebra import *
from exams import LaTeXExamExam

ex = LaTeXExamExam()


@ex.problem
def problema1():
  r"""Trobeu la inversa de la matriu $A = {{MATRIU}}$.

    $$A^{-1}={{INVERSA}}$$
    """
  m = Matriu.invertible()
  inversa = m.inversa()
  return dict(MATRIU=m, INVERSA=inversa)


@ex.problem
def problema2():
  r"""
  Es considera el pla vectorial $P$ de $\mathbb{R}^3$
  $P = \langle{{VECTORS}}\rangle$
  i la base $\mathcal B'={{BASE}}$. Determineu l'equació implícita del pla $P$ respecte de la base $\cal B'$.
  """
  m = Matriu.invertible(maxim=2, mzeros=1, unitaria=True)
  trobat = False
  while not trobat:
    n = Matriu.aleatoria(f=3, c=2, maxim=2, nuls=False)
    if n.rank() < 2:
      continue
    trobat = True
  base = Base.from_matriu(m)
  p = PlaVectorial.from_matriu(n)
  e = p.equacio_implicita(base, prime=1)
  return dict(VECTORS=p, BASE=base, EQUACIO=e)


@ex.problem
def problema3():
  r"""
  Resoleu el següent sistema d'equacions
  $$ {{SISTEMA}}, $$
  """
  x = Vector.aleatori(l=3, maxim=4, nuls=False)
  m = Matriu.invertible(maxim=3, mzeros=1)
  r = m.reordena_aleatoriament_files()
  q = random.randint(0, 3)
  trobat = False
  while not trobat:
    c = Vector.aleatori(l=3, maxim=2, nuls=False)
    c[2] = 0
    c = c * r
    trobat = c.nzeros() == 0
  m = m.inserta_fila(q, c)
  b = m * x
  s = SistemaEquacions(m, b)
  solucio = f"$x={latex(x[0])}$, $y={latex(x[1])}$, $z={latex(x[2])}$."
  return dict(SISTEMA=s, SOLUCIO=solucio)


@ex.problem
def problema4():
  """
  Determineu l'equació implícita del pla perpendicular a la recta
  $${{RECTA}}$$ que passa pel punt ${{PUNT}}$.
  """
  m = Matriu.amb_rang(f=3, c=3, r=3, maxim=5, mzeros=0)
  p0, v, q0 = m.vectors_columna(m)
  p0 = p0.punt()
  q0 = q0.punt()
  r = RectaAfi(p0, v)
  p = PlaAfi.amb_associat(v, q0)
  return dict(RECTA=r.equacio_continua(),
              PUNT=q0,
              SOLUCIO=p.equacio_implicita())


if __name__ == '__main__':
  with open('sample.tex', 'w') as out:
    for _ in range(5):
      ex.shuffle()
      out.write(ex.render(estudiant='Jo mateix'))
      out.write('\n\n')
