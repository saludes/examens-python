from hypothesis.strategies import integers, complex_numbers, composite
from exams import LaTeXExamExam

ex = LaTeXExamExam()

complex = complex_numbers(min_magnitude=1, max_magnitude=2)
smallint = integers(min_value=2,max_value=6)


@composite
def sympy_complex(draw,min_magnitude=1, max_magnitude=2, tolerance=0.01):
    from hypothesis.strategies import fractions
    from hypothesis import assume
    from sympy import I, nsimplify
    z = draw(complex_numbers(min_magnitude=min_magnitude, max_magnitude=max_magnitude))
    z_ = nsimplify(z, rational=True, tolerance=tolerance)
    assume(z_.as_real_imag()[1] != 0)
    return z_




@ex.problem
def potencia(draw):
    """Calculeu la potència {{EXPONENT}} de {{COMPLEX}}"""
    scx = sympy_complex(min_magnitude=1, max_magnitude=2, tolerance=0.1)
    z = draw(scx)
    # z = draw(complex)
    n = draw(smallint)
    zn = z**n
    return dict(COMPLEX=z, EXPONENT=n, SOLUCIO=zn)

@ex.problem
def arrel(draw):
    """Calculeu $\sqrt[{{EXPONENT}}]{ {{COMPLEX}} }$"""
    z = draw(complex)
    n = draw(smallint)
    zn = z**(1/n)
    return dict(COMPLEX=z, EXPONENT=n, SOLUCIO=zn)

