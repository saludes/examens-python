from hypothesis.strategies import integers, complex_numbers
from exams import LaTeXExamExam

ex = LaTeXExamExam()

complex = complex_numbers(min_magnitude=1, max_magnitude=2)
smallint = integers(min_value=2,max_value=6)
    


@ex.problem
def potencia(draw):
    """Calculeu la potència {{EXPONENT}} de {{COMPLEX}}"""
    z = draw(complex)
    n = draw(smallint)
    return dict(COMPLEX=z, EXPONENT=n)

@ex.problem
def arrel(draw):
    """Calculeu $\sqrt[{{EXPONENT}}]{ {{COMPLEX}} }$"""
    return dict(
        COMPLEX=draw(complex),
        EXPONENT=draw(smallint))

