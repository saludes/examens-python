from hypothesis.strategies import integers, complex_numbers, composite
from hypothesis import assume
   
from exams import LaTeXExamExam

ex = LaTeXExamExam()

@composite
def complex(draw,min_magnitude=1, max_magnitude=2,mpart=0.1):
    z = draw(complex_numbers(min_magnitude=min_magnitude, max_magnitude=max_magnitude))
    assume(abs(z.real)>mpart and abs(z.imag)>mpart)
    return z


smallint = integers(min_value=2,max_value=6)


@composite
def sympy_complex(draw,min_magnitude=1, max_magnitude=2, tolerance=0.01):
    from sympy import I, nsimplify
    z_ = draw(complex_numbers(min_magnitude=min_magnitude, max_magnitude=max_magnitude))
    z = nsimplify(z_, rational=True, tolerance=tolerance)
    re,im = z.as_real_imag()
    assume(re!=0 and im!=0)
    return z




@ex.problem
def potencia(draw):
    """Calculeu la potència {{EXPONENT}} de {{COMPLEX}}
    =====
    {% if EXPONENT <= 3 %}
    Podem aplicar la propietat distributiva a {{EXPAND}}, que dóna {{SOLUCIO}},

    O altrament:
    {% endif %}
    Passem a la forma polar, obtenim $|z|={{ZABS}}$ amb argument {{ZARG}}.
    Per tant la potència tindrà $|z^{{EXPONENT}}| = {{WABS}}$ i argument {{WARG}}, que passat a forma binomial és {{SOLUCIO}}.
    """
    from sympy import expand, arg
    scx = sympy_complex(min_magnitude=1, max_magnitude=2, tolerance=0.1)
    z = draw(scx)
    # z = draw(complex)
    n = draw(smallint)
    zn = expand(z**n)
    return dict(
        COMPLEX=z, EXPONENT=n, SOLUCIO=zn,   
        EXPAND= z**n, ZABS = abs(z), ZARG = arg(z),
        WABS= abs(z)**n, WARG = n*arg(z))
    
    
@ex.problem
def arrel(draw):
    """Calculeu $\sqrt[{{EXPONENT}}]{ {{COMPLEX}} }$
    =====
    Passem a la forma polar, amb valor absolut {{ZABS}} i argument {{ZARG}}.
    Fem l'arrel {{EXPONENT}} del valor absolut que dóna {{WABS}}, i dividim l'argument
    per {{EXPONENT}} tenint en compte $2\pi$, que és {{WARG}}, que tornant a la forma binomial dóna {{SOLUCIO}}
    """
    from sympy import Integer, arg
    z = draw(complex())
    n = draw(smallint)
    zn = z**(1/n)
    return dict(
        COMPLEX=z,
        EXPONENT=n,
        SOLUCIO=zn,
        ZABS=abs(z),
        ZARG=arg(z),
        WABS=abs(z)**(1/Integer(n)),
        WARG=arg(z)/n)

