from exams import ClozeExam

ex = ClozeExam()


@ex.problem
def suma():
    "Calculeu: {{A}} + {{B}} = {{SUMA}}"
    from random import randint
    a = randint(0,9)
    b = randint(0,9)
    s = a+b
    dificultat = int(s/8)+1
    return dict(A=str(a),
                B=str(b),
                SUMA=ex.short_answer(
                    (str(s),100),
                    ('*',0,"No és correcte"),
                    value=dificultat))


if __name__ == '__main__':
    for _ in range(5):
        print(ex.render())
