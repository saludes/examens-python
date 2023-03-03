from exams import ClozeExam

ex = ClozeExam()


@ex.problem
def questio1():
    """
Match the following cities with the correct state:
    * San Francisco: {{FRISCO}}
    * Tucson: {{TUCSON}}
    * Los Angeles: {{LA}} 
    * Phoenix: {{PHOENIX}}"""
    states = ('California', 'Arizona')
    isState = lambda i: ex.multichoice(*[(c, c == states[i]) for c in states])
    return dict(
        FRISCO= isState(0),
        TUCSON= isState(1),
        LA= isState(0),
        PHOENIX= isState(1))


@ex.problem
def questio2():
    """
The capital of France is {{FRANCE}}"""
    return dict(
        FRANCE= ex.short_answer(
            ('Paris', 100, "The capital of France is Paris, of course."),
            ('Marseille', 50, "No, that is the second largest city in France (after Paris)"),
            ('*', 0, "Wrong answer. The capital of France is Paris, of course.")))


@ex.problem
def aritmetic():
    """
Digues quan és el mínim comú múltiple de {{N1}} i {{N2}}: {{LCM}}"""
    from random import randint
    from math import lcm
    s1 = randint(10,100)
    s2 = randint(10,100)
    m = lcm(s1, s2)
    return dict(
        N1=s1,
        N2=s2,
        LCM=ex.short_answer(
            (str(m), 100),
            ('*',0, "està malament")))




