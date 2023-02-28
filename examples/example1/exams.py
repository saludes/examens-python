class Problem:
  def __init__(self, name, func, **kw):
    self.name = name
    self.f = func
    self.args = kw

  def __repr__(self):
    return f"<Problem {self.name}>"

  def render(self):
    from jinja2 import Template
    kw = {k: str(v) for k, v in self.f().items()}
    return Template(self.f.__doc__).render(**kw)



class Exam: 
  slots = ['template']

  def __init__(self):
    self.problems = []

  def shuffle(self, seed=None):
    import random
    if seed:
      random.seed(seed)
    random.shuffle(self.problems)

  def __len__(self):
    return len(self.problems)

  def __iter__(self):
    return iter(self.problems)

  def _check_notin(self, name):
    if any(name == p.name for p in self.problems):
      raise KeyError(f"'{name}' ja és a l'examen")

  def problem(self, prob):
    name = prob.__name__
    self._check_notin(name)
    print("Registering", name)
    self.problems.append(Problem(name, prob))
    return prob

  def render(self, **kw):
    from jinja2 import Template
    kw['exam'] = [dict(id=pr.name, r=pr.render()) for pr in self]
    return Template(self.template).render(**kw)


  
class LaTeXExamExam(Exam):
  template = r"""
\documentclass{{doc_options}}{exam}
\usepackage[catalan]{babel} 
\usepackage{amssymb}
\usepackage{graphicx}
\usepackage{amsmath}
\renewcommand{\solutiontitle}{\noindent\textbf{Resolució:}\enspace}
\firstpagefooter{ {{estudiant}} }{}{}

\begin{document}
%% seed: {{seed}}
\begin{questions}
{% for pr in exam %}
  \question ({{pr.id}}) {{pr.r}}
{% endfor %}
\end{questions}
\end{document}
"""

  def __init__(self):
    super().__init__()
    self.size = 12
    self.answers = False

  def configure_docclass(self):
    size_pts = f"{self.size}pt"
    answers = "answers" if self.answers else "noanswers"
    return '[' + ','.join([size_pts, answers]) + ']'

  
  def problems_from_module(self, module):
    from inspect import getmembers, isfunction
    for (name, func) in getmembers(module, isfunction):
      if name.startswith('problem'):
        self._check_notin(name)
        self.problems.append(Problem(name, func))

  def render(self, **kw):
    kw['doc_options'] =self.configure_docclass()
    return super().render(**kw)
  

  if __name__ == '__main__':
    import csv
    from Problemes2 import ex
    with open("estudiants.csv") as ests:
      for r in csv.reader(ests, delimiter=':'):
        nom = r[0] + ' ' + r[1]
        email = r[3]
        with open(f"tex/{nom}.tex", 'w') as tex:
          tex.write(ex.render(estudiant=nom, seed=email))
    