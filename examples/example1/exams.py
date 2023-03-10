class Problem:
  def __init__(self, name, func, **kw):
    from hypothesis.strategies import composite
    self.name = name
    self.content = func.__doc__
    self.f = composite(func)()
    self.args = kw

  def __repr__(self):
    return f"<Problem {self.name}>"
  
  def __call__(self):
     return self.f.example()
  
  def _split(self, text):
    import re
    parts = re.split(r'\n\s*=+\n', text)
    if len(parts) == 2:
       return dict(CONTENT=parts[0], SOLUTION=parts[1])
    elif len(parts)==1:
       return dict(CONTENT=parts[0])
    else:
       raise ValueError("Invalid number of parts.")
     
     
      

  def render(self, template=None):
    from jinja2 import Template
    #kw = {k: str(v) for k, v in self().items()}
    kw = self()
    if template:
       return self.with_parts(template, kw)
    else:
      return Template(self.content).render(**kw)

    
  def with_parts(self, template, params):
    from jinja2 import Template
    d = self._split(self.content)
    tmp = template.render(d)
    return Template(tmp).render(**params)


class Exam: 
  slots = ['template', 'problem_template']

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
    temp = Template(self.problem_template)
    render = lambda pr: pr.render(temp)
    kw['exam'] = [dict(id=pr.name, r=render(pr)) for pr in self]
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
  % {{pr.id}}
  {{pr.r}}
{% endfor %}
\end{questions}
\end{document}
"""

  problem_template=r'''
  \question {{CONTENT}}

  \begin{solution}{{SOLUTION}}\end{solution}
  '''

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
  

class ClozeExam(Exam):
  template = """
    {% for pr in exam %}
      {{pr.id}}: {{pr.r}}
    {% endfor %}
    """

  def multichoice(self, *answers):
      d = {}
      for a in answers:
          if isinstance(a, tuple):
              ans,va = a
              assert ans not in d
              d[ans] = va
          elif isinstance(a, str):
              d[a] = False
          else:
              raise NotImplementedError(a)
      if sum(int(v) for v in d.values()) != 1:
          raise ValueError("Falta resposta correcta o més d'una")
      items = '~'.join([('=' if v else '') + ans for (ans,v) in d.items()])
      #items = items.replace('~=', '=')
      return "{1:MULTICHOICE:" + items + "}"

  def short_answer(self, *answers):
      ans = []
      for a in answers:
          if isinstance(a, tuple):
              fb = '' if len(a) == 2 else a[2]
              v = a[1]
              answer = a[0]
          elif isinstance(a, str):
              fb = ''
              v = 0
              answer = a
          else:
              raise NotImplementedError()
          ans.append((answer, v, fb))
          total = sum(v for (_,v,_) in ans)
      
      def fmt_sa(el):
          ans,v,fb = el
          if v == 0:
              item = ans
          else:
              item =f'%{v}%' + ans
          if fb:
              item += "#" + fb
          return item
      
      ans = "~".join(map(fmt_sa, ans))
      return "{1:SHORTANSWER:" + ans + "}"
  


  if __name__ == '__main__':
    import csv
    from Problemes2 import ex
    with open("estudiants.csv") as ests:
      for r in csv.reader(ests, delimiter=':'):
        nom = r[0] + ' ' + r[1]
        email = r[3]
        fname = f"tex/{nom}.tex"
        with open(fname, 'w') as tex:
          tex.write(ex.render(estudiant=nom, seed=email))
          print("Writing file", fname)
    