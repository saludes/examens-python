#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename:   examen.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2020-2021
Disclaimer: This program is provided "as is", without warranty of any kind,
            either expressed or implied, including, but not linmited to, the
            implied warranties of merchantability and fitness for a particular
            purpose.
            It has been written to generate random models of exams for the
            subject of Linear Algebra at ESEIAAT, Technic University of Catalonia
License:    This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.

 	        See https://www.gnu.org/licenses/
"""

import sys
import os
import re
import unidecode
import random
import glob
import json
import subprocess
from optparse import OptionParser

sys.path.append('.')
try:
  from Problemes import Problemes
except:
  pass


class Examen:

  def __init__(self):
    self.parser = OptionParser()
    self.parser.add_option("--examen", dest="examen", default=None)
    self.parser.add_option("--estudiants", dest="estudiants", default=None)
    self.parser.add_option("--problemes", dest="problemes", default=None)
    self.parser.add_option("--possibles-problemes",
                           dest="possibles",
                           default=None)
    self.parser.add_option("--incompatibles",
                           dest="incompatibles",
                           default=None)
    self.parser.add_option("--grups", dest="grups", default=None)
    self.parser.add_option("--dades", dest="fitxerdades", default=None)
    self.parser.add_option("--tex-engine", dest="engine", default=None)
    self.parser.add_option("--no-solucions",
                           action="store_false",
                           dest="solucions",
                           default=True)
    self.parser.add_option("--aleatori",
                           action="store_true",
                           dest="aleatori",
                           default=False)
    self.parser.add_option("--nombre-examens", dest="nombreexamens")
    self.parser.add_option("--json",
                           action="store_true",
                           dest="json",
                           default=False)
    self.parser.add_option("--ajuda",
                           action="store_true",
                           dest="ajuda",
                           default=False)
    (self.options, self.args) = self.parser.parse_args()
    self.estudiants = []
    self.problemes = []
    self.maxproblema = 0
    self.enunciats = []
    if self.options.ajuda:
      self.ajuda()

  #
  #
  #
  def ajuda(self):
    print(
      "Utilització: examen.py --examen=<fitxer> --estudiants=<fitxer> [--problemes=<enter>] [--dades=<fitxer>] [--no-solucions] [--tex-engine=pdflatex]\n"
    )
    print("Opcions:")
    print(
      "   --examen=<fitxer>                   : Fitxer LaTeX amb el model d'examen"
    )
    print(
      "   --estudiants=<fitxer>               : Fitxer amb nom:cognoms dels estudiants"
    )
    print(
      "   --problemes=<nombre|llista>         : Nombre de problemes o llista de problemes"
    )
    print(
      "   --possibles-problemes=<nombre>      : Nombre de possibles problemes")
    print(
      "                                       : S'escullen aleatòriament \"nombre\" problemes"
    )
    print("   --incompatibles=<incompatibilitats> : Llista d'incompatibiliats")
    print(
      "   --grups=<llista>                    : Llista de grups de problemes")
    print(
      "                                       : Si és possible, sortirà un problema de cad grup"
    )
    print(
      "   --dades=<fitxer>                    : Fitxer amb les dades JSON generades anteriorment"
    )
    print(
      "   --tex-engine=<programa>             : Nom del programa de LaTeX utilitzat"
    )
    print(
      "                                       : Si no s'especifica, no es generen els PDF"
    )
    print(
      "   --aleatori                          : L'ordre dels problemes serà aleatori"
    )
    print(
      "   --nombre-examens=>nombre>           : Identifica els fitxers numèricament i no per nom i cognoms"
    )
    print("                                       : Quantitat d'exàmens a fer")
    print(
      "   --no-solucions                      : No es generen els fitxers amb les solucions"
    )
    print(
      "   --json                              : Es guarden la dades dels enunciats en un fitxer json"
    )
    print("   --ajuda                             : Imprimeix aquesta ajuda")
    sys.exit(0)

  #
  #
  #
  def read_data(self):
    ex = self.options.examen
    est = self.options.estudiants
    prob = self.options.problemes
    possibles = self.options.possibles
    dades = self.options.fitxerdades
    try:
      prob = int(prob)
    except:
      prob = None
    try:
      self.possibles = int(possibles)
    except:
      if isinstance(prob, int):
        self.possibles = prob
      else:
        self.possibles = None
    if prob is None:
      prob = self.options.problemes
      try:
        l = prob.split(",")
        prob = list(map(int, l))
      except:
        prob = None
    if possibles is None and isinstance(prob, int):
      prob = list(range(prob + 1))
    try:
      self.nombreexamens = int(self.options.nombreexamens)
    except:
      self.nombreexamens = None
    regex = re.compile('^\s*#.*$', re.IGNORECASE)
    self.problemes = prob
    #
    # Comprovacions
    #
    if ex is None:
      self.ajuda()
    if prob is not None and dades is not None:
      print(
        "No es poden especificar les opcions --problemes i --dades simultàniament"
      )
      sys.exit(0)
    #
    # Incompatibilitats
    #
    self.incompatibles = None
    if self.options.incompatibles is not None:
      try:
        self.incompatibles = self.options.incompatibles.split(':')
      except:
        print("Incompatibilitats no vàlides")
        sys.exit(0)
      self.incompatibles = [p.split(',') for p in self.incompatibles]
      try:
        self.incompatibles = [list(map(int, p)) for p in self.incompatibles]
      except:
        print("Incompatibilitats no vàlides")
        sys.exit(0)
    #
    # Grups
    #
    self.grups = None
    if self.options.grups is not None:
      try:
        self.grups = self.options.grups.split(':')
      except:
        print("Grups no vàlide")
        sys.exit(0)
      self.grups = [p.split(',') for p in self.grups]
      try:
        self.grups = [list(map(int, p)) for p in self.grups]
        random.shuffle(self.grups)
      except:
        print("Grups no vàlids")
        sys.exit(0)
    if self.grups is not None and isinstance(self.problemes, list):
      print("No es poden especificar llista de problemes i grups de problemes")
      sys.exit(0)
    #
    # Enunciat de l'examen
    #
    try:
      with open(ex) as f:
        self.examen = f.read()
        f.close()
    except:
      print("Error de lectura de l'exàmen")
      sys.exit(0)
    #
    # Dades dels estudiants
    #
    if est is not None:
      try:
        with open(est) as f:
          for line in f:
            line = line.rstrip()
            if regex.match(line):
              continue
            try:
              data = line.split(':')
              self.estudiants.append({
                'nom': data[0].strip(),
                'cognoms': data[1].strip(),
                'email': data[3].strip()
              })
            except:
              continue
        f.close()
      except:
        print("Error de lectura del fitxer d'estudiants")
        sys.exit(0)
    else:
      if self.nombreexamens is None:
        self.ajuda()
    #
    # Enunciats dels problemes
    #
    if isinstance(self.problemes, int):
      self.maxproblema = self.problemes
    elif isinstance(self.problemes, list):
      self.maxproblema = max(self.problemes)
    if self.possibles is not None and self.possibles > self.maxproblema:
      self.maxproblema = self.possibles
    for i in range(1, self.maxproblema + 1):
      try:
        with open(f"p{i}.tex") as f:
          e = f.read()
          f.close()
        self.enunciats.append(e)
      except:
        print("Error en els enunciats dels problemes")
        sys.exit(0)
    #
    # Fitxer JSON en el que guardarem les dades
    #
    dades = glob.glob('examen*.json')
    count = 0
    for f in dades:
      f = f.replace('examen', '')
      f = f.replace('.json', '')
      try:
        v = int(f)
        if v > count:
          count = v
      except:
        pass
    self.count = count + 1

  #
  #
  #
  def generar_examen(self, examen, estudiant, nombre=None):
    engine = self.options.engine
    if self.options.aleatori:
      random.shuffle(examen)
    enunciats = "\n\n".join(examen)
    if self.nombreexamens is not None:
      relacio = {
        'COGNOMS': 'Cognoms',
        'NOM': 'Nom',
        'ENUNCIATS': enunciats,
        'MODEL': f"{nombre}"
      }
      filename = "examen%04d" % nombre
      filename = filename.replace(" ", "0")
    else:
      relacio = {
        'COGNOMS': estudiant['cognoms'],
        'NOM': estudiant['nom'],
        'ENUNCIATS': enunciats,
        'MODEL': f"{nombre}"
      }
      filename = f"{estudiant['cognoms']}-{estudiant['nom']}".lower().replace(
        ' ', '-')
      filename = unidecode.unidecode(filename)
    examen = self.examen
    for k, v in relacio.items():
      examen = examen.replace(k, v)
    with open(f"{filename}.tex", 'w') as f:
      f.write(examen)
      f.close()
    examen = examen.replace('NIC', 'nicsol')
    if self.options.solucions:
      with open(f"{filename}-solucio.tex", 'w') as f:
        f.write(examen)
        f.close()
    if engine is not None:
      comanda = engine.split(' ')
      comanda.extend(["-interaction=nonstopmode", f"{filename}.tex"])
      print(f"S'està executant {engine} {filename}.tex")
      p = subprocess.run(comanda,
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL,
                         shell=False)
      if p.returncode != 0:
        print(f"Hi ha un error en el fitxer {filename}.tex")
        sys.exit(0)
      if self.options.solucions:
        comanda = engine.split(' ')
        comanda.extend(["-interaction=nonstopmode", f"{filename}-solucio.tex"])
        print(f"S'està executant {engine} {filename}-solucio.tex")
        p = subprocess.run(comanda,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL,
                           shell=False)
        if p.returncode != 0:
          print(f"Hi ha un error en el fitxer {filename}-solucio.tex")
          sys.exit(0)

  #
  #
  #
  def crea_carpeta_tex(self):
    dir = os.getcwd()
    if not os.path.exists('tex'):
      try:
        os.mkdir('tex')
      except:
        print("Error en crear la carpeta tex")
        sys.exit(0)
    os.chdir('tex')
    return dir

  #
  #
  #
  def borra_fitxers(self):
    names = [
      '*.log', '*.aux', '*.asy', '*.pre', '*.fls', '*.fdb_*', '*-[0-9].pdf'
    ]
    files = []
    for n in names:
      files += glob.glob(n)
    for f in files:
      os.remove(f)

  #
  #
  #
  def generar_examens(self):
    #try:
    probs = Problemes()
    #except:
    #    print("Error en el fitxer Problemes.py")
    #    sys.exit(0)
    dir = self.crea_carpeta_tex()
    js = {}
    nombre = 1
    if self.nombreexamens is None:
      iterator = self.estudiants
    else:
      iterator = range(self.nombreexamens)
    for e in iterator:
      if isinstance(e, dict):
        js[e['email']] = {}
      else:
        js[nombre] = {}
      examen = []
      problemes = probs.problemes()
      if isinstance(self.problemes, list):
        llista = list(self.problemes)
        random.shuffle(llista)
      else:
        llista = list(range(self.possibles))
        random.shuffle(llista)
        if self.grups is None:
          llista = [x + 1 for x in llista]
          llista = self.comprova_incompatibilitats(llista)
        else:
          llista = self.llista_per_grups()
          if llista is None:
            print("Impossible generar la llista de problemes")
            sys.exit(0)
      for i in range(self.maxproblema):
        if i + 1 not in llista:
          continue
        relacio = problemes[i]()
        p = self.enunciats[i]
        for k, v in relacio.items():
          if v is not None:
            p = p.replace(k, v)
        examen.append(p)
        v = f"problema{i + 1}"
        if self.nombreexamens is None:
          js[e['email']][v] = relacio
        else:
          js[nombre][v] = relacio
      self.generar_examen(examen, e, nombre)
      nombre += 1
    self.borra_fitxers()
    os.chdir(dir)
    jsonfile = self.options.examen.replace('.tex', '')
    t = ("%3d.json" % self.count).replace(' ', '0')
    jsonfile += t
    if self.options.json:
      with open(jsonfile, 'w') as f:
        json.dump(js, f)
      f.close()

  #
  #
  #
  def es_compatible(self, p, llista):
    if self.incompatibles is None:
      return True
    t = set(llista)
    t |= {p}
    for x in self.incompatibles:
      if len(set(x) & t) > 1:
        return False
    return True

  #
  #
  #
  def comprova_incompatibilitats(self, l):
    if self.incompatibles is None:
      return l[0:self.problemes]
    finalitzat = False
    actual = l[0:self.problemes]
    others = l[self.problemes:]
    changes = 0
    while not finalitzat:
      finalitzat = True
      for p in self.incompatibles:
        common = list(set(p) & set(actual))
        if len(common) > 1:
          for x in common[1:]:
            for i, n in enumerate(actual):
              if n == x:
                try:
                  actual[i] = others[changes]
                except:
                  print("Impossible complir les incompatibilitats")
                  sys.exit(0)
                changes += 1
          finalitzat = False
          break
    return actual

  #
  #
  #
  def llista_per_grups(self):
    ng = len(self.grups)
    #
    # Hem especificat els grups --grups=1,2,3:4,5,6,7:8,9,10,11:12,13:14,15,16,17
    # i els nombre de problemes --problemes=7
    #
    llista = []
    k = 0
    count = 0
    while len(llista) < self.problemes:
      if count > ng:
        return None
      k += 1
      g = list(self.grups[k % ng])
      lg = len(g)
      if len(set(g) & set(llista)) == lg:
        count += 1
        continue
      d = list(set(g) - set(llista))
      random.shuffle(d)
      p = [x for x in d if self.es_compatible(x, llista)]
      if len(p) == 0:
        count = +1
      else:
        llista.append(p[0])
        count == 0
    return llista

  #
  #
  #
  def problemes_json(self, js):
    problemes = []
    for k, v in js.items():
      probs = [int(x.replace('problema', '')) for x in v.keys()]
      for p in probs:
        if p not in problemes:
          problemes.append(p)
    return problemes

  #
  #
  #
  def recuperar_examen(self):
    try:
      with open(self.options.fitxerdades) as f:
        js = json.load(f)
      f.close()
    except:
      print(f"Error llegint el fitxer JSON {self.options.fitxerdades}")
      sys.exit(0)
    #
    # Llegim els enunciats dels problemes
    #
    problemes = self.problemes_json(js)
    problemes.sort()
    for i in problemes:
      try:
        with open(f"p{i}.tex") as f:
          e = f.read()
        f.close()
        self.enunciats.append(e)
      except:
        print("Error en els enunciats dels problemes")
        sys.exit(0)

    dir = self.crea_carpeta_tex()
    nombre = 1
    if self.nombreexamens is None:
      iterator = self.estudiants
    else:
      iterator = range(self.nombreexamens)
    for e in iterator:
      examen = []
      if isinstance(e, dict):
        dades = js[e['email']]
      else:
        dades = js[f"{nombre}"]
      probs = [int(x.replace('problema', '')) for x in dades.keys()]
      probs.sort()
      for i in probs:
        relacio = dades[f"problema{i}"]
        p = self.enunciats[i - 1]
        for k, v in relacio.items():
          p = p.replace(k, v)
        examen.append(p)
      self.generar_examen(examen, e, nombre)
      nombre += 1
    self.borra_fitxers()
    os.chdir(dir)

  #
  #
  #
  def main(self):
    self.read_data()
    if self.options.fitxerdades is not None:
      self.recuperar_examen()
    else:
      self.generar_examens()


if __name__ == '__main__':
  main = Examen()
  main.main()
