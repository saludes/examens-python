#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Filename:   examen.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2020
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
from Problemes import Problemes

class Examen:
    def __init__(self):
        self.parser = OptionParser()
        self.parser.add_option("--examen",dest="examen",default=None)
        self.parser.add_option("--estudiants",dest="estudiants",default=None)
        self.parser.add_option("--problemes",dest="problemes",default=None)
        self.parser.add_option("--possibles-problemes",dest="possibles",default=None)
        self.parser.add_option("--tex-engine",dest="engine",default=None)
        self.parser.add_option("--no-solucions",action="store_false",dest="solucions",default=True)
        self.parser.add_option("--aleatori",action="store_true",dest="aleatori",default=False)
        self.parser.add_option("--ajuda",action="store_true",dest="ajuda",default=False)
        (self.options,self.args) = self.parser.parse_args()
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
        print("Utilització: examen.py --examen=<fitxer> --estudiants=<fitxer> --problemes=<enter> [--no-solucions] [--tex-engine=pdflatex]\n")
        print("Opcions:")
        print("   --examen=<fitxer>              : Fitxer LaTeX amb el model d'examen")
        print("   --estudiants=<fitxer>          : Fitxer amb nom:cognoms dels estudiants")
        print("   --problemes=<nombre|llista>    : Nombre de problemes o llista de problemes")
        print("   --possibles-problemes=<nombre> : Nombre de possibles problemes")
        print("                                  : S'escullen aleatòriament \"nombre\" problemes")
        print("   --tex-engine=<programa>        : Nom del programa de LaTeX utilitzat")
        print("                                  : Si no s'especifica, no es generen els PDF")
        print("   --aleatori                     : L'ordre dels problemes serà aleatori")
        print("   --no-solucions                 : No es generen els fitxers amb les solucions")
        sys.exit(0)
    #
    #
    #
    def read_data(self):
        ex = self.options.examen
        est = self.options.estudiants
        prob = self.options.problemes
        possibles = self.options.possibles
        try:
            prob = int(prob)
        except:
            prob = None
        try:
            possibles = int(possibles)
        except:
            possibles = None
        if prob is None:
            prob = self.options.problemes
            try:
                l = prob.split(",")
                prob = list(map(int,l))
            except:
                prob = None
        if possibles is None and isinstance(prob,int):
            prob = list(range(prob + 1))
        regex = re.compile('^\s*#.$',re.IGNORECASE)
        if ex is None or est is None or prob is None:
            self.ajuda()
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
        try:
            with open(est) as f:
                for line in f:
                    line = line.rstrip()
                    if regex.match(line):
                        continue
                    try:
                        data = line.split(':')
                        self.estudiants.append({'nom' : data[0].strip(),
                                                'cognoms' : data[1].strip(),
                                                'email' : data[1].strip()})
                    except:
                        continue
                f.close()
        except:
            print("Error de lectura de l'exàmen")
            sys.exit(0)
        #
        # Enunciats dels problemes
        #
        self.problemes = prob
        if isinstance(self.problemes,int):
            self.maxproblema = self.problemes
        else:
            self.maxproblema = max(self.problemes)
        self.possibles = possibles
        if self.possibles is not None and self.possibles > self.maxproblema:
            self.maxproblema = self.possibles
        for i in range(1,self.maxproblema + 1):
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
            f = f.replace('examen','')
            f = f.replace('.json','')
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
    def generar_examens(self):
        probs = Problemes()
        engine = self.options.engine
        dir = os.getcwd()
        if not os.path.exists('tex'):
            try:
                os.mkdir('tex')
            except:
                print("Error en crear la carpeta tex")
                sys.exit(0)
        os.chdir('tex')

        js = {}
        for e in self.estudiants:
            js[e['email']] = {}
            examen = []
            problemes = probs.problemes()
            if isinstance(self.problemes,list):
                llista = list(self.problemes)
            else:
                llista = random.sample(range(self.possibles),self.problemes)
                llista.sort()
                llista = [x+1 for x in llista]
            for i in range(self.maxproblema):
                if i + 1 not in llista:
                    continue
                relacio = problemes[i]()
                p = self.enunciats[i]
                for k,v in relacio.items():
                    p = p.replace(k,v)
                examen.append(p)
                v = f"problema{i}"
                js[e['email']][v] = relacio
            if self.options.aleatori:
                random.shuffle(examen)
            enunciats = "\n\n".join(examen)
            relacio = {'COGNOMS' : e['cognoms'], 'NOM' : e['nom'], 'ENUNCIATS' : enunciats}
            examen = self.examen
            for k,v in relacio.items():
                examen = examen.replace(k,v)
            filename = f"{e['cognoms']}-{e['nom']}".lower().replace(' ','-')
            filename = unidecode.unidecode(filename)
            with open(f"{filename}.tex",'w') as f:
                f.write(examen)
                f.close()
            examen = examen.replace('NIC','nicsol')
            if self.options.solucions:
                with open(f"{filename}-solucio.tex",'w') as f:
                    f.write(examen)
                    f.close()
            if engine is not None:
                comanda = [f"{engine}","-interaction=nonstopmode", f"{filename}.tex"]
                print (f"S'està executant {engine} {filename}.tex")
                p = subprocess.run(comanda,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                if p.returncode != 0:
                    print (f"Hi ha un error en el fitxer {filename}.tex")
                    sys.exit(0)
                if self.options.solucions:
                    comanda = [f"{engine}","-interaction=nonstopmode", f"{filename}-solucio.tex"]
                    print (f"S'està executant {engine} {filename}-solucio.tex")
                    p = subprocess.run(comanda,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                    if p.returncode != 0:
                        print (f"Hi ha un error en el fitxer {filename}-solucio.tex")
                        sys.exit(0)
        names = ['*.log','*.aux','*.asy','*-1.pdf','*.pre','*.fls','*.fdb_*']
        files = []
        for n in names:
            files += glob.glob(n)
        for f in files:
            os.remove(f)
        os.chdir(dir)
        jsonfile = self.options.examen.replace('.tex','')
        t = ("%3d.json" % self.count).replace(' ','0')
        jsonfile += t
        with open(jsonfile,'w') as f:
            json.dump(js,f)
        f.close()
    #
    #
    #
    def recuperar_examen(self):
        pass
    #
    #
    #
    def main(self):
        self.read_data()
        self.generar_examens()

if __name__ == '__main__':
    main = Examen()
    main.main()
