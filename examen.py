#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Filename:   examen.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2018
Disclaimer: This code is presented "as is" and it has been written to
            generate random models of exams for the subject of Linear
            at ESEIAAT
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
from optparse import OptionParser
sys.path.append('.')
from Problemes import Problemes

class Examen:
    def __init__(self):
        self.parser = OptionParser()
        self.parser.add_option("--examen",dest="examen",default=None)
        self.parser.add_option("--estudiants",dest="estudiants",default=None)
        self.parser.add_option("--problemes",dest="problemes",default=None)
        self.parser.add_option("--tex-engine",dest="engine",default=None)
        self.parser.add_option("--no-solucions",action="store_false",dest="solucions",default=True)
        self.parser.add_option("--aleatori",action="store_false",dest="aleatori",default=True)
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
        print("   --examen=<fitxer>        : Fitxer LaTeX amb el model d'examen")
        print("   --estudiants=<fitxer>    : Fitxer amb nom:cognoms dels estudiants")
        print("   --problemes=<nombre>     : Nombre de problemes")
        print("   --tex-engine=<programa>  : Nom del programa de LaTeX utilitzat")
        print("                            : Si no s'especifica, no es generen els PDF")
        print("   --no-solucions           : L'ordre dels problemes serà aleatori")
        print("   --no-solucions           : No es generen els fitxers amb les solucions")
        sys.exit(0)
    #
    #
    #
    def read_data(self):
        ex = self.options.examen
        est = self.options.estudiants
        prob = self.options.problemes
        if isinstance(prob,int):
            prob = list(range(prob + 1))
        else:
            l = prob.split(",")
            prob = list(map(int,l))
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
                        self.estudiants.append({'nom' : data[0],
                                                'cognoms' : data[1]
                                                })
                    except:
                        continue
                f.close()
        except:
            print("Error de lectura de l'exàmen")
            sys.exit(0)
        #
        # Nombre de problemes
        #
        try:
            self.problemes = prob
        except:
            print("Error en el nombre de problemes")
            sys.exit(0)
        #
        # Enunciats dels problemes
        #
        self.maxproblema = max(self.problemes)
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
        for e in self.estudiants:
            examen = []
            problemes = probs.problemes()
            for i in range(self.maxproblema):
                if i + 1 not in self.problemes:
                    continue
                relacio = problemes[i]()
                p = self.enunciats[i]
                for k,v in relacio.items():
                    p = p.replace(k,v)
                examen.append(p)
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
                comanda = f"{engine} {filename}.tex > /dev/null 2>&1 "
                print (comanda)
                os.system(comanda)
                if self.options.solucions:
                    comanda = f"{engine} {filename}-solucio.tex > /dev/null 2>&1 "
                    print (comanda)
                    os.system(comanda)
        os.system('rm -f *.log *.aux *.asy *-1.pdf *.pre *.fls *.fdb_*')
        os.chdir(dir)
    #
    #
    #
    def main(self):
        self.read_data()
        self.generar_examens()

if __name__ == '__main__':
    main = Examen()
    main.main()
