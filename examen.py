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
        (self.options,self.args) = self.parser.parse_args()
        self.estudiants = []
        self.problemes = 0
        self.enunciats = []
    #
    #
    #
    def read_data(self):
        ex = self.options.examen
        est = self.options.estudiants
        prob = self.options.problemes
        if ex is None or est is None or prob is None:
            print("Error en les dades proporcionades")
            print("Utilització: examen.py --examen=<fitxer> --estudiants=<fitxer> --problemes=<enter>")
            sys.exit('0')
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
                    data = line.split(':')
                    try:
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
            self.problemes = int(prob)
        except:
            print("Error en el nombre de problemes")
            sys.exit(0)
        #
        # Enunciats dels problemes
        #
        for i in range(1,self.problemes + 1):
            try:
                with open(f"p{i}.tex") as f:
                    self.enunciats.append(f.read())
                    f.close()
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
            os.mkdir('tex')
        os.chdir('tex')
        for e in self.estudiants:
            enunciats = ""
            problemes = probs.problemes()
            for i in range(self.problemes):
                relacio = problemes[i]()
                p = self.enunciats[i]
                for k,v in relacio.items():
                    p = p.replace(k,v)
                enunciats += "\n"  + p

            relacio = {'COGNOMS' : e['cognoms'], 'NOM' : e['nom'], 'ENUNCIATS' : enunciats}
            examen = self.examen
            if not os.path.exists('tex'):
                try:
                    os.mkdir('tex')
                except:
                    print("Error en la creació de la carpeta tex")
            for k,v in relacio.items():
                examen = examen.replace(k,v)
            filename = f"{e['cognoms']}-{e['nom']}".lower().replace(' ','-')
            fileneme =
            with open(f"{filename}.tex",'w') as f:
                f.write(examen)
                f.close()
            examen = examen.replace('NIC','sol')
            with open(f"{filename}-solucio.tex",'w') as f:
                f.write(examen)
                f.close()
            if engine is not None:
                comanda = f"{engine} {filename}.tex > /dev/null 2>&1 "
                print (comanda)
                os.system(comanda)
                comanda = f"{engine} {filename}-solucio.tex > /dev/null 2>&1 "
                print (comanda)
                os.system(comanda)
        os.system('rm -f *.log *.aux')
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
