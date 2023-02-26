#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import unidecode
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--estudiants",dest="estudiants",default=None)
parser.add_option("--grup",dest="grup",default=None)
parser.add_option("--carpeta",dest="carpeta",default=None)
parser.add_option("--nomfitxer",dest="nomfitxer",default=None)
parser.add_option("--solucions",action="store_true",dest="solucions",default=False)
(options,args) = parser.parse_args()

regex = re.compile('^\s*#.*$',re.IGNORECASE)
estudiants = []
est = options.estudiants
try:
    grups = options.grup.split(',')
except:
    grups = None
carpeta = options.carpeta
if carpeta is None:
    print("S'ha d'especificar una carpeta")
    sys.exit(0)
nomfitxer = options.nomfitxer
if nomfitxer is None:
    print("S'ha d'especificar un nom de fitxer resultant")
    sys.exit(0)

try:
    f = open(est,encoding='utf8')
    for line in f:
        line = line.rstrip()
        if regex.match(line):
            continue
        try:
            data = line.split(':')
            estudiants.append({'nom' : data[0],'cognoms' : data[1],'grup' : data[4]})
        except:
            continue
    f.close()
except:
    print("Error de lectura del fitxer d'estudiants")
    sys.exit(0)

l = []
for e in estudiants:
    g = e['grup']
    if grups is not None:
        trobat = False
        for grup in grups:
            trobat = trobat or g.find(grup) == 0
        if not trobat:
            continue
    filename = f"{e['cognoms']}-{e['nom']}".lower().replace(' ','-')
    filename = unidecode.unidecode(filename)
    if options.solucions:
        filename += "-solucio.pdf"
    else:
        filename += ".pdf"
    filename = f"{carpeta}/{filename}"
    l.append(filename)
if len(l) == 0:
    print("No hi ha cap estudiants a la llista")
    sys.exit(0)
l.sort()
filenames = " ".join(l)
comanda = f"pdfunite {filenames} {nomfitxer}"
os.system(comanda)
