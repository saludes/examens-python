# Generador d'exàmens amb Python

*examens-python* és una utilitat per a generar exàmens aleatoris a partir de models de preguntes. Per a generar les dades aleatòries es fa servir el llenguatge de programació _Python_ i la seva llibreria _sympy_.

Es fa servir el _LaTeX_ per a generar els diferents models d'examen. Per tant, per fer servir aquesta utilitat necessitem tenir intal·lats els programes següents:

1. **LaTeX**, preferiblement el TeXLive, que està disponible per a Linux, Windows 10 i MacOS (https://www.tug.org/texlive/).
2. **Python**, la versió 3.7 o posterior (https://www.python.org/).
3. **SymPy**, que és una llibreria de _Python_ per a càlcul simbòlic (https://www.sympy.org).

Per fer servir aquesta utilitat des de qualsevol carpeta on tinguem els model d'examen, les diferents preguntes i el fitxer amb les dades dels estudiants, hem de procedir a

* Copiar els fitxers examen.sty i eseiaat-ma.pdf a una carpeta on els TeX pugui trobar-los. Amb Linux poden ser _$HOME/texmf/tex/latex/_ o _/usr/local/share/texmf/tex/latex/_ i en Windows 10 ...
* Copiar el fitxer Algebra.py a una carpeta on el Python pugui trobar-lo. Amb Linux pot ser _/usr/local/lib/python3.8/dist-packages/_ i en Windows 10 ...
* Copiar el fitxer examen.py en una carpeta des d'on es pugui executar des de la línia de comandes. Amb Linux pot ser _/usr/local/bin/_ i en Windpws 10 ...

## Descàrrega

Per a descarregar aquesta utilitat heu d'executar la comanda
```
~$ git clone https://github.com/rafelamer/examens-python.git
```

## Utilització

Un cop tenim en una carpeta els fitxers _examen.tex_, _Problemes.py_, _estudiants.csv_ i _p1.tex_, _p2.tex_, _p3.tex_, etc., executem la comanda
```
~$ examen.py --examen=examen.tex --estudiants=estudiants.csv --problemes=4
```
Les opcions que tenim són
```
Utilització: examen.py --examen=<fitxer> --estudiants=<fitxer> --problemes=<enter> [--no-solucions] [--tex-engine=pdflatex]
   --examen=<fitxer>        : Fitxer LaTeX amb el model d'examen
   --estudiants=<fitxer>    : Fitxer amb nom:cognoms dels estudiants
   --problemes=<nombre>     : Nombre de problemes
   --tex-engine=<programa>  : Nom del programa de LaTeX utilitzat
                            : Si no s'especifica, bo es generen els PDF
   --no-solucions           : No es generen els fitxers amb les solucions

```
