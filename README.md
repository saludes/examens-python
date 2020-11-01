# Generador d'exàmens amb Python

*examens-python* és una utilitat per a generar exàmens aleatoris a partir de models de preguntes. Per a generar les dades aleatòries es fa servir el llenguatge de programació _Python_ i la seva llibreria _sympy_.

Es fa servir el _LaTeX_ per a generar els diferents models d'examen. Per tant, per fer servir aquesta utilitat necessitem tenir intal·lats els programes següents:

1. **LaTeX**, preferiblement el TeXLive, que està disponible per a Linux, Windows 10 i MacOS (https://www.tug.org/texlive/).
2. **Python**, la versió 3.7 o posterior (https://www.python.org/).
3. **SymPy**, que és una llibreria de _Python_ per a càlcul simbòlic (https://www.sympy.org).

Per fer servir aquesta utilitat des de qualsevol carpeta on tinguem els model d'examen, les diferents preguntes i el fitxer amb les dades dels estudiants, hem de procedir a

* Copiar els fitxers examen.sty i eseiaat-ma.pdf a una carpeta on els TeX pugui trobar-los. Amb Linux poden ser _$HOME/texmf/tex/latex/_ o _/usr/local/share/texmf/tex/latex/_ i en Windows 10 ...
* Copiar el fitxer Algebra.py a una carpeta on el Python pugui trobar-lo. Amb Linux pot ser _/usr/local/lib/python3.8/dist-packages/_ i en Windows 10 ...
* Copiar els fitxers _examen.py_ i _credentials.py_  a una carpeta des d'on es pugui executar des de la línia de comandes. Amb Linux pot ser _/usr/local/bin/_ i en Windpws 10 ...

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
                            : Si no s'especifica, no es generen els PDF
   --no-solucions           : No es generen els fitxers amb les solucions

```

## Enviament de correus amb Google API

### Credencials

Per poder enviar a cada estudiant un correu electrònic amb el seu examen adjunt, hem de tenir instal·lada l'API del Google. En Linux la podem instal·lar amb la comanda
```
~# apt install python3-testresources
~# pip3 install --upgrade google-api-python-client google-auth-httplib2 google-auth-oauthlib
```
A continuació, hem d'activar al nostre compte de GMail la utilització d'aquesta API. Per això, accedim al WEB https://developers.google.com/gmail/api/quickstart/python i cliquem al botó _Enable the Gmail API_. Aleshores, ens demanarà el nostre correu electrònic i contrasenya de GMail. Si cal, tornem a clicar a _Enable the Gmail API_ i se'ns obrirà una finestra on ens demana "Enable GMail API: Enter new project name", hi posem, per exemple, "UPC" i cliquem a "NEXT". A continuació a "Configure your OAuth client" hi posem "Desktop App" i cliquem a "CREATE" i després a "DOWNLOAD CLIENT CONFIGURATION".

Com a resultat d'aquesta operació ens haurem descarregat el fitxer _credentials.json_ i el guardem a la carpeta $HOME/credentials/

El segon pas consisteix en executar el programa _credentials.py_ des d'un terminal
```
~$ credentials.py
```
Se'ns obre una finestra del navegador i ens torna a demanar el nostre correu electrònic i contrasenya de GMail. Una vegada completat tindrem el fitxer _token.pickle_ a la carpeta $HOME/credentials/

**Advertència**: Si fem servir un ordinador compartit, hem se tenir present que qualsevol persona que tingui accés als fitxers _credentials.py_ i  _token.pickle_ pot accedir al nostre correu de GMail. És una bona idea que, una vegada utilitzats, els encriptem, per exemple amb les comandes
```
~$ openssl enc -pbkdf2 -aes-256-cbc -in credentials.py -out credentials.py.data
~$ shred -n 5 credentials.py
~$ openssl enc -pbkdf2 -aes-256-cbc -in token.pickle -out token.pickle.data
~$ shred -n 5 token.pickle
```

Quan els vulguem tornar a utilitzar, els hem de desencriptar:
```
~$ openssl enc -aes-256-cbc -d -in credentials.py.data -out credentials.py
~$ openssl enc -aes-256-cbc -d -in token.pickle.data -out token.pickle
```

### Enviament dels correus
