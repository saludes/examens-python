# Generador d'exàmens amb Python

*examens-python* és una utilitat per a generar exàmens aleatoris a partir de models de preguntes. Per a generar les dades aleatòries es fa servir el llenguatge de programació _Python_ i la seva llibreria _sympy_.

Es fa servir el _LaTeX_ per a generar els diferents models d'examen. Per tant, per fer servir aquesta utilitat necessitem tenir intal·lats els programes següents:

1. **LaTeX**, preferiblement el TeXLive, que està disponible per a Linux, Windows 10 i MacOS (https://www.tug.org/texlive/).
2. **Python**, la versió 3.7 o posterior (https://www.python.org/).
3. **SymPy**, que és una llibreria de _Python_ per a càlcul simbòlic (https://www.sympy.org).
4. **Python unidecode**, que serveix, entre altres coses, per treure els accents de qualsevol text. Es fa servir per treure els accents del nom i cognoms dels estudiants ja que poden donar problemes a l'hora d'ajuntar fitxers amb nom que conté aquests nom i cognoms.

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
Cada línia del fitxer _estudiants.csv_ ha de ser de la forma
```
NOM:COGNOMS:CORREU ELECTRÒNIC
```
Les opcions que tenim són
```
Utilització: examen.py --examen=<fitxer> --estudiants=<fitxer> --problemes=<enter> [--no-solucions] [--tex-engine=pdflatex]
   --examen=<fitxer>                 : Fitxer LaTeX amb el model d'examen
   --estudiants=<fitxer>             : Fitxer amb nom:cognoms dels estudiants
   --problemes=<nombre o llista>     : Nombre de problemes o llista de problemes
   --possibles-problemes=<nombre>    : Nombre de possibles problemes
                                     : Si problemes és menor que possibles-problemes, s'escullen a l'atzar per a
                                     : cada estudiants
   --dades=<fitxer>                  : Fitxer amb les dades generades anteriorment                                  
   --tex-engine=<programa>           : Nom del programa de LaTeX utilitzat
                                     : Si no s'especifica, no es generen els PDF
   --aleatori                        : L'ordre dels problemes serà aleatori
   --no-solucions                    : No es generen els fitxers amb les solucions
```

Cada vegada que es fa una col·lecció d'exàmens aleatoris, es guarden les dades aleatòries en un fitxer JSON, que en el cas anterior seria _examen001.json_. Si volem tornar a generar els exàmens amb les mateixes dades, haurem d'executar
```
~$ examen.py --examen=examen.tex --estudiants=estudiants.csv --dades=examen001.json
```

Si en l'examen hi volem incloure gràfics generats amb l'asymptote (https://asymptote.sourceforge.io/), hem d'instal·lar-lo i també convé instal·lar el _latexmk_:
 ```
~# apt install asymptote
~# apt install latexmk
 ```
i copiar el fitxer _LatexMk_ a la carpeta _/etc_. Aleshores, per generar els PDF executarem la comanda  _examen.py_ amb l'opció _--tex-engine='latexmk -pdf'_.

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

El fitxer _token.pickle_ conté una autorització per llegir i enviar correus que té una validesa limitada. Aproximadament al cap d'una hora d'haver-lo descarregat. Per aquest motiu, el programa _enviar-correus.py_, si és necessari actualitza aquest fitxer o torna a descarregar-lo si no el pot actualitzar i per això ha d'obrir un navegador.

Per tant és millor executar el programa _enviar-correus.py_ des d'una sessió en la que es pugui obrir un navegador.

**Advertència**: Si fem servir un ordinador compartit, hem se tenir present que qualsevol persona que tingui accés al fitxer _token.pickle_, pot accedir al nostre correu de GMail. És una bona idea que, una vegada utilitzats, els encriptem, per exemple amb les comandes
```
~$ openssl enc -pbkdf2 -aes-256-cbc -in credentials.json -out credentials.json.data
~$ shred -n 64 credentials.json
~$ openssl enc -pbkdf2 -aes-256-cbc -in token.pickle -out token.pickle.data
~$ shred -n 64 token.pickle
```
Quan els vulguem tornar a utilitzar, els hem de desencriptar:
```
~$ openssl enc -aes-256-cbc -d -in credentials.json.data -out credentials.json
~$ openssl enc -aes-256-cbc -d -in token.pickle.data -out token.pickle
```

### Enviament dels correus

Un cop generat els fitxers PDF amb els enunciats dels exàmens, es poden enviar per correu electrònic als estudiants amb la comanda
```
~$ enviar-examens.py --estudiants=estudiants.csv --subject="Examen Final" --sender=rafel.amer@upc.edu --message=correu.txt
```

En el fitxer _correu.txt_ hi tindrem el cos del missatge que volem enviar a cada estudiant. Per exemple, hi podem escriure les instruccions per a la realització de l'examen

Quan hagi finalitzat l'examen, podem enviar les solucions a cada estudiant afegint l'opció _--solucions_ a la comanda anterior.
```
~$ enviar-examens.py --estudiants=estudiants.csv --subject="Solucions de l'Examen Final" --sender=rafel.amer@upc.edu --message=correu2.txt --solucions
```
