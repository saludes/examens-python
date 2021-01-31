## Segon exemple

La principal diferència entre aquest exemple i l'anterior és que ara, en cada examen s'inclou una figura diferent generada amb el programa _Asymptote_. Per això, en el fitxer _examen.tex_ hem d'incloure unes línies amb comandes del paquet _asymptote_:
```
\begin{asydef}
import eixos;
import graph;
usepackage("amsfonts");
usepackage("amssymb");
usepackage("times");
usepackage("txfonts");
\end{asydef}
```
_eixos.asy_ és un fitxer amb funcions predefinides per generar quadrícules i gràfics de còniques. S'ha de copiar a una carpeta on el programa _asy_ el pugui trobar, per exemple, _$HOME/.asy_ o _/usr/share/asymptote_.

Per generar els fitxers TeX i els PDF amb els exàmens, hem d'executar la comanda
```
~$ examen.py --examen=examen.tex --estudiants=estudiants.csv --problemes=4 --tex-engine='latexmk -pdf'
```
