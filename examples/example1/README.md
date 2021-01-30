## Primer exemple

En aquesta carpeta podem veure l'estructura típica per a elaborar un examen amb 4 exercicis de manera que les dades de cada exercici es generen de manera aleatòria per a cada estudiant. Els fitxers necessaris són els següents:

1. *estudiants.csv*: fitxer de text amb un estudiant a cada línia amb les dades següents separades per dos punts:
```
NOM:GOGNOMS:CORREU ELECTRÒNIC
```  
2. *examen.tex*: conté les dades bàsiques per a generar un examen (fitxer TeX) per a cada estudiant.
3. *p<nombre>.tex*: un fitxer TeX per a cada problema. Ha de contenir l'enunciat del problema i com s'han d'incloure les dades i solucions per a cada estudiant.
4. *Problemes.py*: fitxer en _Python_ que defineix la classe *Problemes* i una funció per a cada un dels problemes. Aquesta funció ha de retornar un diccionari que té com a claus les *paraules* del fitxer TeX que volem substituir i com a *valors*, les cadenes (strings) per les quals s'han de substituir.

Per generar els fitxers TeX i els PDF amb els exàmens, hem d'executar la comanda
```
~$ examen.py --examen=examen.tex --estudiants=estudiants.csv --problemes=4 --tex-engine=pdflatex
```
Aleshores, es crea la carpeta *tex* os s'hi guarden tots els fitxers. 
