Um die Fluideigenschaften von CO2 (Dichte, Viskosit�t, Enthalpie) als 
konstitutive Beziehungen in MUFTE_UG zu erhalten, m�ssen die folgenden
Schritte durchgef�hrt werden:

1. Ermittlung des f�r das gestellte Problem notwenigen Druck- und
   Temperaturbereiches und der n�tigen Aufl�sung der St�tzstellen f�r
   Druck und Temperatur.
2. Anwendung der Beziehung von Span & Wagner (in diesem Verzeichnis bzw. in
   pml_CO2): 
   das Programm "extractproperties" muss entsprechend des p- und T-Bereiches und der 
   Aufl�sung angepasst werden. Die Ausf�hrung von "extractproperties" erzeugt den 
   Datensatz (datensatz.dat) mit Dichte, Enthalpie, und Viskosit�t.  
3. F�r die Erzeugung der c-Funktionen aus dem Datensatz werden PERL-Scripts
   verwendet. Das Skript data.prl nimmt den Datensatz datensatz.dat und 
   erzeugt die Ausgabe. Im Pearl-Skript muss lediglich der Parameter $which
   angepasst werden (1=density, 2=enthalpy, 3=solubility) um die entsprechende 
   Funktion zu schreiben.
4. �berpr�fen ob Datenarrays richtig geschriebn wurden (manchmal werden 2 Zahlen 
   zusammen geschrieben, obwohl ein Komma+Lehrzeichen dazwischen geh�rt ??) 
5. Ein Test der neuen Funktion kann im Verzeichnis ./TEST durchgef�hrt werden.
   Hierzu muss die erzeugte Funktion density.c an die Datei testconstrel.c 
   angeh�ngt werden (cat density.c >> TEST/testconstrel.c). Vorsicht: ist 
   die alte Funktion von density.c noch in testconstrel.c enthalten?
6. Die Ausgabe density.c kann nun an unser constrel-File angeh�ngt werden, z.B.:
	cat density.c >> pml/constrel_CO2.c
7. Die alte Funktion im constrel muss dann nat�rlich gel�scht werden. Falls
   es sich um eine zus�tzliche Funktion handelt, muss diese im constrel_CO2.h
   deklariert werden.


Dieses Vorgehen erm�glicht es, den f�r das gestellte Problem notwendigen 
Druck- und Temperaturbereich vorzugeben, schnell in das Programm zu
implementieren und somit ein aufwendiges "mitschleppen" zus�tzlicher Werte,
die nicht ben�tigt werden, zu vermeiden. Dasselbe gilt auch f�r die Wahl der
St�tzstellen.


Hinweise
--------

- Die Datei datensatz.dat muss das folgende Format haben, da andernfalls
  Probleme bei der Umformatierung im Pearl-Script auftreten k�nnen:
	
	Druck 		 Temperatur	 Dichte 	  Enthalpie	Viskosit�t	  

	1.000000e+05     2.830000e+02    1.881498e+00    -1.370482e+01  1.0E-5
	1.000000e+05     2.840000e+02    1.874743e+00    -1.286816e+01  2.0E-5

  Die Pearl-Scripts geben keine Warnungen oder Fehlermeldungen aus, wenn
  bei der Sortierung etwas schiefgehen sollte, deswegen sollte das Ergebnis
  auch nochmal �berpr�ft werden.


Historie und Fehler
-------------------

- Urspr�nglich wollte ich die Aufl�sung der p-, T-St�tzstellen f�r den
  kritischen Bereich feiner w�hlen als in den Gebieten, in denen nicht 
  so viel passiert. Da bei der Interpolation zwischen zwei Dr�cken und
  zwei Temperaturen gemittelt wird, kam es an dem �bergang zwischen den
  unterschiedlichen Funktionen zu Unstetigkeiten. Deshalb entschied ich
  mich erstmal f�r eine gleichbleibende Aufl�sung �ber den gesamten
  p-, T-Bereich.

- Urspr�nglich habe ich die Felder folgenderma�en deklariert und definiert: 
  
	Deklaration:

	DOUBLE T[41];
	DOUBLE p[110];
	DOUBLE rho[41][110];

	Definition:

	p[0] = 1.000000e+05;
	T[0] = 2.830000e+02;
	rho[0][0] = 1.881498e+00;

  Aus programmiertechnischen Gr�nden habe ich es dann so versucht:
	
	DOUBLE T[41]={273.15,283.15,...};
	DOUBLE p[110]={1.0E5,2.0E5,...};
	DOUBLE rho[41][110]={1.88e+00, 3.78e+00, ..}{...};

