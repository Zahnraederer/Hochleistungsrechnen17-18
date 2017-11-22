#!/bin/bash

#Shell script, dass die Messung der Laufzeit der Programme automatisiert und Die Daten für messung standardmäßig in die Datei data.csv schreibt
# Der Parameter Steht für die Ausgabedatei, die erzeugt werden soll


#Anzahl der Messungen je Parameterwert:
MAXITER=3

#Programm, das ausgeführt wird:
programm='./partdiff-posix'

#standard Ausgabe-Datei:
output="data.csv"

#mindestens 1 Argument übergeben
if (($# <= 1)) ; then
output=$1
fi


echo "beginne Messung"
#setup output:
>$output echo "Number of Threads	Runtime in seconds"
#Datei wird dabei überschrieben!

# für jede anzahl an threads
for (( num_threads=1; num_threads<=12 ; num_threads++ ))
do
#damit man sehen kann, wie weit die Messung ist:
echo "num_threads = $num_threads"

#Argumente für die Messung
args="$num_threads 2 512 2 2 700"

# für jede der MAXITER messungen
for (( i=1; i<=MAXITER ; i++ ))
do
# Nehme Vom Output die 9.Zeile = Berechnungszeit: 999.999999 s
# davon nur das zwishen dem 5. und 6. Leerzeichen = 999.999999
Berechnungszeit=$($programm $args  | sed '9q;d' |cut -d ' ' -f 5 )

# in output-Datei schreiben:
>>$output echo "$num_threads	$Berechnungszeit"

done
# end for i

done
# end for num_threads

echo "fertig"

