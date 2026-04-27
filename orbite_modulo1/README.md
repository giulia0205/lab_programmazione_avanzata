DESCRIZIONE DEL PROGETTO:\
Si vuole studiare il comportamento di un sistema stellare formato da un numero N di pianeti per un tempo pari al periodo di rivoluzione del pianeta più esterno. Si salvano posizioni, velocità e accelerazioni ad ogni step di integrazione all'interno di un file binario per ciascun pianeta. Con uno script apposito i file binari vegnono tradotti in file di testo dai quali poi tramite Gnuplot si sono plottate le orbite del sistema in esame. In questo progetto si solo studiati i seguenti sistemi:
1. Sole + pianeti interni del sistema solare (Mercurio, Venere, Terra, Marte);
2. Trappist-1 + 5 pianeti più interni.

FILE COMPONENTI IL PROGETTO:
1. Input-Output
2. Fisica del sistema
3. Parametri del sistema
4. Makefile
5. Main: creazione e scrittura dei file binari
6. Convertitore binario-testo

ESECUZIONE DEL PROGETTO:
1. make NOME_STELLA (SUN o TRAPPIST1)
2. ./orbitbin2txt NOME_STELLA (SUN o TRAPPIST1) SKIP (rispettivamente 100 e 50)
3. gnuplot \
    3.1 SUN \
        set terminal pngcairo size 900,900  \
        set output "orbite.png" \
        set xlabel "x [AU]" \
        set ylabel "y [AU]" \
        set title "Orbite Sistema Solare Interno" \
        set size ratio -1 \
        set object 1 circle at 0,0 size 0.03 fc rgb "yellow" \
        plot "planet\_0.txt" u 2:3 w l title "Mercurio",\  \
             "planet\_1.txt" u 2:3 w l title "Venere",\  \
             "planet\_2.txt" u 2:3 w l title "Terra",\ \
             "planet\_3.txt" u 2:3 w l title "Marte" \
        exit \
    3.2 TRAPPIST-1 \
        set terminal pngcairo size 900,900 \
        set output "orbite.png" \
        set xlabel "x [AU]" \
        set ylabel "y [AU]" \
        set title "Orbite Sistema Trappist-1 Interno" \
        set size ratio -1 \
        set object 1 circle at 0,0 size 0.03 fc rgb "yellow" \
        plot "planet\_0.txt" u 2:3 w l title "b",\ \
             "planet\_1.txt" u 2:3 w l title "c",\ \
             "planet\_2.txt" u 2:3 w l title "d",\ \
             "planet\_3.txt" u 2:3 w l title "e",\ \
             "planet\_4.txt" u 2:3 w l title "f" \
        exit \
4. make clean \

FILE DI OUTPUT: \
I file contenenti posizioni, velocità e accelerazioni dei pianeti vengono eliminati con "make clean". Rimangono i due file .png raffiguranti le orbite dei sistemi rispettivamente del Sole e di Trappist-1. 