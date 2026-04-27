DESCRIZIONE DEL PROGETTO: \
Si vuole realizzare una simulazione SPH (Smoothed Particle Hydrodynamics) 1D per il caso test dello shock-tube. \
Si è imposta una smoothing lenght adattiva variabile per ogni particella in funzione della densità locale; h viene calcolata dinamicamente affinchè il nucleo di interpolazione contenza un numero pari di vicini pari a 64 particelle. Il calcolo dei vicini viene implelentato secondo due modalità:
1. Ricerca QuickSort (QSORT): le particelle vengono riordinate lungo x;
2. Ricerca LinkedLists (GRID): invece di riordinare tutto l'array di particelle ad ogni passo, il dominio viene suddiviso in celle spaziali e a ciascuna di esse viene associata una lista linkata di particelle. \
L'obiettivo è di paragonare l'efficienza nel tempo di esecuzione delle due modalità. \

FILE COMPONENTI IL PROGETTO:
1. Input-Output
2. Fisica comune alle due modalità
3. Funzioni e strutture per la modalità QSORT
4. Funzioni e strutture per la modalità GRID
5. Parametri del sistema
6. Makefile
7. Simulazione SPH

ESECUZIONE DEL PROGETTO: \
Il progetto viene eseguito lanciando "make QSORT" o "make GRID"; vegnono prodotti i file di testo relativi allo stato delle particelle a ciascun timestep (implelemtazione tramite leapfrog) e viene stampato sul terminale il wall clock time del metodo selezionato. \
Lanciando "make clean" vengono eliminati tutti gli eseguibili e i file di output.

FILE DI OUTPUT: \
I file di testo relativi allo stato delle particelle ad ogni timestep vengono automaticamente eliminati con "make clean". \ 
I seguenti plot sono stati prodotti tramite apposito script python: 
1. Metodo Quicksort \
    I. Plot densità iniziale e finale (shock_tube_p.png); \
    II. item Plot pressione iniziale e finale (shock_tube_rho.png); \
    III. Plot di Velocità, Accelerazione, Densità, Pressione, Energia Interna in funzione della Posizione (5 snapshot) \
    (accelerazione_vs posizione_q.png, densita_vs_posizione_q.png, energia_interna_vs_posizione_q.png, pressione_vs_posizione_q.png, velocita_vs_posizione_q.png).
2. Metodo Linked Lists \
    I. Plot di Velocità, Accelerazione, Densità, Pressione, Energia Interna in funzione della Posizione (5 snapshot) \
    (accelerazione_vs posizione_g.png, densita_vs_posizione_g.png, energia_interna_vs_posizione_g.png, pressione_vs_posizione_g.png, velocita_vs_posizione_g.png).

CONSIDERAZIONI FINALI: \
Si sono valutati separatamente i wall clock time per le modalità Quicksort e LinkedLists, e si è calcolato uno Speedup di circa 0.3. Se ne deduce che in caso di SPH 1D è più efficiente il metodo QuickSort.