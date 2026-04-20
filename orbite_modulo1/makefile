CC = gcc												#compilatore che viene utilizzato
CFLAGS = -O2 -Wall										#flag del compilatore: -O2 ottimizzazione del codice, -Wall attiva quasi tutti i warning

SRC = src/main.c src/physics.c src/io.c src/utils.c		#elenco dei file sorgente
OBJ = $(SRC:.c=.o)										#creazione automatica dei file oggetto .o

CONVERTER = src/orbitbin2txt.c							#definizione del file sorgente del convertitore binario-testo

all: orbit orbitbin2txt									#target principale: quando scrivo nel terminare "make all" vegnono generati orbit e orbitbin2txt


orbit: $(OBJ)											#come viene generato orbit
	$(CC) $(CFLAGS) -o orbit $(OBJ) -lm

orbitbin2txt: $(CONVERTER)								#come viene generato orbitbin2txt
	$(CC) $(CFLAGS) -o orbitbin2txt $(CONVERTER)

%.o: %.c												#regola pattern per creare i file oggetto: $< primo file di dipendenza, $@ nome del target
	$(CC) $(CFLAGS) -c $< -o $@

clean:													#pulisce la directory: scrivendo sul terminale "make clean" vengono cancellati src/*.o, orbit, orbitbin2txt e si può ricompilare da zero
	rm -f src/*.o orbit orbitbin2txt


