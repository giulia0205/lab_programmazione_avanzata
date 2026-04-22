CC = gcc                                                # compilatore che viene utilizzato
CFLAGS = -O2 -Wall                                      # flag del compilatore: -O2 ottimizzazione del codice, -Wall attiva quasi tutti i warning

SRC = src/main.c src/physics.c src/io.c 		        # elenco dei file sorgente per orbit
OBJ = $(SRC:.c=.o)                                      # creazione automatica dei file oggetto .o

CONVERTER_SRC = src/orbitbin2txt.c src/io.c             # sorgenti necessari per il convertitore
CONVERTER_OBJ = $(CONVERTER_SRC:.c=.o)                  # file oggetto del convertitore

TXT_OUT = planet_*.txt                                  # file di testo di output
BIN_OUT = planet_*.bin                                  # file binari di output

all: orbit orbitbin2txt                                 # target principale

orbit: $(OBJ)                                           # come viene generato orbit
	$(CC) $(CFLAGS) -o orbit $(OBJ) -lm

orbitbin2txt: $(CONVERTER_OBJ)                          # come viene generato orbitbin2txt
	$(CC) $(CFLAGS) -o orbitbin2txt $(CONVERTER_OBJ) -lm

%.o: %.c                                                # regola pattern per creare i file oggetto
	$(CC) $(CFLAGS) -c $< -o $@

clean:                                                  # pulizia directory
	rm -f src/*.o orbit orbitbin2txt $(TXT_OUT) $(BIN_OUT)
