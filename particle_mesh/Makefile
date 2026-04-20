
# Compilatore
CC = gcc

# Flags di compilazione
CFLAGS = -Wall -Wextra -O2

# Librerie
LIBS = -lfftw3 -lm

# Schema di assegnazione densità/forza:
# usa make SCHEME=NGP
# oppure make SCHEME=CIC
# oppure make SCHEME=TSC
SCHEME ?= CIC

# Macro passata al compilatore
DEFS = -D$(SCHEME)

# Eseguibili finali
TARGETS = condizioni_iniziali PM1D test_pot

# File di output numerici
OUTPUTS = i_c.txt i_c.bin hist_ic.txt pot_test.txt \
          density_final.txt potential_final.txt force_final.txt final_particles.txt \
          density_step_*.txt \
		  particles_step_*.txt

all: $(TARGETS)

condizioni_iniziali: condizioni_iniziali.o io.o physics.o
	$(CC) $(CFLAGS) $(DEFS) -o $@ $^ $(LIBS)

PM1D: PM1D.o io.o physics.o
	$(CC) $(CFLAGS) $(DEFS) -o $@ $^ $(LIBS)

test_pot: test_pot.o io.o physics.o
	$(CC) $(CFLAGS) $(DEFS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(DEFS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS) $(OUTPUTS)

.PHONY: all clean