CC = gcc
CFLAGS = -Wall -Wextra -O2 
LDFLAGS = -lm

TARGET = sph_sim

SRC = main.c io.c physics.c vicini_qsort.c vicini_grid.c
OBJ = $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

main.o: main.c physics.h vicini_qsort.h vicini_grid.h io.h
	$(CC) $(CFLAGS) -c main.c

io.o: io.c io.h physics.h
	$(CC) $(CFLAGS) -c io.c

physics.o: physics.c physics.h
	$(CC) $(CFLAGS) -c physics.c

vicini_qsort.o: vicini_qsort.c vicini_qsort.h physics.h
	$(CC) $(CFLAGS) -c vicini_qsort.c

vicini_grid.o: vicini_grid.c vicini_grid.h physics.h
	$(CC) $(CFLAGS) -c vicini_grid.c

OUTPUTS = output_qsort_*.txt output_grid_*.txt output_final_qsort.txt output_final_grid.txt

clean:
	rm -f $(OBJ) $(TARGET) $(OUTPUTS)

QSORT: $(TARGET)
	./$(TARGET) QSORT

GRID: $(TARGET)
	./$(TARGET) GRID

.PHONY: all clean QSORT GRID