CC = gcc
CFLAGS = -Wall
LDFLAGS = -lnetcdf -lm

TARGET = MPCCmodel
SRC = MPCCmodel.c ascii.c repwvl_thermal.c repwvl_solar.c
OBJ = $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(TARGET)