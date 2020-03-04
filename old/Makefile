CC = gcc
CFLAGS = -g -O2

TARGET = oplib
OBJS = oplib.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS)

clean:
	$(RM) $(TARGET) $(OBJS)
