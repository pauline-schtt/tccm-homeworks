# Compiler and flags
CC = gcc
CFLAGS = -Wall -lm

# Directories
SRC_DIR = src

# Output
TARGET = MD

# Source files
SRCS = $(SRC_DIR)/main.c $(SRC_DIR)/functions.c

# Rules
all: $(TARGET)

$(TARGET):
	@echo "Building the project..."
	$(CC) $(SRCS) -o $@ $(CFLAGS)
	@echo "Done!"

.PHONY: all clean

# Clean up
clean:
	@echo "Cleaning up..."
	rm -f $(TARGET)
	@echo "Done!"