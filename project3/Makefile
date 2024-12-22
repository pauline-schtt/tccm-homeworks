# Compiler and Flags
CC = gcc
CFLAGS = -Wall -lm

# Directories
SRC_DIR = src

# Output
TARGET = MD

# Source Files
SRCS = $(SRC_DIR)/main.c $(SRC_DIR)/functions.c

# Rules
all: $(TARGET)

$(TARGET):
	@echo "Building the project..."
	$(CC) $(SRCS) -o $@ $(CFLAGS)
	@echo "Done!"

.PHONY: all clean