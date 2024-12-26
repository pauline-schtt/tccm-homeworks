# Compiler and flags
CC = gcc
CFLAGS = -I/usr/local/include -L/usr/local/lib -ltrexio

# Directories
SRC_DIR = src

# Output
TARGET = HF_and_MP2

# Source files
SRCS = $(SRC_DIR)/HF_and_MP2.c

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