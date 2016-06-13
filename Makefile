CC=nvcc
NVCCFLAGS=-arch=sm_50 -rdc=true
CFLAGS=--std=c++11 -O2 -lcufft -lcurand

SOURCE_DIR=src
INCLUDE_DIR=src/headers
BUILD_DIR=debug

SOURCES=$(wildcard $(SOURCE_DIR)/*.cu)
OBJECTS=$(patsubst $(SOURCE_DIR)/%.cu,$(BUILD_DIR)/%.o,$(SOURCES))

vpath %.cu $(SOURCE_DIR)

all: main

main: $(OBJECTS)
	$(CC) $(NVCCFLAGS) $^ $(CFLAGS) -I $(INCLUDE_DIR) -o $@

$(BUILD_DIR)/%.o: %.cu | $(BUILD_DIR)
	$(CC) $(NVCCFLAGS) -c $< $(CFLAGS) -I $(INCLUDE_DIR) -o $@

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

run: main
	./main $(wildcard ./lol_stc_csv/*-AD1-lh.csv)

clean:
	rm -rf $(BUILD_DIR) main


