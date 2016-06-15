CC=nvcc
CFLAGS=--std=c++11 -O2 -lcufft -lcurand
NVCCFLAGS=-arch=sm_50 -rdc=true
NVCCLDFLAGS=-arch=sm_50 -lcufft -lcurand

SOURCE_DIR=src
INCLUDE_DIR=src/headers
BUILD_DIR=debug

SOURCES=$(wildcard $(SOURCE_DIR)/*.cu)
OBJECTS=$(patsubst $(SOURCE_DIR)/%.cu,$(BUILD_DIR)/%.o,$(SOURCES))

vpath %.cu $(SOURCE_DIR)

all: main

main: $(BUILD_DIR)/fmri_corr_coef.o
	$(CC) $(NVCCLDFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: %.cu | $(BUILD_DIR)
	$(CC) $(NVCCFLAGS) -c $< $(CFLAGS) -I $(INCLUDE_DIR) -o $@

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

run: main
	./main

test: main
	./main > /tmp/fmri_corr_coef.m
	matlab -r "run('/tmp/fmri_corr_coef.m');"
	grep "%" /tmp/fmri_corr_coef.m


clean:
	rm -rf $(BUILD_DIR) main


