OS_RAW := $(shell uname -s)

ifneq (,$(findstring MINGW,$(OS_RAW)))
    EXE := .exe
else ifneq (,$(findstring CYGWIN,$(OS_RAW)))
    EXE := .exe
else ifneq (,$(findstring MSYS,$(OS_RAW)))
    EXE := .exe
else
    EXE :=
endif

ifeq ($(OS_RAW), Darwin)
    CCC ?= g++-14
else
    CCC ?= g++
endif

CCFLAGS = -fopenmp -O3 -Wall -std=c++2a -Wno-unused-but-set-variable
SRC = src/Main/main.cpp
HEADERS = $(wildcard src/*.h src/MH/*.h src/Problem/*.h)
TARGET = runTest$(EXE)

$(TARGET): $(SRC) $(HEADERS)
	$(CCC) $(CCFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET)

.PHONY: clean
