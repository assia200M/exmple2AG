# Makefile for the Genetic Algorithm project (az100)

# Compiler and compiler flags
CXX = g++

# CXXFLAGS:
# -std=c++17: Use the C++17 standard. The code uses features from C++11 and later.
# -Wall: Enable all common warnings.
# -Wextra: Enable additional warnings beyond -Wall.
# -O2: Optimization level 2, a good balance between speed and compilation time.
# -I.: Add the current directory to the include search path.
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -I.

# Linker flags (optional)
LDFLAGS =

# Target executable
TARGET = main.out

# Source files
SRCS = main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Remove command
RM = rm -f

# Program-generated files to delete with 'make clean'
PROGRAM_OUTPUT_FILES = log.txt meilleure_solution_finale.txt

# Declare phony targets
.PHONY: all clean run help

# Default target
all: $(TARGET)

# Build the final executable
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

# Compile .cpp to .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build and output files
clean:
	$(RM) $(OBJS) $(TARGET) $(PROGRAM_OUTPUT_FILES)

# Run the program
run: $(TARGET)
	./$(TARGET)

# Help message
help:
	@echo "Makefile for main Genetic Algorithm project"
	@echo ""
	@echo "Available targets:"
	@echo "  all       - Build the project (default target)."
	@echo "  $(TARGET) - Build the project explicitly."
	@echo "  run       - Build (if necessary) and run the project."
	@echo "  clean     - Remove all compiled object files, the executable,"
	@echo "              and program-generated output files (log.txt, meilleure_solution_finale.txt)."
	@echo "  help      - Show this help message."
	@echo ""
	@echo "Prerequisites for running:"
	@echo "  - A C++ compiler (e.g., g++) must be installed and in PATH."
	@echo "  - The 'json.hpp' header file must be in the current directory."
	@echo "  - The 'config1.json' configuration file must be present in the current directory."

