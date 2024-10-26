# Compiler
CC = nvcc

# Compiler flags
COMPUTE_CAPABILITY=61

GENCODE=-gencode arch=compute_$(COMPUTE_CAPABILITY),code=sm_$(COMPUTE_CAPABILITY)
CFLAGS = $(GENCODE) -Xcompiler "-w,-march=native,-DNDEBUG,-Wno-narrowing" -std=c++20 -O3

# Dependency flags (generate dependency files during compilation)
DEPFLAGS = -MMD -MP

# Directories
SRC_DIRS = src ctriepp/ctriepp
OBJ_DIR = obj
BIN_DIR = bin
DOC_DIR = docs
DATA_DIR = data

# Executables names without prefix/suffix (just the target name)
EXEC_NAMES = test benchmark

.PHONY: all directories clean $(EXEC_NAMES) docs

# Generate full path to executables
EXECUTABLES = $(addprefix $(BIN_DIR)/,$(EXEC_NAMES))

# List of .cpp files
CPP_FILES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))

# List of .cu files
CU_FILES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cu))

# Object files corresponding to the .cpp files
BASE_CPP_OBJ_FILES := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(CPP_FILES)))

# Object files corresponding to the .cu files
BASE_CU_OBJ_FILES := $(patsubst %.cu, $(OBJ_DIR)/%.o, $(notdir $(CU_FILES)))

# Combine object files
BASE_OBJ_FILES = $(BASE_CPP_OBJ_FILES) $(BASE_CU_OBJ_FILES)

# Dependency files for each object file
DEP_FILES = $(wildcard $(OBJ_DIR)/*.d)

# Default target
all: directories $(EXECUTABLES)

# Create necessary directories
directories:
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(DOC_DIR) $(DATA_DIR)

# Compile source files into object files
$(OBJ_DIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: %.cu
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

# Include dependency files if they exist
-include $(DEP_FILES)

# Rule to link each executable
define MAKE_EXECUTABLE
$(BIN_DIR)/$(1): $(BASE_OBJ_FILES) $(OBJ_DIR)/$(1).o
	$(CC) $(CFLAGS) -o $$@ $$^
endef

# Generate rules for executables
$(foreach exec,$(EXEC_NAMES),$(eval $(call MAKE_EXECUTABLE,$(exec))))

$(EXEC_NAMES): %: $(BIN_DIR)/%

# Generate documentation
docs:
	doxygen Doxyfile

# Adjust pattern rules for finding source files not just in SRC_DIRS
vpath %.cpp $(SRC_DIRS)
vpath %.cu $(SRC_DIRS)

# Clean
clean:
	rm -f $(BASE_OBJ_FILES) $(DEP_FILES) $(EXECUTABLES)
	rm -rf $(DOC_DIR)/*