# Compiler
CC = g++

# Compiler flags
CFLAGS = -march=native -w -DNDEBUG -Wno-narrowing -std=c++20 -O3
# CFLAGS = -march=native -Wall -Wextra -std=c++20 -O3

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

# Object files corresponding to the .cpp files
BASE_OBJ_FILES = $(filter-out $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(EXEC_NAMES))), $(CPP_FILES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o))

# Dependency files for each object file
# DEP_FILES = $(BASE_OBJ_FILES:.o=.d)

# get all dep files
DEP_FILES = $(wildcard $(OBJ_DIR)/*.d)

# Default target
all: directories $(EXECUTABLES)

# Create necessary directories
directories:
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(DOC_DIR) $(DATA_DIR)

# Compile source files into object files
$(OBJ_DIR)/%.o: %.cpp
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

# Clean
clean:
	rm -f $(BASE_OBJ_FILES) $(DEP_FILES) $(EXECUTABLES)
	rm -rf $(DOC_DIR)/*