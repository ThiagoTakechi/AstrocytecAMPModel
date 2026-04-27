# Compiler
CC = cl

# Compiler flags /W4
CFLAGS = /EHsc /std:c++20

# Target executable
TARGET = RunModel.exe

# Source files
SRC = Astrocyte.cpp StimuliGenerator.cpp Simulator.cpp RunModel.cpp Astrocyte.h StimuliGenerator.h Simulator.h RunModel.h 
OBJ = Astrocyte.obj StimuliGenerator.obj Simulator.obj RunModel.obj

# Default rule: build the target
all: $(TARGET)

# Compile source files into object files
.cpp.obj:
	$(CC) $(CFLAGS) /c $< /Fo$@

# Link object files into the executable
$(TARGET): $(OBJ)
	$(CC) $(OBJ) /Fe$(TARGET)

# Clean build files
clean:
	del /Q $(OBJ) $(TARGET)
