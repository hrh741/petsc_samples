# Include definitions of the compilers, linkers, etc.
include $(PETSC_DIR)/lib/petsc/conf/variables

target = basicOperations.exe
objects = basicOperations.o
CFLAGS = -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

all: $(target)

$(target): $(objects)
	$(CLINKER) -o $@ $^ $(PETSC_LIB)

clean:
	rm -rf $(target) *.o *.dSYM
