static char help[] = "Hello, world!\n\n";

#include <petscsys.h>

/*
   Setting __FUNCT__ is a work-around for C compilers that do not support the
   __FUNCTION__ standard macro. If supported, the __FUNCTION__ macro is
   automatically set by the compiler to the name of the function in which the
   macro is invoked. If defined, this macro is used in the PETSc error
   handlers to provide a complete traceback of routine names.
*/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    PetscMPIInt rank;
    PetscMPIInt size;
    PetscErrorCode ierr;

    PetscInitialize(&argc, &argv, (char*) 0, help);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_SELF, "Hello, world! From process %d of %d\n", rank, size); CHKERRQ(ierr);

    /*
       Notice that the output from PetscPrintf() follows an arbitrary
       sequence. Uncomment the following two lines for a synchronized output
       sequence.
    */
//     ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Hello, world! From process %d of %d\n", rank, size); CHKERRQ(ierr);
//     ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL); CHKERRQ(ierr);

    PetscFinalize();

    return 0;
}
