static char help[] = "Creating and setting vectors.\n\n";

#include <petscvec.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    Vec            x;
    PetscErrorCode ierr;
    PetscInt       i, n = 20;

    PetscInitialize(&argc, &argv, (char*) 0, help);

    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL); CHKERRQ(ierr);

    /*
       In the function VecSetSize(Vec v, int m, int M), m indicates the number
       of components to store on the local process, and M is the total number
       of vector components. Either the local or global dimension, but not
       both can be set to PETSC_DECIDE to indicate that PETSc should determine
       it.
    */
    ierr = VecCreate(PETSC_COMM_WORLD, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, n); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);

    /*
       Set all entries to a constant value.
    */
    // ierr = VecZeroEntries(x); CHKERRQ(ierr);
    // ierr = VecSet(x, 3.14); CHKERRQ(ierr);

    /*
       Set individual elements.
    */
    for (i = 0; i < n; i++) {
        ierr = VecSetValue(x, i, i, INSERT_VALUES); CHKERRQ(ierr);
    }

    /*
       VecSetValue() and VecSetValues() are purely local functions with no
       inter-process communication. Before using the vector, call the assembly
       function pair to exchange values between processors.
    */
    VecAssemblyBegin(x);
    /*
       Optional operations not involving x can be done while MPI messages
       are in transition. This allows overlapping communication and
       computation.
    */
    VecAssemblyEnd(x);

    /*
       Print x to stdout.
    */
    ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    ierr = VecDestroy(&x); CHKERRQ(ierr);

    PetscFinalize();

    return 0;
}
