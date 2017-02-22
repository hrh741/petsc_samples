static char help[] = "Basic operations on vectors.\n\n";

#include <petscvec.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    Vec            x, y, w;
    PetscReal      norm, maxVal, minVal;
    PetscInt       i, n = 20, maxInd, minInd;
    PetscErrorCode ierr;
    PetscScalar    dot;

    PetscInitialize(&argc, &argv, (char*) 0, help);

    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, n); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);

    /*
       Duplicate the vector x
    */
    ierr = VecDuplicate(x, &y); CHKERRQ(ierr);
    ierr = VecDuplicate(x, &w); CHKERRQ(ierr);

    /*
       Set the vector entries
    */
    for (i = 0; i < n; i++) {
        ierr = VecSetValue(x, i, i, INSERT_VALUES); CHKERRQ(ierr);
    }

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    ierr = VecSet(y, 5.0); CHKERRQ(ierr);

    ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = VecView(y, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    /*
       Compute the dot product of x and y.
    */
    ierr = VecDot(x, y, &dot); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "The dot product of x and y equals %g\n", (double) dot); CHKERRQ(ierr);

    /*
       Find the largest and smallest elements of a vector and their indices.
    */
    ierr = VecMax(x, &maxInd, &maxVal); CHKERRQ(ierr);
    ierr = VecMin(x, &minInd, &minVal); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "The largest element of x is %g at index %D\n", (double) maxVal, maxInd); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "The smallest element of x is %g at index %D\n", (double) minVal, minInd); CHKERRQ(ierr);

    /*
       Compute w = a*x + y.
    */
    ierr = VecWAXPY(w, 2.0, x, y); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "w = 2.0 * x + y =\n"); CHKERRQ(ierr);
    ierr = VecView(w, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    /*
       Compute the norm of w.
    */
    ierr = VecNorm(w, NORM_2, &norm); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "The vector norm of w is %g\n", (double) norm); CHKERRQ(ierr);

    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);
    ierr = VecDestroy(&w); CHKERRQ(ierr);

    PetscFinalize();

    return 0;
}
