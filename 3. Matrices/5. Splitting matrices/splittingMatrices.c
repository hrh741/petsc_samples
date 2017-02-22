static char help[] = "Creating and setting vectors.\n\n";

#include <petscmat.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    Mat            A, L, U, D;
    PetscInt       n = 10;
    PetscInt       i, j, Istart, Iend;
    PetscErrorCode ierr;
    MPI_Comm       comm;
    PetscViewer    viewer;

    PetscInitialize(&argc, &argv, (char*) 0, help);
    comm = PETSC_COMM_WORLD;
    viewer = PETSC_VIEWER_STDOUT_WORLD;

    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "(m, n) = (%d, %d)\n", n, n); CHKERRQ(ierr);

    ierr = MatCreate(comm, &A); CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
    ierr = MatSetFromOptions(A); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A, n, PETSC_NULL, n, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A, n, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i < j) {
                ierr = MatSetValue(A, i, j, 3, INSERT_VALUES); CHKERRQ(ierr);
            } else if (i > j) {
                ierr = MatSetValue(A, i, j, 1, INSERT_VALUES); CHKERRQ(ierr);
            } else {
                ierr = MatSetValue(A, i, j, 2, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatView(A, viewer); CHKERRQ(ierr);

    ierr = MatDestroy(&A); CHKERRQ(ierr);

    PetscFinalize();

    return 0;
}
