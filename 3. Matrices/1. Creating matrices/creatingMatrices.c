static char help[] = "Creating and setting vectors.\n\n";

#include <petscmat.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    Vec            b, u;
    Mat            A;
    PetscInt       n = 10;
    PetscInt       i, j, Ii, J, Istart, Iend, columns[3];
    PetscErrorCode ierr;
    PetscScalar    v, values[3];
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
    ierr = MatMPIAIJSetPreallocation(A, 5, PETSC_NULL, 5, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A, 5, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);

    /*
       Implement the following stencil:

       [    -1    ]
       [ -1  4 -1 ]
       [    -1    ]
    */
    values[0] = -1.0;
    values[1] = 4.0;
    values[2] = -1.0;

    for (i = 1; i < n-1; i++) {
        columns[0] = i - 1;
        columns[1] = i;
        columns[2] = i + 1;

        ierr = MatSetValues(A, 1, &i, 3, columns, values, INSERT_VALUES); CHKERRQ(ierr);
    }

    i = n - 1;
    columns[0] = n - 2;
    columns[1] = n - 1;
    ierr = MatSetValues(A, 1, &i, 2, columns, values, INSERT_VALUES); CHKERRQ(ierr);

    i = 0;
    columns[0] = 0;
    columns[1] = 1;
    values[0] = 4.0;
    values[1] = -1.0;
    ierr = MatSetValues(A, 1, &i, 2, columns, values, INSERT_VALUES); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatView(A, viewer); CHKERRQ(ierr);

    ierr = MatCreateVecs(A, &u, &b); CHKERRQ(ierr);
    ierr = VecSetFromOptions(u); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b); CHKERRQ(ierr);

    ierr = VecSet(u, 1.0); CHKERRQ(ierr);
    ierr = VecView(u, viewer); CHKERRQ(ierr);

    ierr = MatMult(A, u, b); CHKERRQ(ierr);
    ierr = VecView(b, viewer); CHKERRQ(ierr);

    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);

    PetscFinalize();

    return 0;
}
