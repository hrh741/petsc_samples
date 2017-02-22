static char help[] = "Creating and setting vectors.\n\n";

#include <petscmat.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    Mat            A, B, C;
    PetscInt       n = 10;
    PetscInt       i, columns[3];
    PetscErrorCode ierr;
    PetscScalar    values[3];
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

    ierr = MatCreate(comm, &B); CHKERRQ(ierr);
    ierr = MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
    ierr = MatSetFromOptions(B); CHKERRQ(ierr);

    ierr = MatCreate(comm, &C); CHKERRQ(ierr);
    ierr = MatSetSizes(C, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
    ierr = MatSetFromOptions(C); CHKERRQ(ierr);

    ierr = MatMPIAIJSetPreallocation(A, 3, PETSC_NULL, 3, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL); CHKERRQ(ierr);

    ierr = MatMPIAIJSetPreallocation(B, 3, PETSC_NULL, 3, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(B, 3, PETSC_NULL); CHKERRQ(ierr);

    //ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &B); CHKERRQ(ierr);
    //ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C); CHKERRQ(ierr);

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
//     ierr = MatView(A, viewer); CHKERRQ(ierr);

    /*
       Implement the following stencil:

       [    -4    ]
       [ 2   3 -4 ]
       [     2    ]
    */
    values[0] = 2.0;
    values[1] = 3.0;
    values[2] = -4.0;

    for (i = 1; i < n-1; i++) {
        columns[0] = i - 1;
        columns[1] = i;
        columns[2] = i + 1;

        ierr = MatSetValues(B, 1, &i, 3, columns, values, INSERT_VALUES); CHKERRQ(ierr);
    }

    i = n - 1;
    columns[0] = n - 2;
    columns[1] = n - 1;
    ierr = MatSetValues(B, 1, &i, 2, columns, values, INSERT_VALUES); CHKERRQ(ierr);

    i = 0;
    columns[0] = 0;
    columns[1] = 1;
    values[0] = 3.0;
    values[1] = -4.0;
    ierr = MatSetValues(B, 1, &i, 2, columns, values, INSERT_VALUES); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//     ierr = MatView(B, viewer); CHKERRQ(ierr);

    ierr = MatMatMult(A, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
//     ierr = MatView(C, viewer); CHKERRQ(ierr);

    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = MatDestroy(&B); CHKERRQ(ierr);
    ierr = MatDestroy(&C); CHKERRQ(ierr);

    PetscFinalize();

    return 0;
}
