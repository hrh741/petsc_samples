static char help[] = "Solves a linear system with KSP.\n\n";

#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    Vec            x, b, u;
    Mat            A;
    KSP            ksp;
    PetscReal      norm;
    PetscInt       m = 4, n = 2;
    PetscInt       i, j, Ii, J, Istart, Iend, its;
    PetscErrorCode ierr;
    PetscScalar    v;
    MPI_Comm       comm;
    PetscViewer    viewer;

    PetscInitialize(&argc, &argv, (char*) 0, help);

    comm = PETSC_COMM_WORLD;
    viewer = PETSC_VIEWER_STDOUT_WORLD;

    ierr = PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "(m, n) = (%d, %d)\n", m, n); CHKERRQ(ierr);

    ierr = MatCreate(comm, &A); CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m*n, m*n); CHKERRQ(ierr);
    ierr = MatSetFromOptions(A); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A, 5, PETSC_NULL, 5, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A, 5, PETSC_NULL); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);

    for (Ii = Istart; Ii < Iend; Ii++) {
        v = -1.0;
        i = Ii / n;
        j = Ii - i * n;

        if (i > 0) {
            J = Ii - n;
            ierr = MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES); CHKERRQ(ierr);
        }

        if (i < m - 1) {
            J = Ii + n;
            ierr = MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES); CHKERRQ(ierr);
        }

        if (j > 0) {
            J = Ii - 1;
            ierr = MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES); CHKERRQ(ierr);
        }

        if (j < n - 1) {
            J = Ii + 1;
            ierr = MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES); CHKERRQ(ierr);
        }

        v = 4.0;

        ierr = MatSetValues(A, 1, &Ii, 1, &Ii, &v, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = MatView(A, viewer); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "----------\n"); CHKERRQ(ierr);

    ierr = MatCreateVecs(A, &u, &b); CHKERRQ(ierr);
    ierr = VecSetFromOptions(u); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b); CHKERRQ(ierr);

    ierr = VecSet(u, 1.0); CHKERRQ(ierr);
    ierr = VecView(u, viewer); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "----------\n"); CHKERRQ(ierr);

    ierr = MatMult(A, u, b); CHKERRQ(ierr);
    ierr = VecView(b, viewer); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "----------\n"); CHKERRQ(ierr);

    ierr = VecDuplicate(u, &x); CHKERRQ(ierr);

    ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-2/((m + 1) * (n + 1)), PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPView(ksp, viewer); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "----------\n"); CHKERRQ(ierr);

    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
    ierr = VecView(x, viewer); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "----------\n"); CHKERRQ(ierr);

    ierr = VecAXPY(x, -1.0, u); CHKERRQ(ierr);
    ierr = VecNorm(x, NORM_2, &norm); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "Norm of error %g\n", norm); CHKERRQ(ierr);

    ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "Iterations %d\n", its); CHKERRQ(ierr);

    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

    ierr = PetscFinalize();

    return 0;
}
