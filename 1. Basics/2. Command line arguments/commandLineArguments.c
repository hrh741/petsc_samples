static char help[] = "Command line arguments.\n\n";

#include <petscsys.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    /*
       myInt will retain its value of -1 unless myInt is provided as a command
       line argument by using the flag "-myInt".
    */
    PetscInt myInt = -1;
    PetscErrorCode ierr;

    PetscInitialize(&argc, &argv, (char*) 0, help);

    ierr = PetscOptionsGetInt(NULL, NULL, "-myInt", &myInt, NULL); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "myInt = %d\n", myInt); CHKERRQ(ierr);

    PetscFinalize();

    return 0;
}
