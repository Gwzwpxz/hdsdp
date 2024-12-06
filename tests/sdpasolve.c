#include <stdio.h>
#include <stdlib.h>

int test_file_io( char *fname );
int test_sdpa_io( char *fname );
int test_bench( char *fname );
int test_primal_primal_dual_bench( char *fname );
int test_mat( char *path );

#include "tests/test_file_io.c"

int main(int argc, const char * argv[]) {
 
    if ( argc > 1 ) {
        return test_file_io((char *) argv[1]);
    } else {
        printf("HDSDP: Software for semidefinite programming\n");
    }
}
