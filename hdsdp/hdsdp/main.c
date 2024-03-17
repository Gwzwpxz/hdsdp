#include <stdio.h>
#include <stdlib.h>

int test_file_io( char *fname );
int test_solver( char *fname );
int test_mat( char *path );

int main(int argc, const char * argv[]) {
 
    if ( argc > 1 ) {
        return test_solver(argv[1]);
    } else {
        
        char *path = "/Users/gaowenzhi/Desktop/gwz/hdsdp/gpu/A/";
        test_mat(path);
        
//        char *fname = "/Users/gaowenzhi/Desktop/gwz/benchmark/sdplib/maxG32.dat-s";
//        return test_solver(fname);
    }
}
