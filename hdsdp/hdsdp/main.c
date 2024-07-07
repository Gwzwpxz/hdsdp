#include <stdio.h>
#include <stdlib.h>

int test_file_io( char *fname );
int test_solver( char *fname );
int test_mat( char *path );

int main(int argc, const char * argv[]) {
 
    if ( argc > 1 ) {
        return test_file_io((char *) argv[1]);
    } else {
        char *fname = "/Users/gaowenzhi/Desktop/potential-reduction/lps/netlib/s_afiro.mps";
        // char *fname = "/Users/gaowenzhi/Desktop/gwz/benchmark/sdplib/mcp100.dat-s";
        return test_file_io(fname);
    }
}
