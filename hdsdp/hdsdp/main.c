#include <stdio.h>
#include <stdlib.h>

int test_file_io( char *fname );
int test_sdpa_io( char *fname );
int test_bench( char *fname );
int test_primal_primal_dual_bench( char *fname );
int test_mat( char *path );


int main(int argc, const char * argv[]) {
 
    if ( argc > 1 ) {
        return test_primal_primal_dual_bench((char *) argv[1]);
    } else {
        char *fname = "/Users/gaowenzhi/Desktop/potential-reduction/lps/miplib/s_bnatt400.mps";
        // char *fname = "/Users/gaowenzhi/Desktop/gwz/benchmark/sdplib/mcp100.dat-s";
        return test_primal_primal_dual_bench(fname);
    }
}
