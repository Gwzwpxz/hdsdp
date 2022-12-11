clear; clc;

src_linsys = ["linsys.c"; "qdldl.c"; "vec_mat.c"];
src_pot = ["pot_constr_mat.c"; "pot_lanczos.c"; "pot_objfunc.c"; "cone_filter.c";...
                 "pot_solver.c"; "pot_utils.c"; "pot_vector.c" ];
src_interface = ["another_lp_solver.c";  "lp_newton.c"; "lp_qmatrix.c"; "qp_avg.c"];

src_linsys = strjoin(fullfile('..', 'linalg', src_linsys));
src_pot = strjoin(fullfile('..', 'potreduce', src_pot));
src_interface = strjoin(fullfile('..', 'interface', src_interface));
src_mex = fullfile('.', 'potlp.c');

src_files = strcat(src_interface, " ", src_pot, " ", src_linsys, " ", src_mex);
delete(fullfile("./", "*." + mexext));

inc = strjoin("-I" + [fullfile("..", "include");
                      fullfile("..", "potreduce");
                      fullfile("..", "linalg");
                      fullfile("..", "interface")]);
no_newton = "-DLINSYSDUMMY";
underblas = "-DUNDERBLAS";

if mexext == "mexw64"
    libblas = '-lut -lmwblas -lmwlapack';
else
    libblas = '-lm -lut -lmwblas -lmwlapack';
end % End if

mex_cmd = sprintf('mex %s %s %s %s %s -output potlp', underblas, no_newton, src_files, inc, libblas);
eval(mex_cmd);

movefile('potlp.mexmaci64', fullfile('..', '..', '..', '..', 'src'));