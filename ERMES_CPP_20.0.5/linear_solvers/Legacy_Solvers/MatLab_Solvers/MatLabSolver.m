%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% - Example of use of an external solver with ERMES 8.0
% 
% - Before using this example:
%   1-) Set the path of this folder inside MATLAB
%   2-) Write the MATLAB.exe path inside GiD (Solving parameters window)
%   3-) Call from ERMES: 
%       Solver path      = C:\\MATLAB7\\bin\\win32\\MATLAB.exe
%		Input parameters = -nosplash -nodesktop -minimize -r MatLabSolverERMESv7p0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear  all;
format long e; 

%%%%%%%%%%%%%%%%%%%%%%% Reading source vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('Vector_b.bin', 'r');
s   = fread(fid,inf,'double');
fclose(fid);

Vector_b=[s(1:2:length(s)) s(2:2:length(s))];

b = Vector_b(:,1) + i*Vector_b(:,2);

clear fid;
clear Vector_b;
clear s;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Reading matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid    = fopen('Matrix_A_int.bin', 'r');
indexs = fread(fid,inf,'int');
fclose(fid);

fid    = fopen('Matrix_A_cmplx.bin', 'r');
compxs = fread(fid, inf, 'double');
fclose(fid);

i_index  = indexs(1:2:length(indexs));
j_index  = indexs(2:2:length(indexs));

c_values = compxs(1:2:length(compxs)) + i*compxs(2:2:length(compxs));

A = sparse(i_index, j_index, c_values);

clear fid;
clear indexs;
clear compxs;
clear i_index;
clear j_index;
clear c_values;

A = A + triu(A,1).';

%%%%%%%%%%%%%%%%%%%%%%%%% Solving system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

%Reordering
p=symamd(A);

rA = A(p,p);
rb = b(p);

clear A;

solver = 'BiConjugate Gradient (BiCG)';

flag   = 0;
relres = 0;
iter   = 0;

tol    = 1e-3;
maxit  = 100;

%Preconditioner
R = cholinc(rA,1e-12);

[x, flag, relres, iter] = bicg(rA, rb, tol, maxit, R.', R);

x(p) = x;

etime = toc;

%%%%%%%%%%%%%%%%%%%%% Writing solution in file %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Resultfile = fopen('Vector_x.dat','w');

for i=1:length(x)
   fprintf(Resultfile,'%12.8f  %12.8f\n',real(x(i)),imag(x(i)));
end

fclose(Resultfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infoFile = fopen('info_solver.dat','w');

fprintf(infoFile,'ENDL\n');

fprintf(infoFile,'MATLAB solver: %s ENDL\n', solver);

fprintf(infoFile,'ENDL\n');

fprintf(infoFile, '- Flag       : %i ENDL\n', flag     );
fprintf(infoFile, '- Size       : %i ENDL\n', length(b));
fprintf(infoFile, '- Iterations : %i ENDL\n', iter     );
fprintf(infoFile, '- Time       : %g ENDL\n', etime    );
fprintf(infoFile, '- Error      : %g ENDL\n', relres   );

fprintf(infoFile,'ENDL\n');

fclose(infoFile);

%%%%%%%%%%%%%%%%%%%%% Info to MATLAB desktop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IsSolvedFile = fopen('IsSolved.dat','w');

if (flag ==  0)
    fprintf(IsSolvedFile,'1\n');
else
    fprintf(IsSolvedFile,'0\n');
end;

fclose(IsSolvedFile);

exit;


