function [numProcs] = check_num_files(basename)

proc = 1;
while exist(sprintf('%s%d.dat', basename, proc), 'file')
    proc = proc+1;
end
numProcs=proc-1;