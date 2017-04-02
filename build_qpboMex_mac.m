function build_qpboMex_mac
% build_qpboMex builds package qpboMex for Mac based on solution from
% http://stackoverflow.com/questions/15298603/linking-errors-from-matlab-mex-library
%
% Anton Osokin (firstname.lastname@gmail.com),  24.09.2014

codePath = 'QPBO-v1.32.src';

srcFiles = { 'qpboMex_mac.cpp' };
allFiles = '';
for iFile = 1 : length(srcFiles)
    allFiles = [allFiles, ' ', srcFiles{iFile}];
end

cmdLine = ['mex ', allFiles, ' -output qpboMex -largeArrayDims ', '-I', codePath];
eval(cmdLine);




