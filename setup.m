% startup script 
%
% Copyright (c) by Tales Imbiriba 2015.

disp(['executing startup script...']);

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located

addpath(mydir(1:end-1))
addpath(genpath([mydir,'Hyperspectral']))
% addpath([mydir,'DataGeneration'])
% addpath([mydir,'Detection'])


stringCommand = ['pctRunOnAll javaaddpath(''',mydir,'Hyperspectral/ParforProgMonv2/java'')'];
fid=fopen('runOnAllJavaMonitorCP.m','wt');
fprintf(fid, '%s',stringCommand);
fclose(fid);


cd gpml-matlab-v3.3-2013-10-19
startup
cd ../sedumi-master
install_sedumi

cd ..
