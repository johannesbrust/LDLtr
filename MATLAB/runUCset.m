function AllStats = runUCset
%        AllStats = runUCset
%
% Runs the BFGS method on UC problems in the CUTEst test set.

jSkipped   = 1;    % Allstats(1) =  average number of skipped updates
jTimes     = 2;	 % Allstats(2) =  total cpu time
jProbs     = 3;	 % Allstats(3) =  total number of problems attempted
jSolved    = 4;	 % Allstats(4) =  total number of problems solved

AllStats   = zeros(1,4);

%
% NOTE: This path has to modified to point to a CUTEst installation
%

%addpath(genpath("/home/jburst/Documents/CUTEST"));
%addpath(genpath("/home/jburst/Dropbox/UCSD/RESEARCH/PGILL/MASTSIF"));
addpath(genpath("/home/jburst/Documents/CUTEST/cutest"));

parms      = bfgs_parms();
write_file = parms.write_file;
solver     = parms.solver;
if write_file
  fid     = fopen(strcat('Summary_', solver, '.csv'), 'w');
end

 [x,status,Stats]=runProb('AKIVA   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ALLINITU');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ARGLINA ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ARGLINB ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ARGLINC ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ARGTRIGLS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ARWHEAD ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BA-L1LS ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BA-L1SPLS');  [AllStats] = collectStats(status,Stats,AllStats); %
 % [x,status,Stats]=runProb('BA-L16LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 % [x,status,Stats]=runProb('BA-L21LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 % [x,status,Stats]=runProb('BA-L49LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 % [x,status,Stats]=runProb('BA-L52LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 % [x,status,Stats]=runProb('BA-L73LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BARD    ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BDQRTIC ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BEALE   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BENNETT5LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 %[x,status,Stats]=runProb('BIGGS6  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BOX     ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BOX3    ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BOXBODLS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BOXPOWER');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BRKMCC  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BROWNAL ');   [AllStats] = collectStats(status,Stats,AllStats); %
%[x,status,Stats]=runProb('BROWNBS ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BROWNDEN');   [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('BROYDN3DLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BROYDN7D');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BROYDNBDLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('BRYBND  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CERI651ALS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CERI651BLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CERI651CLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CERI651DLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CERI651ELS'); [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('CHAINWOO');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CHNROSNB');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CHNRSNBM');   [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('CHWIRUT1LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CHWIRUT2LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CLIFF   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CLUSTERLS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('COATING ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('COOLHANSLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('COSINE  ');   [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('CRAGGLVY');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CUBE    ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CURLY10 ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CURLY20 ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('CURLY30 ');   [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('CYCLOOCFLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DANIWOODLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DANWOODLS');  [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('DECONVU ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DENSCHNA');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DENSCHNB');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DENSCHNC');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DENSCHND');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DENSCHNE');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DENSCHNF');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DEVGLA1 ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DEVGLA2 ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIAMON2DLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIAMON3DLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANA');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANB');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANC');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAAND');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANE');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANF');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANG');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANH');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANI');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANJ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANK');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANL');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANM');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANN');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANO');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXMAANP');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DIXON3DQ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DJTL    ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DMN15102LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DMN15103LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DMN15332LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DMN15333LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DMN37142LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DMN37143LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DQDRTIC ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('DQRTIC  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ECKERLE4LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EDENSCH ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EG2     ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EGGCRATE');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EIGENALS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EIGENBLS');   [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('EIGENCLS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ELATVIDU');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ENGVAL1 ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ENGVAL2 ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ENSOLS  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ERRINROS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ERRINRSM');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EXP2    ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EXPFIT   ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('EXTROSNB ');  [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('FBRAIN3LS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FLETBV3M ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FLETCBV2 ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FLETCBV3 ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FLETCHBV');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FLETCHCR');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FMINSRF2');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FMINSURF');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('FREUROTH');   [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('GAUSS1LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GAUSS2LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GAUSS3LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GAUSSIAN');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GBRAINLS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GENHUMPS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GENROSE ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GROWTHLS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('GULF    ');   [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('HAHN1LS ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HAIRY   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HATFLDD ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HATFLDE ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HATFLDFL');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HATFLDFLS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HATFLDGLS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HEART6LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HEART8LS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HELIX   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HIELOW  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HILBERTA');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HILBERTB');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HIMMELBB ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HIMMELBCLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HIMMELBF ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HIMMELBG ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HIMMELBH ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HUMPS    ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('HYDC20LS ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('INDEF    ');  [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('INDEFM   ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('INTEQNELS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('JENSMP   ');  [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('JIMACK   ');  [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('KIRBY2LS ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('KOWOSB   ');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LANCZOS1LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LANCZOS2LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LANCZOS3LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LIARWHD   '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LOGHAIRY  '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LSC1LS    '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LSC2LS    '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN11LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN12LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN13LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN14LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN15LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN16LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN17LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN21LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('LUKSAN22LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MANCINO   '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MARATOSB  '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MEXHAT    '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MEYER3    '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MGH09LS   '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MGH10LS   '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MGH17LS   '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MISRA1ALS '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MISRA1BLS '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MISRA1CLS '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MISRA1DLS '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MNISTS0LS '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MNISTS5LS '); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MODBEALE');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MOREBV  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MSQRTALS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('MSQRTBLS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NCB20   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NCB20B  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NELSONLS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NONCVXU2');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NONCVXUN');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NONDIA  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NONDQUAR');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('NONMSQRT');   [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('OSBORNEA');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('OSBORNEB');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('OSCIGRAD');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('OSCIPATH');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER1C');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER1D');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER2C');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER3C');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER4C');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER5C');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER6C');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER7C');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PALMER8C');   [AllStats] = collectStats(status,Stats,AllStats); %
  [x,status,Stats]=runProb('PARKCH  ');   [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('PENALTY1');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PENALTY2');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('PENALTY3');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('POWELLBSLS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('POWELLSG');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('POWER   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('QUARTC  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('RAT42LS ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('RAT43LS ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ROSENBR ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ROSENBRTU');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ROSZMAN1LS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('S308    ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SBRYBND ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SCHMVETT');   [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('SCOSINE ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SCURLY10');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SCURLY20');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SCURLY30');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SENSORS ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SINEVAL ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SINQUAD ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SISSER  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SNAIL   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SPARSINE');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SPARSQUR');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SPMSRTLS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SROSENBR');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SSBRYBND');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SSCOSINE');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('SSI'     );   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('STRATEC ');   [AllStats] = collectStats(status,Stats,AllStats); %
   [x,status,Stats]=runProb('TESTQUAD');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('THURBERLS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('TOINTGOR');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('TOINTGSS');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('TOINTPSP');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('TOINTQOR');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('TQUARTIC');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('TRIDIA  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('VARDIM  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('VAREIGVL');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('VESUVIALS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('VESUVIOLS');  [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('VESUVIOULS'); [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('VIBRBEAM');   [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('WATSON  ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('WOODS   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('YATP1LS ');   [AllStats] = collectStats(status,Stats,AllStats); %
[x,status,Stats]=runProb('YATP2LS ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('YFITU   ');   [AllStats] = collectStats(status,Stats,AllStats); %
 [x,status,Stats]=runProb('ZANGWIL2');   [AllStats] = collectStats(status,Stats,AllStats); %

Solved          = AllStats(jSolved);
Time            = AllStats(jTimes);
Probs           = AllStats(jProbs);
iStat           = fix(AllStats);
fprintf('\n STATS    #Problems = %4i, #Solved = %4i,  Skipped = %3i%%, Total time = %8.2f\n\n', ...
                     Probs, Solved, iStat(jSkipped), Time );

return
