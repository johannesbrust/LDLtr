function Allstats = runset(file)
%        Allstats = runset(file)
%
% Runs BFGS method on UC problems listed in the provided file.
%
%  Expected format:
% number name           #  classification   dimension   comments
% 1      AKIVA          #    O-U      UC         2
% 2      ALLINITU       #    O-U      UC         3
% 3      ARGLINA        #    S-U      UC       200
% 4      ARGLINB        #    S-U      UC       200
% ...
%
% Lines beginning with '#' are skipped, only the names are extracted

jSkipped   = 1;    % Allstats(1) =  average number of skipped updates
jTimes     = 2;    % Allstats(2) =  total cpu time
jProbs     = 3;    % Allstats(3) =  total number of problems attempted
jSolved    = 4;    % Allstats(4) =  total number of problems solved

AllStats   = zeros(1,4);

% Open provided probs file and read in the names
id = fopen(file, 'r');
% Setting CommentStyle to '#' prevents reading ignored lines
data   = textscan(id, '%d %10s', 'CommentStyle', '#', 'Delimiter', '');
names  = data{2};
fclose(id);

N = length(names);

for i = 1:N
   [x, status, Stats] = runProb(names{i});
   [AllStats]         = collectStats(status, Stats, AllStats);
end

write_file   = bfgs_parms().write_file;
solver       = bfgs_parms().solver;
if write_file
  pid = fopen('bfgs_parms.m');
  fid = fopen(strcat('./Statistics/', solver, '.stats'), 'a+');
  str = '########################### Current Parameters ###########################';
  fprintf(fid, '\n%s\n', str);
  while ~strcmp(str, 'end')
     str = fgetl(pid);
     fprintf(fid,'%s\n', str);
  end
  str = '############################# End Parameters #############################';
  fprintf(fid, '\n%s\n', str);
  fclose(pid);
  fclose(fid);
end

Solved = AllStats(jSolved);
Time   = AllStats(jTimes);
Probs  = AllStats(jProbs);
iStat  = fix(AllStats);
fprintf('\n STATS    #Problems = %4i, #Solved = %4i,  Skipped = %3i%%, Total time = %8.2f\n\n', ...
                     Probs, Solved, iStat(jSkipped), Time );