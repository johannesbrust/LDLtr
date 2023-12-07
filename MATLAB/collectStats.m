function AllStats = collectStats(status,Stats,AllStats);
%        AllStats = collectStats(status,Stats,AllStats);
% Collects statistics for the quasi-Newton runs

jSkipped   = 1;
jTimes     = 2;
jProbs     = 3;
jSolved    = 4;

AllStats(1:2)    = AllStats(1:2) + Stats(1:2);
AllStats(jProbs) = AllStats(jProbs) + 1;
if Stats(jSkipped) > 0
  AllStats(jSkipped) = AllStats(jSkipped) + 1;
end

if status == 0 || status == 1 || status == 2
  AllStats(jSolved) = AllStats(jSolved) + 1;
end
