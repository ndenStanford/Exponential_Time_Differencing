function [BHP_1_t,BHP_1_p,BPR_1_t,BPR_1_p,FPR_1_t,FPR_1_p,FOPR_1_t,FOPR_1_p,FGPR_1_t,FGPR_1_p]=importEclipse()
filename = 'BHP_case_1';
BHP_1 = xlsread(filename);
BHP_1_t=BHP_1(1:153,1:1);
BHP_1_p=BHP_1(1:153,2:2);

filename = 'BPR_case_1';
BPR_1 = xlsread(filename);
BPR_1_t=BPR_1(1:153,1:1);
BPR_1_p=BPR_1(1:153,2:2);

filename = 'FPR_case_1';
FPR_1 = xlsread(filename);
FPR_1_t=FPR_1(1:153,1:1);
FPR_1_p=FPR_1(1:153,2:2);

filename = 'FOPR_case_1';
FOPR_1 = xlsread(filename);
FOPR_1_t=FOPR_1(1:153,1:1);
FOPR_1_p=FOPR_1(1:153,2:2);

filename = 'FGPR_case_1';
FGPR_1 = xlsread(filename);
FGPR_1_t=FGPR_1(1:153,1:1);
FGPR_1_p=FGPR_1(1:153,2:2);
end