function  [TimeSO,sdinit,sdend] = getsdfTimesSO(task,Serial_order,BehTimesSDF)


 
sdinit = BehTimesSDF(1,task).Taps(Serial_order)+1; %0 is -200

if Serial_order < 4
TimeSO = BehTimesSDF(1,task).Taps(2)+100;
sdend = BehTimesSDF(1,task).Taps(Serial_order+1)+100; 
else
    TimeSO = BehTimesSDF(1,task).Taps(2)+100;
    sdend = BehTimesSDF(1,task).Taps(Serial_order+1)+100; 
end