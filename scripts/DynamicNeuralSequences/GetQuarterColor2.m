function quarterC = GetQuarterColor2(PeakQuarterNum,cellidx)
 

numQ = PeakQuarterNum;

    if cellidx <= numQ(1)
           quarterC =  1;
    elseif cellidx > numQ(1) & cellidx <= numQ(2)
           quarterC =  2;
    elseif cellidx > numQ(2) & cellidx <= numQ(3)
           quarterC =  3;           
    elseif cellidx > numQ(3) & cellidx <= numQ(4)
           quarterC =  4;      
    elseif cellidx > numQ(4) & cellidx <= numQ(5)
           quarterC =  5;   %%%high variability
    elseif cellidx > numQ(5) & cellidx <= numQ(6)
           quarterC =  6;   %%one AP           
    end