function quarterC = GetQuarterColor(PeakQuarterNum,cellidx)
 

psum = 0;
for q = 1:5
    psum = PeakQuarterNum(q)+psum;
    numQ(q) = psum;
end

    if cellidx <= numQ(1)
           quarterC =  1;
    elseif cellidx > numQ(1) & cellidx <= numQ(2)
           quarterC =  2;
    elseif cellidx > numQ(2) & cellidx <= numQ(3)
           quarterC =  3;           
    elseif cellidx > numQ(3) & cellidx <= numQ(4)
           quarterC =  4;      
    elseif cellidx > numQ(4)
           quarterC =  5;                 
    
    end