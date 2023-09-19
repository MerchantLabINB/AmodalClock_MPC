function [movs2,ddwells2,mov_ini,mov_end] = Sort_dwell_movs(movs_aux,dwells_aux,mov_i)
if(mod(length(movs_aux),2)) %impar
    
    movs2 = sum(movs_aux(1:2)) + 3 ;% + 3; %3
    ddwells2 = max(dwells_aux);
    mov_ini = mov_i(1) - movs_aux(1) - 2; %-2
    mov_end =mov_i(2) ;%-1
%     ddwells_test = min(dwells_aux);

else
    
    
    if(length(dwells_aux)==2)
        dwells_aux = max(dwells_aux);
    end
    movs2 = sum(movs_aux(1))+ 3;% + 3 ;%Se agrega  3 para que coincida
    ddwells2 = dwells_aux(1); 
    mov_ini = mov_i(1) - movs2 - 2; %-2
    mov_end = mov_i(1);%-1
end


