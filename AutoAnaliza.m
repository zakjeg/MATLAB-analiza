Ime = "Tabela";
ImenaSpr = [ "dDs_w(m,n,i,p)" "alf1_C" "alf2_C" "alf1_T" "alf2_T" "dDms_w(m,n,i,p)"];
multiplier = 1;
global multip;
multip = [1,1,1,1,1,1];



disp("----------------------- NEW RUN  ------------------------")
disp("------ povecava po 1 --------");;

for inloop=1:11
    multiplier = inloop;
    multip(1)=multiplier
    %disp([ImenaSpr(loop), ' iteration: ', num2str(inloop)]);

    clearvars -except multip ImenaSpr loop inloop multiplier Ime ImeDatoteke
    global multip
    
    run("les.m"); 
    

    ImeDatoteke  = strcat("TAB_",ImenaSpr(1), " multiplier: " , num2str(multiplier), ".xlsx");

    
    disp(ImeDatoteke)

    u = abs(Ru_konstGKSA(8*3+5*3-1, 2:end).');
    
    t = casA(:) / 3600/24 ; 
    
    T = table(t, u, 'VariableNames', {'cas_ure','pomik_abs'});
    
    writetable(T, ImeDatoteke, 'Sheet', 'Pomiki');

end

disp("------ povecava po 0.1 --------");

for loop = 2:6
    multip = [1,1,1,1,1,1];
    for inloop=1:11
        multiplier = 1+((inloop-1)*0.1);
        multip(loop)=multiplier;  %tu se poveca pravi multiplier
        %disp([ImenaSpr(loop), ' iteration: ', num2str(inloop)]);

        
        
        
        

        %KLIC FUNKCIJE

        clearvars -except multip ImenaSpr loop inloop multiplier Ime ImeDatoteke
        global multip
        run("les.m");

        %-----------------------------------------------------------------------------------

        ImeDatoteke  = strcat("TAB_",ImenaSpr(loop), " multiplier: " , num2str(multiplier), ".xlsx");

        disp(ImeDatoteke)

        %USTVARJANJE TABELE

        u = abs(Ru_konstGKSA(8*3+5*3-1, 2:end).');
        
        t = casA(:) / 3600/24 ; 
        
        T = table(t, u, 'VariableNames', {'cas_ure','pomik_abs'});
        
        writetable(T, ImeDatoteke, 'Sheet', 'Pomiki');


        %PLOT TABELE

        %Ru_konstGKSA(8*3+5*3-1,2:end);

        %W5=Ru_konstGKSA(8*3+5*3-1,2:end);

        %plot(casA/60/24,W5)

        %plot(casA/60/24,abs(W5))

        %plot(casA/3600/24,abs(W5))

 


    end
end




