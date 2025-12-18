% ==========================
% PLOTANJE IZRAČUNANIH TABEL (GRAFOV)
% ==========================

for varIdx = 1:6
    figure; hold on; grid on;
    title(ImenaSpr(varIdx));
    xlabel('čas [dni]');
    ylabel('|pomik|');

    for inloop = 1:11

        % multiplier rule differs for varIdx==1
        if varIdx == 1
            mult = inloop;                     % 1,2,...,11
        else
            mult = 1 + (inloop-1)*0.1;         % 1.0,1.1,...,2.0
        end

        fname = strcat("TAB_", ImenaSpr(varIdx), ...
                       " multiplier: ", num2str(mult), ".xlsx");

        if ~isfile(fname)
            warning("Missing file: %s", fname);
            continue
        end

        T = readtable(fname, "Sheet","Pomiki");

        if inloop == 1
            plot(T.cas_ure, T.pomik_abs, 'r', 'LineWidth', 2, ...
                 'DisplayName', 'Prva črta (multip = 1)');
        else
            plot(T.cas_ure, T.pomik_abs, 'k', 'LineWidth', 1, ...
                 'HandleVisibility','off');
        end
    end

    legend('Location','best');
end
