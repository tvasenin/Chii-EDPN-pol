function cuts = VertexCuts(E)
%
n = length(E);
used(1:n) = false; % not from zero!
timer = 0;
tin(1:n) = 0; % not from zero!
fup(1:n) = 0; % not from zero!
cuts(1:n) = false;

for k = n:-1:1
    neis_db{k} = find(E(k,:));
end

dfs(1);
%cuts = find(cuts);

    function dfs(v, p)
        if nargin == 1
            p = -1;
        end
        used(v) = true;
        tin(v) = timer;
        fup(v) = timer;
        timer = timer + 1;
        children = 0;
%        neis = find(E(v,:));
        neis = neis_db{v};
        for i = 1:length(neis)
            to = neis(i);
            if (to == p)
                continue;
            end
            if used(to)
                fup(v) = min (fup(v), tin(to));
            else
                dfs (to, v);
                fup(v) = min (fup(v), fup(to));
                if (fup(to) >= tin(v)) && (p ~= -1)
                    cuts(v) = true;
                end
                children = children + 1;
            end
        end
        if (p == -1) && (children > 1)
            cuts(v) = true;
        end
    end

end