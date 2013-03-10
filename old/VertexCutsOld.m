function cuts = VertexCuts(E)
%
n = length(E);
cuts = false(1,n);
used = false(1,n);
timer = 0;
tin = zeros(1,n);
fup = zeros(1,n);

dfs(1);

    function dfs(v, p)
        if nargin == 1
            p = -1;
        end
        used(v) = true;
        tin(v) = timer;
        fup(v) = timer;
        timer = timer + 1;
        children = 0;
        for to = find(E(v,:))
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