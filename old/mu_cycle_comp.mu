muCycleComp := proc(V:matrix,W:matrix)//,VW:matrix)
    local n, VW, P, mult, mult2, i;//, Q;//, W_n;
begin
    n := nops(V):
    //V := array(1..n,op[V]):
    //W := array(1..n,op[W]):
    //TODO: check if we can use V,W,VW as arrays
    VW := array(1..n,[V[i]*W[i] $ i = 1..n]):

    case n
    of 4 do
        P := (VW[1]+VW[3])*(VW[2]+VW[4]) + VW[1]*VW[3]*(V[2]+V[4]-V[2]*V[4]) + VW[2]*VW[4]*(V[1]+V[3]-V[1]*V[3]):
        break
    otherwise
        //mult = cumprod([V(n) V(1:n-1)]);
    
        P := (1-V[n]) * ECPN_ordered_chain_mu([op(V,1..n-1)],[op(VW,1..n-1)]): //initialization
        mult2:= V[n]:
        V[n] := 1:
        //VW[n]:= VW[n]/mult2:
        //VW[n]:= VW[n]/V[n]:
        VW[n]:= W[n]:
        for i from 1 to n-2 do
            P := P + mult2 * (1-V[i]) * ECPN_ordered_chain_mu([op(V,i+1..n)],[op(VW,i+1..n)]): 
            //            %V(i) = 1; not used anymore anyway
            P := P + mult2 * VW[i] * VW[n]: // contracting 2 reliable nodes
            VW[n] := VW[n] + W[i]:
            mult2 := mult2 * V[i]:
        end_for:
        P := P + mult2 * VW[n-1] * VW[n] // case i = n-1;
    end_case:

    return(P)
end_proc: