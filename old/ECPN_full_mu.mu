ECPN_full_mu := proc(VW:matrix)
    local i, j, P;
begin
    n := nops(VW):
    case n
    of 0 do
        print(Unquoted,"[ERROR] Attempting to call ECPN_full with n=0");
        break
    of 1 do
        print(Unquoted,"[ERROR] Attempting to call ECPN_full with n=1");
        break
    of 2 do
        print(Unquoted,"[WARNING] Attempt to call ECPN_full with n=2");
        P := VW[1]*VW[2];
        break
    
    of 3 do
    print(Unquoted,"[WARNING] Attempt to call ECPN_full with n=3"):
    //no break command needed!
    //P := linalg::scalarProduct(matrix([VW[2],VW[3],VW[1]]), VW):
    otherwise
        P := 0:
        for i from 1 to n-1 do
            P := P + VW[i] * (`+`(VW[j] $ j = i+1..n)):
        end_for:
    
        //  P = sum(sum(triu(VW'*VW,1))); %profiling shows it's the best we can do
//    P = matrix(n,1,[1 $ n]) * triu(VW'*VW,1) * matrix(1,n,[1 $ n]); %well, it's FAR better :)
    end_case:
    return(P)
end_proc:

