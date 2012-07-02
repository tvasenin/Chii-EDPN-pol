ECPN_ordered_chain_mu := proc(V:array,VW:array)
    local P, n; //,iter;
begin
    n := nops(V):
    case n
        of 0 do
            print(Unquoted,"[WARNING] Attemting co call ECPN_ordered_chain with n=0"):
            P := 0:
            break
        of 1 do
            print(Unquoted,"[WARNING] Attemting co call ECPN_ordered_chain with n=1"):
            P := 0:
            break
        of 2 do
            P := VW[1]*VW[2]:
            break
        of 3 do
            //P := VW[2]*(VW[1]+VW[3]) + VW[1]*VW[3]*V[2]:
            P := VW[2]*VW[3] + VW[1]*(VW[2]+V[2]*VW[3]):
            break

            //P :=
            // VW[1]*(VW[2]+V[2]*(VW[3]+V[3]*VW[4])
            //+VW[2]*(VW[3]+V[3]*VW[4])
            //+VW[3]* VW[4]
            //[* 6]
            //[+ 5]

        otherwise
            //iter := VW[n]:
            P := VW[n-1]*VW[n]:
            for i from n-2 downto 1 do
                VW[n] := VW[i+1] + V[i+1] * VW[n]:
                P := P + VW[i]*VW[n]:
            end_for:
    end_case:
    return(P)
end_proc:
