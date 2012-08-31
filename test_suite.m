[VS, ES, ref] = gen_test_S1;            
Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('Test S1       OK!'); else disp('Test S1       FAILED!'); return; end

[VS, ES, ref] = gen_test_S2;            
Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('Test S2       OK!'); else disp('Test S2       FAILED!'); return; end

% S3 is NOT APPLICABLE for POL version

[VS, ES, ref] = gen_test_S4;            
Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('Test S4       OK!'); else disp('Test S4       FAILED!'); return; end

[VS, ES, ref] = gen_test_S5;            
Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('Test S5       OK!'); else disp('Test S5       FAILED!'); return; end


[VS, ES, ref] = gen_test_full(2500);          %~2.5s       %biograph is SLOW
Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('Test FULL     OK!'); else disp('Test FULL     FAILED!'); return; end


[VS, ES, ref] = gen_test_circle(400);
Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('Test CIRCLE   OK!'); else disp('Test CIRCLE   FAILED!'); return; end


[VS, ES, ref] = gen_test_circle_central(400); %~2.9s
Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('Test CIRCLE_C OK!'); else disp('Test CIRCLE_C FAILED!'); return; end



disp('All tests OK!')
