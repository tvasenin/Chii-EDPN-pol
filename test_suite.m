function test_suite

function run_test(VS, ES, ref, text_ok, text_failed)
    Vrel = false(1,length(VS)); Wpol = ones(length(VS),1); ES = sparse(logical(ES));
    ExP_coeffs = ECPN_pol(Vrel, ES, Wpol); ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)
    if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp(text_ok); else disp(text_failed); assert(false,'FAILED TEST!'); end
end

[VS, ES, ref] = gen_test_S1;
run_test(VS, ES, ref, 'Test S1       OK!', 'Test S1       FAILED!');

[VS, ES, ref] = gen_test_S2;
run_test(VS, ES, ref, 'Test S2       OK!', 'Test S2       FAILED!');

% S3 is NOT APPLICABLE for POL version

[VS, ES, ref] = gen_test_S4;
run_test(VS, ES, ref, 'Test S4       OK!', 'Test S4       FAILED!');

[VS, ES, ref] = gen_test_S5;
run_test(VS, ES, ref, 'Test S5       OK!', 'Test S5       FAILED!');

[VS, ES, ref] = gen_test_full(2500);          %~2.5s      %biograph is SLOW
run_test(VS, ES, ref, 'Test FULL     OK!', 'Test FULL     FAILED!');

[VS, ES, ref] = gen_test_circle(400);
run_test(VS, ES, ref, 'Test CIRCLE   OK!', 'Test CIRCLE   FAILED!');

[VS, ES, ref] = gen_test_circle_central(400); %~2.9s
run_test(VS, ES, ref, 'Test CIRCLE_C OK!', 'Test CIRCLE_C FAILED!');

for i = [10 15:26 31]
    [VS, ES, ref] = gen_arpanet_small(i);
    run_test(VS, ES, ref, ['Test ARPA_' int2str(i) '  OK!'], ['Test ARPA_' int2str(i) '  FAILED!']);
end

disp('All tests OK!')

end