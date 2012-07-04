clear all   % Uncomment for release and benchmarks!!
%reset(symengine) % no symengine for POL version
%-read data
clear ES VS p WS ExP_res ExP_coeffs tStart tElapsed BGView
%syms     VS p WS  positive % real
%syms        p      positive % real

%MuPAD versions of critical code
%read(symengine, 'ECPN_full_mu.mu');
%read(symengine, 'ECPN_ordered_chain_mu.mu');
%read(symengine, 'mu_cycle_comp.mu');

%evalin(symengine,'getprop(p)')

clear cnt
global cnt

cnt.NOEDGES   = 0;
cnt.MULTICOMP = 0;
cnt.NUMEL2    = 0;
cnt.NUMEL3    = 0;
cnt.NUMEL4    = 0;
cnt.MAXDEG    = 0;
cnt.FULL      = 0;
cnt.CHAIN     = 0;
cnt.CYCLE     = 0;
cnt.TREE      = 0;
cnt.HNODES    = 0;
cnt.CHAINRED  = 0;
cnt.RELIABLE  = 0;
cnt.BRANCHING = 0;
cnt.TOTAL     = 0;

clear  HIT MISS
global HIT MISS
HIT  = 0;
MISS = 0;

%%
%
disp('---------------------------------------------');
disp('---------------------------------------------');

%[VS, ES] = gen_test_circle_central(60);%
%[VS, ES] = gen_test_circle_central(5);%

%[VS, ES] = gen_test_chain_central(30); %

%[VS, ES] = gen_test_S1;                 %
%[VS, ES] = gen_test_S2;                 %
% S3 is NOT APPLICABLE for POL version
%[VS, ES] = gen_test_S4;                 %
%[VS, ES] = gen_test_S5;                 %


%[VS, ES] = gen_test_chain(200);        %~ 1.5s       %test is deprecated
%%[VS, ES] = gen_test_circle(40);        %~ 1.2s       %test is deprecated
%[VS, ES] = gen_test_circle_central(100);%~ 12s 
%[VS, ES] = gen_test_chain_central(400); %~ 3.5s 
%[VS, ES] = gen_test_tree_balanced(1,1);  %~ ?.?s 
 
%[VS, ES] = gen_test_full(3000);          %~0.6s       %biograph is SLOW
%[VS, ES] = gen_test_big_G1_easy;        %~ 0.3s
%[VS, ES] = gen_test_big_G1;             %~ 0.8s
%[VS, ES] = gen_test_big_G2;             %~ 1.0s

    arpa_time(01:17)=[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; %untested
%arpa_num = 20;                          %~0.3s             %2012-03-07
    arpa_time(18:20)=[ 1 1 1 ];                             %2012-03-07
%arpa_num = 21;                          %~0.9s             %2012-03-07
%arpa_num = 22;                          %~0.9s             %2012-03-07
    arpa_time(21:22)=[ 1 1 ];                               %2012-03-07
%arpa_num = 23;                          %~1s               %2012-07-04
%arpa_num = 31;                          %~1s               %2012-07-04
    arpa_time(23:31)=[ 1 1 1 1 1 1 1 1 1 ];                 %2012-07-04
%arpa_num = 36;                          %~5s               %2012-07-04
    arpa_time(32:36)=[ 5 5 5 5 5 ];                         %2012-07-04
%arpa_num = 37;                          %~23s              %2012-07-04
%arpa_num = 38;                          %
%arpa_num = 39;                          %
    arpa_time(37:39)=[ 23 23 23];                           %2012-07-04
%arpa_num = 42;                          %~98s   1m38s      %2012-07-04
%arpa_num = 44;                          %
%arpa_num = 45;                          %
%arpa_num = 46;                          %
%arpa_num = 47;                          %
%arpa_num = 48;                          %
    arpa_time(40:48)=[ 98 98 98  0  0  0  0  0  0 ];        %2012-07-04
%arpa_num = 49;                          %
%%arpa_num = 50;                          % should be equal to 49
%arpa_num = 51;                          %~902s 15m02s      %2012-03-12
%arpa_num = 51;                          %~758s 12m38s      %2012-03-13
%arpa_num = 51;                          %~722s 12m02s      %2012-07-04
    arpa_time(49:51)=[   0   0  902];
%arpa_num = 52;                          %
%arpa_num = 53;                          %
%arpa_num = 54;                          %
%arpa_num = 55;                          %
%arpa_num = 56;                          %~4617s ~1h17m     %2012-03-15
%arpa_num = 56;                          %~4144s ~1h10m     %2012-07-04
    arpa_time(52:56)=[    0    0    0    0 4617 ];          %2012-03-15
%arpa_num = 57;                          %
    arpa_time(57)=     0;
%arpa_num = 58;                          %
    arpa_time(58)=     0;
%arpa_num = 59;                          %
    arpa_time(59)=     0;

[VS, ES] = gen_arpanet_small(arpa_num);

%[VS, ES] = gen_arpanet;                 %%%EPIC


%plot(arpa_time,'-o');

Vrel = false(1,length(VS));
%WS = ones(1,length(VS));
Wpol = ones(length(VS),1);
ES = sparse(logical(ES));
%ES = single(sparse(ES));
%%ES = int8(ES);

disp(['conncomp = ', int2str(graphconncomp(ES,'Directed',false))])

bio_show = false;
%bio_show = true;

if bio_show
    BGView = biograph(triu(ES));
    set(BGView,'ShowArrows','off');%,'LayoutType','equilibrium');
    view(BGView)
end
%%
%
tStart = tic;
ExP_coeffs = ECPN_pol(Vrel, ES, Wpol);

tElapsed=toc(tStart);
disp('[INFO] ECPN found!')

% ExP_res = poly2sym(ExP_coeffs,p);
% ExP_res = collect(ExP_res); % do we really need this?
%ExP_coeffs = sym2poly(ExP_res);
ExP_coeffs = ExP_coeffs(find(ExP_coeffs,1):end); % fails if all ExP_coeffs are zero :)


%fprintf('\nECP result:   %s\n',char(ExP_res));
fprintf('\nno-edges:  %i\nmulticomp: %i\nnumel2:    %i\nnumel3:    %i\nnumel4:    %i\nmaxdeg:    %i\nfull5+:    %i\nchain5+:   %i\ncycle5+:   %i\ntree5+:    %i\nnhodes:    %i\nchainred:  %i\nreliable:  %i\nbranch:    %i\n----------------\nTOTAL:     %i\n', ...
        [  cnt.NOEDGES    cnt.MULTICOMP  cnt.NUMEL2     cnt.NUMEL3     cnt.NUMEL4     cnt.MAXDEG     cnt.FULL       cnt.CHAIN      cnt.CYCLE      cnt.TREE       cnt.HNODES     cnt.CHAINRED   cnt.RELIABLE   cnt.BRANCHING                    cnt.TOTAL          ]);
%fprintf('\nHIT:  %i\nMISS: %i\nHIT/MISS ratio: %f', [HIT MISS HIT/MISS]);
fprintf('\nElapsed time is %8.6f seconds\n\n',tElapsed);
fprintf('\nECP coeffs:   %s\n',int2str(ExP_coeffs));

[~, ~] = gen_arpanet_small(arpa_num);