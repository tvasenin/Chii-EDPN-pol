clear variables   % Uncomment for release and benchmarks!!
%reset(symengine) % no symengine for POL version
%-read data
clear ES VS p WS ExP_res ExP_coeffs tStart tElapsed BGView

clear global cnt
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

%clear  global HIT MISS
%global HIT MISS
%HIT  = 0;
%MISS = 0;

%%
%
disp('---------------------------------------------');
disp('---------------------------------------------');

%[VS, ES] = gen_test_tree_balanced(1,1);  %~ ?.?s 
 
%[VS, ES] = gen_test_big_G1_easy;        %~ 0.3s
%[VS, ES] = gen_test_big_G1;             %~ 0.8s
%[VS, ES] = gen_test_big_G2;             %~ 1.0s

    arpa_time(01:17)= 0.1;                                  %2013-03-13
    arpa_time(18:20)= 0.1;                                  %2013-03-13
    arpa_time(21:22)= 0.3;                                  %2013-03-13
    arpa_time(23:31)= 0.9;                                  %2013-03-13
%arpa_num = 36;                          %~5s               %2012-07-04
    arpa_time(32:36)=[ 5 5 5 5 5 ];                         %2012-07-04
%arpa_num = 37;                          %~23s              %2012-07-04
%arpa_num = 38;                          %
%arpa_num = 39;                          %
    arpa_time(37:39)=[ 23 23 23];                           %2012-07-04
%arpa_num = 42;                          %~95s              %2012-07-05
%arpa_num = 48;                          %~96s              %2012-07-05
    arpa_time(40:48)=[ 95 95 95  0  0  0  0 96 96 ];        %2012-07-05
%arpa_num = 49;                          %
%%arpa_num = 50;                          % should be equal to 49
%arpa_num = 51;                          %~722s    12m02s   %2012-07-04
%arpa_num = 51;                          %~694s    11m34s   %2012-07-04
%arpa_num = 51;                          %~658s    10m58s   %2012-07-05
%arpa_num = 51;                          %~630s    10m30s   %2012-07-06
    arpa_time(49:51)=[   0   0  630];                       %2012-07-06
%arpa_num = 52;                          %
%arpa_num = 53;                          %
%arpa_num = 54;                          %
%arpa_num = 55;                          %
%arpa_num = 56;                          %~4617s   ~1h17m   %2012-03-15
%arpa_num = 56;                          %~4144s   ~1h10m   %2012-07-04
%arpa_num = 56;                          %~3676s   ~1h02m   %2012-07-07
%arpa_num = 56;                          %~4076s   ~1h08m   %2012-08-27 [S]
%arpa_num = 56;                          %~400Xs   ~1h07m   %2012-09-01 [S]
%arpa_num = 56;                          %~386Xs   ~1h05m   %2013-03-11 [S]
%arpa_num = 56;                          %~379Xs   ~1h04m   %2013-03-12 [S]
    arpa_time(52:56)=[    0    0    0    0 4076 ];          %2012-08-27 [S]
%arpa_num = 57;                          %~87005s ~24h10m   %2012-07-04
%arpa_num = 57;                          %~923XXs ~25h04m   %2012-09-04 [S]
    arpa_time(57)= 87005;
%arpa_num = 58;                          %
    arpa_time(58)=      0;
%arpa_num = 59;                          %
    arpa_time(59)=      0;

[VS, ES, ref] = gen_arpanet(arpa_num);

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

[~, ~, ~] = gen_arpanet(arpa_num);
if isequal(ref, ExP_coeffs') || isequal(ref, ExP_coeffs), disp('OK!'); else disp('FAILED!'); end
