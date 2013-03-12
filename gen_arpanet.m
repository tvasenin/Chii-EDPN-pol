function [V E ref] = gen_arpanet(n)
%gen_arpanet  Generate ArpaNet graph truncated at n edges (no argument means no truncation)
%   Detailed explanation goes here

Kao = [0,3,5,9,10,12,15,17,20,21,24,26,29,33,35,38,40,43,47,48,50,53,55,58,60,62,64,66,69,71,74,77,80,84,86,87,89,92,94,97,100,103,104,107,109,112,115,117,120,122,123,125,127,129,132,134,137,140,142,144];
Fo  = [2,7,26,3,1,4,5,12,2,3,6,3,8,9,5,8,1,10,6,7,6,11,13,8,12,10,18,11,3,59,14,21,10,13,15,16,14,20,17,15,23,18,16,19,17,27,12,18,21,15,22,20,13,23,21,24,22,17,25,23,58,24,57,1,28,18,29,27,32,30,28,37,29,31,30,49,32,33,28,31,38,36,34,32,35,33,34,37,33,40,30,36,39,33,40,38,54,52,39,37,51,59,43,43,45,41,42,46,45,46,44,43,47,44,45,48,46,55,49,47,31,48,53,52,41,40,51,54,50,56,39,53,57,48,58,57,54,26,55,56,25,56,13,41];

V = ones(length(Kao)-1,1)';

E = KaoFo2std(Kao, Fo);

if (nargin < 1)
    n = length(E)
end

V = V(1:n);
E = E(1:n,1:n);

switch n
    case  0, title = 'Correct ECP00:'; ref_str = '';
    case  1, title = 'Correct ECP01:'; ref_str = '';
    case  2, title = 'Correct ECP02:'; ref_str = '';
    case  3, title = 'Correct ECP03:'; ref_str = '';
    case  4, title = 'Correct ECP04:'; ref_str = '';
    case  5, title = 'Correct ECP05:'; ref_str = '';
    case  6, title = 'Correct ECP06:'; ref_str = '';
    case  7, title = 'Correct ECP07:'; ref_str = '';
    case  8, title = 'Correct ECP08:'; ref_str = '';
    case  9, title = 'Correct ECP09:'; ref_str = '';
    case 10, title = 'Correct ECP10:'; ref_str = '-2 -11  -7  14  14  14  13  10   0   0';
    case 11, title = 'Correct ECP11:'; ref_str = '';
    case 12, title = 'Correct ECP12:'; ref_str = '';
    case 13, title = 'Correct ECP13:'; ref_str = '';
    case 14, title = 'Correct ECP14:'; ref_str = '';
    case 15, title = 'Correct ECP15:'; ref_str = '1   6   2  -2   4 -41 -27  13  44  37  29  23  16   0   0';
    case 16, title = 'Correct ECP16:'; ref_str = '1   6   1 -10   0  10 -36 -23  16  46  38  30  24  17   0   0';
    case 17, title = 'Correct ECP17:'; ref_str = '1   6   1 -11  -8   6  15 -32 -20  18  47  39  31  25  18   0   0';
    case 18, title = 'Correct ECP18:'; ref_str = '-1  -17  -47  135   50  -82  -94  -13   13  -46    8   43   65   51   39   29   20    0    0';
    case 19, title = 'Correct ECP19:'; ref_str = '-2  -23  -33  137   45  -91  -91  -15    0  -42   16   50   70   56   42   31   21    0    0';
    case 20, title = 'Correct ECP20:'; ref_str = '-5  -32  -11  153   35 -121 -105   -9    5  -31   25   55   74   58   44   33   22    0    0';
    case 21, title = 'Correct ECP21:'; ref_str = '2   23    7 -229  113  299   31 -208 -133  -14  -42  -13   45   70   84   62   52   37   24    0    0';
    case 22, title = 'Correct ECP22:'; ref_str = '4   24  -23 -194  134  284    6 -204 -133  -33  -34   -3   54   76   87   67   55   39   25    0    0';
    case 23, title = 'Correct ECP23:'; ref_str = '-2  -21   25  322 -525 -367  496  615  -89 -394 -165  -98  -45   21   61   93  112   81   63   43   27    0    0';
    case 24, title = 'Correct ECP24:'; ref_str = '-4  -24   76  242 -554 -294  546  592 -123 -405 -188  -98  -30   26   66  104  119   86   66   45   28    0    0';
    case 25, title = 'Correct ECP25:'; ref_str = '-2   -7   27   -4  213 -481 -244  523  558 -134 -428 -188  -83  -25   31   77  111  124   89   68   46   29    0    0';
    case 26, title = 'Correct ECP26:'; ref_str = '-4   -1   16   59   82 -461 -158  589  500 -191 -461 -170  -67  -26   40   85  120  131   94   70   48   30    0    0';
    case 27, title = 'Correct ECP27:'; ref_str = '';
    case 28, title = 'Correct ECP28:'; ref_str = '';
    case 29, title = 'Correct ECP29:'; ref_str = '';
    case 30, title = 'Correct ECP30:'; ref_str = '';
%    case 31, title = 'Correct?ECP31:'; ref_str = '-4   -5   24   17  -62  -44  102  342 -478 -721  777  949 -385 -598 -209  -27  -61  119  141  154  154  110   80   55   35    0    0';
    case 31, title = 'Correct ECP31:'; ref_str = '-2   -1   22   -9  -39   -6   98  101 -454 -125  613  453 -265 -488 -174  -64   -6   91  132  154  154  110   80   55   35    0    0';
    case 32, title = 'Correct?ECP32:'; ref_str = '';
    case 33, title = 'Correct?ECP33:'; ref_str = '';
    case 34, title = 'Correct?ECP34:'; ref_str = '';
    case 35, title = 'Correct?ECP35:'; ref_str = '';
%    case 36, title = 'Correct?ECP36:'; ref_str = '4  -15  -45  108  132 -246 -217  165  631 -443 -840  741  915 -479 -645 -175   42   -4  164  170  175  170  121   95   65   41    0    0';
    case 36, title = 'Correct?ECP36:'; ref_str = '2   -9  -25  118  -10 -183 -119  276  287 -496 -226  584  453 -358 -537 -146    4   51  136  161  175  170  121   95   65   41    0    0';
%    case 37, title = 'Correct?ECP37:'; ref_str = '2  -15   -2  131 -135 -191  246  208 -122 -748  571  386  283 -744 -628  817  737 -493 -599 -181   77   17  182  183  162  174  134  102   69   43    0    0';
    case 37, title = 'Correct?ECP37:'; ref_str = '2   -7   -6   85 -211  114  202   79 -494 -181  334  342  121 -563  -60  494  304 -393 -496 -118   34   72  154  174  162  174  134  102   69   43    0    0';
    case 38, title = 'Correct?ECP38:'; ref_str = '';
    case 39, title = 'Correct?ECP39:'; ref_str = '';
%    case 40, title = 'Correct?ECP40:'; ref_str = '-2    23   -20  -162   336  -174   -61   145    83    12  -656   293  -584  1166   155   385  -904  -880   730   713  -464  -524  -159   132   102   166   159   180   194   156   116    77    47     0     0';
%    case 41, title = 'Correct?ECP41:'; ref_str = '-2    23   -20  -162   336  -174   -61   145    83    12  -656   293  -584  1166   155   385  -904  -880   730   713  -464  -524  -159   132   102   166   159   180   194   156   116    77    47     0     0';
%    case 42, title = 'Correct?ECP42:'; ref_str = '-2    23   -20  -162   336  -174   -61   145    83    12  -656   293  -584  1166   155   385  -904  -880   730   713  -464  -524  -159   132   102   166   159   180   194   156   116    77    47     0     0';
    case 40, title = 'Correct?ECP40:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   77   47    0    0';
    case 41, title = 'Correct ECP41:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   77   47    0    0';
    case 42, title = 'Correct ECP42:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   77   47    0    0';
    case 43, title = 'Correct?ECP43:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   78   49    0    0';
    case 44, title = 'Correct ECP44:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   78   49    0    0';
    case 45, title = 'Correct ECP45:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  118   81   51    0    0';
    case 46, title = 'Correct ECP46:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  120   82   53    0    0';
    case 47, title = 'Correct ECP47:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  158  121   84   54    0    0';
    case 48, title = 'Correct ECP48:'; ref_str = '-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  196  159  123   85   55    0    0';
    case 49, title = 'Correct?ECP49:'; ref_str = '';
    case 50, title = 'Correct?ECP50:'; ref_str = '';
    case 51, title = 'Correct?ECP51:'; ref_str = '2   -1  -17  -10   13  438 -784 -144  514  969 -785 -818  172  382  265  -66  411 -555 -287  302  264 -321 -352   21  150  209  176  189  210  219  178  133   90   58    0    0';
%    case 52, title = 'Correct?ECP52:'; ref_str = '6     83   -186  -2026   4319   6885 -13816 -19068  23794  38816 -25693 -48086  13102  37207  -1246 -15393   -378    677  -2924   3851   4221  -1466  -3588  -1783   2218    810   -596   -178    416   -163   -485    -46    182    398    248    251    225    194    195    131     99     66     -1      3';
    case 52, title = 'Correct?ECP52:'; ref_str = '';
    case 53, title = 'Correct?ECP53:'; ref_str = '';
    case 54, title = 'Correct?ECP54:'; ref_str = '';
    case 55, title = 'Correct?ECP55:'; ref_str = '';
%    case 56, title = 'Correct?ECP56:'; ref_str = '2    -5   -21    81   -71  -147   646  -387 -1261   328  2428   351 -2206 -1474  1020   430    -1   341   879  -286  -415   150  -131  -102  -141   251  -812  -323   219   256  -161  -220   209   240   269   255   273   281   269   208   151   102    65     0     0';
    case 56, title = 'Correct?ECP56:'; ref_str = '2    -5   -21    81   -71  -147   646  -387 -1260   328  2422   370 -2272 -1386  1056   311   -18   408   944  -385  -377   152  -152   -91  -143   258  -816  -323   219   256  -161  -220   209   240   269   255   273   281   269   208   151   102    65     0     0';
    case 57, title = 'CorrectXECP57:'; ref_str = '1       6      -9    -104    -469    3112   -9484   45853 -123343   79444  153809 -143354 -145458   36604  191689   40851 -114941  -75156    5651   60757   19477   -9179  -26415   -9337   10514   13195    5874   -2337     611   -2580   -3608    -867   -2883   -1462     349     571     531      61     525     511     537     435     389     349     307     230     163     108      68       0       0';
    case 58, title = 'Correct?ECP58:'; ref_str = '';
    case 59, title = 'Correct?ECP59:'; ref_str = '';
    case 60, title = 'Correct?ECP60:'; ref_str = '';
end

ref = sscanf(ref_str,'%d');

disp([title ref_str]);

end
