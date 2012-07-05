function [V E] = gen_arpanet_small(n)
%gen_arpanet_small  Generates SMALL ArpaNet graph
%   Detailed explanation goes here

Kao = [0,3,5,9,10,12,15,17,20,21,24,26,29,33,35,38,40,43,47,48,50,53,55,58,60,62,64,66,69,71,74,77,80,84,86,87,89,92,94,97,100,103,104,107,109,112,115,117,120,122,123,125,127,129,132,134,137,140,142,144];
Fo  = [2,7,26,3,1,4,5,12,2,3,6,3,8,9,5,8,1,10,6,7,6,11,13,8,12,10,18,11,3,59,14,21,10,13,15,16,14,20,17,15,23,18,16,19,17,27,12,18,21,15,22,20,13,23,21,24,22,17,25,23,58,24,57,1,28,18,29,27,32,30,28,37,29,31,30,49,32,33,28,31,38,36,34,32,35,33,34,37,33,40,30,36,39,33,40,38,54,52,39,37,51,59,43,43,45,41,42,46,45,46,44,43,47,44,45,48,46,55,49,47,31,48,53,52,41,40,51,54,50,56,39,53,57,48,58,57,54,26,55,56,25,56,13,41];

V = ones(length(Kao)-1,1)';

E = KaoFo2std(Kao, Fo);

V = V(1:n);
E = E(1:n,1:n);

switch n
    case  0, disp('Correct ECP00:');
    case  1, disp('Correct ECP01:');
    case  2, disp('Correct ECP02:');
    case  3, disp('Correct ECP03:');
    case  4, disp('Correct ECP04:');
    case  5, disp('Correct ECP05:');
    case  6, disp('Correct ECP06:');
    case  7, disp('Correct ECP07:');
    case  8, disp('Correct ECP08:');
    case  9, disp('Correct ECP09:');
    case 10, disp('Correct ECP10:-2 -11  -7  14  14  14  13  10   0   0');
    case 11, disp('Correct ECP11:');
    case 12, disp('Correct ECP12:');
    case 13, disp('Correct ECP13:');
    case 14, disp('Correct ECP14:');
    case 15, disp('Correct ECP15:1   6   2  -2   4 -41 -27  13  44  37  29  23  16   0   0');
    case 16, disp('Correct ECP16:1   6   1 -10   0  10 -36 -23  16  46  38  30  24  17   0   0');
    case 17, disp('Correct ECP17:1   6   1 -11  -8   6  15 -32 -20  18  47  39  31  25  18   0   0');
    case 18, disp('Correct ECP18:-1  -17  -47  135   50  -82  -94  -13   13  -46    8   43   65   51   39   29   20    0    0');
    case 19, disp('Correct ECP19:-2  -23  -33  137   45  -91  -91  -15    0  -42   16   50   70   56   42   31   21    0    0');
    case 20, disp('Correct ECP20:-5  -32  -11  153   35 -121 -105   -9    5  -31   25   55   74   58   44   33   22    0    0');
    case 21, disp('Correct ECP21:2   23    7 -229  113  299   31 -208 -133  -14  -42  -13   45   70   84   62   52   37   24    0    0');
    case 22, disp('Correct ECP22:4   24  -23 -194  134  284    6 -204 -133  -33  -34   -3   54   76   87   67   55   39   25    0    0');
    case 23, disp('Correct ECP23:-2  -21   25  322 -525 -367  496  615  -89 -394 -165  -98  -45   21   61   93  112   81   63   43   27    0    0');
    case 24, disp('Correct ECP24:-4  -24   76  242 -554 -294  546  592 -123 -405 -188  -98  -30   26   66  104  119   86   66   45   28    0    0');
    case 25, disp('Correct ECP25:-2   -7   27   -4  213 -481 -244  523  558 -134 -428 -188  -83  -25   31   77  111  124   89   68   46   29    0    0');
    case 26, disp('Correct ECP26:-4   -1   16   59   82 -461 -158  589  500 -191 -461 -170  -67  -26   40   85  120  131   94   70   48   30    0    0');
    case 27, disp('Correct ECP27:');
    case 28, disp('Correct ECP28:');
    case 29, disp('Correct ECP29:');
    case 30, disp('Correct ECP30:');
%    case 31, disp('Correct?ECP31:-4   -5   24   17  -62  -44  102  342 -478 -721  777  949 -385 -598 -209  -27  -61  119  141  154  154  110   80   55   35    0    0');
    case 31, disp('Correct ECP31:-2   -1   22   -9  -39   -6   98  101 -454 -125  613  453 -265 -488 -174  -64   -6   91  132  154  154  110   80   55   35    0    0');
    case 32, disp('Correct?ECP32:');
    case 33, disp('Correct?ECP33:');
    case 34, disp('Correct?ECP34:');
    case 35, disp('Correct?ECP35:');
%    case 36, disp('Correct?ECP36:4  -15  -45  108  132 -246 -217  165  631 -443 -840  741  915 -479 -645 -175   42   -4  164  170  175  170  121   95   65   41    0    0');
    case 36, disp('Correct?ECP36:2   -9  -25  118  -10 -183 -119  276  287 -496 -226  584  453 -358 -537 -146    4   51  136  161  175  170  121   95   65   41    0    0');
%    case 37, disp('Correct?ECP37:2  -15   -2  131 -135 -191  246  208 -122 -748  571  386  283 -744 -628  817  737 -493 -599 -181   77   17  182  183  162  174  134  102   69   43    0    0');
    case 37, disp('Correct?ECP37:2   -7   -6   85 -211  114  202   79 -494 -181  334  342  121 -563  -60  494  304 -393 -496 -118   34   72  154  174  162  174  134  102   69   43    0    0');
    case 38, disp('Correct?ECP38:');
    case 39, disp('Correct?ECP39:');
%    case 40, disp('Correct?ECP40:-2    23   -20  -162   336  -174   -61   145    83    12  -656   293  -584  1166   155   385  -904  -880   730   713  -464  -524  -159   132   102   166   159   180   194   156   116    77    47     0     0');
%    case 41, disp('Correct?ECP41:-2    23   -20  -162   336  -174   -61   145    83    12  -656   293  -584  1166   155   385  -904  -880   730   713  -464  -524  -159   132   102   166   159   180   194   156   116    77    47     0     0');
%    case 42, disp('Correct?ECP42:-2    23   -20  -162   336  -174   -61   145    83    12  -656   293  -584  1166   155   385  -904  -880   730   713  -464  -524  -159   132   102   166   159   180   194   156   116    77    47     0     0');
    case 40, disp('Correct?ECP40:-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   77   47    0    0');
    case 41, disp('Correct ECP41:-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   77   47    0    0');
    case 42, disp('Correct ECP42:-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   77   47    0    0');
    case 43, disp('Correct?ECP43:-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   78   49    0    0');
    case 44, disp('Correct ECP44:-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  116   78   49    0    0');
    case 45, disp('Correct ECP45:-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  118   81   51    0    0');
    case 46, disp('Correct ECP46:-2   11  -14  -90  452 -661  -65  440  897 -851 -849  225  411  248  -44  533 -492 -326  262  266 -382 -426  -62   84  157  138  150  180  194  156  120   82   53    0    0');
    case 47, disp('Correct?ECP47:');
    case 48, disp('Correct?ECP48:');
    case 49, disp('Correct?ECP49:');
    case 50, disp('Correct?ECP50:');
    case 51, disp('Correct?ECP51:2   -1  -17  -10   13  438 -784 -144  514  969 -785 -818  172  382  265  -66  411 -555 -287  302  264 -321 -352   21  150  209  176  189  210  219  178  133   90   58    0    0');
%    case 52, disp('Correct?ECP52:6     83   -186  -2026   4319   6885 -13816 -19068  23794  38816 -25693 -48086  13102  37207  -1246 -15393   -378    677  -2924   3851   4221  -1466  -3588  -1783   2218    810   -596   -178    416   -163   -485    -46    182    398    248    251    225    194    195    131     99     66     -1      3');
    case 52, disp('Correct?ECP52:');
    case 53, disp('Correct?ECP53:');
    case 54, disp('Correct?ECP54:');
    case 55, disp('Correct?ECP55:');
    case 56, disp('Correct?ECP56:2    -5   -21    81   -71  -147   646  -387 -1261   328  2428   351 -2206 -1474  1020   430    -1   341   879  -286  -415   150  -131  -102  -141   251  -812  -323   219   256  -161  -220   209   240   269   255   273   281   269   208   151   102    65     0     0');
    case 57, disp('Correct?ECP57:1       6      -9    -104    -469    3112   -9484   45853 -123343   79444  153809 -143354 -145458   36604  191689   40851 -114941  -75156    5651   60757   19477   -9179  -26415   -9337   10514   13195    5874   -2337     611   -2580   -3608    -867   -2883   -1462     349     571     531      61     525     511     537     435     389     349     307     230     163     108      68       0       0');
    case 58, disp('Correct?ECP58:');
    case 59, disp('Correct?ECP59:');
end

end
