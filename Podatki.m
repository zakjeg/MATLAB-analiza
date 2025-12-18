
%  ***********************************
%   PODATKI O GEOMETRIJI KONSTRUKCIJE
%  ***********************************



%  ********************
%   Koordinate vozlišè
%  ********************
% ------------------------------------------------------------------------------
% voz=[1  x-1 z-1 
%      :   :   :
%      i  z-i z-i
%      :   :   :
%      N  x-N z-N ]
%
% x-i: x koordinata i-tega vozlišèa (v metrih) v globalnem koordinatnem sistemu 
% z-i: z koordinata i-tega vozlišèa (v metrih) v glob. koord. sistemu
% N: število vseh vozlišè
% ------------------------------------------------------------------------------

voz=[
1 0.00 0.
2 0.080 0.
3 0.3535 0.
4 0.6270 0.
5 0.900 0.
6 1.1730 0.
7 1.4465 0.
8 1.7200 0.
9 1.80 0.
];



%   *****************
%    Potek elementov
%   *****************
% ------------------------------------------------------------------------------
%  pq=[ 1  p-1 q-1 
%       :   :   :
%       i  p-i q-i
%       :   :   :
%       m  p-m q-m ]
%
% i-ti element poteka od "p-i"-tega do "q-i"-tega vozlišèa 
% ------------------------------------------------------------------------------
pq=[
1 1 2
2 2 3
3 3 4
4 4 5
5 5 6
6 6 7
7 7 8
8 8 9
];
 
 
%   **************************
%    Preèni prerezi elementov
%   **************************
% ------------------------------------------------------------------------------
%  pr=[ e-1 pr-1 
%       e-2 pr-2
%        :   :
%       p-i pr-j
%        :   :
%       e-m pr-p ]
%
% i-ti element "e-i" ima j-ti preèni prerez "pr-j"; 
% število preènih prerezov "p" je lahko enako ali pa manjše številu elementov "m"
% V primeru, da so vsi preèni prerezi enaki, ni potrebno nièesar podajati.
% ------------------------------------------------------------------------------
pr=[ 
];


%   *************************************
%    Temperaturna polja preènih prerezov
%   *************************************
% ------------------------------------------------------------------------------
%  tpr=[ pr-1 tpr-1 
%        pr-2 tpr-2
%         :   :
%        pr-i tpr-j
%        :   :
%        pr-p tpr-t ]
%
% i-temu preènemu prerezu "pr-i" pripada j-to temperaturno polje "tpr-j"; 
% število temperaturnih polj "t" je lahko enako ali pa manjše številu preènih prerezov "p"
% V primeru, da imajo vsi prerezi enako temperaturno polje, ni potrebno nièesar podajati.
% ------------------------------------------------------------------------------
tpr=[ 
];

% 1- ali 2-dimenzionalno temperaturno polje
razp_T=2;


%   **********************
%    Sprostitve elementov
%   **********************
%  Sprostitev elementa pomeni, da je doloèena komponenta sile ali momenta 
% v krajišèih elementa enaka niè. Sprošèene komponente se nanašajo na LKS elementa.
% ------------------------------------------------------------------------------
%  spr=[ i s-i(1) ... s-i(6)];
%
%      s-i(j)=1 -> ni sprostitve
%      s-i(j)=0 -> sprostitev
%
% i-ti element ima 6 komponent krajišènih sil ali momentov s-i
%    s-i(1): N(0)
%    s-i(2): Q(0)
%    s-i(3): M(0)
%    s-i(4): N(L)
%    s-i(5): Q(L)
%    s-i(6): M(L)
%     
% ------------------------------------------------------------------------------
spr=[
   %1 1 1 0 1 1 1
   %2 1 1 0 1 1 0
   %3 1 1 1 1 1 0
];

%   **********************
%    Podpore konstrukcije
%   **********************
% ------------------------------------------------------------------------------
%  fix=[ fix(1) ... fix(n)];
%
% Konstrukcija je podprta v n-prostostnih stopnjah glede na GKS
%    i-to vozlišèe konstrukcije 
%		ima sledeèe prostostne stopnje: U(i),W(i)in FI(i)
%  U(i) je pomik v smeri X-osi
%  W(i) je pomik v smeri Z-osi
%  FI(i) je zasuk okrog Y-osi
% ------------------------------------------------------------------------------

fix=[4 5 23]; 



% ****************************************
%  Naèin in stopnja numeriène integracije 
% ****************************************
%   nacin=1 -> Gauss
%   nacin=2 -> Lobatto
nacin=2;
stopnja=5;


% *****************************************************
%  Stopnja interpolacijskega polinoma za EPSILON (N-1) 
% *****************************************************
% interpolacija specifiène spremembe dolžine (epsilona) 
%    z interpolacijskim polinomom stopnje N-1
%             N
% epsilon(x)=SUM[epsilon_n * P_n(x)]
%            n=1

polinom(1)=4;


% ***************************************************
%  Stopnja interpolacijskega polinoma za KAPPO (M-1) 
% ***************************************************
% interpolacija psevdoukrivljenosti (kappe) 
%    z interpolacijskim polinomom stopnje M-1
%           M
% kappa(x)=SUM[kappa_m * P_m(x)]
%          m=1

polinom(2)=4;


% elementi s konstantnim potekom psevdoukrivljenosti (kappe)
kap_const=[];


% *********************************************************
%  PODATKI O RAÈUNU ODZIVA KONSTRUKCIJE NA POŽARNO OBTEŽBO 
% *********************************************************
% zaèetni èasovni korak (izhodišèni inkrement = 0)
zacetek=0;

%ORIGINAL
% število èasovnih korakov
% st_korakov=24*3600*192/(3600*4)+1;

%ZMANJSANO ST KORAKOV ZA TEST
st_korakov=20;

% èasovni prirastek
ink_cas=4*3600;  %v sekundah
% avt=1 -> èasovni inkrement se izbira samodejno glede na podane èase
% avt=0 -> èasovni inkrement izbiramo sami
avt=0;

% projektni obtežni faktor
lambda_mej=1e800;
odmik=5e-3;

% èasi posameznih èasovnih korakov (v minutah)
cas=[0:3600*4:192*24*3600];




% **************
%  Krèenje lesa
% **************
% krcenje_w=0 ->  raèun brez upoštevanja krèenja
% krcenje_w=1 ->  v raèunu se upošteva krèenje
krcenje_w=1;


% **************
%  Lezenje lesa 
% **************
% lezenje_w=0 ->  raèun brez lezenja lesa
% lezenje_w=1 ->  v raèunu se upošteva lezenje lesa
lezenje_w=1;
% upoštevan model
%mod_cw=1 -> model A
%mod_cw=2 -> model B
mod_cw=1;

% ******************************
%  Mehano-sorptivna deformacija
% ******************************
% mehso_w=0 ->  raèun brez upoštevanja mehanosorptivne def.
% mehso_w=1 ->  v raèunu se upošteva mehanosorptivna def.
mehso_w=1;


% ********************************
%  Èasovno odvisni obtežni faktor 
% ********************************
% faktor velikosti spremenljive obtežbe za posamezen èasovni inkrement
%   glede na zaèetno vrednost
obt_faktor=ones(1,size(cas,2));
% obt_faktor(1:200)=1:1:200;




