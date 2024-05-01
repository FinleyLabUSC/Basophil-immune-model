function dydt = core_TregDeath(t,y,params)

% Specify state variables
T           = y(1);
T_dead      = y(2);
Treg        = y(3);
Treg_dead   = y(4);
CD8         = y(5);
CD8_sup     = y(6);


% Specify parameters
kpT             = params(1,1); % kpT
kdT             = params(2,1); % kdT
kpTr            = params(3,1); % kpTr
kdTr            = params(4,1); % kdTr
str_TregDeath   = params(5,1); % str_TregDeath
kpCD8           = params(6,1); % kpCD8
ksupCD8         = params(7,1); % ksupCD8
str_CD8supp     = params(8,1); % str_CD8supp
kdCD8           = params(9,1); % kdCD8
Baso            = params(10,1); % number of basophils


% Specify ODEs
dydt(1,1) = kpT*T - kdT*T*CD8;
dydt(2,1) = kdT*T*CD8;
dydt(3,1) = kpTr*Treg - kdTr*Treg*(1+Baso*str_TregDeath);
dydt(4,1) = kdTr*Treg*(1+Baso*str_TregDeath);
dydt(5,1) = kpCD8*CD8 - ksupCD8*CD8*(Treg) - kdCD8*CD8;
dydt(6,1) = ksupCD8*CD8*(Treg);


return
