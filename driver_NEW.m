%% driver file for Basophil-T cell-cancer cell model


clear all;
%-- Specify parameters
% cancer cells
    params(1,1) = 8.2e-6;           % kpT -- fixed;  double once in 24 hours
    params(2,1) = 3.85e-10;         % kdT -- fixed;  +CD8 alone: 21.85% death (achieves 21.7%)

% Treg
    params(3,1) = 0;                % kpTr
    params(4,1) = 3e-6;             % kdTr -- fixed; basal 41% death when cultured alone for 72h (achieves 40%)
    params(5,1) = 5.6e-5;           % str_TregDeath: strength of basophil effect on Treg death
                                    %           -- fixed; +Baso: 71% death at 72h (achieves 71%)

% CD8 cells
    params(6,1) = 1e-9;             % kpCD8 -- fixed to balance death 
    params(7,1) = 3.71e-10;         % ksupCD8 -- fixed; +Treg: 0.21% cancer cell death in 72h (achieves 0.23%)
    params(8,1) = 5.24e-2;          % str_CD8prolif strength of basophil effect on promoting CD8 proliferation 
                                    %           -- tuned to achieve 19.2% cancer cell death when considering both basophil-mediated mechanisms
    params(9,1) = 1e-9;             % kdCD8 -- fixed to exactly balance production

% Basophils
    params(10,1) = 25000;           % number of basophils



%-- set details
numDays = 2; % 48h beyond D1
tstep = 1000;
stopTime = 24*numDays*3600; %time
maxY = 1e5;
n_diffeqn = 6;
options = odeset('RelTol',1e-12);


%-- initial values
init_cancer = 12500; %cancer cells are added at D1
init_Treg = 6000;
init_CD8 = 25000;
init_Basophil = params(10,1);
initvalue = zeros(n_diffeqn,1);
initvalue(1,1) = init_cancer;  %T


%-- solve ODEs
tspan = [0:tstep:stopTime];
[tsim, results_onlycancer] = ode15s(@core_base,tspan,initvalue,options,params);

initvalue(5,1) = init_CD8;  %CD8
[tsim, results_cancerCD8] = ode15s(@core_base,tspan,initvalue,options,params);

initvalue(5,1) = init_CD8;  %CD8
[tsim, results_cancerCD8Baso] = ode15s(@core_base,tspan,initvalue,options,params);

initvalue(3,1) = init_Treg;  %Treg
initvalue(5,1) = init_CD8;  %CD8
[tsim, results_noB] = ode15s(@core_base,tspan,initvalue,options,params);
[tsim, results_TregDeath] = ode15s(@core_TregDeath,tspan,initvalue,options,params);
[tsim, results_CD8prolif] = ode15s(@core_CD8prolif,tspan,initvalue,options,params);
[tsim, results_TregDeath_CD8prolif] = ode15s(@core_TregDeath_CD8prolif,tspan,initvalue,options,params);


%-- calculate percent cancer cell death
percentages = (1-([results_cancerCD8(end,1);results_cancerCD8Baso(end,1);results_noB(end,1);...
        results_TregDeath(end,1);results_CD8prolif(end,1);...
        results_TregDeath_CD8prolif(end,1)]/init_cancer))*100;

cancerCellCount(:,1:6) = [results_cancerCD8(:,1) results_cancerCD8Baso(:,1) results_noB(:,1)...
        results_TregDeath(:,1) results_CD8prolif(:,1) ...
        results_TregDeath_CD8prolif(:,1)];


%-- plot results
w = 6; % Linewidth for plotting
tsim = tsim/3600;

figure(1)
plot(tsim,results_onlycancer(:,1),'LineWidth',w,'Color','k')
xlim([0 numDays*24])
xlabel('time (hrs)')
ylabel('number of cancer cells')
title('cancer cells')
set(gca,'FontSize',14,'LineWidth',2);



figure(2);
subplot(1,2,1)
whichSpecies = 1; % cancer cells
plot(tsim,results_onlycancer(:,whichSpecies),'LineWidth',w,'Color','k');
hold on; plot(tsim,results_cancerCD8(:,whichSpecies),'LineWidth',w,'Color',[0.64,0.08,0.18]);
hold on; plot(tsim,results_cancerCD8Baso(:,whichSpecies),'LineWidth',w,'Color',[0 0.8 0.8]);
hold on; plot(tsim,results_noB(:,whichSpecies),'LineWidth',w,'Color',[0.5,0.5,0.5]);
hold on; plot(tsim,results_TregDeath(:,whichSpecies),'LineWidth',w,'Color',[0.49,0.18,0.56]);
hold on; plot(tsim,results_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.07,0.62,1.00]);
hold on; plot(tsim,results_TregDeath_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.47,0.67,0.19]);
hold off;
xlim([0 numDays*24])
xlabel('time (hrs)')
ylabel('number of cancer cells')
legend('cancer cells only', 'cancer,CD8','cancer,CD8,Baso(no effect)','cancer,CD8,Treg',...
    'cancer,CD8,Treg,Baso(+Treg death)','cancer,CD8,Treg,Baso(+CD8 prolif)',...
    'cancer,CD8,Treg,Baso(+Treg death,+CD8 prolif)')
title('cancer cells')
set(gca,'FontSize',14,'LineWidth',2);



subplot(1,2,2)
whichSpecies = 5; % CD8 T cells
plot(tsim,results_onlycancer(:,whichSpecies),'LineWidth',w,'Color','k');
hold on; plot(tsim,results_cancerCD8(:,whichSpecies),'LineWidth',w,'Color',[0.64,0.08,0.18]);
hold on; plot(tsim,results_cancerCD8Baso(:,whichSpecies),'LineWidth',w,'Color',[0 0.8 0.8]);
hold on; plot(tsim,results_noB(:,whichSpecies),'LineWidth',w,'Color',[0.5,0.5,0.5]);
hold on; plot(tsim,results_TregDeath(:,whichSpecies),'LineWidth',w,'Color',[0.49,0.18,0.56]);
hold on; plot(tsim,results_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.07,0.62,1.00]);
hold on; plot(tsim,results_TregDeath_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.47,0.67,0.19]);
hold off;
xlim([0 numDays*24])
xlabel('time (hrs)')
ylabel('number of CD8 T cells')
title('CD8')
set(gca,'FontSize',14,'LineWidth',2);


figure(3)
subplot(2,2,1)
whichSpecies = 3; % Treg
plot(tsim,results_onlycancer(:,whichSpecies),'LineWidth',w,'Color','k');
hold on; plot(tsim,results_cancerCD8(:,whichSpecies),'LineWidth',w,'Color',[0.64,0.08,0.18]);
hold on; plot(tsim,results_cancerCD8Baso(:,whichSpecies),'LineWidth',w,'Color',[0 0.8 0.8]);
hold on; plot(tsim,results_noB(:,whichSpecies),'LineWidth',w,'Color',[0.5,0.5,0.5]);
hold on; plot(tsim,results_TregDeath(:,whichSpecies),'LineWidth',w,'Color',[0.49,0.18,0.56]);
hold on; plot(tsim,results_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.07,0.62,1.00]);
hold on; plot(tsim,results_TregDeath_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.47,0.67,0.19]);
hold off;
xlim([0 numDays*24])
xlabel('time (hrs)')
ylabel('number of cells')
title('Treg')
ylim([-500 7500])
set(gca,'FontSize',14,'LineWidth',2);


subplot(2,2,2)
whichSpecies = 5; % CD8 T cells
plot(tsim,results_onlycancer(:,whichSpecies),'LineWidth',w,'Color','k');
hold on; plot(tsim,results_cancerCD8(:,whichSpecies),'LineWidth',w,'Color',[0.64,0.08,0.18]);
hold on; plot(tsim,results_cancerCD8Baso(:,whichSpecies),'LineWidth',w,'Color',[0 0.8 0.8]);
hold on; plot(tsim,results_noB(:,whichSpecies),'LineWidth',w,'Color',[0.5,0.5,0.5]);
hold on; plot(tsim,results_TregDeath(:,whichSpecies),'LineWidth',w,'Color',[0.49,0.18,0.56]);
hold on; plot(tsim,results_TregDeath_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.47,0.67,0.19]);
hold off;
xlim([0 numDays*24])
xlabel('time (hrs)')
ylabel('number of cells')
title('CD8')
set(gca,'FontSize',14,'LineWidth',2);



subplot(2,2,3)
whichSpecies = 4; % dead Treg
plot(tsim,results_onlycancer(:,whichSpecies),'LineWidth',w,'Color','k');
hold on; plot(tsim,results_cancerCD8(:,whichSpecies),'LineWidth',w,'Color',[0.64,0.08,0.18]);
hold on; plot(tsim,results_cancerCD8Baso(:,whichSpecies),'LineWidth',w,'Color',[0 0.8 0.8]);
hold on; plot(tsim,results_noB(:,whichSpecies),'LineWidth',w,'Color',[0.5,0.5,0.5]);
hold on; plot(tsim,results_TregDeath(:,whichSpecies),'LineWidth',w,'Color',[0.49,0.18,0.56]);
hold on; plot(tsim,results_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.07,0.62,1.00]);
hold on; plot(tsim,results_TregDeath_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.47,0.67,0.19]);
hold off;
xlim([0 numDays*24])
xlabel('time (hrs)')
ylabel('number of cells')
title('dead Treg')
set(gca,'FontSize',14,'LineWidth',2);


subplot(2,2,4)
whichSpecies = 6; % dead CD8 T cells
plot(tsim,results_onlycancer(:,whichSpecies),'LineWidth',w,'Color','k');
hold on; plot(tsim,results_cancerCD8(:,whichSpecies),'LineWidth',w,'Color',[0.64,0.08,0.18]);
hold on; plot(tsim,results_cancerCD8Baso(:,whichSpecies),'LineWidth',w,'Color',[0 0.8 0.8]);
hold on; plot(tsim,results_noB(:,whichSpecies),'LineWidth',w,'Color',[0.5,0.5,0.5]);
hold on; plot(tsim,results_TregDeath(:,whichSpecies),'LineWidth',w,'Color',[0.49,0.18,0.56]);
hold on; plot(tsim,results_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.07,0.62,1.00]);
hold on; plot(tsim,results_TregDeath_CD8prolif(:,whichSpecies),'LineWidth',w,'Color',[0.47,0.67,0.19]);
hold off;
xlim([0 numDays*24])
xlabel('time (hrs)')
ylabel('number of cells')
title('dead CD8')
set(gca,'FontSize',14,'LineWidth',2);




