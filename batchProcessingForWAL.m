%% Ne vs temperature
figure
scatter(temForCheck(1:number/2),ne,'o');
xlim([0 1.1*max(temForCheck)])
% upLim=max(ne)+0.1*(max(ne)-min(ne));
% downLim=min(ne)-0.1*(max(ne)-min(ne));
% ylim([downLim upLim])
% clearvars upLim downLim

xlabel {T (K)}
ylabel {n_e (cm^{-3})}
grid on
box on

%% D vs temperature
figure
scatter(temForCheck(1:number/2),dif,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {D (cm^2/s)}
grid on
box on

%% L_SO vs temperature
figure
scatter(temForCheck(1:number/2),lso,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {L_{SO} (m)}
grid on
box on

%% L_Phi vs temperature
figure
scatter(temForCheck(1:number/2),lphi,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {L(\phi) (m)}
grid on
box on

%% L_SO vs D
figure
scatter(dif,lso,'o');
% xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {D (m)}
ylabel {L(\phi) (m^2/s)}
grid on
box on

%% Tau_SO vs Tau_p
figure
scatter(taup,tauso,'o');
% xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {\tau_p (s)}
ylabel {\tau_{SO} (s)}
grid on
box on