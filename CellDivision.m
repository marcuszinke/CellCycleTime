%% cdc2-33 mutant cell - cycle time modeling %%
% Code written by Cell Division team to model the cycle time of cdc2-33 
% mutant cells based on cell birth length.
% Modeling approaches include an analytical solution to the growth rate
% ODE, a two-term pertubation approximation, and a 2nd order polynomial
% best fit equation. 

%% DECLARE CONSTANTS %%
T0 = 105;           %min % Timer phase constant time
lc = 12.75;          %μm % %Critical length
g0 = 0.006;        %Constant expression value given by article
B = 0.0046;        %Constant expression value given by article
lb = 0:0.01:20;    %Birth length
e_lb = [5.2 6.5 7.5 8.5 9.8 9.8 11.8 16.2];  %experimental birth lengths
e_T = [217 149 137 197 161 115 131 109];     %experimental Cycle times


%% ANALYTICAL CYCLE TIME MODEL %%
Tc = T0 + (lc-lb)./(g0+B.*lb);      % Calculate analytical results

Tc_fit = polyfit(e_lb,e_T,2);       % Employ MATLAB ployfit function

%% PERTURBATION MODEL %%

% set constant expressions
a=(g0+B.*lb)./B;
b=g0/B;

% solve for lambertW values
for k2 = 0:0.0001:0.001
    lamW_fxn = (exp(b./a).*(-b.*B+lc.*B+b.*k2))./(a.*k2); % Calculate interior of lambertW()
end
lamW_2 = lambertw(lamW_fxn); %run lambertw() function      
num_lamW = double(lamW_2);   %convert lambertW values from symbolic to numeric values

% Calculate perturbation approximation curves
iter = 1;
for k2 = 0:0.00001:0.0001
    Tc_p(iter, :) = T0 + (a.*k2.*num_lamW-a.*B-b.*k2)./(a.*B.*k2); % General case
    Tc_p_2(iter, :) = T0 + (lc-lb)./(g0+B.*lb+k2.*lb);             % Growth rate fixed by birth length
    iter = iter+1;
end


%% PLOTTING %%

%------%Plot 1: Analytical vs Real Data / Polyfit%------%
% figure; hold all;
% xlim([6 22]);
% ylim([80 230]);
% 
% %Plot experimental data
% x1 = linspace(0,20);
% f1 = polyval(Tc_fit,x1);
% plot(x1,f1,'b--')
% plot(e_lb, e_T,'o', 'LineWidth', 2);
% 
% %Plot analytical curve
% % labelText = '-Analytical Model';
% % labelX = 5;   % Adjust the x-coordinate of the label
% % labelY = 190; % Adjust the y-coordinate of the label
% % text(labelX, labelY, labelText, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'b');
% yline(95,'--', 'LineWidth', 2);
% plot(lb, Tc,'-.','LineWidth', 2);
% 
% legend('Experimental Data - 2nd Order Polynomial Fit','Experimental Data',sprintf('Constant Timer Period, T_0 = %d', T0),'Analytical Approach', 'FontSize', 16);
% xlabel('Birth length (μm)', 'FontSize', 16);
% ylabel('Cycle time (min)', 'FontSize', 16);
% subtitle('Cell cycle time in cdc-33 wee1-6 mutant fission yeast cells', 'FontSize', 16);
% title('Comparision of analytical model and experimental data w/ polyfit', 'FontSize', 20);

%------%Plot 2: Analytical vs Perturbation vs Real%------%
figure; hold all;
xlim([6 22]);
ylim([80 230]);

%Plot experimental data
x1 = linspace(0,20);
f1 = polyval(Tc_fit,x1);
plot(e_lb, e_T, 'o', 'LineWidth', 2);
plot(x1,f1,'b--')

%Plot analytical curve
labelText = '-Analytical model';
labelX = 8.65;   % Adjust the x-coordinate of the label
labelY = 200; % Adjust the y-coordinate of the label
text(labelX, labelY, labelText, 'FontSize', 20,'Color', '#FF8800');
yline(95,'--', 'LineWidth', 2);
plot(lb, Tc,'-.','LineWidth', 3);

%Plot Perturbation curves
labelText = 'Perturbation curves-';
labelX = 6.1;   % Adjust the x-coordinate of the label
labelY = 187.5; % Adjust the y-coordinate of the label
text(labelX, labelY, labelText, 'FontSize', 20, 'Color', '#006400');
% plot(lb, Tc_p, '-', 'LineWidth', 1);
plot(lb, Tc_p_2, '-', 'LineWidth', 0.5);

legend('Experimental Data','Experimental Data - 2nd Order Polynomial Fit',sprintf('Constant Timer Period, T_0 = %d', T0),'Analytical Approach', 'FontSize', 25);
xlabel('Birth length (μm)', 'FontSize', 25);
ylabel('Cycle time (min)', 'FontSize', 25);
title('Cell cycle time vs birth length in cdc-33 wee1-6 mutant fission yeast cells', 'FontSize', 30);
subtitle('Comparision of Analytical model, Perturbation approximation, and Experimental data', 'FontSize', 25);

%------%Plot 3: Perturbation vs Real Data / Polyfit%------%
% figure; hold all;
% xlim([6 22]);
% ylim([80 230]);
% 
% %Plot experimental data
% x1 = linspace(0,20);
% f1 = polyval(Tc_fit,x1);
% plot(x1,f1,'r--')
% plot(e_lb, e_T, 'o');
% 
% %Plot Perturbation curves
% labelText = 'Perturbation curves-';
% labelX = 7;   % Adjust the x-coordinate of the label
% labelY = 185; % Adjust the y-coordinate of the label
% text(labelX, labelY, labelText, 'FontSize', 14, 'Color', '#006400');
% % plot(lb, Tc_p, '-', 'LineWidth', 1);
% plot(lb, Tc_p_2, '-', 'LineWidth', 0.5);
% 
% yline(95,'--', 'LineWidth', 2);
% 
% legend('Experimental Data - 2nd Order Polynomial Fit','Experimental Data',sprintf('Constant Timer Period, T_0 = %d', T0), 'FontSize', 16);
% xlabel('Birth length (μm)', 'FontSize', 16);
% ylabel('Cycle time (min)', 'FontSize', 16);
% subtitle('Cell cycle time in cdc-33 wee1-6 mutant fission yeast cells', 'FontSize', 16);
% title('Comparision of perturbation approx and experimental data', 'FontSize', 20);
