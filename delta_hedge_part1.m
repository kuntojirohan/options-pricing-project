clear;
K =1120;
TC = 1e-3;
dt = 1/252;
% Time to expiry on each date to expiry
timeToExp =[21;20;19;18;17;16;15;14;13;12;11;10;9;8;7;6;5;4;3;2;1]/252;
%Index on each date to expiry

St = [1134.609985;1128.839966;1128.170044;1129.439941;...
      1145.199951;1139.319946;1140.530029;1148.160034;...
      1150.569946;1141.810059;1132.170044;1126.209961;...
      1127.000000;1122.469971;1108.060059;1109.189941;...
      1091.329956;1093.949951;1095.400024;1109.780029;...
      1122.319946];

% St(6:end) = 1503; 
subplot(2,1,1),plot(St); hold all; plot(0*St+K);title('Index Level');
axis([-inf,inf,-inf,inf]); grid on;
r = 0.0110362;
y=0.01632286;
Ca = 22.92;
IV = 0.17058;
N = numel(St)-1;
deltah = zeros(N-1,1);
fprintf('------------------------------------------------------------------\n');
fprintf('Delta Hedging \n');
fprintf('------------------------------------------------------------------\n');

% Call premium

CallPremium = Ca * (1-TC);

BankBalance = CallPremium;
stockCosts=0;
% Get the vector of delta positions
for t =1:N
    S0 = St(t);
    t2e = timeToExp(t);
    [ cprice, deltah(t), ~ ] = blackScholesCallPrice( K, t2e, S0, r, y, IV ); 
    if t>1
        oldStockPosition = deltah(t-1);
    else
        oldStockPosition = 0;
    end
    amtBuy = (deltah(t)-oldStockPosition)*St(t);
    stockCosts = amtBuy +abs(amtBuy)*TC;
    BankBalance = BankBalance*exp(r*dt) - stockCosts;

    fprintf(' t= %i; delta = %.2f; bought $ %.2f of the index; Bank $ %.2f \n', t,deltah(t),amtBuy,BankBalance);
   
end

subplot(2,1,2),plot(deltah(1:end));
axis([-inf,inf,-inf,inf]); grid on;
title('Delta Position');
fprintf('------------------------------------------------------------------\n');

CallPayoff = max(0,St(end)-K);
profit = deltah(end)*St(end) +BankBalance - CallPayoff;

fprintf('Total Hedge Profit: $ %.2f \n', profit);
fprintf('------------------------------------------------------------------\n');



%% Black Scholes Merton Formula in Matlab - Call Option

function [ cprice, delta, gamma ] = blackScholesCallPrice( K, T, S0, r, y, sigma )
numerator = log(S0./K) + (r-y+0.5*sigma.^2).*T;
denominator = sigma.*sqrt(T);
d1 = numerator./denominator;
d2 = d1 - denominator;
cprice = S0 *exp(-y*T).* normcdf(d1) - exp(-r.*T).*K.*normcdf(d2);
delta = exp(-y*T)*normcdf(d1);
gamma = exp(-y*T)*normpdf(d1) ./ (S0.*denominator);

end