%% Preparing Data

optionsdata = readtable('/Users/shraddha/Documents/UCD/Trimester1/Block2/Derivative Securities/Assignment/OptionData Group6.xlsx');

%Converting 'Put = 1 and Call = 0' and Delete (Y/N) as categorical
%variables
optionsdata.PutCall = categorical(optionsdata.PutCall);
optionadata.Delete = categorical(optionsdata.Delete);

% Deleting Rows with two strikes
optionsdata = optionsdata(optionsdata.Delete == "N",:);

putoptions = optionsdata( optionsdata.PutCall == '1',:);
calloptions = optionsdata(optionsdata.PutCall == '0',:);

%% Using blackScholesCallPrice function to get the Price for Call Options
% K: Strike Price = data.Strike
% T: expiry time (in years) = 1/12
% S0: Current Stock Price/Spot Rate = 1122.319946 (Adjusted Closing Price)
% r: annualised risk-free interest rate = 1.10362% (34 days)
% y: Dividend Yield = 1.632286% 
% sigma: Volatitlity of the stock =  data.ImpliedVol

calloptions.CallOptionsPrice = blackScholesCallPrice(calloptions.Strike, 1/12 , 1122.319946 , 0.0110362 , 0.01632286, calloptions.ImpliedVol );

%% Calculating Mid-point of Bid-Ask Spread

s = height(calloptions);

calloptions.BidAskMid = (calloptions.BidPrice + calloptions.AskPrice)/2;
format short
calloptions.DerivedImpliedVol = zeros(s, 1); 
DerivedImpliedVol = matrix(s,6)

%% Computing ImpliedVol from Black-Scholes Model


DerivedImpliedVol = fmincon(@(sigma) blackScholesCallPrice(calloptions.Strike, 1/12 , 1122.319946 , 0.0110362 , 0.01632286, calloptions.ImpliedVol),[calloptions.Strike, 1/12 , 1122.319946 , 0.0110362 , 0.01632286 , calloptions.ImpliedVol],[],[],[],[],[],[],@constraint);
calloptions.DerivedImpliedVol = DerivedImpliedVol(1,6);

function[c,ceq]= constraint (sigma)
c = calloptions.BidAskMid - blackScholesCallPrice(calloptions.Strike, 1/12 , 1122.319946 , 0.0110362 , 0.01632286, sigma );
ceq=[];
end


%% Functions

% Black Scholes Merton Formula in Matlab

function [ cprice, delta, gamma ] = blackScholesCallPrice( K, T, S0, r, y, sigma )
numerator = log(S0./K) + (r-y+0.5*sigma.^2).*T;
denominator = sigma.*sqrt(T);
d1 = numerator./denominator;
d2 = d1 - denominator;
cprice = S0 *exp(-y*T).* normcdf(d1) - exp(-r.*T).*K.*normcdf(d2);
delta = normcdf(d1);
gamma = normpdf(d1) ./ (S0.*denominator);

end

