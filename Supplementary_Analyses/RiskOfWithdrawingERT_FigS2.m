clc
clear all
close all

%% This code generates the results shown in Figure S2 of the manuscript:
%% Thompson et al. Using real-time modelling to optimise an outbreak response: Insights from the 2017 Ebola outbreak in the Democratic Republic of the Congo.

%% All code was written in MATLAB, compatible with version R2022a.

%% ©2024 Robin Thompson <robin.thompson@maths.ox.ac.uk>


% Set up and plot the incidence data
incidence_data = zeros(1000,1);
incidence_data(1) = 1;
incidence_data(daysact('27-mar-2017',  '18-apr-2017') + 1) = incidence_data(daysact('27-mar-2017',  '18-apr-2017') + 1) + 1;
incidence_data(daysact('27-mar-2017',  '24-apr-2017') + 1) = incidence_data(daysact('27-mar-2017',  '24-apr-2017') + 1) + 1;
incidence_data(daysact('27-mar-2017',  '26-apr-2017') + 1) = incidence_data(daysact('27-mar-2017',  '26-apr-2017') + 1) + 1;
incidence_data(daysact('27-mar-2017',  '11-may-2017') + 1) = incidence_data(daysact('27-mar-2017',  '11-may-2017') + 1) + 1;
incidence_data(daysact('27-mar-2017',  '24-apr-2017') + 1) = incidence_data(daysact('27-mar-2017',  '24-apr-2017') + 1) + 1;
incidence_data(daysact('27-mar-2017',  '30-apr-2017') + 1) = incidence_data(daysact('27-mar-2017',  '30-apr-2017') + 1) + 1;
incidence_data(daysact('27-mar-2017',  '02-may-2017') + 1) = incidence_data(daysact('27-mar-2017',  '02-may-2017') + 1) + 1;

ERTArrivalDay = daysact('27-mar-2017',  '15-may-2017') + 1;


% Discretise a gamma distributed serial interval distribution using the method of
% Cori et al.;
SI_mean = 15.3;
SI_sd = 9.3;

SI_scale = SI_sd^2/SI_mean;
SI_shape = SI_mean/SI_scale;


SI_discrete = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape,SI_scale);
    SI_discrete(k) = trapz(intVals, funcVals);
end




%Now calculate the likelihood of R using the data up until the day before the arrival of the ERT (ERT arrival: day 50; calculation up to day 49) assuming a constant R across this time period
%We assume that the number of cases each day is drawn from a negative binomial distribution, and repeat all of our results for two values of the dispersion parameter, k

for dispersion_k = [0.2 10]

RVals = [0.001:0.001:20];
logL = zeros(length(RVals),1);

for rIt = 1:length(RVals)
    R = RVals(rIt);
    for t = 2:(ERTArrivalDay - 1)
        I_t = incidence_data(t);
        if I_t > 1
        for j = 1:I_t
            logL(rIt) = logL(rIt) - log(j);
        end
        end
        sumV = 0;
        for s = 1:(t - 1)
            sumV = sumV + incidence_data(t - s)*SI_discrete(s);
        end
        p0 = R*sumV/(dispersion_k + R*sumV);
        logL(rIt) = logL(rIt) + log(gamma(dispersion_k + I_t)) + I_t*log(p0) + dispersion_k*log(1 - p0);
    end
end

logL = logL - max(logL);
likelihoodVals = exp(logL);

normfac = trapz(RVals,likelihoodVals);
likelihoodVals = likelihoodVals./normfac; % Normalising the likelihood so that it is a valid pdf

% Now calculate the risk of withdrawing the ERT each day from the day of
% the arrival of the ERT onwards - using the entire likelihood
% function for R (rather than a single value of R). Note that EOO (the end
% of outbreak probability) is 1 - Risk of withdrawing ERT.

logEOO = zeros(ERTArrivalDay + 100,length(RVals));
EOO = zeros(ERTArrivalDay + 100,length(RVals));
for RIt = 1:length(RVals)
    R = RVals(RIt)
    for t = ERTArrivalDay:(ERTArrivalDay + 100)

        logEOO(t,RIt) = 0;
        for j = t:(t + length(SI_discrete))
            sumV = 0;
            for s = 1:(j-1)
                 if s < (length(SI_discrete) + 0.5)
                sumV = sumV + incidence_data(j-s)*SI_discrete(s);
                end
            end
            p0 = R*sumV/(dispersion_k + R*sumV);
            logEOO(t,RIt) = logEOO(t,RIt) + dispersion_k*log(1 - p0);
        end

        EOO(t,RIt) = exp(logEOO(t, RIt));

    end
end
combinedEOO = zeros(ERTArrivalDay + 100,1);
for t = ERTArrivalDay:(ERTArrivalDay + 100)
    funcV = zeros(length(RVals),1);
    for i = 1:length(RVals)
        funcV(i) = EOO(t,i)*likelihoodVals(i);
    end
    combinedEOO(t) = trapz(RVals, funcV);
end


%%% Now plot the figures

serial_first_case = datenum('27-mar-2017');
serials = [serial_first_case:(serial_first_case+61)];
datesHere = datestr(serials);

if dispersion_k == 0.2
figure(1)
bar([1:61], [incidence_data(1:60); 0], 'BarWidth', 1)
hold on
plot([ERTArrivalDay ERTArrivalDay], [0 2], 'k--')
xticks([1:10:61])
xticklabels(datesHere(1:10:61,:))
yticks([0 1 2])
box off
xlabel('Date')
ylabel('Number of cases')

figure(2)
bar([1:length(SI_discrete)], SI_discrete, 'BarWidth', 1)
xlim([0.5 50.5])
xticks([1 [5:5:50]])
box off

figure(3)
plot(RVals, likelihoodVals)
xlim([0 5])
box off

figure(4)
plot([0:56], 1 - combinedEOO(ERTArrivalDay:(ERTArrivalDay + 56)))
serial_first_case = datenum('15-may-2017');
serials = [serial_first_case:(serial_first_case+51)];
datesHere = datestr(serials);
xticks([0:10:50])
xticklabels(datesHere(1:10:51,:))
xlim([0 55])
box off
hold on
plot([48 48], [0 1], 'k--')
end



if dispersion_k == 10

figure(3)
hold on
plot(RVals, likelihoodVals, 'r')
xlim([0 10])
box off

figure(4)
hold on
plot([0:56], 1 - combinedEOO(ERTArrivalDay:(ERTArrivalDay + 56)), 'r')
end

end