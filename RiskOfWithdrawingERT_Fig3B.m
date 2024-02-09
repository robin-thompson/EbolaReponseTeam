clc
clear all
close all

%% This code generates Figure 3B of the manuscript:
%% Thompson et al. Using real-time modelling to optimise an outbreak response: Insights from the 2017 Ebola outbreak in the Democratic Republic of the Congo.

%% All code was written in MATLAB, compatible with version R2022a.

%% ©2024 Robin Thompson <robin.thompson@maths.ox.ac.uk>


% Set up the incidence data
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


%% Now calculate the likelihood of R using the data up until the day before the arrival of the ERT (ERT arrival: day 50; calculation up to day 49) assuming a constant R across this time period
RVals = [0.01:0.001:20];
logL = zeros(length(RVals),1);

for rIt = 1:length(RVals)
    R = RVals(rIt);
    for t = 2:(ERTArrivalDay - 1)
        k = incidence_data(t);
        lambda = 0;
        for s = 1:(t - 1)
            lambda = lambda + R*incidence_data(t - s)*SI_discrete(s);
        end
        logL(rIt) = logL(rIt) + k*log(lambda) - lambda - log(factorial(k));
    end
end

logL = logL - max(logL);
likelihoodVals = exp(logL);

normfac = trapz(RVals,likelihoodVals);
likelihoodVals = likelihoodVals./normfac; % Normalising the likelihood so that it is a valid pdf


% Now find the percentile values of R from the likelihood function, and
% calculate the corresponding risk of withdrawing the ERT each day for each
% individual value of R. Note that the Risk of withdrawing the ERT = 1 - EOO.
percentileVals = [0.5 0.6 0.7 0.8 0.9 0.95 0.99];
hold on
plot([0.5 (length(percentileVals) + 0.5)], [98 98], 'k--')

for pIt = 1:length(percentileVals)

    RVals = [0.01:0.001:20];

    percentileVal = percentileVals(pIt);
    val = 0;
    RIt = 2;
    while val < percentileVal
        val = trapz(RVals(1:RIt), likelihoodVals(1:RIt));
        RIt = RIt + 1;
    end
    RIt = RIt - 1;
    REstimate = (RVals(RIt) + RVals(RIt - 1))/2;

    RVals = REstimate;

    logEOO = zeros(ERTArrivalDay + 100,length(RVals));
    EOO = zeros(ERTArrivalDay + 100,length(RVals));
    for RIt = 1:length(RVals)
        R = RVals(RIt)
        for t = ERTArrivalDay:(ERTArrivalDay + 100)

            logEOO(t,RIt) = 0;
            for j = t:(t + length(SI_discrete))
                for s = 1:(j-1)
                    if s < length(SI_discrete + 0.5)
                        logEOO(t,RIt) = logEOO(t,RIt) - R*incidence_data(j-s)*SI_discrete(s);
                    end
                end
            end

            EOO(t,RIt) = exp(logEOO(t, RIt));

        end
    end

    % Now find the days on which these end of outbreak probability thresholds
    % are exceeded. Note again that the Risk of withdrawing the ERT = 1 - EOO.
    riskThresholds = [0.9 0.95 0.99];
    dayExceeded = zeros(length(riskThresholds),1);
    for riskIt = 1:length(riskThresholds)

        foundIt = 0;
        for i = 1:length(EOO)
            if foundIt == 0 && EOO(i) > riskThresholds(riskIt)

                foundIt = 1;
                dayExceeded(riskIt) = i;

            end
        end

        % Plot the results
        hold on
        if riskIt == 1
            plot([(pIt - 0.25) (pIt + 0.25)], [dayExceeded(riskIt) dayExceeded(riskIt)], 'b')
        elseif riskIt == 2
            plot([(pIt - 0.25) (pIt + 0.25)], [dayExceeded(riskIt) dayExceeded(riskIt)], 'g')
        else
            plot([(pIt - 0.25) (pIt + 0.25)], [dayExceeded(riskIt) dayExceeded(riskIt)], 'r')
        end
        serial_first_case = datenum('27-mar-2017');
        serials = [serial_first_case:(serial_first_case+200)];
        datesHere = datestr(serials);
        ylim([80 105])
        yticks([80:5:105])
        yticklabels(datesHere(80:5:105,:))
        box off
        xlim([0.5 (length(percentileVals) + 0.5)])
    end

end
