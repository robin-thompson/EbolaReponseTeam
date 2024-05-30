clc
clear all
close all

%% This code generates Figure 3A of the manuscript:
%% Thompson et al. Using real-time modelling to inform the 2017 Ebola outbreak response in DR Congo. Nature Communications, 2024.

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


%Now calculate the likelihood of R using the data up until the day before the arrival of the ERT (ERT arrival: day 50; calculation up to day 49) assuming a constant R across this time period
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

for pIt = 1:length(percentileVals)

    RVals = [0.01:0.001:20];

    percentileVal = percentileVals(pIt);
    val = 0;
    RIt = 2;
    while val < percentileVal
        val = trapz(RVals(1:RIt), likelihoodVals(1:RIt))
        RIt = RIt + 1;
    end
    RIt = RIt - 1;
    REstimate = RVals(RIt)

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


    % Now plot the figure
    serial_first_case = datenum('27-mar-2017');
    serials = [serial_first_case:(serial_first_case+61)];
    datesHere = datestr(serials);

    hold on

    plot([0:56], 1 - EOO(ERTArrivalDay:(ERTArrivalDay + 56)))
    serial_first_case = datenum('15-may-2017');
    serials = [serial_first_case:(serial_first_case+51)];
    datesHere = datestr(serials);
    xticks([0:10:50])
    xticklabels(datesHere(1:10:51,:))
    xlim([0 55])
    box off
    hold on
    plot([48 48], [0 1], 'k--')
    yticks([0:0.2:1])

end
