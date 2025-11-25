

%period_length: Length of the actual test: total number of days

period_length = 1000
evaluation_period = 1000; %when evaluating, you may not use the entire 1000 days depending on the number of timestamps.
maxEvaluationPeriod = 1000; 
%-----------------------------------------------------------------



%For the 4th paper, we use notransform, NaN, adj = 1

%correl_opts: how we transform the correlation matrix.
% 
correl_opts = {'notransform', 'subplus', 'absolute value', 'subplus negative elements'};
thresholdIndices = {NaN, 0.3,0.5};
shrinkage_constants = {0.1,NaN,0.1,0.3,0.5};
adjs = {'1','2','3','4','5','6'};
alphas = {0.5,0.7,0.9};
measures = {'deg'}   ;
number_runs = {5, 10, 15};
no_stocks = {20};

stocks = {'c', 'p'}; %We always compute the central/peripheral portfolios

% %initalise the column names for the results table
VarNames =   ["Measure",    "runs",     "adj" ,     "cp" ,     "tv" ,     "alp" ,     "No_stock" ,     "ShrinkCon" ,     "CorrOpts" ,     "PREqualWeights" ,     "PRminvar_long_all" ,     "PRminvar_short_all" ,     "EREqualWeights" ,     "ERminvar_long_all" ,     "ERminvar_short_all" ,     "SDEqualWeights" ,     "SDminvar_long_all" ,     "SDminvar_short_all" ,     "VaR005EqualWeights" ,     "VaR005minvar_long_all" ,     "VaR005minvar_short_all" ,     "CVaR005EqualWeights" ,     "CVaR005minvar_long_all" ,     "CVaR005minvar_short_all" ,     "SREqualWeights" ,     "SRminvar_long_all" ,     "SRminvar_short_all" ,     "InvestedStocks" ,     "meanIndexEqual" ,      "meanIndexminvar_long_all" ,     "meanIndexminvar_short_all" ,     "Index" ] ;



%data folder path
root_path = split(cd,'\');
data_path = join([join(root_path(1:end), "\"), "\Data"], "");
addpath(genpath(data_path))


%Load additional helper methods
methods_path = join([join(root_path(1:end), '\'), "\Methods"], "");
addpath(genpath(methods_path))

%Add a save path for the results
save_path =  join([join(root_path(1:end), '\'), "\Results\sp500"], "");

%Load the table of stock prices
T = readtable('Daily S&P.xlsx',  'VariableNamingRule', 'preserve' );
full_data_length = 6131;

%Extract the test and evaluation periods
if evaluation_period~= 0
    T = T(end-period_length-evaluation_period+1:end- evaluation_period,:);
else
    T = T(end-maxEvaluationPeriod:end-maxEvaluationPeriod+period_length,:);
end

stockNames =  T.Properties.VariableNames(2:end);
numberStocks = size(T,2)-1;
Date = T.Date;



%Make the table of returns 
changeTable = table(Date);
for stock = 1:numberStocks
    %logarithmic returns ie. ln(p_t/p_t-1) = ln(p_t) - ln(p_t-1)
    changeTable.(stockNames{stock}) = [0;100*diff(log(T.(stockNames{stock})))];
end




marketAvg = [];


arrayChangeTable = table2array(changeTable(:,2:end));





numberExperimentsTotal = size(correl_opts,2)*size(thresholdIndices,2)*size(shrinkage_constants,2)*size(adjs,2)*size(measures,2)*size(alphas,2)*size(no_stocks, 2)*size(stocks,2)*size(number_runs,2);




%Replace the nans and infs that arise in the changeTable and then compute the average for the market average..
for i = 1:size(changeTable,1)

    returnRow = arrayChangeTable(i,:);
    infs = ~isinf(returnRow);
    returnRow = returnRow(infs);
    nans = ~isnan(returnRow);
    returnRow = returnRow(nans);
    
    marketAvg = [marketAvg; mean(returnRow)];
end
% % %----------------------------------------------------------------------

%To identify trading periods
tempTable = table2array(T(:,2:end));
colsums = sum(abs(tempTable));
entirely_forbidden = find(colsums==0);
partially_allowed = find(colsums~=0);
stockStarts  = zeros(1,size(tempTable,2));
stockStarts(entirely_forbidden) = inf;

for i = 1:size(partially_allowed,2)
    if partially_allowed(i)== 206
            partially_allowed(i)
    end
    
    nzs = find(tempTable(:,partially_allowed(i)) == 0);
    if size(nzs,1) == 0
        stockStarts(partially_allowed(i)) =1;
    else
        stockStarts(partially_allowed(i)) = max(nzs)+2; 
    end
end



graph_options = {number_runs ,thresholdIndices ,correl_opts,shrinkage_constants ,adjs} ;
data_options = {changeTable, evaluation_period, stockStarts};
analysis_options = {no_stocks, alphas, measures};
tic
numberExperimentsRan= 0;
numberExperiments = size(graph_options{1},2)*size(graph_options{2},2)*size(graph_options{3},2)*size(graph_options{4},2)*size(graph_options{5},2);      



% numberExperimentsTotal
TMeasure = cell(numberExperimentsTotal,1);
Truns = zeros(numberExperimentsTotal,1);
Tadj = cell(numberExperimentsTotal,1);
Tcp = cell(numberExperimentsTotal,1);
Ttv = zeros(numberExperimentsTotal,1);
Talp = zeros(numberExperimentsTotal,1);
TNo_stock = zeros(numberExperimentsTotal,1);
TShrinkCon = zeros(numberExperimentsTotal,1);
TCorrOpts = cell(numberExperimentsTotal,1);
TPREqualWeights =  zeros(numberExperimentsTotal,1);
TPRminvar_long_all =  zeros(numberExperimentsTotal,1);
TPRminvar_short_all =  zeros(numberExperimentsTotal,1);
TEREqualWeights =  zeros(numberExperimentsTotal,1);
TERminvar_long_all =  zeros(numberExperimentsTotal,1);
TERminvar_short_all =  zeros(numberExperimentsTotal,1);
TSDEqualWeights =  zeros(numberExperimentsTotal,1);
TSDminvar_long_all =  zeros(numberExperimentsTotal,1);
TSDminvar_short_all =  zeros(numberExperimentsTotal,1);
TVaR005EqualWeights =  zeros(numberExperimentsTotal,1);
TVaR005minvar_long_all =  zeros(numberExperimentsTotal,1);
TVaR005minvar_short_all =  zeros(numberExperimentsTotal,1);
TCVaR005EqualWeights =  zeros(numberExperimentsTotal,1);
TCVaR005minvar_long_all =  zeros(numberExperimentsTotal,1);
TCVaR005minvar_short_all =  zeros(numberExperimentsTotal,1);
TSREqualWeights =  zeros(numberExperimentsTotal,1);
TSRminvar_long_all = zeros(numberExperimentsTotal,1);
TSRminvar_short_all =  zeros(numberExperimentsTotal,1);
TAvg_equal_weights =  zeros(numberExperimentsTotal,1);
TAvg_minvar_long_all =  zeros(numberExperimentsTotal,1);
TAvg_minvar_short_all =  zeros(numberExperimentsTotal,1);
TInvestedStocks = cell(numberExperimentsTotal,1);
TmeanIndexEqual =  zeros(numberExperimentsTotal,1);
TmeanIndexminvar_long_all =  zeros(numberExperimentsTotal,1);
TmeanIndexminvar_short_all  =  zeros(numberExperimentsTotal,1);
TIndex = cell(numberExperimentsTotal, 1);
marketSR = (mean(marketAvg)*252)/(std(marketAvg)*sqrt(252))
marketER = mean(marketAvg)*252
marketSD = std(marketAvg)*sqrt(252)

for i_1 = 1:size(graph_options{1},2)
    number_run = graph_options{1}{i_1};
    for i_2 = 1:size(graph_options{2},2)  
        thresholdIndex = graph_options{2}{i_2}; 
        for i_3 = 1:size(graph_options{3},2)
            correl_opt = graph_options{3}{i_3};
            for i_4 = 1:size(graph_options{4},2)
                shrinkage_constant = graph_options{4}{i_4};
                for i_5 = 1:size(graph_options{5},2)
                    adj = graph_options{5}{i_5};


                                                                                         [cum_sum,  invested_in, ER, SD, SR, VaR005,CVaR005,Weights_minvarlong,Weights_minvarshort, indices,indicesavg] =simulation(number_run, thresholdIndex, correl_opt, shrinkage_constant, adj, data_options, analysis_options);

                                                                                       %Now we loop over the remaining options

                                                                                        for m = 1:size(measures,2)
                                                                                            for a = 1:size(alphas, 2)
                                                                                                for No = 1:size(no_stocks,2)
                                                                                                    for stock = 1:2
                                                                                                        %stock = 1: central
                                                                                                        %stock = 2: periphery
                                                                                                        nthrow = numberExperimentsRan+1;
                                                                                                        TMeasure{nthrow}=  string(measures{m});
                                                                                                        Truns(nthrow) = number_run;
                                                                                                        Tadj{nthrow} = adj;
                                                                                                        Tcp{nthrow} = string(stocks{stock});
                                                                                                        Ttv(nthrow) = thresholdIndex;
                                                                                                        Talp(nthrow) = alphas{a};
                                                                                                        TNo_stock(nthrow) = no_stocks{No};
                                                                                                        TShrinkCon(nthrow) = shrinkage_constant;
                                                                                                        TCorrOpts{nthrow} = string(correl_opt);
                                                                                                        TPREqualWeights(nthrow) = cum_sum{m,a,No,stock}(1);
                                                                                                        TPRminvar_long_all(nthrow) = cum_sum{m,a,No,stock}(2);
                                                                                                        TPRminvar_short_all(nthrow) = cum_sum{m,a,No,stock}(3);
                                                                                                        TEREqualWeights(nthrow) = ER{m,a,No,stock}(1);
                                                                                                        TERminvar_long_all(nthrow) = ER{m,a,No,stock}(2);
                                                                                                        TERminvar_short_all(nthrow) = ER{m,a,No,stock}(3);
                                                                                                        TSDEqualWeights(nthrow) = SD{m,a,No,stock}(1);
                                                                                                        TSDminvar_long_all(nthrow) = SD{m,a,No,stock}(2);
                                                                                                        TSDminvar_short_all(nthrow) = SD{m,a,No,stock}(3);
                                                                                                        TVaR005EqualWeights(nthrow) = VaR005{m,a,No,stock}(1);
                                                                                                        TVaR005minvar_long_all(nthrow) = VaR005{m,a,No,stock}(2);
                                                                                                        TVaR005minvar_short_all(nthrow) = VaR005{m,a,No,stock}(3);
                                                                                                        TCVaR005EqualWeights(nthrow) = CVaR005{m,a,No,stock}(1);
                                                                                                        TCVaR005minvar_long_all(nthrow) = CVaR005{m,a,No,stock}(2);
                                                                                                        TCVaR005minvar_short_all(nthrow) = CVaR005{m,a,No,stock}(3);

                                                                                                        TSREqualWeights(nthrow) = SR{m,a,No,stock}(1);

                                                                                                        TSRminvar_long_all(nthrow) = SR{m,a,No,stock}(2);
                                                                                                        TSRminvar_short_all(nthrow) = SR{m,a,No,stock}(3);
                                                                                                        TInvestedStocks{nthrow} = invested_in{m,a,No,stock};
                                                                                                        TmeanIndexEqual(nthrow) = indicesavg{m,a,No,stock}(1);
                                                                                                        TmeanIndexminvar_long_all(nthrow) =  indicesavg{m,a,No,stock}(2);

                                                                                                        TmeanIndexminvar_short_all(nthrow)  =  indicesavg{m,a,No,stock}(3);
                                                                                                        TIndex{nthrow} =  indices{m,a,No,stock};
% 
%                                                                                                         
%                                                                                                         
                                                                                                        
                                                                                                        
           
                                                                                                        numberExperimentsRan = (numberExperimentsRan+1);

                                                                                                                                                                                clc
                                                                                        percentage_complete = (numberExperimentsRan/numberExperimentsTotal) * 100
                                                                                        remaining_time = ((100-percentage_complete)/percentage_complete * toc)/60
                                                                                                    end
                                                                                                end
                                                                                            end
                                                                                        end
                                                                                        

                                                                                        
                end
            end
        end
    end
end
results_table = table([TMeasure{:}].', Truns,[Tadj{:}].', [Tcp{:}].', Ttv, Talp, TNo_stock, TShrinkCon, [TCorrOpts{:}].', TPREqualWeights, TPRminvar_long_all, TPRminvar_short_all, TEREqualWeights, TERminvar_long_all, TERminvar_short_all, TSDEqualWeights, TSDminvar_long_all, TSDminvar_short_all,TVaR005EqualWeights,TVaR005minvar_long_all, TVaR005minvar_short_all, TCVaR005EqualWeights, TCVaR005minvar_long_all, TCVaR005minvar_short_all, TSREqualWeights, TSRminvar_long_all, TSRminvar_short_all, TInvestedStocks, TmeanIndexEqual, TmeanIndexminvar_long_all,   TmeanIndexminvar_short_all, TIndex, 'VariableNames', VarNames  )

                                                                                                    
TotalTime = toc

%252 trading days in a year
marketSR = (mean(marketAvg)*252)/(std(marketAvg)*sqrt(252));
marketER = mean(marketAvg)*252;
marketSD = std(marketAvg)*sqrt(252);
marketVaR005 = quantile(marketAvg,0.05);
Index005 = marketAvg <= marketVaR005;
marketCVaR005 = mean(marketAvg(Index005));
marketAvgs = {marketAvg};
startDate = Date(1);
endDate = Date(end);
marketMetrics = table(marketSR, marketER, marketSD, marketVaR005, marketCVaR005,marketAvgs, startDate, endDate);

save_name =  split(string(datetime), " ");
save_name = join(save_name, "-");
save_name = strrep(save_name, ':', '-');
% save(join([save_path, save_name], '\'), 'results_table');
% save(join([save_path, "marketMetrics"], '\'), 'marketMetrics');
sortrows(results_table,    "PRminvar_long_all", 'descend')


