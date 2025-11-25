function  [cum_sum,  invested_in, ER, SD, SR, VaR005,CVaR005,Weights_minvarlong,Weights_minvarshort, indices, indicesavg]= simulation(number_runs, thresholdIndex, correl_opt, shrinkage_constant, adj, data_options, analysis_options)


%Graph options:
%number_runs,thresholdValue,correl_opt,shrinkage_constant, adj

root_path = split(cd,'\');
methods_path = join([join(root_path(1:end-1), '\'), "\Methods"], "");
addpath(genpath(methods_path))

%Data options:
% changeTable, evaluation_period, stockStarts
changeTable = data_options{1};
evaluation_period = data_options{2};
stockStarts = data_options{3};
%

%Analysis options
% investedStocks, alpha, measure
no_stocks = analysis_options{1};
alphas = analysis_options{2};
measures = analysis_options{3};
%construct the graph


%Note that the results will be measure, alpha, nodes



%---------------------------Information about Data--------------------------------------
stockNames =  changeTable.Properties.VariableNames(2:end);
Date = changeTable.Date;
numberStocks = size(changeTable,2)-1;
numberDays = size(changeTable,1);

%----------------------------------------------------------------------------------------





%--------------------------------Initalize variables--------------------------------------
PREqualWeightsAll = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
PRminvarlongAll =  cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
PRminvarshortAll =  cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);

cum_sum = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
indicesavg = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
ER = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
SR = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
SD = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
VaR005 = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
Index005 = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
CVaR005 = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
M = [];
WindowA = {};

Weights_minvarlong = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
Weights_minvarshort = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
invested_in = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
indices = cell(size(measures, 2), size(alphas, 2),size(no_stocks,2),2);
%----------------------------------------------------------------------------------------


%This will run over the whole period, leaving no evaluation period.
% wrongdefaultWindowSize =floor((numberDays-evaluation_period)/number_runs)

%^Wrong! We've already removed the evaluation period.
defaultWindowSize = floor(numberDays/number_runs);



for nthWindow = 1:number_runs-1
    nthWindow;
    %Indices
    windowStart = (nthWindow-1)*(defaultWindowSize) + 1;
    windowEnd =  (nthWindow)*(defaultWindowSize) ;
    nextWindowStart = (nthWindow)*defaultWindowSize + 1;
    nextWindowEnd = (nthWindow+1)*defaultWindowSize;
    
   
    %Which stocks are available to be traded in this window
     existing_stocks = stockStarts <=  windowStart;
 
    
    %get the windows
    windowStocks = changeTable(windowStart:windowEnd, 2:end);
    nextWindowStocks = changeTable(nextWindowStart:nextWindowEnd, 2:end);

    %convert to array
    windowStocks = windowStocks{:,:};
    nextWindowStocks = nextWindowStocks{:,:};
    
    %corrcoef matrix:
    corrStocks = corrcoef(windowStocks);
    
   
    
    existing_stocks = existing_stocks.*~isnan(corrStocks(1,:));
    nonNan_indices = find(existing_stocks == 1);
    %Tranform corr matrix
     %remove NaNs 
    corrStocks = corrStocks(nonNan_indices,nonNan_indices) ;
    
    
    
    %possibly shrink the corr matrix
    if ~isnan(shrinkage_constant)
        corrStocks = Shrinkage(windowStocks(:,nonNan_indices), corrStocks, shrinkage_constant);
    else
        %
    end
    switch correl_opt
        case 'notransform'
        corrStocks = corrStocks;
        case 'subplus'
        corrStocks = subplus(corrStocks);
        case 'subplus negative elements'
        corrStocks = subplus(-corrStocks);
        case 'absolute value'
        corrStocks = abs(corrStocks);
    end

    
  %corrs sent to nan by shrinkage are set to 0
    corrStocks( find(isnan(corrStocks))) = 0;

    if isnan(thresholdIndex)
        A = Adjacency(corrStocks,adj);
    else
        A = Adjacency(corrStocks,adj,thresholdIndex);
    end
%To ensure the A's have the correct size we add back zeros to NaN slots
    newA = zeros(numberStocks,numberStocks);
    newA(nonNan_indices,nonNan_indices) = A;
    A = newA;
        WindowA{nthWindow} = A;
       
        M = buildM(M, A);

     %we loop over the analysis options
     for m = 1:size(measures,2)
         
         measure = measures{m};

         for a = 1:size(alphas, 2)

            alpha = alphas{a};
             
             
            centrality_vector = Centrality_Measure(WindowA,M,measure,alpha);
                          
            [~, order] = sort(centrality_vector,'descend'); 

            order = nonzeros((order).*existing_stocks(order).');

        for No = 1:size(no_stocks,2)


           %No = 1 denotes investing in the central stocks, 2 denotes investing in the periphery 
            investedStocks = no_stocks{No};
            invested_in{m,a,No,1} = [invested_in{m,a,No,1}; string(stockNames(order(1:investedStocks))), string(Date(nextWindowStart)), string(Date(nextWindowEnd)) ];
            invested_in{m,a,No,2} = [invested_in{m,a,No,2}; string(stockNames(order(end-investedStocks+1:end))), string(Date(nextWindowStart)), string(Date(nextWindowEnd)) ];

            Y_oldc = windowStocks(:, order(1:investedStocks));
            Yc = nextWindowStocks(:,order(1:investedStocks));
            w_equal = weight_portfolio(Y_oldc,'equal');
            w_minvarlong = weight_portfolio(Y_oldc,'min_var_withoutshortsel');
            w_minvarshort = weight_portfolio(Y_oldc,'min_var_withshortsel');
            PREqualWeights= Yc*w_equal;
            PRminvarlong = Yc*w_minvarlong;
            PRminvarshort = Yc*w_minvarshort;
            PREqualWeightsAll{m,a, No,1} = [PREqualWeightsAll{m,a, No,1}; PREqualWeights];
            PRminvarlongAll{m,a, No,1} =  [PRminvarlongAll{m,a, No,1}; PRminvarlong];
            PRminvarshortAll{m,a, No,1} = [PRminvarshortAll{m,a, No,1};   PRminvarshort];
            Weights_minvarlong{m,a, No,1}= [Weights_minvarlong{m,a, No,1}; w_minvarlong.'];
            Weights_minvarshort{m,a, No,1}= [Weights_minvarshort{m,a, No,1}; w_minvarshort.'];
            
            
            
             Y_oldp = windowStocks(:, order(end-investedStocks+1:end));
            Yp = nextWindowStocks(:,order(end-investedStocks+1:end));
            w_equal = weight_portfolio(Y_oldp,'equal');
            w_minvarlong = weight_portfolio(Y_oldp,'min_var_withoutshortsel');
            w_minvarshort = weight_portfolio(Y_oldp,'min_var_withshortsel');
            PREqualWeights= Yp*w_equal;
            PRminvarlong = Yp*w_minvarlong;
            PRminvarshort = Yp*w_minvarshort;
            PREqualWeightsAll{m,a, No,2} = [PREqualWeightsAll{m,a, No,2}; PREqualWeights];
            PRminvarlongAll{m,a, No,2} =  [PRminvarlongAll{m,a, No,2}; PRminvarlong];
            PRminvarshortAll{m,a, No,2} = [PRminvarshortAll{m,a, No,2};   PRminvarshort];
            Weights_minvarlong{m,a, No,2}= [Weights_minvarlong{m,a, No,2}; w_minvarlong.'];
            Weights_minvarshort{m,a, No,2}= [Weights_minvarshort{m,a, No,2}; w_minvarshort.'];
            
            
            
            
             if nthWindow == 1
   
             indexminvarlongallc = norm(Weights_minvarlong{m,a,No,1},1);
             indexminvarshortallc = norm(Weights_minvarshort{m,a,No,1},1);
             indexeqc = norm(w_equal,1); %should always be 1
             
             indexminvarlongallp = norm(Weights_minvarlong{m,a,No,2},1);
             indexminvarshortallp = norm(Weights_minvarshort{m,a,No,2},1);
             indexeqp = norm(w_equal,1); %should always be 1
             else
            oldPortfolioc = invested_in{m,a,No,1}(nthWindow-1,1:end-2);
            currentPortfolioc = invested_in{m,a,No,1}(nthWindow,1:end-2);
            %CENTRAL
           %Seeing how the actual stocks we invested in change:
           
           %Stocks that remain the same, but whose weight potentially
           %changes A n B
           [~, iioc]  = intersect(oldPortfolioc,currentPortfolioc, 'sorted');
           %Entirely new stocks A \ B
           [~, dioc]  = setdiff(oldPortfolioc,currentPortfolioc, 'sorted');
           
           %A n B
           [~, iinc]  = intersect(currentPortfolioc,oldPortfolioc, 'sorted');
           % B \ A 
           [~, dinc]  = setdiff(currentPortfolioc,oldPortfolioc, 'sorted');
            
            
            
            oldWeights_minvarlongc =Weights_minvarlong{m,a,No,1}(nthWindow-1, :);
            newWeights_minvarlongc = Weights_minvarlong{m,a,No,1}(nthWindow, :);
                        
            oldWeights_minvarshortc =Weights_minvarshort{m,a,No,1}(nthWindow-1, :);
            newWeights_minvarshortc = Weights_minvarshort{m,a,No,1}(nthWindow, :);
            
            %PERIPHERAL
            oldPortfoliop = invested_in{m,a,No,2}(nthWindow-1,1:end-2);
            currentPortfoliop = invested_in{m,a,No,2}(nthWindow,1:end-2);
            
           %Seeing how the actual stocks we invested in change:
           
           %Stocks that remain the same, but whose weight potentially
           %changes An B
           [~, iiop]  = intersect(oldPortfoliop,currentPortfoliop, 'sorted');
           %Entirely new stocks A \ B
           [~, diop]  = setdiff(oldPortfoliop,currentPortfoliop, 'sorted');
           
           %A n B
           [~, iinp]  = intersect(currentPortfoliop,oldPortfoliop, 'sorted');
           % B \ A 
           [~, dinp]  = setdiff(currentPortfoliop,oldPortfoliop, 'sorted');
            
            
            
            oldWeights_minvarlongp =Weights_minvarlong{m,a,No,2}(nthWindow-1, :);
            newWeights_minvarlongp = Weights_minvarlong{m,a,No,2}(nthWindow, :);
                        
            oldWeights_minvarshortp =Weights_minvarshort{m,a,No,2}(nthWindow-1, :);
            newWeights_minvarshortp = Weights_minvarshort{m,a,No,2}(nthWindow, :);
            



%CENTRAL
   differenceSameStocksminvarlongc = abs(Weights_minvarlong{m,a,No,1}(nthWindow-1,iioc) - Weights_minvarlong{m,a,No,1}(nthWindow,iinc));
   differenceDifferentStocksminvarlongc = abs(Weights_minvarlong{m,a,No,1}(nthWindow-1,dioc) + Weights_minvarlong{m,a,No,1}(nthWindow,dinc));
   
   
   indexminvarlongallc = sum(differenceSameStocksminvarlongc) + sum(differenceDifferentStocksminvarlongc);
   indexminvarlongallc = indexminvarlongallc/2;
   
   
   differenceSameStocksminvarshortc = abs(Weights_minvarshort{m,a,No,1}(nthWindow-1,iioc) - Weights_minvarshort{m,a,No,1}(nthWindow,iinc));
   differenceDifferentStocksminvarshortc = abs(Weights_minvarshort{m,a,No,1}(nthWindow-1,dioc) + Weights_minvarshort{m,a,No,1}(nthWindow,dinc));
   indexminvarshortallc = sum(differenceSameStocksminvarshortc) + sum(differenceDifferentStocksminvarshortc);
   indexminvarshortallc = indexminvarlongallc/2;
   
    indexeqc = size(oldPortfolioc(dioc),2) + size(currentPortfolioc(dinc),2);
    indexeqc = indexeqc/(2*size(oldPortfolioc,2));
    
    
    %PERIPHERAL
       differenceSameStocksminvarlongp= abs(Weights_minvarlong{m,a,No,2}(nthWindow-1,iiop) - Weights_minvarlong{m,a,No,2}(nthWindow,iinp));
   differenceDifferentStocksminvarlongp = abs(Weights_minvarlong{m,a,No,2}(nthWindow-1,diop) + Weights_minvarlong{m,a,No,2}(nthWindow,dinp));
   
   
   indexminvarlongallp = sum(differenceSameStocksminvarlongp) + sum(differenceDifferentStocksminvarlongp);
   indexminvarlongallp = indexminvarlongallp/2;
   
   
   differenceSameStocksminvarshortp = abs(Weights_minvarshort{m,a,No,2}(nthWindow-1,iiop) - Weights_minvarshort{m,a,No,2}(nthWindow,iinp));
   differenceDifferentStocksminvarshortp = abs(Weights_minvarshort{m,a,No,2}(nthWindow-1,diop) + Weights_minvarshort{m,a,No,2}(nthWindow,dinp));
   indexminvarshortallp = sum(differenceSameStocksminvarshortp) + sum(differenceDifferentStocksminvarshortp);
   indexminvarshortallp = indexminvarlongallp/2;
   
    indexeqp = size(oldPortfoliop(diop),2) + size(currentPortfoliop(dinp),2);
    indexeqp = indexeqp/(2*size(oldPortfoliop,2));
    
   
    end
   indices{m,a,No,1} = [indices{m,a,No,1};indexeqc, indexminvarlongallc, indexminvarshortallc];
   indices{m,a,No,2} = [indices{m,a,No,2};indexeqp, indexminvarlongallp, indexminvarshortallp];
        end



                
         end
     end
     
     
     
end

%Now we have all the data, so we just compute cumulative returns, looping
%all combinations of parameters we had

    for m = 1:size(measures,2)
        for a  = 1:size(alphas, 2)
            for No = 1:size(no_stocks,2)
                %CENTRAL
                cum_sum{m,a,No,1} = [ sum(PREqualWeightsAll{m,a,No,1}), sum(PRminvarlongAll{m,a,No,1}),  sum(PRminvarshortAll{m,a,No,1})];
                indicesavg{m,a,No,1} =  mean(indices{m,a,No,1});
                ER{m,a,No,1} = mean([ PREqualWeightsAll{m,a,No,1}, PRminvarlongAll{m,a,No,1},  PRminvarshortAll{m,a,No,1}])*252;
                SD{m,a,No,1} = std([ PREqualWeightsAll{m,a,No,1}, PRminvarlongAll{m,a,No,1},  PRminvarshortAll{m,a,No,1}])*sqrt(252);
                SR{m,a,No,1} = ER{m,a,No,1}./SD{m,a,No,1};
                VaR005{m,a,No,1} = quantile( [PREqualWeightsAll{m,a,No,1}, PRminvarlongAll{m,a,No,1},  PRminvarshortAll{m,a,No,1}], 0.05);
                Index005{m,a,No,1}  = [PREqualWeightsAll{m,a,No,1} <= VaR005{m,a,No,1}(1),  PRminvarlongAll{m,a,No,1} <= VaR005{m,a,No,1}(2),   PRminvarshortAll{m,a,No,1} <= VaR005{m,a,No,1}(3)];
                CVaR005{m,a,No,1} = mean([PREqualWeightsAll{m,a,No,1}(Index005{m,a,No,1}(:,1)), PRminvarlongAll{m,a,No,1}(Index005{m,a,No,1}(:,2)), PRminvarshortAll{m,a,No,1}(Index005{m,a,No,1}(:,3))]);
 %PERIPHERAL
                cum_sum{m,a,No,2} = [ sum(PREqualWeightsAll{m,a,No,2}), sum(PRminvarlongAll{m,a,No,2}),  sum(PRminvarshortAll{m,a,No,2})];
                indicesavg{m,a,No,2} =  mean(indices{m,a,No,2});
                ER{m,a,No,2} = mean([ PREqualWeightsAll{m,a,No,2}, PRminvarlongAll{m,a,No,2},  PRminvarshortAll{m,a,No,2}])*252;
                SD{m,a,No,2} = std([ PREqualWeightsAll{m,a,No,2}, PRminvarlongAll{m,a,No,2},  PRminvarshortAll{m,a,No,2}])*sqrt(252);
                SR{m,a,No,2} = ER{m,a,No,2}./SD{m,a,No,2};
                SR{m,a,No,2} = ER{m,a,No,2}./SD{m,a,No,2};
                VaR005{m,a,No,2} = quantile( [PREqualWeightsAll{m,a,No,2}, PRminvarlongAll{m,a,No,2},  PRminvarshortAll{m,a,No,2}], 0.05);
                Index005{m,a,No,2}  = [PREqualWeightsAll{m,a,No,2} <= VaR005{m,a,No,2}(1),  PRminvarlongAll{m,a,No,2} <= VaR005{m,a,No,2}(2),   PRminvarshortAll{m,a,No,2} <= VaR005{m,a,No,2}(3)];
                CVaR005{m,a,No,2} = mean([PREqualWeightsAll{m,a,No,2}(Index005{m,a,No,2}(:,1)), PRminvarlongAll{m,a,No,2}(Index005{m,a,No,2}(:,2)), PRminvarshortAll{m,a,No,2}(Index005{m,a,No,2}(:,3))]);
            end
        end
    end
        


end