function weg = weight_portfolio(R,weight_opt)

% WEIGHT_PORTFOLIO determines the proportion
% of m assets to be invested.
% 'weight_opt'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'equal'                    : Equally weighted portfolio
% 'min_var_withshortsel'     : Minimum variance with shortselling
% 'min_var_withoutshortsel'  : Minimum variance without shortselling
% 'mean_var_withshortsel'    : Markowitz mean variance with shortselling
% 'mean_var_withoutshortsel' : Markowitz mean variance without shortselling

[~, m] = size(R);

switch weight_opt
    
    case 'equal'
       % The investent is made equally to each portfolio. 
        weg = repmat(1/m,m,1);
        
    case 'min_var_withshortsel'
        % Covariance matrix
        covmat = cov(R);
        zerosvec = zeros(m,1);
        AeqMat = ones(1,m);
        BeqMat = 1;
        % Upper and lower bound for the weights.
        LBMat = repmat(-0.25,m,1);
        UBMat = repmat( 0.25,m,1);
        initweights = repmat(1/m,m,1);
        options = optimset('display','off','largescale','off');
        weg = quadprog(covmat,zerosvec,[],[],AeqMat,BeqMat,LBMat,UBMat,initweights,options);
        %[weg, fval, exitflag] = quadprog(covmat,zerosvec,[],[],AeqMat,BeqMat,LBMat,UBMat,initweights,options)
    
	case 'min_var_withoutshortsel'
        
        % Covariance matrix
        covmat = cov(R);
        zerosvec = zeros(m,1);
        AeqMat = ones(1,m);
        BeqMat = 1;
        % Upper and lower bound for the weights.
        LBMat = zeros(m,1);
        UBMat = repmat(0.25,m,1);
        initweights = repmat(1/m,m,1);
        options = optimset('display','off','largescale','off');
        weg = quadprog(covmat,zerosvec,[],[],AeqMat,BeqMat,LBMat,UBMat,initweights,options);
        
    case 'mean_var_withshortsel'
		r_target = 0.00100;
        [weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);

		if exitflag ~= 1
			r_target = 0.00090;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00080;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00070;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;
		
		if exitflag ~= 1
			r_target = 0.00060;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00050;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;
		
		if exitflag ~= 1
			r_target = 0.00045;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;
			
		if exitflag ~= 1
			r_target = 0.00040;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00035;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00030;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00025;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00020;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00015;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00010;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;
		
		if exitflag ~= 1
			r_target = 0.00005;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00002;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;
		
		if exitflag ~= 1
			r_target = 0.00001;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			weg = zeros(m,1);
		end;
		
		
		
 case 'mean_var_withoutshortsel'
 	    r_target = 0.00100;
        [weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		
		if exitflag ~= 1
			r_target = 0.00090;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00080;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00070;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00060;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00050;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00045;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;
			
		if exitflag ~= 1
			r_target = 0.00040;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00035;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00030;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00025;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00020;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00015;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00010;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;
		
		if exitflag ~= 1
			r_target = 0.00005;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00002;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithShortSelling(m, R, r_target);
		end;

		if exitflag ~= 1
			r_target = 0.00001;
			[weg,fval,exitflag,output,lambdaqp] = MeanVarianceWithoutShortSelling(m, R, r_target);
		end;
		
		if exitflag ~= 1
			weg = zeros(m,1);
		end;
		
end


