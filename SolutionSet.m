classdef SolutionSet < handle
    %R. Dan, Aug 2017
    %Processed density measurements for cold cryoprotectants.
    %   Detailed explanation goes here
    
    properties
        % initial values
        name               % string
        mol_weight         % double
        e_molecule         % double
        eDensity           % vector of double
        concww             % vector of double
        concwv             % vector of double
        floats             % vector of double
        sinks              % vector of double
        density            % vector of double
        upuncertain        % vector of double
        downuncertain      % vector of double
        index              % vector of int, generally a temp place holder
        rho                % struct 'cryo' struct: 'cold': double
        %                       'rt'  : double
        %        'h2o'  struct: 'cold': double
        %                       'rt'  : double
        temp               % string
        degtemp            % double
        label              % string
        
        % derived values
        VE                 % vector of double
        dVE                % struct 'up', 'down'   uncertainty of VE
        modeledVE          % vector of double, length: 101
        modeledDensity     % struct 'fitted'       vector of double
        %        'fit'          cell array: NonlinearModel
        %                                   vector of double
        %         'interpolated' vector of double, length: 101 
        modeledElecDensity % struct 'fitted'       vector of double
        %        'fit'          cell array: NonlinearModel
        %                                   vector of double
        %         'interpolated' vector of double, length: 101
        eDensityDiff % struct 'fitted' vector of double
        %         'interpolated' vector of double, length: 101
        modeledSpecificVol % vector of double, length: 101
        concInterpolated   % vector of double, length: 101
        SpecificVol        % vector of double
        predictedDensity   % vector of double
        Vapph2o            % struct 'measured':     vector of double
        %        'interpolated': vector of double, length 101
        %        'forced':       vector of double
        % ID double
        Vapp               % struct 'measured':     vector of double
        %        'interpolated': vector of double, length 101
        %        'forced':       vector of double
        % ID double
    end
    
    methods
        %% Data Preparation
        function Solu = SolutionSet(name, truncate)
            % Constructor:  construct a SolutionSet object.  Name is a string
            % that gives the type of cryoprotectant
            
            % classify data
            if nargin==0
                error('enter a cryoprotectant name or temperature')
            else
                if nargin~=2
                    truncate = 1; % truncate
                end
                Solu.name = name;
            end
            
            % Load cryoprotectant data
            if strcmp(Solu.name, 'PEG 200')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'PEG 200');
                rtrho = 1.1194;
                Solu.mol_weight = 200.00;
                Solu.e_molecule = 106.00;
            elseif strcmp(Solu.name, 'PEG 200 RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'PEG 200 RT');
                rtrho = 1.1094;
                Solu.mol_weight = 200.00;
                Solu.e_molecule = 106.00;
            elseif strcmp(Solu.name, 'PEG 200 HT') %HT is Hot
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'PEG 200 hot');
                rtrho = 1.096;
                Solu.mol_weight = 200.00;
                Solu.e_molecule = 106.00;
            elseif strcmp(Solu.name, 'Ethanol')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Ethanol');
                rtrho = 0.78508;
                Solu.mol_weight = 46.07;
                Solu.e_molecule = 26.00;
            elseif strcmp(Solu.name, 'Ethanol RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Ethanol RT');
                rtrho = 0.78508;
                Solu.mol_weight = 46.07;
                Solu.e_molecule = 26.00;
            elseif strcmp(Solu.name, '2-Propanol')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', '2-Propanol');
                rtrho = 0.78498;
                Solu.mol_weight = 60.10;
                Solu.e_molecule = 34.00;
            elseif strcmp(Solu.name, '2-Propanol RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', '2-Propanol RT');
                rtrho = 0.78498;
                Solu.mol_weight = 60.10;
                Solu.e_molecule = 34.00;
            elseif strcmp(Solu.name, 'Ethylene Glycol')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Ethylene Glycol');
                rtrho = 1.1099;
                Solu.mol_weight = 62.07;
                Solu.e_molecule = 34.00;
            elseif strcmp(Solu.name, 'Glycerol')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Glycerol');
                rtrho = 1.2580;
                Solu.mol_weight = 92.09;
                Solu.e_molecule = 50;
            elseif strcmp(Solu.name, 'Ethylene Glycol RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Ethylene Glycol RT');
                rtrho = 1.1099;
                Solu.mol_weight = 62.07;
                Solu.e_molecule = 34.00;
            elseif strcmp(Solu.name, 'Glycerol RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Glycerol RT');
                rtrho = 1.2580;
                Solu.mol_weight = 92.09;
                Solu.e_molecule = 50.00;
            elseif strcmp(Solu.name, 'MPD RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'MPD RT');
                rtrho = 0.9181;
                Solu.mol_weight = 118.18;
                Solu.e_molecule = 66.00;
            elseif strcmp(Solu.name, 'MPD')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'MPD');
                rtrho = 0.9181;
                Solu.mol_weight = 118.18;
                Solu.e_molecule = 66.00;
            elseif strcmp(Solu.name, 'PPG')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'PPG');
                rtrho = 1.00311;
                Solu.mol_weight = 425.00;
                Solu.e_molecule = 234.00;
            elseif strcmp(Solu.name, 'PPG RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'PPG RT');
                rtrho = 1.00311;
                Solu.mol_weight = 425.00;
                Solu.e_molecule = 234.00;
            elseif strcmp(Solu.name, 'Methanol RT')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Methanol RT');
                rtrho = 0.7869;
                Solu.mol_weight = 32.04;
                Solu.e_molecule = 18.00;
            elseif strcmp(Solu.name, 'Methanol RT 2')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Methanol RT 2');
                rtrho = 0.7876;
                Solu.mol_weight = 32.04;
                Solu.e_molecule = 18.00;
            elseif strcmp(Solu.name, 'Methanol')
                [~, ~, raw] = xlsread('/Users/Ritwik Dan/My Documents/DensityPaperDraftFile/DensityPaperDraft/MATLAB/labdata.xlsx', 'Methanol');
                rtrho = 0.7869;
                Solu.mol_weight = 32.04;
                Solu.e_molecule = 18.00;
            else
                error('enter a valid cryoprotectant name')
            end
            raw = raw(2:end,:);
            
            % Replace non-numeric cells with -9999
            R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
            raw(R) = {-9999}; % Replace non-numeric cells
            
            % Create output variable
            in = cell2mat(raw);
            
            % Clear temporary variables
            clearvars raw R;
            
            %truncate incomplete data
            if truncate == 1
                tmp   = [];
                for i = 1:length(in(:,1))
                    k = 1;
                    if in(i,3) == -9999
                        k = 0;
                    end
                    if k
                        tmp = [tmp ; in(i,:)];
                    end
                end
                in = tmp;
            end
            Solu.degtemp = in(1,8);
            if Solu.degtemp > 80
                Solu.temp = 'rt';
            else
                Solu.temp = 'cold';
            end
            
            if strcmp(Solu.temp, 'rt') % if room temperature
                Solu.concwv = in(:,1);
                Solu.concww = in(:,7);
                Solu.density = in(:,3)/1000;
                tmp = Solu.name(1:length(Solu.name)-3);
                if length(tmp)==6 && strcmp(tmp(1:6), 'PEG 200')
                    tmp = 'PEG';
                end
                %Solu.label = strcat(tmp,{', '},num2str(floor(Solu.degtemp)), 'K');
                Solu.label = [tmp,', ',num2str(floor(Solu.degtemp)), 'K'];
            else % else 77K
                Solu.concwv = in(:,1);
                Solu.concww = in(:,7);
                Solu.floats = in(:,2)./1000;
                Solu.density = in(:,3)./1000;
                Solu.sinks = in(:,4)./1000;
                Solu.label = strcat(Solu.name,', 77K');
                Solu.upuncertain   = Solu.floats - Solu.density;
                Solu.downuncertain = Solu.density - Solu.sinks;
                Solu.upuncertain(1) = 0;
                Solu.downuncertain(1) = 0;
            end
            
            
            Solu.rho =struct('cryo',[],'h2o',[]);
            Solu.rho.h2o = struct('cold', 0.94, 'rt',0.997044);
            if strcmp(Solu.temp,'rt');
                Solu.rho.h2o.rt = Solu.density(1);
            end
            Solu.rho.cryo = struct('cold', Solu.density(end), 'rt', rtrho);
            
            % Calculate Excess Volume
            x = 0:100; x = x';
            Solu.concInterpolated = x;
            Solu.VE = makeExcessVolume(Solu, Solu.temp);
            Solu.modeledDensity = struct('fitted',[], 'fit',[], 'interpolated', []);
            %Solu.modeledVEonly = struct('fitted',[], 'fit',[], 'interpolated', []);
            
             % Calculate uncertainty in Excess Volume
             Solu.dVE = struct('up','down');
             if strcmp(Solu.temp,'cold')
                 Solu.dVE.up = makeExcessVolume(Solu, Solu.temp, 'up') - Solu.VE;  Solu.dVE.up(1) = 0;
                 Solu.dVE.down = -1*(makeExcessVolume(Solu, Solu.temp, 'down') - Solu.VE);  Solu.dVE.down(1) = 0;
             end
            
            % Interpolate Excess Volume
                %%%%%%% LEAVE DESIRED MODEL UNCOMMENTED %%%%%%%
                %[Solu.modeledVEonly.fitted, Solu.modeledVEonly.fit] = Solu.gaussianReg();
                %Solu.modeledVEonly.interpolated = predict(Solu.modeledVEonly.fit{1}, x);
                %Solu.modeledVE.interpolated = predict(Solu.modeledVE.fit{1}, x) + polyval(Solu.modeledVE.fit{2}, x);
                %[Solu.modeledVE.fitted Solu.modeledVE.fit] = Solu.quarticReg();  
                %Solu.modeledVE.interpolated = polyval(Solu.modeledVE.fit, x-50);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            h2o_M = 18.0; % molecular weight of water
            h2o_elecs = 10.0; % electrons/mol pure water
            avog_N = 6.022e23; % avo gadro's number
            cc_A3 = 1.00e-24; % cm^3/A^3

            % Calculate Electron Density
            Solu.eDensity = Solu.makeElecDensity();
                
            % Calculate Specific Volume
            Solu.SpecificVol = 1./Solu.density;
            
            % Model density and interpolated electron density
            
            Solu.SpecificVol = Solu.idealVolume(Solu.temp, Solu.concww)+Solu.VE;
            %Solu.modeledSpecificVol = Solu.idealVolume(Solu.temp, Solu.concInterpolated)+Solu.modeledVE.interpolated;
            [Solu.modeledDensity.fitted, Solu.modeledDensity.fit] = Solu.quarticRegDens();  
            Solu.modeledDensity.interpolated = polyval(Solu.modeledDensity.fit, Solu.concInterpolated/100);
            Solu.modeledVE = 1./(Solu.modeledDensity.interpolated) - Solu.idealVolume(Solu.temp, Solu.concInterpolated);
            
            
            Solu.modeledElecDensity = Solu.modeledDensity.interpolated.*(Solu.concInterpolated/100.0*(Solu.e_molecule/Solu.mol_weight)+(100.0-Solu.concInterpolated)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
            %[Solu.modeledElecDensity.fitted, Solu.modeledElecDensity.fit] = Solu.quarticRegElecDens();
            %Solu.modeledElecDensity.interpolated = polyval(Solu.modeledElecDensity.fit, Solu.concInterpolated/100);
            
            Solu.eDensityDiff = struct('fitted', [], 'interpolated', []);
            Solu.eDensityDiff = Solu.makeElecDensityDiff();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%CHANGE MADE HERE FOR POLYINTERP %%%%
            % make predictions for lab
            %Solu.predictedDensity = 1./(Solu.idealVolume(Solu.temp, Solu.concww)...
            %    +predict(Solu.modeledVE.fit{1}, Solu.concww) + polyval(Solu.modeledVE.fit{2}, Solu.concww));
            
            % calculate apparent volume of cryoprotectant and water
            Solu.Vapp = struct('measured',[],'interpolated',[], 'forced', [], 'ID', []);
            Solu.Vapph2o = struct('measured',[],'interpolated',[], 'forced', [], 'ID', []);
            if strcmp(Solu.temp, 'rt')
                
                Solu.Vapp.measured = 1./Solu.density +(100./Solu.concww - 1).*(1./Solu.density - 1/Solu.rho.h2o.rt);
                Solu.Vapph2o.measured = 1./Solu.density +Solu.concww/100.*(1./(1-Solu.concww./100)).*(1./Solu.density - 1/Solu.rho.cryo.rt);
                
                Solu.Vapp.interpolated = 1./Solu.modeledDensity.interpolated +(100./Solu.concInterpolated - 1).*(1./Solu.modeledDensity.interpolated - 1/Solu.rho.h2o.rt);
                Solu.Vapph2o.interpolated = 1./Solu.modeledDensity.interpolated +Solu.concInterpolated/100.*(1./(1-Solu.concInterpolated./100)).*(1./Solu.modeledDensity.interpolated - 1/Solu.rho.cryo.rt);
                
            else
                
                Solu.Vapp.measured = 1./Solu.density +(100./Solu.concww - 1).*(1./Solu.density - 1/Solu.rho.h2o.cold);
                Solu.Vapph2o.measured = 1./Solu.density +Solu.concww/100.*(1./(1-Solu.concww./100)).*(1./Solu.density - 1/Solu.rho.cryo.cold);
                
                Solu.Vapp.interpolated = 1./Solu.modeledDensity.interpolated +(100./Solu.concInterpolated - 1).*(1./Solu.modeledDensity.interpolated - 1/Solu.rho.h2o.cold);
                Solu.Vapph2o.interpolated = 1./Solu.modeledDensity.interpolated +Solu.concInterpolated/100.*(1./(1-Solu.concInterpolated./100)).*(1./Solu.modeledDensity.interpolated - 1/Solu.rho.cryo.cold);
                
            end
            
            % force a fit on the high concentration region
            y = Solu.Vapp.interpolated;
            M = min(y(2:end));
            T=0;
            for i = 2:length(y)
                if T == 0
                    if y(i) == M && T==0
                        T = 1;
                        Solu.Vapp.ID=i;
                    end
                end
            end
            Solu.Vapp.forced = y(Solu.Vapp.ID).*ones(1,Solu.Vapp.ID);
            
            
            % force a fit on the low concentration region
            y = Solu.Vapph2o.interpolated;
            M = min(y(2:end));
            T=0;
            for i = 2:length(y)
                if T == 0
                    if y(i) == M && T==0
                        T = 1;
                        Solu.Vapph2o.ID=i;
                    end
                end
            end
            Solu.Vapph2o.forced = y(Solu.Vapph2o.ID).*ones(1,length(y) - Solu.Vapph2o.ID);
        end
        
        %% Data Processing
        
        function VE  = makeExcessVolume(self, temp, dVE)
            % calculates the excess specific volume of solution from the measured
            % solution density, and the density of either pure substance at the
            % temperature, T.  Assume concentrations are in %w/w.
            if nargin == 2
                dVE = '';
                % calculate actual volume
                volAct = 1./self.density;
            elseif strcmp(dVE, 'up')
                volAct = 1./self.floats;
            elseif strcmp(dVE, 'down')
                volAct = 1./self.sinks;
            end
            
            % calculate ideal volume
            volId = self.idealVolume(temp);
            % calculate excess specific volume
            VE = volAct - volId;
        end
        
        function idealVol = idealVolume(self, temp, x)
            % calculates the ideal specific volume of a mixture (under the simplistic
            % assumption)
            if nargin ==3
                conc = x;
            else
                conc = self.concww;
            end
            if strcmp(temp, 'cold')
                cryorho = self.rho.cryo.cold;
                h2orho = self.rho.h2o.cold;
            elseif (strcmp(temp, 'rt'))
                cryorho = self.rho.cryo.rt;
                h2orho = self.rho.h2o.rt;
            else
                error('enter valid temperature')
            end
            idealVol = conc/100/cryorho + (1-conc/100)/h2orho;
        end
        
        function [Y_hat FitArray] = gaussianReg(self)
            % This function performs a nonlinear regression on the functional form of a gaussian curve
            % double reg is 0 for a normal regressian, and is 1 if the variance will
            %have different values on either side of the gaussian.  ve=1 if calculating
            %the fit for excess specific volume. weightwater = 1 if the density of
            %amorphous ice should be weighted (not set up to work
            %with doublereg). FitArray{1} gives gaussian fit.  FitArray{2} gives its quartic ECM.
            
            %choose sign convention for VE and choose to weight 0%
            %concentration.
            X = self.concww; Y = self.VE;
            if strcmp(self.temp,'cold')
                DY = 1./abs(self.dVE.up + self.dVE.down);
            else
                DY = ones(length(X),1);
            end
                   
            weightwater = 1; % weightwater = 1 if conc = 0 and conc = 100 should be weighted
            ve = 1; %%sign convention = +1 if concave gaussian expected
            if strcmp(self.temp,'rt')
                if (strcmp(self.name,'Ethanol') || strcmp(self.name,'2-Propanol'))
                    ve = -1;
                end
            end
                       
            
            % weight by 1/uncertainty = 1/DY
            mw = max(DY(2:end));
            M = max(X);
            L = length(X);
            if weightwater==1
                Weight = ones(L,1);
                for i = 1:L            
                    if X(i)==0
                        Weight(i) = 10*mw;
                    elseif X(i) ==M
                        Weight(i) = 10*mw;
                    else
                        Weight(i) = DY(i);
                    end
                end
            else
                Weight = ones(L,1);
            end
            
            %  Approximate an initial fit           
            if ve==1
                [B0 i] = min(Y); % Amplitude
            else
                [B0 i] = max(Y); % Amplitude
            end
            B1 = X(i);  % Horizontal shift
            
                        

      
            
            B3 = 0; % Vertical shift
            
            % if cold monoalcohol model
            if 0 == 1% strcmp(self.temp,'cold') && (strcmp(self.name,'Ethanol') || strcmp(self.name,'2-Propanol'))
                %2 sided gaussian model, piecewise given by heaviside(xmax)

                % sort data
                [X, I] = sort(X);
                Ytemp = Y;
                DYtemp = DY;
                for j = 1:length(I)
                    Y(j) = Ytemp(I(j));
                    DY(j) = DYtemp(I(j));
                end
                i = I(i);
            
            
                %approximate LHS width
                j = 1; k = 0;
                    while (j < i && k ==0)
                        if Y(j) > B0/2
                            B2LHS = abs(X(j) - B1);  % Width of curve. ~ Standard deviation or width-of-half-max.
                            k = 1;
                        end
                        j = j +1;
                    end
                
                %approximate RHS width
                j = length(Y); k = 0;
                    while (j > i && k == 0)
                        if Y(j) > B0/2
                            B2RHS = abs(X(j) - B1);  % Width of curve. ~ Standard deviation or width-of-half-max.
                            k = 1;
                        end
                        j = j -1;
                    end
                    
                myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*heaviside(B1-x).*exp(-((x-B1)/100).^2/(2*b(2)^2)) + ...
                    b(1)*heaviside(x-B1).*exp(-((x-B1)/100).^2/(2*b(3)^2))+heaviside(x-B1).*((x-B1)/100).*(b(4)+b(5)*x/100)+b(6)),...
                    [B0,  B2LHS/100, B2RHS/100, 0,0, B3],'Weight', Weight);    
                
                %next to last model
                %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*(heaviside(B1-x).*exp(-(x-B1).^2/(2*b(3)^2)) + ...
                %    heaviside(x-B1).*(exp(-(x-B1).^2/(2*b(4)^2))+(x-B1).*(b(5)+b(6)*x))+b(2))),...
                %    [B0, B3, B2LHS, B2RHS, 0,0],'Weight', Weight);    
                    
                %hermite polynomial decomposition
                %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*(heaviside(b(2)-x).*exp(-(x-b(2)).^2/(2*b(3)^2)) + ...
                %    heaviside(x-b(2)).*(1+b(5).*(1-(x-b(2)).^2/(2*b(4)^2))).*exp(-(x-b(2)).^2/(2*b(4)^2)))), ... 
                %    [B0, B1, B2LHS, B2RHS, 0],'Weight', Weight);
                    
                %quartic poly linearly on RHS . "mean" held fixed
                %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*(heaviside(B1-x).*exp(-(x-B1).^2/(2*b(3)^2)) + ...
                %    heaviside(x-B1).*(exp(-(x-B1).^2/(2*b(4)^2))+(x-B1).*(b(5)))+b(2))),...
                %    [B0, B3, B2LHS, B2RHS, 0],'Weight', Weight);
                
                %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*(heaviside(B1-x).*exp(-(x-B1).^2/(2*b(3)^2)) + ...
                %    heaviside(x-B1).*(exp(-(x-B1).^2/(2*b(4)^2))+(x-B1).*(b(5)+b(6)*x+b(7)*x.^2+b(8)*x.^3))+b(2))),...
                %    [B0, B3, B2LHS, B2RHS, 0, 0,0,0],'Weight', Weight);
                
                %2sided gaussian
                %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*(heaviside(b(2)-x).*exp(-(x-b(2)).^2/(2*b(4)^2)) +  heaviside(x-b(2)).*exp(-(x-b(2)).^2/(2*b(5)^2)))+b(3)), [B0, B1, B3, B2LHS, B2RHS],'Weight', Weight);
                
                %Bp = polyfit(X(i:end),Y(i:end),4);
                %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*(heaviside(b(2)-x).*exp(-(x-b(2)).^2/(2*b(3)^2)) + ...
                %    heaviside(x-b(2)).*(b(5)*x.^4+b(6)*x.^3+b(7)*x.^2+b(8)*x+b(9)))), [B0, B1, B2LHS, B2RHS, Bp],'Weight', Weight);
                
                % first term written is LHS Gaussian, second term is RHS
%                 Y_hat = myFit.Fitted;
%                 pFit = polyfit(X, Y - myFit.Fitted, 4);
%                 
%                 %rerun nonlinear fit with pFit as initial values for the ECM
%                 myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*(heaviside(b(2)-x).*exp(-(x-b(2)).^2/(2*b(3)^2)) ...
%                     +  heaviside(x-b(2)).*exp(-(x-b(2)).^2/(2*b(4)^2)))) + polyval([b(5),b(6),b(7),b(8),b(9)],x) ...
%                     , [B0, B1, B2LHS, B2RHS, pFit],'Weight', Weight);
%                 
            else
                
%                     [A I] = sort(Y);
%                     j = 1; k = 0;
%                     while j < length(Y) && k ==0
%                         if A(j) > B0/2
%                             B2 = abs(X(I(j)) - B1);  % Width of curve. ~ Standard deviation or width-of-half-max.
%                             k = 1;
%                         end
%                         j = j +1;
%                     end
%                     
%                     %estimate amplitude of first hermite excitation
%                     B3 = Y(1) - B0*exp(-(X(1)-B1).^2/(2*B2)^2);
%                     B4 = max(abs(Y - (B0*exp(-(X-B1).^2/(2*B2)^2)+B3)));
                    
                    % quartic model
                    %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1) + b(2)*x + b(3)*x.^2+ b(4)*x.^3+ b(5)*x.^4), [0,0,0,0,0],'Weight', Weight);
                    
                    % quartic model with B.C.'s
                    myFit = fitnlm(X,Y, @(b,x)(b(1)*x/100.*(1-x/100).*(1 + b(2)*x/100 + b(3)*(x/100).^2)), [B0,0,0],'Weight', Weight);
                    
                    
                    % gaussian model  
                    %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*exp(-(x-b(2)).^2/(2*b(3)^2))...
                    %    +b(4) +b(5)*(x-b(2))+ b(6)*(x-b(2)).^2), [B0,B1, B2, B3,0,0],'Weight', Weight);
                    
                    %oldmodel
                    %myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*exp(-(x-b(2)).^2/(2*b(3)^2))...
                    %    +b(4) +b(5)*(x-b(2))+ b(6)*(x-b(2)).^2+ b(7)*(x-b(2)).^3 + b(8)*(x-b(2)).^4), [B0,B1, B2, B3,0,0,0,0],'Weight', Weight);
                    
                    
                    %myFit = NonLinearModel.fit(X,Y, @(b,x)((b(1)-b(5).*(x-b(2))/(sqrt(2)*b(3))).*exp(-(x-b(2)).^2/(2*b(3)^2))+b(4)), [B0, B1, B2, B3, B4],'Weight', Weight);

                    
            end


                    Y_hat = myFit.Fitted;
                    FitArray = {myFit, polyfit(X, Y - myFit.Fitted, 4)};

        end
        
        function [Y_hat fit] = quarticReg(self)
            % quartic regression
            X = self.concww-50; Y = self.VE;
            degree = 4;
            
            %increase the degree of fit for monoalcohols
            if strcmp(self.temp,'cold')
                if (strcmp(self.name,'Ethanol') || strcmp(self.name,'2-Propanol'))
                    degree = 4;
                end
            end
            
            
            fit = polyfit(X, Y, degree);
            Y_hat = polyval(fit, X);
        end
        
        function [Y_hat, fit] = quarticRegDens(self)
            % quartic regression
            X = self.concww/100; Y = self.density;
            degree = 4;
            fit = polyfit(X, Y, degree);
            Y_hat = polyval(fit, X);
        end
        
        function [Y_hat, fit] = quarticRegElecDens(self)
            % quartic regression
            X = self.concww/100; Y = self.eDensity;
            degree = 5;
            fit = polyfit(X, Y, degree);
            Y_hat = polyval(fit, X);
        end
        
        function [resid, rsq, rsq_adj, stdslope] = residCheck(self, p, y, y_hat, x)
            % this function analyzes at the residuals of a given regression.  the
            % inputs and resid are vectors and rsq and stdslope are double.
            
            resid = y-y_hat;
            ssr = sum(resid.^2);
            sst = (length(y) -1)*var(y);
            rsq = 1-ssr/sst;
            
            n = length(x);
            dfe = n-p-1;
            dft = n-1;
            rsq_adj = 1.0-(ssr/dfe)/(sst/dft);
            
            ssx = sum((x - mean(x)).^2);
            stdslope = sqrt(var(resid)/ssx);
            
        end
        
        function [volfrac, expfrac] = calcContraction(coldObj, warmObj, dvolfrac)
            %if dvolfrac is 'up' or 'down', volfrac will be the error
            %estimates for the data in percent change of volume
            %%%%%%%%%%%%%%%%%%%%%%%%CHANGE MADE HERE FOR POLYINTERP %%%%
            %     if nargin == 2
            %         warm = idealVolume(coldObj,'rt')+predict(warmObj.modeledVE.fit{1}, coldObj.concww) + polyval(warmObj.modeledVE.fit{2}, coldObj.concww);
            %         volfrac = 100*(coldObj.density.^-1 - warm)./warm;
            %         expfrac  = 100*(coldObj.modeledDensity.^-1 - warmObj.modeledDensity.^-1)./warmObj.modeledDensity.^-1;
            %     elseif strcmp(dvolfrac, 'up')
            %         warm = idealVolume(coldObj,'rt')+predict(warmObj.modeledVE.fit{1}, coldObj.concww) + polyval(warmObj.modeledVE.fit{2}, coldObj.concww);
            %         volfrac = 100*(coldObj.floats.^-1 - warm)./warm;
            %         % needs continuous dDensity model: expfrac  = 100*(coldObj.modeledDensity.^-1 - warmObj.modeledDensity.^-1)./warmObj.modeledDensity.^-1;
            %         expfrac = [];
            %     elseif strcmp(dvolfrac, 'down')
            %         warm = idealVolume(coldObj,'rt')+predict(warmObj.modeledVE.fit{1}, coldObj.concww) + polyval(warmObj.modeledVE.fit{2}, coldObj.concww);
            %         volfrac = 100*(coldObj.sinks.^-1 - warm)./warm;
            %         % needs continuous dDensity model: expfrac  = 100*(coldObj.modeledDensity.^-1 - warmObj.modeledDensity.^-1)./warmObj.modeledDensity.^-1;
            %         expfrac = [];
            %     end
            
            if nargin == 2
                %warm = idealVolume(coldObj,'rt') + polyval(warmObj.modeledVE.fit, coldObj.concww);
                %warm = idealVolume(coldObj,'rt') + predict(warmObj.modeledVEonly.fit{1}, coldObj.concww);
                warm_dens = polyval(warmObj.modeledDensity.fit, coldObj.concww/100);
                warm_VE = 1./(warm_dens) - warmObj.idealVolume(warmObj.temp, coldObj.concww);
                warm = idealVolume(coldObj, 'rt') + warm_VE;
                volfrac = 100*(coldObj.density.^-1 - warm)./warm;
                expfrac  = 100*(coldObj.modeledDensity.interpolated.^-1 - warmObj.modeledDensity.interpolated.^-1)./warmObj.modeledDensity.interpolated.^-1;
            elseif strcmp(dvolfrac, 'up')
                warm_dens = polyval(warmObj.modeledDensity.fit, coldObj.concww/100);
                warm_VE = 1./(warm_dens) - warmObj.idealVolume(warmObj.temp, coldObj.concww);
                warm = idealVolume(coldObj, 'rt') + warm_VE;
                volfrac = 100*(coldObj.floats.^-1 - warm)./warm;
                % needs continuous dDensity model: expfrac  = 100*(coldObj.modeledDensity.^-1 - warmObj.modeledDensity.^-1)./warmObj.modeledDensity.^-1;
                expfrac = [];
            elseif strcmp(dvolfrac, 'down')
                warm_dens = polyval(warmObj.modeledDensity.fit, coldObj.concww/100);
                warm_VE = 1./(warm_dens) - warmObj.idealVolume(warmObj.temp, coldObj.concww);
                warm = idealVolume(coldObj, 'rt') + warm_VE;
                volfrac = 100*(coldObj.sinks.^-1 - warm)./warm;
                % needs continuous dDensity model: expfrac  = 100*(coldObj.modeledDensity.^-1 - warmObj.modeledDensity.^-1)./warmObj.modeledDensity.^-1;
                expfrac = [];
            end
            
            
        end
        
        function [concentration tempRange] = findZero(coldObj, warmObj)
            % Calculates the concentration at which rapidly freezing to the
            % temperature of coldObj from the that of warmObj will yield no
            % change in density. coldObj and warmObj are assumed to be of
            % the cryoprotectant.
            
            % get volumetric contraction model
            conc = coldObj.concInterpolated;
            %warm = idealVolume(coldObj,'rt',conc)+predict(warmObj.modeledVE.fit{1}, conc) + polyval(warmObj.modeledVE.fit{2}, conc);
            warm_dens = warmObj.modeledDensity.interpolated;
            warm_VE = 1./(warm_dens) - warmObj.idealVolume(warmObj.temp, conc);
            warm = idealVolume(coldObj, 'rt', conc) + warm_VE;
            cold = coldObj.modeledDensity.interpolated;
            volfrac = 100*(cold.^-1 - warm)./warm;
            
            % find two concentrations surrounding zero
            len = length(volfrac);
            if sign(volfrac(1))==1 && sign(volfrac(end))==-1
                i=1;
                while sign(volfrac(i))==1
                    i=i+1;
                end
                if sign(volfrac(i))==0
                    concentration = i;
                    highc = 'stop';
                else
                    highc = i; lowc = i-1;
                end
            else
                error('Make sure the volume fraction has a zero point if you want to find a zero point.')
            end
            % weighted average of those two concentrations
            if ~ischar(highc)
                highv = volfrac(highc); lowv = volfrac(lowc);
                rise = highv-lowv;
                run = highc-lowc; % always 1
                slope  = rise/run; % units of volfrac/conc
                
                if abs(volfrac(highc)) > abs(volfrac(lowc))
                    concentration = lowc + lowv/slope;
                else
                    concentration = highc + highv/slope;
                end
            end
            
            % build tempRange, in the form [warmTemp, coldTemp]
            
            tempRange = strcat(int2str(warmObj.degtemp), 'K-', int2str(coldObj.degtemp),'K');
            
        end
        
        function [vec, rsq, rsq_adj, rmse] = getParameters(self)
            % returns a vector of parameters in a quartic regression, (a4,
            % a3, a2, a1, a0) of the excess volume in self
            vec = self.modeledVE.fit;
            [resid, rsq, rsq_adj, ~] = self.residCheck(4, self.VE, self.modeledVE.fitted, self.concww);
            sum = 0; len = length(resid);
            for i=1:len
                sum = sum + resid(i).^2;
            end
            rmse = (sum./len).^0.5;
        end
        
        function [vec, rsq, rsq_adj, rmse] = getParametersDens(self)
            vec = self.modeledDensity.fit;
            [resid, rsq, rsq_adj, ~] = self.residCheck(4, self.density, self.modeledDensity.fitted, self.concww/100);
            sum = 0; len = length(resid);
            for i=1:len
                sum = sum + resid(i).^2;
            end
            rmse = (sum./len).^0.5;
        end
        
        function printParametersDens(self)
            [vec, ~, rsq_adj, rmse] = getParametersDens(self);
            str = '';
            for i = 1:length(vec)
                str = strcat(str, ' &  ', ' ', num2str(round(vec(i), 4)));
            end
            str = strcat(str, ' &  ', num2str(round(rmse, 4,'significant')), ' &  ', num2str(round(rsq_adj, 4, 'significant')), '\n');
            fprintf(str)
        end
        
        function [vec, rsq, rsq_adj, rmse] = getParametersElecDens(self)
            vec = self.modeledElecDensity.fit;
            [resid, rsq, rsq_adj, ~] = self.residCheck(5, self.eDensity, self.modeledElecDensity.fitted, self.concww/100);
            sum = 0; len = length(resid);
            for i=1:len
                sum = sum + resid(i).^2;
            end
            rmse = (sum./len).^0.5;
        end
        
        function [conc, Dconc] = density2conc(self, density, Ddensity)
            %maps density (double) to a %w/w concentration in self.
            %Linearly searches through interpolated fit, and then linearly
            %interpolates the two nearest neighbors.  Expected that
            %self.modeledDensity monotonically decreases with i (i.e. only for Ethanol RT
            %or 2-Propanol RT). Includes unertainty
                conc = self.inversedensity(density);
                up = self.inversedensity(density+Ddensity);
                down = self.inversedensity(density-Ddensity);
                Dconc = max(abs(conc-up), abs(conc-down));           
        end
        
        function conc = inversedensity(self, density)
            %maps density (double) to a %w/w concentration in self.
            %Linearly searches through interpolated fit, and then linearly
            %interpolates the two nearest neighbors.  Expected that
            %self.modeledDensity monotonically decreases with i (i.e. only for Ethanol RT
            %or 2-Propanol RT)
            y = self.modeledDensity;
            x = self.concInterpolated;
            found = 0;
            try
                for i=1:length(x)
                    if found==0 && y(i)<density
                        a = [x(i-1) y(i-1)];
                        b = [x(i) y(i)];
                        found = 1;
                    end
                end
            catch
                sprintf('Error: Density out of Range');
            end
            
            
            dxdy = (b(1)-a(1))/(b(2)-a(2));
            
            conc = dxdy*(density-a(2)) + a(1);
        end
        
        function density = inverseconc(self, conc)
            %maps density (double) to a %w/w concentration in self.
            %Linearly searches through interpolated fit, and then linearly
            %interpolates the two nearest neighbors.  Expected that
            %self.modeledDensity monotonically decreases with i (i.e. only for Ethanol RT
            %or 2-Propanol RT)
            y = self.modeledDensity;
            x = self.concInterpolated;
            
            i = ceil(conc);
            a = [x(i-1) y(i-1)];
            b = [x(i) y(i)];
                      
            dxdy = (b(1)-a(1))/(b(2)-a(2));
            
            density = dxdy^-1*(conc-a(1)) + a(2);
        end
        
        function densityArray = inverseconcArray(self, concArray)
            %maps density (double) to a %w/w concentration in self.
            %Linearly searches through interpolated fit, and then linearly
            %interpolates the two nearest neighbors.  Expected that
            %self.modeledDensity monotonically decreases with i (i.e. only for Ethanol RT
            %or 2-Propanol RT)
            densityArray = zeros(length(concArray),1);
            for i = 1:length(concArray)
                densityArray(i) = self.inverseconc(concArray(i));
            end
        end
        
        function eDensity = makeElecDensity(self)
            h2o_M = 18.0; % molecular weight of water
            h2o_elecs = 10.0; % electrons/mol pure water
            avog_N = 6.022e23;
            cc_A3 = 1.00e-24;
            eDensity = self.density.*(self.concww/100.0*(self.e_molecule/self.mol_weight)+(100.0-self.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
        end
        
        function eDensityDiff = makeElecDensityDiff(self)
            if strcmp(self.temp, 'rt')
                h2o_ed = 0.3344;
                protein_ed = 0.43;
            else
                h2o_ed = 0.3144;
                protein_ed = 0.443;
            end
            norm = (protein_ed - h2o_ed)^2;
            eDensityDiff = struct('fitted', [], 'interpolated', []);
            eDensityDiff.fitted = (self.eDensity - protein_ed).^2/norm;
            eDensityDiff.interpolated = (self.modeledElecDensity - protein_ed).^2/norm;
        end
        
        %% Data Presentation
        
        function makePlots(self)
            %makes all figures, this should be modifiable to work in any basis
            %of concentration units (but isn't yet). (Same goes for error
            %bars.)
            
            close all
            
            %% Plot Density
            x = self.concInterpolated;
            y = self.modeledDensity;
            figure('Color','white');
            hold on
            box ON
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Density (g/mL)','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            if strcmp(self.temp, 'rt') % if room temperature
                A=plot(self.concww, self.density, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            else
                A = errorbar(self.concww, self.density, self.downuncertain, self.upuncertain,  'go','MarkerSize', 12,'MarkerFaceColor','g');%,'LineWidth',2);
                %A=plot(self.concww, self.density, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            end
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            
            plot(x,y,'g-','LineWidth',3.5)
            
            %% Plot Excess Specific Volume
            y = self.modeledVE.interpolated;
            figure('Color','white');
            hold on
            box ON
            
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Excess Specific Vol. (mL/g)','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            A=plot(self.concww, self.VE, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
            %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            
            plot(x,y,'g-','LineWidth',3.5)
            
            
            %% Plot of Specific Volume with asymptotes
            
            y = self.modeledSpecificVol;
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            
            plot(x,y,'g-','LineWidth',3.5)
            
            
            % calculate linear pure cryo asymptote
            fit = polyfit(x(91:101), y(91:101), 1);
            Vpurecryo = polyval(fit, x);
            
            % calculate specific volume of ideal mixture
            idealV = self.idealVolume(self.temp, x);
            
            % calculate pure water asymptotes
            fit = polyfit(x(1:11), y(1:11), 1);
            Vpurewater = polyval(fit, x);
            
            figure('Color','white');
            hold on
            hold on
            box ON
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Specific Vol. (mL/g)','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            A=plot(self.concww, self.SpecificVol, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            
            plot(x, idealV, 'k:', 'LineWidth', 2)
            %plot(x, Vpurecryo, 'g--', 'LineWidth', 2)
            %plot(x, Vpurewater, 'g--', 'LineWidth', 2)
            plot(x,self.modeledSpecificVol, 'g-', 'LineWidth',3.5)
            
            uistack(A, 'top')
            
            
            %% Plot apparent specific volumes
            
            y = self.Vapp.interpolated;
            figure('Color','white');
            hold on
            box ON
            
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('\upsilon^{cryo}_{app} (mL/g)','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            A=plot(self.concww, self.Vapp.measured, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
            %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            
            Id = self.Vapp.ID; yForced = self.Vapp.forced;
            plot(x(1:Id),yForced,'g--','LineWidth',3.5)
            plot(x(Id:end),y(Id:end),'g-','LineWidth',3.5)
            %plot(x,y,'g-','LineWidth',3.5)
            
            
            
            %% Plot apparent volume of water
            
            y = self.Vapph2o.interpolated;
            figure('Color','white');
            hold on
            box ON
            
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('\upsilon^{water}_{app} (mL/g)','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            A=plot(self.concww, self.Vapph2o.measured, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
            %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            
            
            Id = self.Vapph2o.ID; yForced = self.Vapph2o.forced;
            plot(x(1:Id),y(1:Id),'g-','LineWidth',3.5)
            plot(x(Id+1:end),yForced,'g--','LineWidth',3.5)
            %plot(x,y,'g-','LineWidth',3.5)
            uistack(A, 'top')
            
            
        end
        
        function plotDensity(self)
            %% Plot Density
            x = self.concInterpolated;
            y = self.modeledDensity.interpolated;
            figure('Color','white');
            hold on
            box ON
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Density (g/mL)','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y)+0.025;  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            if strcmp(self.temp, 'rt') % if room temperature
                A=plot(self.concww, self.density, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            else
                A = errorbar(self.concww, self.density, self.downuncertain, self.upuncertain,  'go','MarkerSize', 12,'MarkerFaceColor','g','LineWidth',2);
                %A=plot(self.concww, self.density, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            end
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            
            plot(x,y,'g-','LineWidth',3.5)
        end
        
        function comparePlots(coldObj, rtObj)
            % vol(77K)-vol(rt)/vol(rt) vs conc
            [volfrac, expfrac] = calcContraction(coldObj, rtObj);
            y = expfrac;
            x = coldObj.concInterpolated;
            
            figure('Color','white');
            hold on
            box ON
            
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Volume Fraction','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            
            %plotmax=6;%max([volfracG;volfracEG]);
            %plotmin=-10;%min([volfracG;volfracEG]);
            axis([xmin, xmax, ymin, ymax]);
            
            A=plot(coldObj.concww, volfrac , 'bo','MarkerSize', 12,'MarkerFaceColor','b');
            
            h = legend(coldObj.name);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
            %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
            %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            plot(x,y,'b-','LineWidth',3.5)
            
            
            %             % calculate errorbars for volume fraction
            %                 up = 100*(G.floats.^-1 - warmG)./warmG; up(1) = volfracG(1);
            %                 down = 100*(G.sinks.^-1 - warmG)./warmG; down(1)=volfracG(1);
            %             dvolfracG = struct('up',up-volfracG,'down',volfracG-down);
            %                 up = 100*(EG.floats.^-1 - warmEG)./warmEG; up(1) = volfracEG(1);
            %                 down = 100*(EG.sinks.^-1 - warmEG)./warmEG; down(1)=volfracEG(1);
            %             dvolfracEG = struct('up',up-volfracEG,'down',volfracEG-down);
            
            %A=plot(G.conc, volfracG, 'ro','MarkerSize', 12,'MarkerFaceColor','r');
            %A=errorbar(G.conc, volfracG,dvolfracG.up, dvolfracG.down, 'ro','MarkerSize', 12,'MarkerFaceColor','r','LineWidth', 2);
            %B=plot(EG.conc, volfracEG, 'bv','MarkerSize', 12,'MarkerFaceColor','b');
            %B=errorbar(EG.conc, volfracEG,dvolfracEG.up, dvolfracEG.down, 'bv','MarkerSize', 12,'MarkerFaceColor','b','LineWidth',2);
            %plot(GRT.conc, devGRT, 'ro','MarkerSize', 12,'MarkerFaceColor',[1 1 1])
            %plot(EGRT.conc, devEGRT, 'bv','MarkerSize', 12,'MarkerFaceColor',[1 1 1])
        end
        
        function makeTicks(self, xmin, xmax, dX,dx, ymin, ymax, dY, dy)
            %%% ADD CODE TO MAKE MINOR TICKS FOR X AND Y AXES
            
            %makes the major and minor x and y axis ticks for a plot. dX and dY are the
            %major axis tick distance and dx and dy are the minor tick distance
            
            % specify the minor grid vectors
            xminor = xmin:dx:xmax;
            xmajor = xmin:dX:xmax;
            yminor = ymin:dy:ymax;
            ymajor = ymin:dY:ymax;
            
            % specify the tick heights
            Dx = ymax-ymin;
            Dy = xmax-xmin;
            hX = Dx/75;  hY = Dy/75;
            hx = hX*3/4;  hy = hY*3/4;
            
            % make and plot ticks
            
            %X major
            [a,b] = self.prepareTicks(xmajor, ymin, ymin+hX);
            plot(a,b, 'k', 'LineWidth', 3.5)
            [a,b] = self.prepareTicks(xmajor, ymax, ymax-hX);
            plot(a,b, 'k', 'LineWidth', 3.5)
            
            %X minor
            [a,b] = self.prepareTicks(xminor, ymin, ymin+hx);
            plot(a,b, 'k', 'LineWidth', 3)
            [a,b] = self.prepareTicks(xminor, ymax, ymax-hx);
            plot(a,b, 'k', 'LineWidth', 3)
            
            %Y major
            [a,b] = self.prepareTicks(ymajor, xmin, xmin+hY);
            plot(b,a, 'k', 'LineWidth', 3.5)
            [a,b] = self.prepareTicks(ymajor, xmax, xmax-hY);
            plot(b,a, 'k', 'LineWidth', 3.5)
            
            %Y minor
            [a,b] = self.prepareTicks(yminor, xmin, xmin+hy);
            plot(b,a, 'k', 'LineWidth', 3)
            [a,b] = self.prepareTicks(yminor, xmax, xmax-hy);
            plot(b,a, 'k', 'LineWidth', 3)
        end
        
        function [a,b] = prepareTicks(self, ticks, bound, height)
            % specify the Y-position and the height of minor grids
            a = reshape([ticks;ticks;NaN(1,length(ticks))],1,length(ticks)*3);
            b = repmat([bound height NaN],1,length(ticks));
        end
        
        function printTable(coldObj, warmObj)
            % prints a table in the command window that displays density, excess volume, and contraction upon cooling as a function of temperature
            % coldObj is a cold solution set, and warmObj is a warmer
            % solution set from that same cryoprotectant agent
            %
            %                         & name                  &                            &                                 \\
            % \hline
            % Concentration (wt. frac.)    & Density (g/mL)        & Excess Volume (mL/g)       & Contraction upon Cooling %     \\
            % \hline
            %  entry      & entry      & entry      & entry      \\
            %  entry      & entry      & entry      & entry      \\
            %  entry      & entry      & entry      & entry      \\
            
            %data prep
            name = coldObj.name;
            conc = coldObj.concww;
            density = coldObj.density;
            VE = coldObj.VE;
            [volfrac, ~] = calcContraction(coldObj, warmObj);
            len  = length(conc);
            
            if len ~= length(density)
                if len ~= length(VE)
                    if len ~= length(volfrac)
                        error('Not all data vectors are of same length')
                    end
                end
            end
            
            %data print
            
            fprintf('\\begin{table}\n')
            fprintf('\\caption{}\n')
            fprintf('\\begin{tabular}{lccc}      % Alignment for each cell: l=left, c=center, r=right\n')
            
            
            fprintf('%s \t &  &  \t& \t\t\\\\\n', name);
            fprintf('\\hline \n')
            fprintf('$W_1$ &  $\\rho/(g/ml)$   &  $\\upsilon^E/(ml/g)$   &  $\\Delta V (\\%%)$  \\\\\n')
            fprintf('\\hline\n')
            for i = 1:len
                fprintf('%.4f \t\t& %0.4f \t& %.4f \t& %.2f \\\\\n',conc(i)/100, density(i), VE(i), volfrac(i))
            end
            
            fprintf('\\end{tabular}\n')
            fprintf('\\end{table}\n\n')
            
        end
        
        function [Ddens,Pdens, DVE, Dc] = printUncertainty(coldObj, warmObj)
            % prints the mean uncertainties in density, excess volume, and contraction upon cooling
            % coldObj is a cold solution set, and warmObj is a warmer
            % solution set from that same cryoprotectant agent
            %
            
            name = coldObj.name;
            
            l = length(coldObj.concww);
            sumdens = 0;    sumdVE = 0; sumpVE = 0;
            sumdc = 0;     sumpdens = 0;
            [dcu, ~] = coldObj.calcContraction(warmObj, 'up');
            [c, ~] = coldObj.calcContraction(warmObj);
            [dcd, ~] = coldObj.calcContraction(warmObj, 'down');
            for i = 2:l
                if sign(coldObj.upuncertain(i))==-1 || sign(coldObj.downuncertain(i))==-1
                    if sign(coldObj.upuncertain(i))==-1
                        coldObj.upuncertain(i) = -coldObj.upuncertain(i);
                    end
                    if sign(coldObj.downuncertain(i))==-1
                        coldObj.downuncertain(i) = -coldObj.downuncertain(i);
                    end
                    %error('sign error in raw a float or sink value');
                end
                dens = (coldObj.upuncertain(i) + coldObj.downuncertain(i))./2;
                ve = (abs(coldObj.dVE.up(i)) + abs(coldObj.dVE.down(i)))./2;
                sumdens = sumdens + dens.^2;
                sumpdens = sumpdens + (dens./coldObj.density(i)).^2;
                sumdVE = sumdVE + ve.^2;
                if i ~= l
                    sumpVE = sumpVE + (ve/coldObj.VE(i)).^2;
                end
                sumdc = sumdc + (0.5*(dcu(i)-c(i)) + 0.5*(dcd(i)-c(i))).^2;
            end
            % no uncertainty added for pure water, so divide by l-1
            Ddens = sumdens.^0.5/(l-1); %mean uncertainty of cold density
            DVE = sumdVE.^0.5/(l-1); %mean uncertainty of cold VE
            Pdens = 100*sumpdens.^0.5/(l-1); %mean percent uncertainty in cold density
            PVE = 100*sumpVE.^0.5/(l-2); %mean uncertainty of cold V
            Dc = sumdc.^0.5/(l-1); %mean uncertainty of contraction
            
            [~, rsq, ~] = coldObj.residCheck(coldObj.VE, coldObj.modeledVE.fitted, coldObj.concww);
            [~, rsqrt, ~] = warmObj.residCheck(warmObj.VE, warmObj.modeledVE.fitted, warmObj.concww);
            
            
            fprintf('\n for %s, we have: \n', name);
            fprintf('rms uncertainty (g/ml) in density is %f \n', Ddens);
            fprintf('rms percent uncertainty in density  is %f \n', Pdens);
            fprintf('rms uncertainty (ml/g) in VE  is %f \n', DVE);
            fprintf('rms percent uncertainty in VE  is %f \n', PVE);
            fprintf('rms uncertainty (percent contraction) in contraction is %f .\n', Dc);
            fprintf('r-squared (cold) is %f .\n', rsq);
            fprintf('r-squared (warm) is %f .\n', rsqrt);
        end
        
        function plotElecDensity(self)
            x = self.concInterpolated;
            y = self.modeledElecDensity;
            figure('Color','white');
            hold on
            box ON
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Electron Density (e-/A^3)','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y)+0.01;  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            h2o_M = 18.0; % molecular weight of water
            h2o_elecs = 10.0; % electrons/mol pure water
            avog_N = 6.022e23; % avo gadro's number
            cc_A3 = 1.00e-24; % cm^3/A^3
            if strcmp(self.temp, 'rt') % if room temperature
                A=plot(self.concww, self.eDensity, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            else
                upuncertainElec = self.upuncertain.*(self.concww/100.0*(self.e_molecule/self.mol_weight)+(100.0-self.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                downuncertainElec = self.downuncertain.*(self.concww/100.0*(self.e_molecule/self.mol_weight)+(100.0-self.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                A = errorbar(self.concww, self.eDensity, downuncertainElec, upuncertainElec,  'go','MarkerSize', 12,'MarkerFaceColor','g','LineWidth',2);
            end
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            
            plot(x,y,'g-','LineWidth',3.5)
        end
        
        function plotElecDensityDiff(self)
            x = self.concInterpolated;
            y = self.eDensityDiff.interpolated;
            figure('Color','white');
            hold on
            box ON
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Squared Electron Density Difference','FontSize', 30,'FontName','Arial');
            ymin=min(y);    ymax=max(y)+0.01;  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
            h2o_M = 18.0; % molecular weight of water
            h2o_elecs = 10.0; % electrons/mol pure water
            avog_N = 6.022e23; % avo gadro's number
            cc_A3 = 1.00e-24; % cm^3/A^3
            if strcmp(self.temp, 'rt') % if room temperature
                A=plot(self.concww, self.eDensityDiff.fitted, 'go','MarkerSize', 12,'MarkerFaceColor','g');
            else
                upuncertainElec = self.upuncertain.*(self.concww/100.0*(self.e_molecule/self.mol_weight)+(100.0-self.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                downuncertainElec = self.downuncertain.*(self.concww/100.0*(self.e_molecule/self.mol_weight)+(100.0-self.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                h2o_ed = 0.3144; protein_ed = 0.443; norm = (protein_ed - h2o_ed)^2;
                upuncertainElecDiff = (upuncertainElec).^2/norm; % is this correct?
                downuncertainElecDiff = (downuncertainElec).^2/norm; % is this correct?
                A = errorbar(self.concww, self.eDensityDiff.fitted, downuncertainElecDiff, upuncertainElecDiff,  'go','MarkerSize', 12,'MarkerFaceColor','g','LineWidth',2);
            end
            
            h = legend(self.label);
            set(h,'Location', 'Best', 'Orientation', ...
                'vertical', 'Fontname', 'Arial','FontSize', 30,...
                'EdgeColor',[1 1 1]);
            % 'YColor',[1 1 1],'XColor',[1 1 1]
            self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            uistack(A, 'top')
            
            plot(x,y,'g-','LineWidth',3.5)
            
        end
        
        
        %% Miscillaneous
        function disp(self)
            % Display self, if not empty, in this format: (left,right)
            % If empty, display 'Empty <classname>'
            if isempty(self)
                fprintf('Empty %s\n', class(self))
            elseif length(self)>1
                disp@handle(self)
            else
                %fprintf('(%f,%f)\n', self.left, self.right)
            end
        end
        
        %doesn't work yet
        
        function newconc = unitConvert(conc, bestModel, u, molarmass, rho, t)
            % convert conc in units of % w/v to newconc in units of u.  bestModel is
            % the best model (assumed to be polynomial) for the density of
            % a given mixture. type is 'g' or 'eg'. t is the target temperature, rt is
            % 298K and 20 is 293K.
            % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % % concentration units are                         %
            % % u = 1: mole fraction                            %
            % % u = 2: % w/w                                    %
            % % u = 3: % w/v                                    %
            % % % % % % % % % % % % % % % % % % % % % % % % % % %
            
            if u~=3  % if units aren't %w/v, convert to %w/w
                % %w/w = %w/v divided by average density of solution
                % determine density of data before cooling and convert to %w/w
                %dens = polyval(bestModel, conc, 1);
                if nargin==5
                    t='twenty';
                end
                
                
                if molarmass > 70
                    type = 'g';
                else
                    type = 'eg';
                end
                if strcmp(type, 'g')
                    if strcmp(t, 'rt')
                        dens = transTheory(conc, rho, type, 'rt') + predict(bestModel{1}{1},conc)+ polyval(bestModel{2}, conc);
                    else
                        dens = transTheory(conc, rho, type, 'twenty') + predict(bestModel{1}{1},conc)+ polyval(bestModel{2}, conc);
                    end
                else
                    dens = transTheory(conc, rho, type, 'rt') + predict(bestModel{1}{1},conc)+ polyval(bestModel{2}, conc);
                end
                newconc = conc./dens;
                
                % if u==1, convert to mole fraction, x
                if u==1
                    % initialize molar mass of water (g/mol)
                    Mh2o = 18.01528;
                    % convert to mole fraction
                    newconc = (molarmass/Mh2o*(100./newconc-1)+1).^-1;
                end
            else
                newconc = conc;
            end
        end
        
        %doesn't work yet
        
        function [Y_hat rsq bestecm kitty] = bestECM(X,Y)
            % finds the best ECM of a gaussian, a quartic, and a quadratic regression.
            [Y1_hat ecm1] = quarticReg(X,Y);
            [a rsq1 ~ ]  = residCheck(Y1_hat, Y, X);
            [Y2_hat ecm2]= quadraticReg(X,Y);
            [a rsq2 ~ ]  = residCheck(Y2_hat, Y, X);
            [Y3_hat ecm3] = gaussianReg(X,Y);
            [a rsq3 ~ ]  = residCheck(Y3_hat, Y, X);
            
            if  rsq1>rsq2 && rsq1>rsq3
                Y_hat = Y1_hat; rsq = rsq1;
                kitty = 'quartic';
                bestecm = ecm1;
            elseif  rsq2 > rsq3
                Y_hat = Y2_hat; rsq = rsq2;
                kitty = 'quadratic';
                bestecm = ecm2;
            else % it is rsq3
                Y_hat = Y3_hat; rsq = rsq3;
                kitty = 'gaussian';
                bestecm = ecm3;
            end
        end
        
        %obsolete
        
        
        function [Y_hat fit] = quadraticReg(X,Y)
            % quadratic regression
            fit = polyfit(X, Y, 2);
            Y_hat = polyval(fit, X);
        end
        
        %obsolete
        
        function density = transTheory(conc, rho, type, t)
            % calculates the predicted density as a function of concentration
            % based on simplistic theoretical values. rho is a structure of relevant densities,
            % and concentration is a vector. Type is 'g' or 'eg'. Supposes conc unit
            % is %w/v. t is a string, either 'cold' or 'rt'.
            
            
            % if the data is cold
            if strcmp(t, 'cold')
                
                if strcmp(type, 'g')
                    
                    % calculate mass fraction
                    x = rho.h2o.rt*(100./conc - 1/(rho.g.rt));
                    
                    % calculate predicted density
                    
                    density = rho.g.cold*rho.h2o.cold*(1+x)./(rho.h2o.cold + rho.g.cold*x);
                    
                elseif strcmp(type, 'eg')
                    
                    % calculate mass fraction
                    x = rho.h2o.rt*(100./conc - 1/(rho.eg.rt));
                    
                    % calculate predicted density
                    density = rho.eg.cold*rho.h2o.cold*(1+x)./(rho.h2o.cold + rho.eg.cold*x);
                    
                else
                    error('input type as either `g` or  `eg`');
                end
                
                % set pure water to the density of pure water
                l = length(density);
                for i  = 1:l
                    if isnan(density(i))
                        density(i) = rho.h2o.cold;
                    end
                end
                
                
            elseif strcmp(t, 'rt')
                if strcmp(type, 'g')
                    
                    rho.h2o.rt=1;
                    rho.g.rt = 1.25802;
                    
                    % calculate mass fraction
                    x = rho.h2o.rt*(100./conc - 1/(rho.g.rt));
                    
                    % calculate predicted density
                    
                    density = rho.g.rt*rho.h2o.rt*(1+x)./(rho.h2o.rt + rho.g.rt*x);
                    
                elseif strcmp(type, 'eg')
                    
                    % calculate mass fraction
                    
                    x = rho.h2o.rt*(100./conc - 1/(rho.eg.rt));
                    
                    % calculate predicted density
                    density = rho.eg.rt*rho.h2o.rt*(1+x)./(rho.h2o.rt + rho.eg.rt*x);
                    
                else
                    error('input type as either `g` or  `eg`');
                end
                
                % set pure water to the density of pure water
                l = length(density);
                for i  = 1:l
                    if isnan(density(i))
                        density(i) = rho.h2o.rt;
                    end
                end
            elseif strcmp(t, 'twenty')
                if strcmp(type, 'g')
                    
                    rho.h2o.rt=1;
                    rho.g.rt = 1.26108;
                    
                    % calculate mass fraction
                    x = rho.h2o.rt*(100./conc - 1/(rho.g.rt));
                    
                    % calculate predicted density
                    
                    density = rho.g.rt*rho.h2o.rt*(1+x)./(rho.h2o.rt + rho.g.rt*x);
                elseif strcmp(type, 'eg')
                    error('I couldn`t find that data yet. :-(');
                end
            else
                error('input temperature of data set');
            end
        end
        
        %obsolete
        
        function [lnx, index] = logtrans(x)
            %transforms data of vector, x, but throws out sad data (because of 0 %w/v
            %conc).
            
            a = length(x);
            lnx = log(x);
            
            index = [];
            for i  = 1:a
                if ~isinf(lnx(i))
                    index = [index i];
                end
            end
        end
        
        %obsolete
        
        function out = extrapolateTrans2zero(conc, rho, type, t)
            % does a nonlinear fit on high concentration data to provide a stochastic
            % estimate for the y-intercept of a simple mixture model.  This will be
            % applied to find the molecular density of pure water.
            %
            % % select relevent parameters
            %
            %
            % % calculate input
            %
            % % define model
            %
            %
            % % regress with a nonlinear fit
            % myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*exp(-(x-b(2)).^2/(2*b(3)^2))+b(4)), [B0, B1, B2, B3]);
            %
            % density = rho.g.cold*rho.h2o.cold*(1+x)./(rho.h2o.cold + rho.g.cold*x);
            %
            % % extract relevant values
            %
            %
            %         % calculate mass fraction
            %         x = rho.h2o.rt*(100./conc - 1/(rho.g.rt));
            %
            %         % calculate predicted density
            %
            %       [B0 i] = max(Y); % Amplitude
            % B1 = X(i);  % Vertical shift
            %
            % [A I]  = sort(Y);
            % j = 1; k = 0;
            % while j < length(Y) && k ==0
            %     if A(j) > B0/2
            %         B2 = abs(X(I(j)) - B1);  % Width of curve. ~ Standard deviation or width-of-half-max.
            %         k = 1;
            %     end
            %     j = j +1;
            % end
            % B3 = 0; % vertical shift
            % myFit = NonLinearModel.fit(X,Y, @(b,x)(b(1)*exp(-(x-b(2)).^2/(2*b(3)^2))+b(4)), [B0, B1, B2, B3]);
            %  %set initial values
            %     amp = myFit.Coefficients{1,1}; %height of gaussian
            %     mean = myFit.Coefficients{2,1}; %middle of gaussian
            %     B0 = myFit.Coefficients{3,1}; %width of curve
            %     shift = myFit.Coefficients{4,1}; %vertical shift
            %
            %
        end
        
        %doesn't work yet.  Doesn't need to work yet.
    end
    
    
end

