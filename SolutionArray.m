classdef SolutionArray < handle
    %R. Dan Aug 2017
    %Processed density measurements for cold cryoprotectants.
    %   enter "myArray = SolutionArray('cold');" to generate an array of
    %   SolutionSet objects, which can operate the given "methods".
    %
    %   List of ID's:
    %       'cold': all 77K data sets
    %       'warm': all ~298K data sets 
   properties
      Sol
   end
   methods
      function obj = SolutionArray(ID)
         if nargin ~= 0
            if strcmp(ID,'cold')
                N = 7;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('PPG');
                 obj(2,1).Sol = SolutionSet('PEG 200');
                 obj(3,1).Sol = SolutionSet('MPD');
                 obj(4,1).Sol = SolutionSet('Ethanol');
                 obj(5,1).Sol = SolutionSet('2-Propanol');
                 obj(6,1).Sol = SolutionSet('Glycerol');
                 obj(7,1).Sol = SolutionSet('Ethylene Glycol');
                 % obj(8,1).Sol = SolutionSet('Methanol');
            
            elseif strcmp(ID,'warm')
                N = 7;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('PPG RT');
                 obj(2,1).Sol = SolutionSet('PEG 200 RT');
                 obj(3,1).Sol = SolutionSet('MPD RT');
                 obj(4,1).Sol = SolutionSet('Ethanol RT');
                 obj(5,1).Sol = SolutionSet('2-Propanol RT');
                 obj(6,1).Sol = SolutionSet('Glycerol RT');
                 obj(7,1).Sol = SolutionSet('Ethylene Glycol RT');
                 obj(8,1).Sol = SolutionSet('Methanol RT');
                 
            elseif strcmp(ID,'warm_e')
                N = 5;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('MPD RT');
                 obj(2,1).Sol = SolutionSet('Ethanol RT');
                 obj(3,1).Sol = SolutionSet('2-Propanol RT');
                 obj(4,1).Sol = SolutionSet('Glycerol RT');
                 obj(5,1).Sol = SolutionSet('Ethylene Glycol RT');
                 
            elseif strcmp(ID,'cold_e')
                N = 5;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('MPD');
                 obj(2,1).Sol = SolutionSet('Ethanol');
                 obj(3,1).Sol = SolutionSet('2-Propanol');
                 obj(4,1).Sol = SolutionSet('Glycerol');
                 obj(5,1).Sol = SolutionSet('Ethylene Glycol');
           
            elseif strcmp(ID,'cold new')
                N = 5;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('PEG 200');
                 obj(2,1).Sol = SolutionSet('PPG');
                 obj(3,1).Sol = SolutionSet('MPD');
                 obj(4,1).Sol = SolutionSet('2-Propanol');
                 obj(5,1).Sol = SolutionSet('Ethanol');

            elseif strcmp(ID,'warm new')
                N = 5;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('PEG 200 RT');
                 obj(2,1).Sol = SolutionSet('PPG RT');
                 obj(3,1).Sol = SolutionSet('MPD RT');
                 obj(4,1).Sol = SolutionSet('2-Propanol RT');
                 obj(5,1).Sol = SolutionSet('Ethanol RT');
                 
            elseif strcmp(ID,'cold mono')
                N = 2;
                 obj(1,1).Sol = SolutionSet('2-Propanol');
                 obj(2,1).Sol = SolutionSet('Ethanol');
                 
            elseif strcmp(ID,'Ethanol')
                N = 2;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('Ethanol');
                 obj(2,1).Sol = SolutionSet('Ethanol RT');        
            elseif strcmp(ID,'PEG 200')
                N = 3;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('PEG 200');
                 obj(2,1).Sol = SolutionSet('PEG 200 RT'); 
                 obj(3,1).Sol = SolutionSet('PEG 200 HT'); 
            elseif strcmp(ID,'2-Propanol')
                N = 2;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('2-Propanol');
                 obj(2,1).Sol = SolutionSet('2-Propanol RT'); 
            elseif strcmp(ID,'Glycerol')
                N = 2;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('Glycerol');
                 obj(2,1).Sol = SolutionSet('Glycerol RT');      
            elseif strcmp(ID,'Ethylene Glycol')
                N = 2;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('Ethylene Glycol');
                 obj(2,1).Sol = SolutionSet('Ethylene Glycol RT');      
            elseif strcmp(ID,'MPD')
                N = 2;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('MPD');
                 obj(2,1).Sol = SolutionSet('MPD RT');      
            elseif strcmp(ID,'Methanol')
                N = 2;
                 obj(N,1) = SolutionArray; 
                 obj(1,1).Sol = SolutionSet('Methanol');
                 obj(2,1).Sol = SolutionSet('Methanol RT'); 
            end
         end
      end
      
      function overlayVApp(self)
          close all
          
            N = length(self);
            x = self(1,1).Sol.concInterpolated;
            % build label vector
            vlabels = [];
            for i = 1:N
                vlabels{i} = strcat(self(i,1).Sol.label(1:end-1),' K');
            end
            % build color vector
            colors = get(0,'DefaultAxesColorOrder');
            colors = cat(1, colors, [0 0 0]);
            
            % plot apparent specific volumes of cryoprotectants
            
            figure('Color','white');
            hold on
            box ON 
            
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('\upsilon^{cryo}_{app} (mL/g)','FontSize', 30,'FontName','Arial');
            for i = 1:N

                A(i)=plot(self(i,1).Sol.concww, self(i,1).Sol.Vapp.measured, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            end
            %ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, 0.5, 1.3]);
            h = legend(vlabels);
                set(h,'Location', 'Best', 'Orientation', ...
                    'vertical', 'Fontname', 'Arial','FontSize', 30,...
                    'EdgeColor',[1 1 1]);
                    % 'YColor',[1 1 1],'XColor',[1 1 1]
            for i = 1:N
                y = self(i,1).Sol.Vapp.interpolated;
                %plot(x,y,'-', 'Color', colors(i,:),'LineWidth',3.5)
                Id = self(i,1).Sol.Vapp.ID; yForced = self(i,1).Sol.Vapp.forced;
                plot(x(1:Id),yForced,'--','LineWidth',3.5,'Color', colors(i,:))
                plot(x(Id:end),y(Id:end),'-','LineWidth',3.5,'Color', colors(i,:))
                uistack(A, 'top')
            end

                
               
            % plot apparent specific volumes of water
            
            figure('Color','white');
            hold on
            box ON 
            
            set(gca,'FontSize', 30)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('\upsilon^{water}_{app} (mL/g)','FontSize', 30,'FontName','Arial');
            for i = 1:N

                A(i)=plot(self(i,1).Sol.concww, self(i,1).Sol.Vapph2o.measured, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
            end
            %ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, 0.7, 1.1]);
            h = legend(vlabels);
                set(h,'Location', 'Best', 'Orientation', ...
                    'vertical', 'Fontname', 'Arial','FontSize', 30,...
                    'EdgeColor',[1 1 1]);    
                    % 'YColor',[1 1 1],'XColor',[1 1 1]
            for i = 1:N
                y = self(i,1).Sol.Vapph2o.interpolated;
                %plot(x,y,'-', 'Color', colors(i,:),'LineWidth',3.5)
                Id = self(i,1).Sol.Vapph2o.ID; yForced = self(i,1).Sol.Vapph2o.forced;
                plot(x(1:Id),y(1:Id),'-','LineWidth',3.5, 'Color', colors(i,:))
                plot(x(Id+1:end),yForced,'--','LineWidth',3.5,'Color', colors(i,:))
                uistack(A, 'top')
            end
            
                
      end
      
      function overlayVE(self)
        close all
          
            N = length(self);
            x = self(1,1).Sol.concInterpolated;
            % build label vector
            vlabels = [];
            for i = 1:N
                vlabels{i} = strcat(self(i,1).Sol.label(1:end-1),' K');
            end
            
            % build color vector
            %colors=jet(N);  % n is the number of different items you have
            colors = get(0,'DefaultAxesColorOrder');
            colors = cat(1, colors, [0 0 0]);
            
            % plot excess volumes of cryoprotectants
            
            figure('Color','white');
            hold on
            box ON 
            
            set(gca,'FontSize', 25)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('\upsilon^{E} (mL/g)','FontSize', 30,'FontName','Arial');
            y = [];
            bonus = 0;
            for i = 1:N
                if strcmp(self(i,1).Sol.temp,'cold')
                    A(i)=errorbar(self(i,1).Sol.concww, self(i,1).Sol.VE,self(i,1).Sol.dVE.down,self(i,1).Sol.dVE.up,...
                        'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:),'LineWidth',2);
                    set(A(i),'Color',colors(i,:));
                else
                    A(i)=plot(self(i,1).Sol.concww, self(i,1).Sol.VE, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                    if strcmp(self(i,1).Sol.name,'2-Propanol RT')
                        bonus = 0.005;
                    end
                end
                
                %A(i)=plot(self(i,1).Sol.concww, self(i,1).Sol.VE, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
                y = [y, self(i,1).Sol.VE'];
            end
            
            ymin=min(y);    ymax=max(y)+bonus;  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
%             self(1,1).Sol.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
%                 ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20));
            ax = gca; ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
            clear y
            h = legend(vlabels);
                set(h,'Location', 'Best', 'Orientation', ...
                    'vertical', 'Fontname', 'Arial','FontSize', 25,...
                    'EdgeColor',[1 1 1]);
                    % 'YColor',[1 1 1],'XColor',[1 1 1]
            for i = 1:N
                y = self(i,1).Sol.modeledVE;
                %plot(x,y,'-', 'Color', colors(i,:),'LineWidth',3.5)
                %Id = self(i,1).Sol.VE.ID; yForced = self(i,1).Sol.Vapp.forced;
                %plot(x(1:Id),yForced,'--','LineWidth',3.5,'Color', colors(i,:))
                plot(x,y,'-','LineWidth',3.5,'Color', colors(i,:))
                uistack(A, 'top')
            end                        
      end
      
      function overlayDensity(self)
            close all
          
            N = length(self);
            x = self(1,1).Sol.concInterpolated;
            % build label vector
            vlabels = [];
            for i = 1:N
                vlabels{i} = strcat(self(i,1).Sol.label(1:end-1),' K');
            end
            
            % build color vector
            %colors=jet(N);  % n is the number of different items you have
            colors = get(0,'DefaultAxesColorOrder');
            colors = cat(1, colors, [0 0 0]);
            
            % plot excess volumes of cryoprotectants
            
            figure('Color','white');
            hold on
            box ON 
            
            set(gca,'FontSize', 25)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Density (g/mL)','FontSize', 30,'FontName','Arial');
            y = [];
            for i = 1:N
                if strcmp(self(i,1).Sol.temp,'cold')
                    A(i)=errorbar(self(i,1).Sol.concww, self(i,1).Sol.density,self(i,1).Sol.downuncertain,self(i,1).Sol.upuncertain,...
                        'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:),'LineWidth',2);
                    set(A(i),'Color',colors(i,:));
                else
                    A(i)=plot(self(i,1).Sol.concww, self(i,1).Sol.density, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                end
                   %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
                y = [y, self(i,1).Sol.density'];
                uistack(A(i), 'top')
            end
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            %self(1,1).Sol.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
            %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20));
            %set(gca,'YMinorTick','on')
            ax = gca; ax.YMinorTick = 'on'; ax.XMinorTick = 'on';

            %h = legend(vlabels);
            %h = legend(strcat(self(1,1).Sol.label(end-3:end-1),' K')); 
            %set(h,'Location', 'Best', 'Orientation', ...
            %        'vertical', 'Fontname', 'Arial','FontSize', 25,...
            %        'EdgeColor',[1 1 1]);
                    % 'YColor',[1 1 1],'XColor',[1 1 1]
            for i = 1:N
                y = self(i,1).Sol.modeledDensity.interpolated;
                %plot(x,y,'-', 'Color', colors(i,:),'LineWidth',3.5)
                %Id = self(i,1).Sol.Vapp.ID; yForced = self(i,1).Sol.Vapp.forced;
                %plot(x(1:Id),yForced,'--','LineWidth',3.5,'Color', colors(i,:))
                plot(x,y,'-','LineWidth',3.5,'Color', colors(i,:))
                uistack(A, 'top')
            end

                
                
      end
      
      function overlayElecDensity(self)
          close all
          
            N = length(self);
            x = self(1,1).Sol.concInterpolated;
            % build label vector
            vlabels = [];
            for i = 1:N
                vlabels{i} = strcat(self(i,1).Sol.label(1:end-1),' K');
            end
            
            % build color vector
            %colors=jet(N);  % n is the number of different items you have
            colors = get(0,'DefaultAxesColorOrder');
            colors = cat(1, colors, [0 0 0]);
            
            % plot excess volumes of cryoprotectants
            
            figure('Color','white');
            hold on
            box ON 
            
            set(gca,'FontSize', 25)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Electron Density (e-/A^3)','FontSize', 30,'FontName','Arial');
            y = [];
            h2o_M = 18.0; % molecular weight of water
            h2o_elecs = 10.0; % electrons/mol pure water
            avog_N = 6.022e23; % avo gadro's number
            cc_A3 = 1.00e-24; % cm^3/A^3
            for i = 1:N
                if strcmp(self(i,1).Sol.temp,'cold')
                    upuncertainElec = self(i,1).Sol.upuncertain.*(self(i,1).Sol.concww/100.0*(self(i,1).Sol.e_molecule/self(i,1).Sol.mol_weight)+(100.0-self(i,1).Sol.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                    downuncertainElec = self(i,1).Sol.downuncertain.*(self(i,1).Sol.concww/100.0*(self(i,1).Sol.e_molecule/self(i,1).Sol.mol_weight)+(100.0-self(i,1).Sol.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                    A(i)=errorbar(self(i,1).Sol.concww, self(i,1).Sol.eDensity,downuncertainElec,upuncertainElec,...
                        'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:),'LineWidth',2);
                    set(A(i),'Color',colors(i,:));
                else
                    A(i)=plot(self(i,1).Sol.concww, self(i,1).Sol.eDensity, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                end
                   %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
                y = [y, self(i,1).Sol.eDensity'];
                uistack(A(i), 'top')
            end
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
%             self(1,1).Sol.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
%                 ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20));
            ax = gca; ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
            h = legend(vlabels);
                set(h,'Location', 'Best', 'Orientation', ...
                    'vertical', 'Fontname', 'Arial','FontSize', 25,...
                    'EdgeColor',[1 1 1]);
                    % 'YColor',[1 1 1],'XColor',[1 1 1]
            for i = 1:N
                y = self(i,1).Sol.modeledElecDensity;
                %plot(x,y,'-', 'Color', colors(i,:),'LineWidth',3.5)
                %Id = self(i,1).Sol.Vapp.ID; yForced = self(i,1).Sol.Vapp.forced;
                %plot(x(1:Id),yForced,'--','LineWidth',3.5,'Color', colors(i,:))
                plot(x,y,'-','LineWidth',3.5,'Color', colors(i,:))
                uistack(A, 'top')
            end
          
          
          
      end
      
      function overlayContrast(self)
          close all
          
            N = length(self);
            x = self(1,1).Sol.concInterpolated;
            % build label vector
            vlabels = [];
            for i = 1:N
                vlabels{i} = strcat(self(i,1).Sol.label(1:end-1),' K');
            end
            
            % build color vector
            %colors=jet(N);  % n is the number of different items you have
            colors = get(0,'DefaultAxesColorOrder');
            colors = cat(1, colors, [0 0 0]);
            
            % plot excess volumes of cryoprotectants
            
            figure('Color','white');
            hold on
            box ON 
            
            set(gca,'FontSize', 25)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('Normalized Squared Contrast','FontSize', 30,'FontName','Arial');
            y = [];
            h2o_M = 18.0; % molecular weight of water
            h2o_elecs = 10.0; % electrons/mol pure water
            avog_N = 6.022e23; % avo gadro's number
            cc_A3 = 1.00e-24; % cm^3/A^3
            for i = 1:N
                if strcmp(self(i,1).Sol.temp,'cold')
                    upuncertainElec = self(i,1).Sol.upuncertain.*(self(i,1).Sol.concww/100.0*(self(i,1).Sol.e_molecule/self(i,1).Sol.mol_weight)+(100.0-self(i,1).Sol.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                    downuncertainElec = self(i,1).Sol.downuncertain.*(self(i,1).Sol.concww/100.0*(self(i,1).Sol.e_molecule/self(i,1).Sol.mol_weight)+(100.0-self(i,1).Sol.concww)/100.0*(h2o_elecs/h2o_M))*avog_N*cc_A3;
                    h2o_ed = 0.3144; protein_ed = 0.443; norm = (protein_ed - h2o_ed)^2;
                    upuncertainElecDiff = (upuncertainElec).^2/norm; % is this correct?
                    downuncertainElecDiff = (downuncertainElec).^2/norm;
                    A(i)=errorbar(self(i,1).Sol.concww, self(i,1).Sol.eDensityDiff.fitted,downuncertainElecDiff,upuncertainElecDiff,...
                        'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:),'LineWidth',2);
                    set(A(i),'Color',colors(i,:));
                else
                    A(i)=plot(self(i,1).Sol.concww, self(i,1).Sol.eDensityDiff.fitted, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                end
                   %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
                y = [y, self(i,1).Sol.eDensityDiff.fitted'];
                uistack(A(i), 'top')
            end
            ymin=min(y);    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
%             self(1,1).Sol.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
%                 ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20));

            ax = gca; ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
            h = legend(vlabels);
                set(h,'Location', 'Best', 'Orientation', ...
                    'vertical', 'Fontname', 'Arial','FontSize', 25,...
                    'EdgeColor',[1 1 1]);
                    % 'YColor',[1 1 1],'XColor',[1 1 1]
            for i = 1:N
                y = self(i,1).Sol.eDensityDiff.interpolated;
                %plot(x,y,'-', 'Color', colors(i,:),'LineWidth',3.5)
                %Id = self(i,1).Sol.Vapp.ID; yForced = self(i,1).Sol.Vapp.forced;
                %plot(x(1:Id),yForced,'--','LineWidth',3.5,'Color', colors(i,:))
                plot(x,y,'-','LineWidth',3.5,'Color', colors(i,:))
                uistack(A, 'top')
            end
          
      end
      
      function overlayContraction(self,other)
          close all
          % self is cold, other is warm
          
            N = length(self);
            x = self(1,1).Sol.concInterpolated;
            % build label vector
            vlabels = [];
            for i = 1:N               
                vlabels{i} = self(i,1).Sol.label;
                [volfrac{i}, expfrac{i}] = calcContraction(self(i,1).Sol, other(i,1).Sol);
                
                if strcmp(self(i,1).Sol.temp,'cold')
                    up{i} = abs(calcContraction(self(i,1).Sol, other(i,1).Sol,'up') - volfrac{i}); up{i}(1)=0;%volfracG(1);
                    down{i} = abs(calcContraction(self(i,1).Sol, other(i,1).Sol,'down') - volfrac{i}); down{i}(1)=0;%volfracG(1);
                end
            end
            
            % build color vector
            %colors=jet(N);  % n is the number of different items you have
            colors = get(0,'DefaultAxesColorOrder');
            colors = cat(1, colors, [0 0 0]);
            
            % plot excess volumes of cryoprotectants
            
            figure('Color','white');
            hold on
            box ON 
            
            set(gca,'FontSize', 25)
            set(gca, 'LineWidth', 3)
            xlabel('Concentration (% w/w)','FontSize', 30,'FontName','Arial');
            ylabel('% Contraction','FontSize', 30,'FontName','Arial');
            y = [];
            for i = 1:N
                
                if strcmp(self(i,1).Sol.temp,'cold')
                    A(i)=errorbar(self(i,1).Sol.concww, volfrac{i},down{i},up{i},...
                        'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:),'LineWidth',2);
                    set(A(i),'Color',colors(i,:));
                else
                    A(i)=plot(self(i,1).Sol.concww, volfrac{i}, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                end
                                                              
               % A(i)=plot(self(i,1).Sol.concww, volfrac{i}, 'o','MarkerEdgeColor',colors(i,:),'MarkerSize', 12,'MarkerFaceColor',colors(i,:));
                %self.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
                %    ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20))
                y = [y, volfrac{i}'];
            end
            ymin=min(y) - 2;    ymax=max(y);  %set axis limits
            xmin=min(x);    xmax=max(x);
            axis([xmin, xmax, ymin, ymax]);
            
%             self(1,1).Sol.makeTicks(xmin, xmax, floor((xmax-xmin)/5), floor((xmax-xmin)/20),...
%                 ymin, ymax, floor((ymax-ymin)/5), floor((ymax-ymin)/20));
            
            ax = gca; ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
            h = legend(vlabels);
                set(h,'Location', 'Best', 'Orientation', ...
                    'vertical', 'Fontname', 'Arial','FontSize', 25,...
                    'EdgeColor',[1 1 1]);
                    %,'YColor',[1 1 1],'XColor',[1 1 1]
            for i = 1:N
                y = expfrac{i};
                %plot(x,y,'-', 'Color', colors(i,:),'LineWidth',3.5)
                %Id = self(i,1).Sol.Vapp.ID; yForced = self(i,1).Sol.Vapp.forced;
                %plot(x(1:Id),yForced,'--','LineWidth',3.5,'Color', colors(i,:))
                plot(x,y,'-','LineWidth',3.5,'Color', colors(i,:))
                uistack(A, 'top')
            end
            y = zeros(length(x),1);
            plot(x,y,'-','LineWidth',3.5,'Color', 'k')
            
      end
                  
      function printZeroPoint(coldArray, warmArray)
            % prints a table in the command window that displays [name],
            % temperature cooling range and zero estimated contraction
            %concentration
            %           
            % & & Calculated Isochoric Point & &               \\
            % \hline
            %  cryoprotectant     & entry      & entry  [...]  \\
            %  concentration      & entry      & entry  [...]  \\

            len = length(coldArray);
            if len ~= length(warmArray)
                error('Warm/Cold array lengths not equal')
            end
            fprintf('\\hline \n')
            fprintf('Cryoprotectant & Conc. (\\%% w/w) \\\\\n')
            fprintf('\\hline \n')
            for i = 1:len
                fprintf('%s & %.1f \\\\\n',coldArray(i,1).Sol.name,findZero(coldArray(i,1).Sol, warmArray(i,1).Sol))
            end

      end
      
      function printContractionTable(self, other)
            % prints a table in the command window that displays
            % contraction upon cooling at -4%, -3% -2%, -1%, 0%, and 1%
            % coldObj is a cold solution set, and warmObj is a warmer
            % solution set from that same cryoprotectant agent
            % self is a cold array, and other is a warm array.
            %
            %  name       & T0/Tf (K)  &  2%        & 1%         & 0%         & -1%         & -2%        & -3%         & -4%         \\
            % \hline     
            %  entry      & entry      & entry      & entry      & entry      & entry       & entry      & entry       & entry       \\
            %  entry      & entry      & entry      & entry      & entry      & entry       & entry      & entry       & entry       \\
            %  entry      & entry      & entry      & entry      & entry      & entryy      & entry      & entry       & entry       \\
            
            
            
            % header
            fprintf('&  $T_0,T_f (K)$ & 2\\%% & 1\\%% & 0\\%% & -1\\%% & -2\\%% & -3\\%% & -4\\%% \\\\\n')
            fprintf('\\hline\n')
            
            %for each cryoprotectant
            for n = 1:length(self)
                %data prep
                name = self(n,1).Sol.name;
                Tf = floor(self(n,1).Sol.degtemp);
                T0 = floor(other(n,1).Sol.degtemp);
                [~, expfrac]= calcContraction(self(n,1).Sol, other(n,1).Sol);
                len  = length(expfrac);
                target = 2;
                c = [];
                for i = 1:len        
                    if expfrac(i)<=target && target>=-4
                       a = expfrac(i); b = expfrac(i-1);
                       adj = (a-target)/(a-b);
                       c = [c, i+adj];
                       target = target -1;
                    end
                end
                if length(c) ~= 7       
                    error('Wrong number of concentractions extracted.')    
                end
                %print data
                fprintf('%s  \t & %.0f - %.0f \t &', name, T0, Tf);
                fprintf('%.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\\n', c(1), c(2), c(3), c(4), c(5), c(6), c(7));
            end

            
            


      end
      
      function printParameterTable(self)
         %needs updated to nonlinear model used
          % prints a table in the command window that displays
            % contraction upon cooling at -4%, -3% -2%, -1%, 0%, and 1%
            % coldObj is a cold solution set, and warmObj is a warmer
            % solution set from that same cryoprotectant agent
            % self is a cold array, and other is a warm array. 
            % a_i is scaled by 10^{i*2} so the dependent variable is weight
            % fraction.
            %
            %         & $a_0$  &  $a_1$        & $a_2$         & $a_3$         & $a_4$      & RMSE     & $R^2$      \\
            % \hline     
            %  name      & entry      & entry      & entry      & entry      & entry       & entry  & entry          \\
            %  name      & entry      & entry      & entry      & entry      & entry       & entry  & entry           \\
            %  name      & entry      & entry      & entry      & entry      & entryy      & entry  & entry           \\
            
            
            
            % header
            fprintf('& $a_0$  &  $a_1$        & $a_2$         & $a_3$         & $a_4$     & RMSE    & $R^2$       \\\\\n')
            fprintf('\\hline\n')
            
            %for each cryoprotectant
            for n = 1:length(self)
                %data prep
                name = self(n,1).Sol.name;
                [vec, rsq, rmse] = self(n).Sol.getParameters();
                %a4 = vec(1); a3 = 10^(0)*vec(2); a2 = 10^(0)*vec(3); a1 =
                %10^(0)*vec(4); a0 = vec(5); %unscaled values
                a4 = 10^(8)*vec{1}; a3 = 10^(6)*vec{2}; a2 = 10^(4)*vec{3}; a1 = 10^(2)*vec{4}; a0 = vec{5};
                fprintf('%s &  %.3e  &  %.3f  & %.3f & %.3f  & %.3f   & %.3f  & %.3f  \\\\\n', name,a0, a1, a2, a3, a4, 1000*rmse, rsq)
            end
            fprintf('\n')
      end      
      
      function printParameters(self)
          %prints parameters
          for n = 1:length(self)
                %data prep
                name = self(n,1).Sol.name;
                fprintf(name)
                self(n,1).Sol.modeledVE.fit{1}
                fprintf('\n\n')
          end
            fprintf('\n')
          
      end
      
      function printTables(self,other)
          %applies the printTable funciton to all elements of the Solution
          %Array      
          for i  = 1:length(self)
              self(i,1).Sol.printTable(other(i,1).Sol)
          end
      end
      
      function printError(self, other)
          %prints the average values and specific values for the
          %uncertainties output from printUncertainties() in SolutionSet.
          % self is a cold Array and other is the corresponding warm array.
          l = length(self);
          sum = [0,0, 0, 0];
          for i=1:l
              [Ddens,Pdens, DVE, Dc] =  printUncertainty(self(i,1).Sol(1),other(i,1).Sol(1));
              sum = sum + [Ddens,Pdens, DVE, Dc];
          end
          avg = sum./(l);
            fprintf('\n and the mean of these are:');
            fprintf('average uncertainty (g/ml) in density is %f \n', avg(1));
            fprintf('average percent uncertainty in density  is %f \n', avg(2));
            fprintf('average uncertainty (ml/g) in VE  is %f \n', avg(3));
            fprintf('average uncertainty (percent contraction) in contraction is %f .\n', avg(4));
      end
      
      function printMeasUncertContract(self, other)
                      
          
          N = length(self);
            x = self(1,1).Sol.concInterpolated;
            l = 0; s = 0;
            % build label vector
            vlabels = [];
            for i = 1:N               
                vlabels{i} = self(i,1).Sol.label;
                [volfrac{i}, expfrac{i}] = calcContraction(self(i,1).Sol, other(i,1).Sol);
                
                if strcmp(self(i,1).Sol.temp,'cold')
                    up{i} = abs(calcContraction(self(i,1).Sol, other(i,1).Sol,'up') - volfrac{i}); up{i}(1)=0;%volfracG(1);
                    down{i} = abs(calcContraction(self(i,1).Sol, other(i,1).Sol,'down') - volfrac{i}); down{i}(1)=0;%volfracG(1);
                    
                    s = s + sum((up{i}-down{i})/2);
                    l = l + length(self(i,1).Sol);
                end
            end
          
            fprintf('Mean measurement uncertainty in contraction is %f', s./l);
      end
      
   end 
end