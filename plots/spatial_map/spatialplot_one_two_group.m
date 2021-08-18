%% Spatial map of Kenya - which counties have one-group and which two-group

S = shaperead('./shapefiles/ken_adm_iebc_20191031_shp/ken_admbnda_adm1_iebc_20191031.shp');
[~,index] = sortrows({S.ADM1_EN}.'); S = S(index); clear index;
for i = 1:47
      S(i).ID = i;
end
%%
T = readtable('../../forecasts/parameter_posterior_means.csv');
assortative_county = T.epsilon > 0.8;
upper_SES = T.P_eff < 0.85;

%%
n = 10001;
% map =parula(n);
map = cool(n);


%%
% greyshade = [0.7 0.7 0.7];

I = assortative_county(1).*upper_SES(1).*(n-1) + 1;

ColoredConstituencies = makesymbolspec('Polygon',{ 'ID',1,'FaceColor',map(I,:)} );
    
    
 for i = 2:47
    I =  assortative_county(i).*upper_SES(i).*(n-1) + 1;
    ColoredConstituencies.FaceColor{i,1} = 'ID';
    ColoredConstituencies.FaceColor{i,2} = i;
    ColoredConstituencies.FaceColor{i,3} = map(I,:);
 end
 

%  
 %%
colormap(map);
figure(1);
clf;
%Create three axes
h_kenya = axes;
h_kenya.Position =  [0.100 0.100 0.55 0.85];
h_kenya.FontSize = 18;
h_kenya.XAxis.Visible = 'off';
h_kenya.YAxis.Visible = 'off';


h_Nairobi = axes;
h_Nairobi.Position =  [0.7 0.600 0.25 0.35];
h_Nairobi.FontSize = 18;
h_Nairobi.XAxis.Visible = 'off';
h_Nairobi.YAxis.Visible = 'off';

h_mombasa = axes;
h_mombasa.Position =  [0.7 0.100 0.25 0.35];
h_mombasa.FontSize = 18;
h_mombasa.XAxis.Visible = 'off';
h_mombasa.YAxis.Visible = 'off';

axes(h_kenya);
% mapshow(S,'LineStyle','none','SymbolSpec',ColoredConstituencies);
% h = mapshow(S,'SymbolSpec',ColoredConstituencies,'LineStyle','none');
title('Kenyan counties')
h = mapshow(S,'SymbolSpec',ColoredConstituencies);



% cb = colorbar;
% cb.Location = 'westoutside';
% cb.Ticks = (0:0.1:0.5)/max_value;
% 
% cb.TickLabels = {'0%','10%','20%','30%','40%','50%'} ;
% cb.Label.String = 'Percentage infected';
% cb.FontSize = 28;



axes(h_Nairobi);
    mapshow(S(30),'SymbolSpec',ColoredConstituencies);
    title('Nairobi')
    
   axes(h_mombasa);
    mapshow(S(28),'SymbolSpec',ColoredConstituencies);
    title('Mombasa')  

%%

%  dim = [.496 .227 .025 .025];
annotation('ellipse',[0.49365625,0.227008149010477,0.0211875,0.026155995343423]);
annotation('ellipse',[0.356345911949686,0.415263748597082,0.032018867924528,0.032547699214368]);


%%

annotation('arrow',[0.5140625,0.73046875],[0.238813736903376,0.2782305005]);
annotation('arrow',[0.38671875,0.7203125],[0.437882421420256,0.707799767171129]);  

