clear all;
clc;close all;
% colors={'b','r','c','m','g','y',[0 0.7 0],'k'};
colors=get(groot,'DefaultAxesColorOrder');
markers={'s','o','v','>','*','s','o'};
set(0, 'DefaultAxesFontSize',22)
set(0,'defaultlinelinewidth',1)
set(0,'defaultlinemarkersize',12)

% === Options === %
% ofDir='/home/thijsgillebaart/foam/thijsgillebaart-3.1/';
ofDir='/media/thijsgillebaart/THIJSDATA/foam/thijsgillebaart-3.1/';
% ofDir='/Users/thijsgillebaart/foam/thijsgillebaart-3.1';
% ofDir='/Volumes/THIJSDATA/foam/thijsgillebaart-3.1/';

% ==== Rigid body ==== %
% mainDir='/run/papers/consistentFsiRBF/consistent/rigidBody/';
% cases={'3DoF-RK4-1_forcing','3DoF-RK4-2_forcing','3DoF-RK4-4_forcing','3DoF-BDF2_forcing','3DoF-BDF3_forcing'};
% schemes={'bdf3','bdf3','bdf3','bdf3','bdf3'};
% legendNames={'3DoF-RK4-1','3DoF-RK4-2','3DoF-RK4-4','3DoF-BDF2','3DoF-BDF3'};

% ==== Circle In Circle ==== %
mainDir='/run/papers/consistentFsiRBF/consistent/circleInCircle/U10/';

% cases={'fsiBodyMesh4_RK4-1','fsiBodyMesh4_RK4-2','fsiBodyMesh4_RK4-4','fsiBodyMesh4_BDF3'};
% cases={'fsiBodyMesh_RK4-1','fsiBodyMesh_RK4-2','fsiBodyMesh_RK4-4','fsiBodyMesh_BDF3'};
% schemes={'bdf2','bdf2','bdf2','bdf3'};
% legendNames={'RK4 F(t^{n+1})','RK4 F(t^{n+1/2})','RK4 F(t^{n}+k\Delta t)','BDF3'};

% cases={'fsiBodyMesh_RK4-4','fsiBodyMesh_RK4-4','fsiBodyMesh_BDF2','fsiBodyMesh_RK4-4','fsiBodyMesh_BDF3'};
% schemes={'bdf1','bdf2','bdf2','bdf3','bdf3'};
% legendNames={'RK4 BDF1','RK4 BDF2','BDF2','RK4 BDF3','BDF3'};

% cases={'fsiBodyMesh_RK4-1','fsiBodyMesh_RK4-2','fsiBodyMesh_RK4-4','fsiBodyMesh_BDF2','fsiBodyMesh_BDF3'};
% schemes={'bdf3','bdf3','bdf3','bdf3','bdf3'};
% legendNames={'RK4 F(t^{n+1})','RK4 F(t^{n+1/2})','RK4 F(t^{n}+k\Deltat)','BDF2','BDF3'};

% cases={'fsiBodyMesh_BDF1','fsiBodyMesh_BDF2','fsiBodyMesh_BDF3'};
% schemes={'bdf1','bdf2','bdf3'};
% legendNames={'bdf 1','bdf 2','bdf 3'};

% cases={'fsiBody3Mesh_RK4-4','fsiBody3Mesh_RK4-4','fsiBody3Mesh_RK4-4'};
% schemes={'bdf1','bdf2','bdf3'};
% legendNames={'bdf 1','bdf 2','bdf 3'};

% cases={'fsiBody3Mesh_BDF3','fsiBody3Mesh_BDF3','fsiBody3Mesh_BDF3'};
% schemes={'bdf1','bdf2','bdf3'};
% legendNames={'bdf 1','bdf 2','bdf 3'};

cases={'fsiBody3Mesh4_RK4-1','fsiBody3Mesh4_RK4-2','fsiBody3Mesh4_RK4-4','fsiBody3Mesh4_BDF2','fsiBody3Mesh4_BDF3'};
schemes={'bdf3','bdf3','bdf3','bdf3','bdf3'};
legendNames={'RK4 F(t^{n+1})','RK4 F(t^{n+1/2})','RK4 F(t^{n}+k\Deltat)','bdf2','bdf3'};

%% Options
dataName='body-block-state.dat';

orderDim = 4;
normN=2;

lastSchemeRef=0;
orderlines=[1,2,3];
writeFigure=1;
xlimits=[3e-4,5e-2];

figureName=strrep(['ORDER_fsiBodyScheme_theta_consistent'],'/','');
%% Read data
nCases=length(cases);
for i=1:nCases
    runDir{i}=[ofDir,mainDir,cases{i},'/',schemes{i}];
end

% data to be saved
Op=cell(nCases,1);
Ep=cell(nCases,1);
dataSol=cell(nCases,1);
dts=cell(nCases,1);
for q=1:nCases
    disp(['|| === ',cases{q},' === ||']);
    % === Set time steps === %
    files=dir(runDir{q});%Get all dirs in mainDir
    dt=[];counter=0;
    for i=1:length(files)
        timeStringIndexStart = strfind(files(i).name,'dt') + length('dt');
        if(~isempty(timeStringIndexStart) && timeStringIndexStart >= 0 && (files(i).isdir == 1))
            counter=counter+1;
            tmptime=str2double(files(i).name(timeStringIndexStart:end));
            dt(counter)=tmptime;
            dtNames{counter}=num2str(dt(counter));
            dirs{counter}=[runDir{q},'/',files(i).name];
        end
    end
    % sort time steps from big to small
    [dt,a]=sort(dt,'descend');
    dtNames=dtNames(a);
    dirs=dirs(a);
    nDts=length(dt);
    dts{q}=dt;

    % === Set times === %
    testData=dlmread([dirs{1},'/',dataName],'\t',1,0);
    endTime=testData(end,1);
    times=testData(:,1);
    nTimes=length(times);
  
    % === Read in all data === %
    dataSol{q}=cell(nDts,1);
    for i=1:nDts
        data=dlmread([dirs{i},'/',dataName],'\t',1,0);
        dataSol{q}{i}=data(:,orderDim);
        clear data;
    end 
end

% === Error === %
for q=1:nCases  
    nDts=length(dts{q});
    % === Calculate temporal and model error === %
    Op{q}=zeros(1,nDts-2);
    Ep{q}=zeros(1,nDts-1);
    for i=1:nDts-1
        steps=1;
        
        %Takes last scheme last time step as reference solution
        if(lastSchemeRef)
            stepsRef=dts{q}(i)/dts{end}(end);
            Ep{q}(i)=norm(dataSol{q}{i}(steps:steps:end)-dataSol{end}{end}(stepsRef:stepsRef:end),normN)/norm(dataSol{end}{end}(stepsRef:stepsRef:end),normN);
        else
            stepsRef=dts{q}(i)/dts{q}(end);
            Ep{q}(i)=norm(dataSol{q}{i}(steps:steps:end)-dataSol{q}{end}(stepsRef:stepsRef:end),normN)/norm(dataSol{q}{end}(stepsRef:stepsRef:end),normN);
        end
    end

    % === Error&Order based on vector with ALL time data ===%
    for i=1:nDts-2
        Op{q}(i)=log(Ep{q}(i+1)/Ep{q}(i))/log(dt(i+1)/dt(i));
    end
    disp(['Order = ',num2str(Op{q})]);
end

%% Plotting

ylabelText='|\epsilon_{\theta}|_{2}';

% order lines
Epv=[];
EpvEnd=[];
EpvEnd1=[];
for i=1:length(Ep)
    Epv=[Epv,Ep{i}(1)];
    Epv=sort(Epv,'descend');
    EpvEnd=[EpvEnd,Ep{i}(end)];
    EpvEnd=sort(EpvEnd,'descend');
    EpvEnd1=[EpvEnd1,Ep{i}(end-1)];
    EpvEnd1=sort(EpvEnd1,'descend');
end

o1ErrorStart = 4/3*Epv(1);
o2ErrorStart = sqrt(Epv(1)*Epv(end));
o3ErrorStart = 2/3*Epv(end);
o1Error=4/3*EpvEnd1(1);
o2Error=4/3*EpvEnd1(end-1);
o3Error=3/4*min(EpvEnd);
thirdOrderLine=[o3Error*((dt(end-2)/dt(end-1))^3),o3Error];
secondOrderLine=[o2Error,o2Error/((dt(end-2)/dt(end-1))^2)];
firstOrderLine=[o1Error,o1Error/((dt(end-2)/dt(end-1))^1)];
o1Line=[o1ErrorStart,o1ErrorStart/((dt(1)/dt(end-1))^1)];
o2Line=[o2ErrorStart,o2ErrorStart/((dt(1)/dt(end-1))^2)];
o3Line=[o3ErrorStart,o3ErrorStart/((dt(1)/dt(end-1))^3)];
oLineDt=[dt(1),dt(end-1)];

%%
h=figure('Position',[0,0,1200,900]);
loglog(dts{1}(1:end-1),Ep{1},'linestyle','-','color',colors(1,:),'marker',markers{1});
hold on;
for i=2:length(runDir)
    loglog(dts{i}(1:end-1),Ep{i},'linestyle','-','color',colors(i,:),'marker',markers{i});
end
grid on;
xlabel('\Delta t')
ylabel(ylabelText)
xlim(xlimits)

hl=legend(legendNames,'Location','NorthOutside','Orientation','Horizontal');
hc = get(hl, 'children');
% for i=1:nCases
%     hcpos=get(hc(i*3),'position');
%     set(hc(i*3), 'position', [hcpos(1),hcpos(2)*0.7,hcpos(3)]) % Centering the third    
% end

% % order triangles/lines
% if(maxOrder>=1)
%     loglog([dt(end-2),dt(end-1)],firstOrderLine,'-k');
%     loglog([dt(end-2),dt(end-1)],[firstOrderLine(1),firstOrderLine(1)],'-k');
%     loglog([dt(end-1),dt(end-1)],[firstOrderLine(1),firstOrderLine(2)],'-k');
%     t1=text((dt(end-2)+dt(end-1))/2.5,sum(firstOrderLine)/1.9,'1','Fontsize',22);
% end
% if(maxOrder>=2)
%     loglog([dt(end-2),dt(end-1)],secondOrderLine,'-k');
%     loglog([dt(end-2),dt(end-1)],[secondOrderLine(1),secondOrderLine(1)],'-k');
%     loglog([dt(end-1),dt(end-1)],[secondOrderLine(1),secondOrderLine(2)],'-k');
%     t2=text((dt(end-2)+dt(end-1))/2.5,sum(secondOrderLine)/1.9,'2','Fontsize',22);
% end
% if(maxOrder>=3)
%     loglog([dt(end-2),dt(end-1)],thirdOrderLine,'-k');
%     loglog([dt(end-2),dt(end-1)],[thirdOrderLine(2),thirdOrderLine(2)],'-k');
%     loglog([dt(end-2),dt(end-2)],[thirdOrderLine(1),thirdOrderLine(2)],'-k');
%     t3=text((dt(end-2)+dt(end-1))/2,sum(thirdOrderLine)/4,'3','Fontsize',22);
% end
if(~isempty(find(orderlines==1, 1)))
    l1=loglog(oLineDt,o1Line,'-k','linewidth',1);
    uistack(l1,'bottom')
    t1=text(sqrt(dt(end-1)*dt(1)),sqrt(o1Line(1)*o1Line(2))*4/3,'1','Fontsize',16);
end
if(~isempty(find(orderlines==2, 1)))
    l2=loglog(oLineDt,o2Line,'-k','linewidth',1);
    uistack(l2,'bottom')
    t2=text(sqrt(dt(end-1)*dt(1)),sqrt(o2Line(1)*o2Line(2))*4/3,'2','Fontsize',16);
end
if(~isempty(find(orderlines==3, 1)))
    l3=loglog(oLineDt,o3Line,'-k','linewidth',1);
    uistack(l3,'bottom')
    t3=text(sqrt(dt(end-1)*dt(1)),sqrt(o3Line(1)*o3Line(2))*2/3,'3','Fontsize',16);
end
set(gca,'Xdir','reverse');

if(writeFigure)
    hgexport(gcf,['figures/',figureName], hgexport('factorystyle'), 'Format','png');
    hgexport(gcf,['figures/',figureName], hgexport('factorystyle'), 'Format','eps');
end