clear all;clc;
close all;
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

% ==== Circle In Circle ==== %
mainDir='/run/papers/consistentFsiRBF/consistent/circleInCircle/U10/';

%% All runs

% ==== bdf 1, bdf 2, bdf 3 ==== %
% nomovingMesh4
% cases={'nomovingMesh4','nomovingMesh4','nomovingMesh4'};
% schemes={'bdf1','bdf2','bdf3'};
% legendNames={'bdf1','bdf2','bdf3'};

% ExpandingRotatingMesh
% cases={'expandingRotatingMesh4_20','expandingRotatingMesh4_20','expandingRotatingMesh4_20'};
% schemes={'bdf1','bdf2','bdf3'};
% legendNames={'bdf1','bdf2','bdf3'};

% movingBodyMesh
% cases={'movingBodyMesh4','movingBodyMesh4','movingBodyMesh4'};
% schemes={'bdf1','bdf2','bdf3'};
% legendNames={'bdf1','bdf2','bdf3'};

% FSI Order standard
% cases={'fsiBody3Mesh4_BDF1','fsiBody3Mesh4_RK4-4','fsiBody3Mesh4_BDF2','fsiBody3Mesh4_BDF3'};
% schemes={'bdf1','bdf2','bdf2','bdf3'};
% legendNames={'bdf1','RK4 F(t^{n}+k\Deltat)','bdf2','bdf3'};

% ==== RK FSI ==== %
% rotatingMesh
% cases={'rotatingMesh4_20','rotatingMesh4_20','rotatingMesh4_20'};
% schemes={'bdf2nor','bdf2area','bdf2'};
% legendNames={'bdf2_{nor}','bdf2_{area}','bdf2'};

% expandingMesh
% cases={'expandingMesh4','expandingMesh4','expandingMesh4'};
% schemes={'bdf2nor','bdf2area','bdf2'};
% legendNames={'bdf2_{nor}','bdf2_{area}','bdf2'};

% RK interpolation for rigid body
% cases={'fsiBody3Mesh4_RK4-1','fsiBody3Mesh4_RK4-2','fsiBody3Mesh4_RK4-4','fsiBody3Mesh4_BDF2'};
% schemes={'bdf2','bdf2','bdf2','bdf2'};
% legendNames={'RK4 F(t^{n+1})','RK4 F(t^{n+1/2})','RK4 F(t^{n}+k\Deltat)','BDF2'};

% ==== Relaxation & Unstructured ==== %
% cases={'fsiBody3Mesh4_BDF1_rel0.9-n','fsiBody3Mesh4_BDF2_rel0.9-n','fsiBody3Mesh4_BDF3_rel0.9-n','fsiBody3Mesh4_BDF3'};
% schemes={'bdf1','bdf2','bdf3','bdf3'};
% legendNames={'bdf1 relax 0.9','bdf2 relax 0.9','bdf3 relax 0.9','bdf3'};

cases={'fsiBody3MeshUn4_BDF3_rel0.99-n2','fsiBody3Mesh4_BDF3_rel0.9-n2'};
schemes={'bdf3','bdf3'};
legendNames={'bdf3','bdf3 unstructured'};

%% Extra
% cases={'nomovingMesh_rel0.9-n','nomovingMesh_rel0.9-n','nomovingMesh_rel0.9-n','nomovingMesh'};
% schemes={'bdf1','bdf2','bdf3','bdf3'};
% legendNames={'bdf 1','bdf 2','bdf 3','bdf 3'};

% cases={'nomovingMesh_rel0.9','nomovingMesh_rel0.9','nomovingMesh_rel0.9','nomovingMesh'};
% schemes={'bdf1','bdf2','bdf3','bdf3'};
% legendNames={'bdf 1','bdf 2','bdf 3','bdf 3'};

% cases={'fsiBodyMeshUn4_BDF3','fsiBodyMeshUn4_BDF3','fsiBodyMeshUn4_BDF3'};
% cases={'fsiBodyMesh4_BDF3_rel0.9-n','fsiBodyMesh4_BDF3_rel0.9-n','fsiBodyMesh4_BDF3_rel0.9-n'};
% schemes={'bdf1','bdf2','bdf3'};
% legendNames={'bdf 1','bdf 2','bdf 3'};

% cases={'fsiBody3Mesh_BDF1','fsiBody3Mesh_BDF2','fsiBody3Mesh_RK4-4','fsiBody3Mesh_BDF3'};
% schemes={'bdf1','bdf2','bdf2','bdf3'};
% legendNames={'bdf^1','bdf^2','RK4 F(t^{n}+k\Deltat)','bdf^3'};

% cases={'fsiBody3Mesh4_BDF3_rel0.9-n'};
% schemes={'bdf3'};
% legendNames={'bdf^3_{rel0.9}'};

%% Options
% Set characteristics
setName='cellSet_convergenceSet';
lastSchemeRef=0;

velocity=1;
pressure=1;
normN=Inf;
writeFigure=1;
orderlines=[2,3];

xlimits=[3e-4,5e-2];

%% Show settings
disp(['setName = ',setName]);
disp(['last solution as reference = ',num2str(lastSchemeRef)]);
disp(['order for pressure: ',num2str(pressure),', velocity: ',num2str(velocity)]);
disp(['norm of error is ',num2str(normN)]);

%% Creating variables
nCases=length(cases);
for i=1:nCases
    runDir{i}=[ofDir,mainDir,cases{i},'/',schemes{i}];
end

figureName=strrep(['ORDER_fsiBodyScheme_pU_',cases{1}],'/','');
%% Read data & Calculate Errors

% data to be saved
OpP=cell(nCases,1);
OpU=cell(nCases,1);
EpP=cell(nCases,1);
EpU=cell(nCases,1);
dataSolP=cell(nCases,1);
dataSolU=cell(nCases,1);
dts=cell(nCases,1);
for q=1:nCases
    disp(['|| === Reading ',cases{q},': ',schemes{q},' === ||']);
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
    files=dir([dirs{1},'/sets/']);
    times=[];counter=0;
    for i=1:length(files)
        tmptime=sscanf(files(i).name,['%f']);
        if(~isempty(tmptime) && length(tmptime)==1)
            counter=counter+1;
            times(counter)=tmptime;
        end
    end
    times=sort(times,'ascend');
    time=times(end);
    
    %get lenght of data
    nDatas=size(dlmread([dirs{1},'/sets/',num2str(time),'/',setName,'_p.xy'],'\t',1,0),1);
    
    % === Read in all data === %
    dataSolP{q}=zeros(nDatas,nDts);
    dataSolU{q}=zeros(nDatas,nDts);
    for i=1:nDts
        data=dlmread([dirs{i},'/sets/',num2str(time),'/',setName,'_p.xy'],'\t',1,0);
        dataSolP{q}(:,i)=data(:,1);
        data=dlmread([dirs{i},'/sets/',num2str(time),'/',setName,'_U.xy'],'\t',1,0);
        dataSolU{q}(:,i)=sqrt(data(:,1).^2 + data(:,2).^2);
    end
    clear data;
end

% === Error === %
for q=1:nCases
    disp(['|| === Order ',cases{q},': ',schemes{q},' === ||']);
    nDts=length(dts{q});
    % === Calculate temporal and model error === %
    OpP{q}=zeros(1,nDts-2);
    OpU{q}=zeros(1,nDts-2);
    EpP{q}=zeros(1,nDts-1);
    EpU{q}=zeros(1,nDts-1);
    for i=1:nDts-1
        %Takes last scheme last time step as reference solution
        if(lastSchemeRef)
            EpP{q}(i)=norm(dataSolP{q}(:,i)-dataSolP{end}(:,end),normN)/norm(dataSolP{end}(:,end),normN);
            EpU{q}(i)=norm(dataSolU{q}(:,i)-dataSolU{end}(:,end),normN)/norm(dataSolU{end}(:,end),normN);
        else
            EpP{q}(i)=norm(dataSolP{q}(:,i)-dataSolP{q}(:,end),normN)/norm(dataSolP{q}(:,end),normN);
            EpU{q}(i)=norm(dataSolU{q}(:,i)-dataSolU{q}(:,end),normN)/norm(dataSolU{q}(:,end),normN);
        end
    end

    % === Error&Order based on vector with ALL time data ===%
    for i=1:nDts-2
        OpP{q}(i)=log(EpP{q}(i+1)/EpP{q}(i))/log(dt(i+1)/dt(i));
        OpU{q}(i)=log(EpU{q}(i+1)/EpU{q}(i))/log(dt(i+1)/dt(i));
    end
    disp(['Order pressure = ',num2str(OpP{q})]);
    disp(['Order velocity = ',num2str(OpU{q})]);
end

%% Prepare Plotting
ylabelText='|\epsilon|_{\infty}';

% order lines
Epv=[];
EpvEnd=[];
EpvEnd1=[];
for i=1:length(EpP)
    Epv=[Epv,EpP{i}(1),EpU{i}(1)];
    Epv=sort(Epv,'descend');
    EpvEnd=[EpvEnd,EpP{i}(end),EpU{i}(end)];
    EpvEnd=sort(EpvEnd,'descend');
    EpvEnd1=[EpvEnd1,EpP{i}(end-1),EpU{i}(end-1)];
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


%% Plotting
h=figure('Position',[0,0,1200,900]);
if(pressure)
    loglog(dts{1}(1:end-1),EpP{1},'linestyle','-','color',colors(1,:),'marker',markers{1});
    hold on;
    for i=2:length(runDir)
        loglog(dts{i}(1:end-1),EpP{i},'linestyle','-','color',colors(i,:),'marker',markers{i});
    end
end
if(velocity)
    loglog(dts{1}(1:end-1),EpU{1},'linestyle','--','color',colors(1,:),'marker',markers{1});
    hold on;
    for i=2:length(runDir)
        loglog(dts{i}(1:end-1),EpU{i},'linestyle','--','color',colors(i,:),'marker',markers{i});
    end
end
grid on;
legend(legendNames,'Location','NorthOutside','Orientation','Horizontal');
xlabel('\Delta t')
ylabel(ylabelText)
xlim(xlimits)

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
%     hgexport(gcf,['figures/',figureName], hgexport('factorystyle'), 'Format','png');
    hgexport(gcf,['figures/',figureName], hgexport('factorystyle'), 'Format','eps');
end