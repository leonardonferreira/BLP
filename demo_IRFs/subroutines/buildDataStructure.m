
function modelSpec=buildDataStructure(dataSetSpec,modelSpec)
%packs together data from different sources


%-main endogenous---------------------------------------------------------%

dataSetSpec.sourceFile         ='FRED-MD_2015SM';
% available sourceFile
% FRED-MD_2015m1: 134 monthly series McCracken & Ng(2014)


%data transformations
dataSetSpec.plotData           =false; %plots individual(!) charts 
dataSetSpec.interpolateMissing =false;


%load base set using St. Louis FRED data
dataStructure=loadFREDmonthly(dataSetSpec);


%-policy variable---------------------------------------------------------%

%replace average rates w\ end of month for psvar identification
if ismember('GS1',dataSetSpec.dataList) && strcmp(modelSpec.identification,'PSVAR')

    load endOfMonthGS1data
    
    isOneYearRate =ismember(dataStructure.varname,'GS1');
    oneyearDates  =datenum(year(eomGS1.dates),month(eomGS1.dates),1);
    
    %add to main data
    dataStructure.data(:,isOneYearRate)=eomGS1.data(ismember(oneyearDates,dataStructure.dates),:);
    
    %add to presample
    if isfield(dataStructure,'preSdates')
        dataStructure.preSdata(:,isOneYearRate)=eomGS1.data(ismember(oneyearDates,dataStructure.preSdates),:);
    end
end


%-add expectations to endogenous set--------------------------------------%

if isfield(dataSetSpec,'eDataList') && ~isempty(dataSetSpec.eDataList)
    
    eDataLabel={'mean','median','dispersion','skewness','iqr'};
    eDataShortLabel={'M','Md','D','S','I'};
    
    
    %load raw data
    load ConsensusForecastsUS
    
    %find variables
    eDataList=dataSetSpec.eDataList;
    [~,eDataSelect]=ismember(eDataList,CFUS1Y.labels);
    
    %define common time bounds to align data and expectations
    consensusDates    =datenum(year(CFUS1Y.dates),month(CFUS1Y.dates),1);

    lowerT=max(dataStructure.dates(1),consensusDates(1));
    upperT=min(dataStructure.dates(end),consensusDates(end-1));
    
    commonTimeLine=lowerT;
    while commonTimeLine(end)<upperT
        commonTimeLine=[commonTimeLine; addtodate(commonTimeLine(end),1,'month')];
    end
    
    %find relevant expectation series
    if length(dataSetSpec.eDataType)==1 %only one moment of distribution
        
        %1Year ahead forecasts approximated from calendar year forecasts
        eval(['eDataTemp=CFUS1Y.' dataSetSpec.eDataType{1} ';'])
        
        eData=eDataTemp(:,eDataSelect);
        
    else
    
        %more than one moment
        eData=NaN(length(consensusDates),length(dataSetSpec.eDataType));
        for j=1:length(dataSetSpec.eDataType)

            eDataType=dataSetSpec.eDataType{j};

            %1Year ahead forecasts approximated from calendar year forecasts
            eval(['eDataTemp=CFUS1Y.' eDataType ';'])

            eData(:,j)=eDataTemp(:,eDataSelect);
        end
        
    end
        
        
        
    %trim unmatching data points in original set
    dataStructure.data(~ismember(dataStructure.dates,commonTimeLine),:)=[];

    %add expectations (shift one month back to align with information set)
    dataStructure.data=[dataStructure.data eData(find(ismember(consensusDates,commonTimeLine))+1,:)];

    
    %reset timeline
    dataStructure.dates=commonTimeLine;
    
    %add labels
    cflab=eDataShortLabel(ismember(eDataLabel,dataSetSpec.eDataType));
    dataStructure.varname=[dataStructure.varname strcat(strcat('Cf',cflab),eDataList)];
    
    %add varname
    tempName=dataSetSpec.eDataType; tempName=upper(tempName);
    for name=tempName'
        
        dataStructure.varLongName=[dataStructure.varLongName strcat(strcat({'Consensus '},{name{1}},{' '}),eDataList)];
    end
    
    %add expectations to presample
    if modelSpec.presample
    
        dataStructure.preSdata=[dataStructure.preSdata eData(find(ismember(consensusDates,dataStructure.preSdates))+1,:)];
    end

        
end



%-load external proxies for shock identification--------------------------%

%load external instruments if PSVAR identification
if strcmp(modelSpec.identification,'PSVAR')
    
    
    if strfind(modelSpec.selectedInstrument,'MPN')
        
        %load narrative Romer&Romer series (updated)
        load RRnarrative    
        
    elseif strfind(modelSpec.selectedInstrument,'FF4')
        
        %load FF4-based instruments
        load FF4_IV_variants
        
    elseif strfind(modelSpec.selectedInstrument,'MM_IV')
        
        %load instruments corrected for Fed info
        load MM_IVinstruments
        
    end
    
    %add instrument details to model specification
    modelSpec.instrument =externalInstrument;

end        
        

%load all into final structure
modelSpec.dataStructure=dataStructure;








%children functions
%-------------------------------------------------------------------------%

function res=loadFREDmonthly(dataSetSpec)
%
% loads up raw data and applies required transformations
% 
% inputs (structure)
% sourceFile         = mat file name (see content below)
% seriesList         = if useStoredList is false, allows to define a new 
%                      one, comes in the form of a cell string {'';'';''}
% beginSet           = date
% endSet             = date
% takeFirstDiff      = true/false, applies to all
% plotData           = true/false if true plots original and transformed
% interpolateMissing = true/false, fills up missing values with centered MA
%                      does not apply to beginning and end of series
%
% output (structure)
% data        = [Txn] matrix of transformed data
% dates       = [Tx1] vector of reference dates
% varname     = {1xn} cell of data identifiers
% varLongName = {1xn} cell of variables names
% 
% miranda 2015 smirandaagrippino@london.edu


%load raw data from file
load(dataSetSpec.sourceFile);

%sourceFile is in mat format; contains:
% data            = [TxN] matrix of raw data
% dates           = [Tx1] vector of dates
% logTransform    = [1xN] logical for log transformations
% dataName        = {1xN} cell of data identifiers
% dataDescription = {1xN} cell of variables names



%detect data list
dataList=dataSetSpec.dataList;


[~,dataSelect]=ismember(dataList,dataName); dataSelect(dataSelect==0)=[];

%trim relevant items
dataRaw         = data(:,dataSelect);
datesRaw        = dates;
dataName        = dataName(dataSelect);
dataDescription = dataDescription(dataSelect);
logTransform    = logTransform(dataSelect);


data  =dataRaw; 
dates =datesRaw;


%take logs
logTransform(and(any(data<0),logTransform))=false; %remove transform if <0
data(:,logTransform)=log(data(:,logTransform))*100;

if isfield(dataSetSpec,'eDataList') && ~isempty(dataSetSpec.eDataList)
    %transform to yoy growth rates (for comparison with expectations)
    if any(dataSetSpec.yoyTransform)

        if dataSetSpec.yoyTransform(ismember(dataSetSpec.dataList,'UNRATE'))

            data(:,ismember(dataSetSpec.dataList,'UNRATE'))=...
                filter(ones(12,1)./12,1,data(:,ismember(dataSetSpec.dataList,'UNRATE')));


            data(1:12,ismember(dataSetSpec.dataList,'UNRATE')) =nan;
            dataSetSpec.yoyTransform(ismember(dataSetSpec.dataList,'UNRATE'))=false;
        end

        data(2:end,dataSetSpec.yoyTransform)=filter(ones(12,1),1,diff(data(:,dataSetSpec.yoyTransform)));
        data(1:12,dataSetSpec.yoyTransform) =nan;


        if all(dataSetSpec.yoyTransform)

            dates(1:12)=[];
        end

    end
end


%fill up NaNs
if dataSetSpec.interpolateMissing
    
    nanOpt.method  = 2;
    nanOpt.winsize = 1;
    
    data=removeNaNs(data,nanOpt);
    
end

%select relevant time span
timeSelect= dates >= dataSetSpec.beginSet & dates <= dataSetSpec.endSet;

dataM  = data(timeSelect,:);
datesM = dates(timeSelect);

%remove rigged edges
edges = any(isnan(dataM),2);

dataM  = dataM(~edges,:);
datesM = datesM(~edges);

%plot data
if dataSetSpec.plotData
    
    for j=1:size(dataM,2)
        
        figure;
        subplot(2,1,1)
        plot(datesRaw,dataRaw(:,j)); axis tight; grid on; dateaxis('x',10);
        set(gca,'FontSize',9)
        title([dataDescription{j}, ' all obs'],'FontSize',11)
        %
        subplot(2,1,2)
        plot(datesM,dataM(:,j)); axis tight; grid on; dateaxis('x',10);
        set(gca,'FontSize',9)
        title('series over selected sample')
        %
        pause;
        close(gcf)

        
    end
    
end

if isfield(dataSetSpec,'beginPreSample') && isfield(dataSetSpec,'endPreSample')

    %-build presample----------------------------------------------------------
    timeSelect= dates >= dataSetSpec.beginPreSample & dates <= dataSetSpec.endPreSample;

    preSdates = dates(timeSelect);
    preSdata  = data(timeSelect,:);


    %only for initialization
    nanOpt.method  = 1;
    nanOpt.winsize = 1;

    %preSdata=removeNaNs(preSdata,nanOpt);

    
    %load output    
    res.preSdata    = preSdata;
    res.preSdates   = preSdates;
    
end

%load output
res.data        = dataM;
res.dates       = datesM;
res.varname     = dataName;
res.varLongName = dataDescription;
