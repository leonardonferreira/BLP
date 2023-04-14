
function modelSpec=buildDataStructure(fileName)
%packs together data from different sources


%read FRED-MD data sheet 
[tempData,tempText]=xlsread([pwd '/datafolder/' fileName]);

%unpack
data   =tempData(2:end,2:end);
dates  =x2mdate(tempData(2:end,1));
tcodes =tempData(1,2:end);
labels =tempText(1,2:end);


%remove reserves from data set
remove   =[67 68 69];
data(:,remove)   =[];
tcodes(:,remove) =[];
labels(remove)   =[];


%store original numbers
leveldata     =data;

%transform for stationarity (St. Louis Fed code)
data          =prepare_missing(data,tcodes);


%remove outliers (replace with centered MA(3))
n =size(data,2);

for j=1:n
   
    data(:,j) =remove_outliers(data(:,j));
end



%load all into final structure
dataStructure.data      =data;
dataStructure.dates     =dates;
dataStructure.varname   =labels;
dataStructure.leveldata =leveldata;
dataStructure.tcodes    =tcodes;
modelSpec.dataStructure =dataStructure;


