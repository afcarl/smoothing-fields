% function [outVol param]=ReadXml(inVolName)
%
% Input inVolName is a .xml file path (full) containing dimension, type etc 
% information of the input volume. The input .raw file should be at the same 
% location of input xml file. param.dim=dimension, param.type=type of data,
% param.endian=endianness and param.res=resolution of the data,

function [outVol param]=ReadXml(inVolName)

xDoc=xmlread(inVolName);
t1=xDoc.getElementsByTagName('Data-type');
type=t1.item(0).getFirstChild.getData.toCharArray;
type=type';
type=lower(type);
if strcmp(type,'unsigned byte')==1
    type='uint8';
elseif strcmp(type,'float')==1
    type='float32';
elseif strcmp(type,'unsigned short')==1
    type='ushort';
elseif strcmp(type,'short')==1
    type='short';
elseif strcmp(type,'byte')==1
    type='int8';
elseif strcmp(type,'integer')==1
    type='int32';
end   
flag=0;
fp1=fopen(inVolName,'r');
while feof(fp1)~=1
    L=fgetl(fp1);
    k=findstr(L,'nDimensions');
    if k
        flag=1;
        dim=str2num(L(k+length('nDimensions')+2));
        break;
    end
end
if flag==0
    dim=3 % default dimension
end

fclose(fp1);
    
t2=xDoc.getElementsByTagName('Endianess');
e=t2.item(0).getFirstChild.getData.toCharArray;
e=e(1);

e=lower(e);
t3=xDoc.getElementsByTagName('Extents');
for i=0:dim-1
    S=t3.item(i).getFirstChild.getData.toCharArray;
    N1(i+1)=str2num(S');
end


t4=xDoc.getElementsByTagName('Resolution');
for i=0:dim-1
    S=t4.item(i).getFirstChild.getData.toCharArray;
    res(i+1)=str2num(S');
end

% fprintf('Size of input volume=[%d %d %d],type=%s\n',N1(1),N1(2),N1(3),type);
volname=inVolName(1:end-4);
Involname=strcat(volname,'.raw');
fp1=fopen(Involname,'r');

if strcmp(e,'l')
    outVol=fread(fp1,type,'l');
    param.endian='little';
elseif strcmp(e,'b')
    outVol=fread(fp1,type,'b');
    param.endian='big';
else
    fprintf('unknown type, check xml file for Endianness\n');
    return;
end

outVol=reshape(outVol,N1);
param.dim=N1;
param.type=type;
param.res=res;

fclose(fp1);




