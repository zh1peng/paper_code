function [output,rows]=select_group(input,column,match,option)
% input is a cell or matrix
%column is a number
%match is a string
%option: 1.select 2.exclude (in this mode rows are left rows)
%Example:
% [male_cov,male_rows]=select_group(cov,8,'0');
% male_ids=subid(male_rows);
% [female_cov,female_rows]=select_group(cov,8,'1');
% female_ids=subid(female_rows);

%%%%%%%not test on cell. it works well with matrix 2016.06.10;

if option==1 %select
   if iscell(input)&&ischar(match)
    rows=find(strcmp(input(:,column),match)==1);
elseif ismatrix(input)&&ischar(match)
    tmp=str2num(match);
    rows=find(input(:,column)==tmp);
end
output=input(rows,:);

elseif option==2 %exclude
  
    if iscell(input)&&ischar(match)
  rows=cellfun('isequal', input(:,column),match)
  
 elseif ismatrix(input)&&ischar(match)
    tmp=str2num(match);
    rows=(input(:,column)==tmp);
  end
  output=input(~rows,:);
  rows=~rows;
    
end
end
