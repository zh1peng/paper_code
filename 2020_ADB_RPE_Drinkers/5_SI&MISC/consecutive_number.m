function [idx, consecutive_n]=consecutive_number(vector_in,filter)
% vector_in = [34 35 36 78 79 80 81 82 84 85 86 102 103 104 105 106 107 201 202 203 204]
% consecutive_n:  number consecutive numbers
% idx: first number idx
% find consecutive numbers in a vector. Use it to find same button press
% markers.
% author zhipeng 2018/07/24
k = [true;diff(vector_in(:))~=1 ];
s = cumsum(k);
tmp_consecutive_n =  histc(s,1:s(end));
tmp_idx = find(k);
idx=tmp_idx(tmp_consecutive_n>filter);
consecutive_n=tmp_consecutive_n(tmp_consecutive_n>filter);
end
