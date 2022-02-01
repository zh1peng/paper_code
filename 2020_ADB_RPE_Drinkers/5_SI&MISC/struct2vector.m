function varargout=struct2vector(str_in)
%example: [type, lat, urevent]=struct2vector(EEG.event);
field_names=fieldnames(str_in);
for i=1:length(field_names)
for n=1:length(str_in)
    eval(sprintf('%s(n)=str_in(n).%s;' ,field_names{i},field_names{i}));
end
end
for ii=1:nargout
varargout{ii}=eval(sprintf('transpose(%s);',field_names{ii}));
end
end