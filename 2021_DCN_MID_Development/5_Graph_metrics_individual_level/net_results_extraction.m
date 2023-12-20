function [output]=net_results_extraction(subfiles)
% extract sw reuslts from parfor version of sw_batch
% Outputs:
%            'Cp:
%                   Clustering coefficients of real networks.
%            'Lp:
%                   Shortest path lengths of real networks.
%          'locE:
%                   Local efficiencies of real networks.
%            'gE:
%                   Global efficiencies of real networks.
%           'deg:
%                   Degrees of real networks.
%            'bw:
%                   Betweennesses of real networks.
%           'mod:
%                   The modularities of real networks.
%
%
%       'Cpratio':
%                   The ratio' of clustering coefficients between real
%                   networks and comparable random networks (i.e., gamma).
%       'Lpratio':
%                   The ratio' of shortest path lengths between real
%                   networks and comparable random networks (i.e., lambda).
%     'locEratio':
%                   The ratio' of local efficiencies between real
%                   networks and comparable random networks (i.e., gamma).
%       'gEratio':
%                   The ratio' of global efficiencies between real
%                   networks and comparable random networks (i.e., lambda).
%      'degratio':
%                   The ratio' of degrees between real networks and comparable
%                   random networks.
%       'bwratio':
%                   The ratio' of Betweennesses between real networks and
%                   comparable random networks.
%      'modratio':
%                   The ratio' of modularities between real networks and
%                   comparable random networks.
%
%           'aCp:
%                   The area under curve of clustering coefficients.
%           'aLp:
%                   The area under curve of shortest path lengths.
%         'alocE:
%                   The area under curve of local efficiencies.
%           'agE:
%                   The area under curve of global efficiencies.
%          'adeg:
%                   The area under curve of degrees.
%           'abw:
%                   The area under curve of betweennesses.
%          'amod:
%                   The area under curve of modularities.
%
%      'aCpratio':
%                   The area under curve of ratio's of clustering coefficients.
%      'aLpratio':
%                   The area under curve of ratio's of shortest path lengths.
%    'alocEratio':
%                   The area under curve of ratio's of local efficiencies.
%      'agEratio':
%                   The area under curve of ratio's of global efficiencies.
%     'adegratio':
%                   The area under curve of ratio's of degrees.
%      'abwratio':
%                   The area under curve of ratio's of betweennesses.
%     'amodratio':
%                   The area under curve of ratio's of modularities.

var2ext1={'Cp'  % densityx1
'Lp'
'locE'
'gE'
'deg'
'bw'
'mod'
'Cpratio'
'Lpratio'
'locEratio'
'gEratio'
'degratio'
'bwratio'
'modratio'}
var2ext2={    %density independent results one value
'aCp'  % densityx1
'aLp'
'alocE'
'agE'
'adeg'
'abw'
'amod'
'aCpratio'
'aLpratio'
'alocEratio'
'agEratio'
'adegratio'
'abwratio'
'amodratio'}


for subi=1:length(subfiles)
    disp(['subject========', num2str(subi)])
     disp(subfiles{subi})
    tmp_load=load(subfiles{subi},'net_sum');
    sub_net_results=tmp_load.net_sum;

    for var_i=1:length(var2ext1)
        varname=var2ext1{var_i};
        eval(['output.',varname,'(:,subi)=sub_net_results.',varname,';']);
    end
    
      for var_i=1:length(var2ext2)
        varname=var2ext2{var_i};
        eval(['output.',varname,'(subi)=sub_net_results.',varname,';']);
      end
    
end

clear tmp_load sub_net_results
end