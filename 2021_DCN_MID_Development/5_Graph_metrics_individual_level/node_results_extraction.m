function [node_output]=node_results_extraction(subfiles)
% node.Cp:
%                   Clustering coefficients of each node in real networks.
%           node.Lp:
%                   Shortest path lengths of each node in real networks.
%         node.locE:
%                   Local efficiencies of each node in real networks.
%           node.gE:
%                   Global efficiencies of each node in real networks.
%          node.deg:
%                   Degrees of each node in real networks.
%           node.bw:
%                   Betweennesses of each node in real networks.
%
%       node.Cprand:
%                   Clustering coefficients of each node in random networks.
%       node.Lprand:
%                   Shortest path lengths of each node in random networks.
%     node.locErand:
%                   Local efficiencies of each node in random networks.
%       node.gErand:
%                   Global efficiencies of each node in random networks.
%      node.degrand:
%                   Degrees of each node in random networks.
%       node.bwrand:
%                   Betweennesses of each node in random networks.
%
%      node.Cpratio:
%                   The ratio of nodal clustering coefficients between real
%                   and comparable random networks.
%      node.Lpratio:
%                   The ratio of nodal shortest path lengths between real
%                   and comparable random networks.
%    node.locEratio:
%                   The ratio of nodal local efficiencies between real and
%                   comparable random networks.
%      node.gEratio:
%                   The ratio of nodal global efficiencies between real and
%                   comparable random networks.
%     node.degratio:
%                   The ratio of nodal degrees between real and comparable
%                   random networks.
%      node.bwratio:
%                   The ratio of nodal betweennesses between real and
%                   comparable random networks.
%
%          node.aCp:
%                   The area under curve of clustering coefficients.
%          node.aLp:
%                   The area under curve of shortest path lengths.
%        node.alocE:
%                   The area under curve of local efficiencies.
%          node.agE:
%                   The area under curve of global efficiencies.
%         node.adeg:
%                   The area under curve of degrees.
%          node.abw:
%                   The area under curve of betweennesses.
%
%     node.aCpratio:
%                   The area under curve of ratios of clustering coefficients.
%     node.aLpratio:
%                   The area under curve of ratios of shortest path lengths.
%   node.alocEratio:
%                   The area under curve of ratios of local efficiencies.
%     node.agEratio:
%                   The area under curve of ratios of global efficiencies.
%    node.adegratio:
%                   The area under curve of ratios of degrees.
%     node.abwratio:
%                   The area under curve of ratios of betweennesses.

var2ext1={'Cp'  % densityx1
'Lp'
'locE'
'gE'
'deg'
'bw'
'Cpratio'
'Lpratio'
'locEratio'
'gEratio'
'degratio'
'bwratio'}
var2ext2={
'aCp'  % densityx1
'aLp'
'alocE'
'agE'
'adeg'
'abw'%density independent results one value
'aCpratio'
'aLpratio'
'alocEratio'
'agEratio'
'adegratio'
'abwratio'
}

for subi=1:length(subfiles)
    disp(['subject========', num2str(subi)])
    disp(subfiles{subi})
    tmp_load=load(subfiles{subi},'node_sum');
    sub_node_results=tmp_load.node_sum;
    for var_i=1:length(var2ext1)
        varname=var2ext1{var_i};
        eval(['node_output.',varname,'(:,subi,:)=sub_node_results.',varname,';']);
    end
      for var_i=1:length(var2ext2)
        varname=var2ext2{var_i};
        eval(['node_output.',varname,'(:,subi)=sub_node_results.',varname,';']);
      end
      clear tmp_load sub_node_results
end





