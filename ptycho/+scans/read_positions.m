%READ_POSITIONS load positions
function p = read_positions(p)
import utils.verbose

% check for continuous scan
p.scan.is_cont = false;

% Variables for all scans
p.positions_real = [];
p.positions_orig = [];


% read position data
if ~isfield(p,'src_positions') || isempty(p.src_positions)
    error(' p.src_positions is not set')
end

% load the positions from a provided function if possible 
position_func = str2func(sprintf('scans.positions.%s', p.src_positions));
p = position_func(p);


scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
for ii = 1:p.numscans
    p.scanindexrange(ii,:) = [scanfirstindex(ii) scanfirstindex(ii+1)-1];
    p.scanidxs{ii} = p.scanindexrange(ii,1):p.scanindexrange(ii,end); 
end

p.positions_orig = p.positions_real;
%size(p.positions_orig)
%scatter(p.positions_orig(:,1),p.positions_orig(:,2),'x')

end