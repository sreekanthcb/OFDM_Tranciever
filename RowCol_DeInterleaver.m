function out = RowCol_DeInterleaver(inp,k)

out = reshape(inp,[],k)';
out = out(:)';

end
