function out = RowCol_Interleaver(inp,k)

out = reshape(inp,[],length(inp)/k)';
out = out(:)';

end
