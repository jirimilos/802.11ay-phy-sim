function output = fftPerSymbol(input)
% The input could be either row or column vector, the output is always
% column vector
[nRows, nCols] = size(input);
if (nRows ~= 1) && (nCols ~= 1)
    error('The input is not a vector')
end

if nRows == 1 % the input is a row vector
    output = fft(input,[],1); 
    % the output should be column vector
    output = output.';
elseif nCols == 1 % the input is a column vector
    output = fft(input,[],2); 
else
    error('The input is not a vector')
end

end

