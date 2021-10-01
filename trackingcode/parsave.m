function [ ] = parsave( fname, var1, varname, append )
% inputs: fname   - filename you wish to save to
%         varname - name of variable you wish to save
%         append  - a flag (1 or 0) telling function to append to existing
%         .mat file or to save it as a new file entirely. It will overwrite
%         an existing file with the same name if not appended.

% This code 'parsave.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

eval([varname '=var1;']);

if append
    save(fname, '-append', varname);
else
    save(fname, varname);
end

end

