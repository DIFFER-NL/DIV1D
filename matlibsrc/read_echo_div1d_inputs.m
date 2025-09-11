function input = read_echo_div1d_inputs(filepath)
% input = read_echo_div1d_inputs(filepath) reads echo_div1d_inputs.txt from
% specified filepath
% returns input struct used for further matlab analysis functionality

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% Jan 2025

try % open the output file
    disp(['reading:  ',filepath])
    fid = fopen(filepath,'r');
catch
    input = struct; %output = struct;
    disp('echo_div1d_input.txt is not available');
    return
end
Text1 = textscan(fid,'git tag: %s','delimiter','\n');
version = Text1{1};
input = struct;
input.version = version;
input.numerics = read_namelist(fid, 'DIV1D_NUMERICS');
input.physics  = read_namelist(fid, 'DIV1D_PHYSICS' );
input.grid = []; %read_namelist(fid, 'DIV1D_GRID'); 
input.physics.num_impurities = length(input.physics.impurity_concentration);
fclose(fid);

end
