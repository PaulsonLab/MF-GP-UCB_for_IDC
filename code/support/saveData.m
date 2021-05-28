% function [] = saveData(data, file_save_name, struct_name, overwrite =0)
% saves the [data] (must be flat) to a file named [file_save_name]. It will attempt to add
% to an existing file, or create a new file if one does not yet exist.
% [struct_name] is the structure name.
%
% a fourth optional variable, overwrite, can be passed if the file_save_name
% should overwrite an existing file.
%
% ie. file_save_name.struct_name(i) will return the i-th data
%
function [] = saveData(data, file_save_name, struct_name, overwrite)
if nargin<4
    overwrite =0;
end

try
    if ~overwrite   %try to save to existing file
        p1 = load(file_save_name);
        varname = eval(['p1.',struct_name]);
        ps = [varname;data];
        swap(struct_name,ps)
    elseif overwrite    %overwrite existing file
        swap(struct_name,data)
    end
catch   % the file does not exist so make a new one
    swap(struct_name,data)
end

save(file_save_name,struct_name)

end

% A support function to change variable names
function swap(struct_name, data)
    assignin('caller',struct_name,data)
end
