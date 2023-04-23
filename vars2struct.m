function s = vars2struct(varnames)

s = struct();
for varname = string(varnames(:)')
    s.(varname) = evalin("caller", varname);
end