import json

f_json = '/shares/departments/AUG/users/git/a8//ncdf_out/aug34954test-1.json'
with open(f_json, 'r') as fjson:
    var_d = json.load(fjson)

print(var_d['ABC'])
print(var_d['TIX'])
print(var_d['SLATX'])
