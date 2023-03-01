import os
import sys
import glob
import numpy as np
import pandas as pd
from datetime import datetime
from gemmi import cif

path_to_output = sys.argv[1]

outputs = glob.glob(os.path.join(path_to_output, "System_0", "*.data"))

# parse outputs to list of dictionaries
data = []
for output in outputs:
    with open(output) as f:
        output_split = os.path.splitext(os.path.basename(output))[0].split('_')
        temperature = output_split[-2]
        pressure = output_split[-1]
        code = ""
        time = ""
        operator = ""
        framework_name = ""
        adsorbate_name = ""
        adsorbate_definition = ""
        ff_definition = ""
        loading = ""
        loading_error = ""
        excess_loading = ""
        excess_loading_error = ""
        enthalpy = ""
        enthalpy_error = ""

        for index, line in enumerate(f):
            if index == 2:
                code = line.strip('\n').replace(' ','-')
            elif index == 7:
                time = datetime.strptime(line.strip('\n'), "%a %b %d %H:%M:%S %Y")
            elif "Hostname:" in line:
                operator = line.split()[-1]
            elif "Framework name:" in line:
                framework_name = line.split()[-1]
            elif "Forcefield: " in line:
                ff_definition = line.split()[-1]
            elif "(Adsorbate molecule)" in line:
                adsorbate_name = line.split()[2].strip("[]")
            elif "MoleculeDefinitions: " in line:
                adsorbate_definition = line.split()[-1]
            elif "Average loading absolute [mol/kg framework]" in line:
                loading = float(line.split()[5])
                loading_error = float(line.split()[7])
            elif "Average loading excess [mol/kg framework]" in line:
                excess_loading = float(line.split()[5])
                excess_loading_error = float(line.split()[7])
            elif "[KJ/MOL]" in line:
                enthalpy = float(line.split()[0])
                enthalpy_error = float(line.split()[2])

        data.append({
            "temperature": temperature,
            "pressure": pressure,
            "code": code,
            "time": time,
            "operator": operator,
            "framework_name": framework_name,
            "adsorbate_name": adsorbate_name,
            "adsorbate_definition": adsorbate_definition,
            "ff_definition": ff_definition,
            "loading": loading,
            "loading_error": loading_error,
            "excess_loading": excess_loading,
            "excess_loading_error": excess_loading_error,
            "enthalpy": enthalpy,
            "enthalpy_error": enthalpy_error
        })


# covert to dataframe
df = pd.DataFrame(data)
df=df.sort_values(by=['pressure'])

def checkifunique(dataframe, column):
    try:
        nunique = dataframe[column].nunique()
        assert nunique == 1
    except:
        print(column+' not unique')


d = cif.Document()
d.add_new_block('raspa2aif')

block = d.sole_block()
block.set_pair('_audit_aif_version', '6acf6ef')

#label metadata

checkifunique(df,'operator')
block.set_pair('_exptl_operator',  str(df['operator'][0]))
#using mean of the simulation dates
block.set_pair('_simltn_date', str(df['time'].mean().isoformat()))

checkifunique(df,'code')
block.set_pair('_simltn_code', str(df['code'][0]))

block.set_pair('_exptl_method', 'GCMC')

checkifunique(df,'adsorbate_name')
block.set_pair('_exptl_adsorptive', str(df['adsorbate_name'][0]))
checkifunique(df,'temperature')
block.set_pair('_exptl_temperature', str(df['temperature'][0]))

checkifunique(df,'adsorbate_definition')
block.set_pair('_simltn_forcefield_adsorbent', str(df['adsorbate_definition'][0]))
checkifunique(df,'ff_definition')
block.set_pair('_simltn_forcefield_adsorbent', str(df['ff_definition'][0]))

checkifunique(df,'framework_name')
block.set_pair('_adsnt_material_id', str(df['framework_name'][0]))

#currently assuming units in raspa outputs don't differ
block.set_pair('_units_temperature', 'K')
block.set_pair('_units_energy', 'kJ/mol')
block.set_pair('_units_loading','mol/kg')
block.set_pair('_units_pressure','Pa')

#format adsorption

loop_ads = block.init_loop('_adsorp_', ['pressure', 'amount_absolute', 'amount_absolute_error', 'amount_excess', 'amount_excess_error', 'enthalpy', 'enthalpy_error'])
loop_ads.set_all_values([
    ['%.5E' % item for item in df['pressure'].astype(float)],
    ['%.5E' % item for item in df['loading'].astype(float)],
    ['%.5E' % item for item in df['loading_error'].astype(float)],
    ['%.5E' % item for item in df['excess_loading'].astype(float)],
    ['%.5E' % item for item in df['excess_loading_error'].astype(float)],
    ['%.5E' % item for item in df['enthalpy'].astype(float)],
    ['%.5E' % item for item in df['enthalpy_error'].astype(float)],
])

d.write_file('test.aif')