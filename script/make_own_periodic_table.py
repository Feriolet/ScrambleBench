from rdkit import Chem
import pickle

# Rationale:
# The Periodic Table instance is not pickable in the GenBench3D Analysis, so we will
# pre-generated a dictionary with the necessary atom property in pkl format
# This will enable us to use GenBench with multiprocessing
# Granted, this mean that we will edit the function used to call periodic table in
# the GenBench3D source code

pt = Chem.GetPeriodicTable()

# Extract relevant information into a picklable format
periodic_table_data = {}
for i in range(1, pt.GetMaxAtomicNumber() + 1):
    try:
        symbol = pt.GetElementSymbol(i)
        atomic_weight = pt.GetAtomicNumber(symbol)
        vdw = pt.GetRvdw(i)
        periodic_table_data[i] = {'GetElementSymbol': symbol, 'GetAtomicNumber': atomic_weight, 'GetRvdw': vdw}
        periodic_table_data[symbol] = {'GetElementSymbol': symbol, 'GetAtomicNumber': atomic_weight, 'GetRvdw': vdw}
    except RuntimeError:  # Handle cases where an atomic number might not exist
        pass


with open('periodic_table_info.pkl', 'wb') as f:
    pickle.dump(periodic_table_data, f)

# Test to see it works:
with open('periodic_table_info.pkl', 'rb') as f:
    loaded_data = pickle.load(f)
