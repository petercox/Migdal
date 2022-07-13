# Migdal

*Precise Predictions and New Insights for Atomic Ionisation from the Migdal Effect*  
P. Cox, M. J. Dolan, C. McCabe, H. M. Quiney
arXiv:2208.xxxxx

## Probability tables

Exclusive and semi-inclusive single-ionisation probabilities are available for  
**He, C, F, Ne, Si, Ar, Ge, Kr, Xe**

Exclusive double ionisation probabilities are available for  
**He, C, F, Ne, Ar, Xe**

Single & Double bound excitation probabilities are available for  
**He, C, F, Ne, Ar, Xe**

Excitation+ionisation probailities are available for  
**He, C, F, Ne**

## Example usage

For more detailed examples see examples.py, examples.nb.

### Python 

```Python
# Initialisation
Ne = Migdal('Ne')

# Orbital binding energies
for (orbital, energy) in Ne.orbitals:
    print('{:4}: {:.3e} keV'.format(orbital, energy))

# Semi-inclusive differential probabilities for neutron scattering
Ne.load_probabilities(inclusive=True)
dPdE = Ne.dpI1(np.array([[log(energy), log(velocity)]]), '1s')

# Exclusive differential probabilities for dark matter
Ne.load_orbitals(dark_matter=True, double=False)
dPdE = Ne.dpI1(np.array([[log(energy), log(velocity)]]), '1s')
```

### Mathematica 

```Mathematica
(* Orbital binding energies *)
$orbitals["Ne"]

(* Semi-inclusive differential probabilities for neutron scattering *)
LoadProbabilities("Ne", Inclusive->True)
dPdE = dpI1["1s"][Log[energy],Log[velocity]]

(* Exclusive differential probabilities for dark matter *)
LoadProbabilities("Ne", DarkMatter->True, Double->False)
dPdE = dpI1["1s"][Log[energy],Log[velocity]]
```

For support/questions please contact peter.cox@unimelb.edu.au