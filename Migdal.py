##################################################
"""Migdal: Read/interpolate Migdal ionisation
        probability tables."""

print("""------------------------------------------
Migdal ionisation probabilities
P. Cox, M. Dolan, C. McCabe, H. Quiney (2022)
arXiv:2207.xxxxx
------------------------------------------""")

__author__ = "Peter Cox"
__email__ = "peter.cox@unimelb.edu.au"
__version__ = 1.0

#*******************************************************************************

import numpy as np
import os.path
from scipy.interpolate import CubicSpline, RectBivariateSpline, RegularGridInterpolator

#*******************************************************************************

class Migdal:
    """Class to store/access differential Migdal probabilities."""

    alpha = 1/137.0359895

    # Occupied orbitals and energies (keV)
    elements = ['He','C','F','Ne','Si','Ar','Ge','Kr','Xe']
    orbitals = {}
    orbitals['He'] = [('1s',2.497980e-02)]
    orbitals['C']  = [('1s',3.083175e-01),('2s',1.921393e-02),('2p-',1.178626e-02),('2p',1.179453e-02)]
    orbitals['F']  = [('1s',7.186886e-01),('2s',4.288152e-02),('2p-',1.966527e-02),('2p',1.997926e-02)]
    orbitals['Ne'] = [('1s',8.930084e-01),('2s',5.267702e-02),('2p-',2.320668e-02),('2p',2.308252e-02)]
    orbitals['Si'] = [('1s',1.877873e+00),('2s',1.687017e-01),('2p-',1.167549e-01),('2p',1.159407e-01),('3s',1.500285e-02),('3p-',7.254801e-03),('3p',6.468341e-03)]
    orbitals['Ar'] = [('1s',3.241600e+00),('2s',3.377363e-01),('2p-',2.620990e-01),('2p',2.597887e-01),('3s',3.500978e-02),('3p-',1.620128e-02),('3p',1.599535e-02)]
    orbitals['Ge'] = [('1s',1.118596e+01),('2s',1.454912e+00),('2p-',1.288286e+00),('2p',1.256025e+00),('3s',2.019073e-01),('3p-',1.452096e-01),('3p',1.405931e-01),('3d-',4.426870e-02),('3d',4.359360e-02),('4s',1.566679e-02),('4p-',7.047754e-03),('4p',6.185060e-03)]
    orbitals['Kr'] = [('1s',1.441346e+01),('2s',1.961391e+00),('2p-',1.765333e+00),('2p',1.711030e+00),('3s',3.054330e-01),('3p-',2.345592e-01),('3p',2.262024e-01),('3d-',1.027950e-01),('3d',1.014111e-01),('4s',3.232023e-02),('4p-',1.473540e-02),('4p',1.399614e-02)]
    orbitals['Xe'] = [('1s',3.475594e+01),('2s',5.509354e+00),('2p-',5.161449e+00),('2p',4.835587e+00),('3s',1.170374e+00),('3p-',1.024780e+00),('3p',9.612494e-01),('3d-',7.081319e-01),('3d',6.948998e-01),('4s',2.293898e-01),('4p-',1.755814e-01),('4p',1.628001e-01),('4d-',7.377911e-02),('4d',7.166829e-02),('5s',2.748725e-02),('5p-',1.340357e-02),('5p',1.196770e-02)]

    ##################################################

    def __init__(self, element):

        self._dpI1 = None
        self._dpIE = None
        self._dpI2 = None

        self._dpI1orbital = None
        self._dpI2orbital = None

        self._p0 = None
        self._pE1 = None
        self._pE2 = None
        self._pI1 = None
        self._pIE = None
        self._pI21 = None
        self._pI2 = None
        self._ptotal = None

        self._pI1orbital = None
        self._pI21orbital = None
        self._pI2orbital = None

        self.element = element
        self.orbitals = Migdal.orbitals[element]
        self.emin = 1.0E-4 # keV
        self.emax = 20.0 # keV
        self.vmin = None
        self.vmax = None

    ##################################################

    def dpI1(self, points, orbital=None):
        """Returns differential probability for single ionisation without excitation."""

        if orbital is None:
            return Migdal._return_probability(self._dpI1, (points[:,0],points[:,1]), ev=True)
        else:
            return Migdal._return_probability(self._dpI1orbital, (points[:,0],points[:,1]), orbital, ev=True)

    ##################################################

    def dpIE(self, points):
        """Returns differential probability for single ionisation with excitation."""

        return Migdal._return_probability(self._dpIE, (points[:,0],points[:,1]))

    ##################################################

    def dpI2(self, points, orbital=None):
        """Returns differential probability for double ionisation."""

        if orbital is None:
            return Migdal._return_probability(self._dpI2, points)
        else:
            return Migdal._return_probability(self._dpI2orbital, points, orbital)

    ##################################################

    def p0(self, points):
        """Returns probability for no transition."""

        return Migdal._return_probability(self._p0, points)

    ##################################################

    def pE1(self, points):
        """Returns probability for excitation."""

        return Migdal._return_probability(self._pE1, points)

    ##################################################

    def pE2(self, points):
        """Returns probability for excitation."""

        return Migdal._return_probability(self._pE2, points)

    ##################################################

    def pI1(self, points, orbital=None):
        """Returns integrated probability for single ionisation without excitation."""

        if orbital is None:
            return Migdal._return_probability(self._pI1, points)
        else:
            return Migdal._return_probability(self._pI1orbital, points, orbital)

    ##################################################

    def pIE(self, points):
        """Returns integrated probability for single ionisation with excitation."""

        return Migdal._return_probability(self._pIE, points)

    ##################################################

    def pI2(self, points, orbital=None):
        """Returns integrated probability for double ionisation with both electrons above threshold."""

        if orbital is None:
            return Migdal._return_probability(self._pI2, points)
        else:
            return Migdal._return_probability(self._pI2orbital, points, orbital)

    ##################################################

    def pI21(self, points, orbital=None):
        """Returns integrated probability for double ionisation with one electron above threshold."""

        if orbital is None:
            return Migdal._return_probability(self._pI21, points)
        else:
            return Migdal._return_probability(self._pI21orbital, points, orbital)

    ##################################################

    def Ptotal(self, points):
        """Returns total integrated probability."""

        if self._pI1 is None:
            print('Call Load() method to initialise.')
            return None

        Ptotal = np.zeros(len(points))

        if self._p0 is not None:
            Ptotal = Ptotal + Migdal._return_probability(self._p0, points)
        if self._pE1 is not None:
            Ptotal = Ptotal + Migdal._return_probability(self._pE1, points)
        if self._pE2 is not None:
            Ptotal = Ptotal + Migdal._return_probability(self._pE2, points)
        if self._pI1 is not None:
            Ptotal = Ptotal + Migdal._return_probability(self._pI1, points)
        if self._pIE is not None:
            Ptotal = Ptotal + Migdal._return_probability(self._pIE, points)
        if self._pI2 is not None:
            Ptotal = Ptotal + Migdal._return_probability(self._pI2, points)

        return Ptotal

    ##################################################

    def load_total_probabilities(self, dark_matter=False, double=True, e_threshold=None, excitations=False, inclusive=False, integrated=False, velocity_grid='linear'):
        """Initialise probabilities."""

        if inclusive:
            excitations = False
            double = False

        # No transition
        self._p0 = Migdal._interpolate('{0}/{0}_no-transition.txt'.format(self.element), 'N')

        # Excitations
        if excitations:
            # Single excitation
            self._pE1 = Migdal._interpolate('{0}/excitations/{0}_single-excitation.txt'.format(self.element), 'E')
            if double:
                # Double excitation
                self._pE2 = Migdal._interpolate('{0}/excitations/{0}_double-excitation.txt'.format(self.element), 'E')

                # Single ionisation + excitation
                pIE = Migdal._interpolate('{0}/excitations/{0}_excitation+ionisation.txt'.format(self.element), 'EI', integrated=integrated, e_threshold=e_threshold)
                if integrated:
                    self._pIE = pIE
                else:
                    self._dpIE = pIE

        # Single ionisation
        if inclusive:
            pI1 = Migdal._interpolate('{0}/semi-inclusive/{0}_semi-inclusive.txt'.format(self.element), 'SI', integrated=integrated, e_threshold=e_threshold)
        else:
            # Logarithmic velocity grid
            if dark_matter or velocity_grid == 'log':
                pI1 = Migdal._interpolate('{0}/single-ionisation/{0}_single-ionisation_DM.txt'.format(self.element), 'SI', integrated=integrated, e_threshold=e_threshold)
            # Linear velocity grid
            else:
                pI1 = Migdal._interpolate('{0}/single-ionisation/{0}_single-ionisation.txt'.format(self.element), 'SI', integrated=integrated, e_threshold=e_threshold)

        if integrated:
            self._pI1 = pI1
        else:
            self._dpI1 = pI1

        # Double ionisation
        if double:
            pI2 = Migdal._interpolate('{0}/double-ionisation/{0}_double-ionisation.txt'.format(self.element), 'DI', integrated=integrated, e_threshold=e_threshold)
            if integrated:
                self._pI21 = pI2[0]
                self._pI2 = pI2[1]
            else:
                self._dpI2 = pI2

    ##################################################

    def load_probabilities(self, dark_matter=False, double=True, e_threshold=None, inclusive=False, integrated=False, velocity_grid='linear'):
        """Initialise differential probabilities for individual orbitals."""

        if inclusive:
            double = False

        # Single ionisation
        pI1orbital = {}
        for orbital in Migdal.orbitals[self.element]:
            if inclusive:
                pI1orbital[orbital[0]] = self._interpolate('{0}/semi-inclusive/{0}_{1}_semi-inclusive.txt'.format(self.element,orbital[0]), 'SI', integrated=integrated, e_threshold=e_threshold)
            else:
                # Logarithmic velocity grid
                if dark_matter or velocity_grid == 'log':
                    pI1orbital[orbital[0]] = self._interpolate('{0}/single-ionisation/{0}_{1}_DM.txt'.format(self.element,orbital[0]), 'SI', integrated=integrated, e_threshold=e_threshold)
                # Linear velocity grid
                else:
                    pI1orbital[orbital[0]] = self._interpolate('{0}/single-ionisation/{0}_{1}.txt'.format(self.element,orbital[0]), 'SI', integrated=integrated, e_threshold=e_threshold)

        if integrated:
            self._pI1orbital = pI1orbital
        else:
            self._dpI1orbital = pI1orbital

        # Double ionisation
        if double:
            if integrated:
                self._pI21orbital = {}
                self._pI2orbital = {}
            self._dpI2orbital = {}

            for i in range(len(Migdal.orbitals[self.element])):
                for j in range(i,len(Migdal.orbitals[self.element])):
                    pair = Migdal.orbitals[self.element][i][0] + Migdal.orbitals[self.element][j][0]
                    if self.element == 'C' and pair == '2p-2p':
                        continue

                    pI2orbital = self._interpolate('{0}/double-ionisation/{0}_{1}.txt'.format(self.element,pair), 'DI', integrated=integrated, e_threshold=e_threshold)

                    if integrated:
                        self._pI21orbital[pair] = pI2orbital[0]
                        self._pI2orbital[pair] = pI2orbital[1]
                    else:
                        self._dpI2orbital[pair] = pI2orbital

    ##################################################

    @staticmethod
    def _interpolate(filename, type, integrated=False, e_threshold=None):
        """Read Migdal table and return interpolating function."""

        if not os.path.exists(filename):
            print('Migdal::Interpolate: Missing file {}'.format(filename))
            return None

        ne = 0
        nv = 0
        
        # Read Migdal data from file
        with open (filename) as f:
            for line in f:
                line = line.split()

                try:
                    # Read data points
                    if nv > 0:
                        # Energies
                        if (type =='SI' or type == 'EI') and ipt % nv == 0:
                            energies[ien] = float(line[0])
                            ien += 1
                        elif type == 'DI' and ipt % (ne*nv) == 0:
                            energies[ien] = float(line[0])
                            ien += 1

                        # Velocities
                        if ipt < nv:
                            if type == 'N' or type == 'E':
                                ln_velocities[ipt] = np.log(float(line[0]))
                            elif type == 'SI' or type == 'EI':
                                ln_velocities[ipt] = np.log(float(line[1]))
                            elif type == 'DI':
                                ln_velocities[ipt] = np.log(float(line[2]))

                        # Probabilities
                        if type == 'N' or type == 'E':
                            dPdE[ipt] = float(line[1])
                        elif type == 'SI' or type == 'EI':
                            dPdE[ipt] = float(line[2])
                        elif type == 'DI':
                            dPdE[ipt] = float(line[3])
                        ipt += 1

                    # Set-up arrays
                    elif 'energies' in line[2]:
                        ne = int(line[4])
                        energies = np.zeros(ne)
                    
                    elif 'velocities' in line[2]:
                        nv = int(line[4])
                        ln_velocities = np.zeros(nv)

                        if type == 'N' or type == 'E':
                            dPdE = np.zeros(nv)
                        elif type == 'SI':
                            dPdE = np.zeros(ne * nv)
                        elif type == 'EI':
                            dPdE = np.zeros(ne * nv)
                        elif type == 'DI':
                            dPdE = np.zeros(ne * ne * nv)
                        ipt = 0
                        ien = 0

                except (IndexError, ValueError):
                    continue

        # Interpolating functions
        # No transition, excitation
        if type == 'N' or type == 'E':
            dPdE = np.log(dPdE)
            return CubicSpline(ln_velocities, dPdE)

        ln_energies = np.log(energies)

        # Interpolating functions for differential probabilities
        if not integrated:
            # Single ionisation
            if type == 'SI' or type == 'EI':
                dPdE = np.log(dPdE)
                dPdE = dPdE.reshape((ne,nv))
                return RectBivariateSpline(ln_energies, ln_velocities, dPdE, kx=3, ky=3)

            # double ionisation
            elif type == 'DI':
                dPdE = np.log(dPdE)
                dPdE = dPdE.reshape((ne,ne,nv))
                return RegularGridInterpolator(np.array([ln_energies,ln_energies,ln_velocities]), dPdE, method='linear')

        # Interpolating functions for integrated probabilities
        else:
            # Integral limits
            ln_emin = ln_energies[0]
            ln_emax = ln_energies[-1]
            if e_threshold is None:
                ln_ecut = ln_emin
            elif np.log(e_threshold) < ln_emin:
                e_threshold = None
                ln_ecut = ln_emin
            else:
                ln_ecut = np.log(e_threshold)

            # Single ionisation
            if type == 'SI' or type == 'EI':
                dPdE = dPdE.reshape((ne,nv))
                integrated = np.zeros(nv)

                for i in range(len(ln_velocities)):
                    integrand = energies * dPdE[:,i]
                    integrated[i] = np.log(CubicSpline(ln_energies, integrand, extrapolate=False).integrate(ln_ecut, ln_emax))
                    
                return CubicSpline(ln_velocities, integrated)

            # double ionisation
            elif type == 'DI':
                dPdE = dPdE.reshape((ne,ne,nv))
                integrated = np.zeros((2,nv))

                for i in range(len(ln_velocities)):
                    integrand = dPdE[:,:,i] 
                    for j in range(len(energies)):
                        for k in range(len(energies)):
                            integrand[j,k] = integrand[j,k] * energies[j] * energies[k]

                    # One electron above threshold
                    if e_threshold is not None:
                        integrated[0,i] = np.log(2 * RectBivariateSpline(ln_energies, ln_energies, dPdE[:,:,i], kx=3, ky=3).integral(ln_emin, ln_ecut, ln_ecut, ln_emax))
                    else:
                        integrated[0,i] = -1e10

                    # Both electrons above threshold
                    integrated[1,i] = np.log(RectBivariateSpline(ln_energies, ln_energies, dPdE[:,:,i], kx=3, ky=3).integral(ln_ecut, ln_emax, ln_ecut, ln_emax))

                return (CubicSpline(ln_velocities, integrated[0]), CubicSpline(ln_velocities, integrated[1]))

    ##################################################

    @staticmethod
    def _return_probability(fn, points, orbital=None, ev=False):

        if orbital is not None:
            if fn is None:
                print('Call LoadOrbitals() method to initialise.')
                return None

            try:
                if ev:
                    return np.exp(fn[orbital].ev(points[0], points[1]))
                else:
                    return np.exp(fn[orbital](points))
            except KeyError:
                print("'{}' is not a valid orbital.".format(orbital))
                return None

        else:
            if fn is None:
                print('Call Load() method to initialise.')
                return None

            if ev:
                return np.exp(fn.ev(points[0], points[1]))
            else:
                return np.exp(fn(points))

    ##################################################

#*******************************************************************************