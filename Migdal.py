##################################################
"""Migdal: Read/interpolate Migdal ionisation
        probability tables."""

print("""------------------------------------------
Migdal ionisation probabilities
P. Cox, M. Dolan, C. McCabe, H. Quiney (2022)
arXiv:2208.12222
------------------------------------------""")

__author__ = "Peter Cox"
__email__ = "peter.cox@unimelb.edu.au"
__version__ = '1.1.0'

#*******************************************************************************

import numpy as np
import os.path
from scipy.interpolate import CubicSpline, NdBSpline, RectBivariateSpline, make_interp_spline

#*******************************************************************************

class Migdal:
    """Class to store/access differential Migdal probabilities."""

    # Directory containing the probability tables
    data_dir = os.path.dirname(os.path.abspath(__file__))

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
        self._dpI1_dipole = None
        self._dpIE = None
        self._dpI21 = None
        self._dpI2 = None

        self._dpI1_orbital = None
        self._dpI1_dipole_orbital = None
        self._dpI21_orbital = None
        self._dpI2_orbital = None

        self._p0 = None
        self._pE1 = None
        self._pE2 = None
        self._pI1 = None
        self._pI1_dipole = None
        self._pIE = None
        self._pI21 = None
        self._pI2 = None
        self._ptotal = None

        self._pE1_orbital = None
        self._pI1_orbital = None
        self._pI1_dipole_orbital = None
        self._pI21_orbital = None
        self._pI2_orbital = None

        self.element = element
        self.orbitals = Migdal.orbitals[element]
        self.emin = 1.0E-4 # keV
        self.emax = 20.0 # keV

    ##################################################

    def dpI1(self, points, orbital=None):
        """Returns differential probability for single ionisation without excitation."""

        if orbital is None:
            return Migdal._return_probability(self._dpI1, points)
        else:
            return Migdal._return_probability(self._dpI1_orbital, points, orbital)

    ##################################################

    def dpI1dipole(self, points, orbital=None):
        """Returns differential probability for single ionisation using dipole approximation."""

        if orbital is None:
            return Migdal._return_probability(self._dpI1_dipole, points)
        else:
            return Migdal._return_probability(self._dpI1_dipole_orbital, points, orbital)

    ##################################################

    def dpIE(self, points):
        """Returns differential probability for single ionisation with excitation."""

        return Migdal._return_probability(self._dpIE, points)

    ##################################################

    def dpI2(self, points, orbital=None):
        """Returns differential probability for double ionisation."""

        if orbital is None:
            return Migdal._return_probability(self._dpI2, points)
        else:
            return Migdal._return_probability(self._dpI2_orbital, points, orbital)

    ##################################################

    def dpI21(self, points, orbital=None):
        """Returns differential probability for double ionisation with one electron above threshold."""

        if orbital is None:
            return Migdal._return_probability(self._dpI21, points)
        else:
            return Migdal._return_probability(self._dpI21_orbital, points, orbital)

    ##################################################

    def p0(self, points):
        """Returns probability for no transition."""

        return Migdal._return_probability(self._p0, points)

    ##################################################

    def pE1(self, points, orbital=None):
        """Returns probability for excitation."""

        if orbital is None:
            return Migdal._return_probability(self._pE1, points)
        else:
            return Migdal._return_probability(self._pE1_orbital, points, orbital)

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
            return Migdal._return_probability(self._pI1_orbital, points, orbital)

    ##################################################

    def pI1dipole(self, points, orbital=None):
        """Returns probability for single ionisation using dipole approximation."""

        if orbital is None:
            return Migdal._return_probability(self._pI1_dipole, points)
        else:
            return Migdal._return_probability(self._pI1_dipole_orbital, points, orbital)

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
            return Migdal._return_probability(self._pI2_orbital, points, orbital)

    ##################################################

    def pI21(self, points, orbital=None):
        """Returns integrated probability for double ionisation with one electron above threshold."""

        if orbital is None:
            return Migdal._return_probability(self._pI21, points)
        else:
            return Migdal._return_probability(self._pI21_orbital, points, orbital)

    ##################################################

    def ptotal(self, points):
        """Returns total integrated probability."""

        if self._pI1 is None:
            print('Call load_total_probabilities() method to initialise.')
            return None

        ptotal = np.zeros(len(points))

        if self._p0 is not None:
            ptotal = ptotal + Migdal._return_probability(self._p0, points)
        if self._pE1 is not None:
            ptotal = ptotal + Migdal._return_probability(self._pE1, points)
        if self._pE2 is not None:
            ptotal = ptotal + Migdal._return_probability(self._pE2, points)
        if self._pI1 is not None:
            ptotal = ptotal + Migdal._return_probability(self._pI1, points)
        if self._pIE is not None:
            ptotal = ptotal + Migdal._return_probability(self._pIE, points)
        if self._pI2 is not None:
            ptotal = ptotal + Migdal._return_probability(self._pI2, points)

        return ptotal

    ##################################################

    def load_probabilities(self, dark_matter=False, dipole=False, double=False, e_threshold=None, excitations=False, inclusive=False, integrated=False, velocity_grid='linear'):
        """Initialise probabilities for individual orbitals."""

        if inclusive and dipole:
            print('migdal::load_probabilities: dipole and inclusive are incompatible options.')
            return False

        if inclusive:
            if (e_threshold is None or e_threshold < 0.56):
                print('migdal::load_probabilities: e_threshold must be at least 0.56 keV for semi-inclusive probabilities.')
                return False
            excitations = False
            double = False

        # Single excitation
        if excitations:
            self._pE1_orbital = {}
            for orbital in self.orbitals:
                self._pE1_orbital[orbital[0]] = Migdal._interpolate('{0}/excitations/{0}_{1}_single-excitation.txt'.format(self.element,orbital[0]), 'E', required=False)

        # Dipole single ionisation
        if dipole:
            pI1_orbital = {}
            for orbital in self.orbitals:
                pI1_orbital[orbital[0]] = self._interpolate('{0}/single-ionisation/{0}_{1}_dipole.txt'.format(self.element,orbital[0]), 'SID', integrated=integrated, e_threshold=e_threshold)
            if integrated:
                self._pI1_dipole_orbital = pI1_orbital
            else:
                self._dpI1_dipole_orbital = pI1_orbital
            return True

        # Single ionisation
        pI1_orbital = {}
        for orbital in self.orbitals:
            if inclusive:
                pI1_orbital[orbital[0]] = self._interpolate('{0}/semi-inclusive/{0}_{1}_semi-inclusive.txt'.format(self.element,orbital[0]), 'SI', integrated=integrated, e_threshold=e_threshold)
            else:
                # Logarithmic velocity grid
                if dark_matter or velocity_grid == 'log':
                    pI1_orbital[orbital[0]] = self._interpolate('{0}/single-ionisation/{0}_{1}_DM.txt'.format(self.element,orbital[0]), 'SI', integrated=integrated, e_threshold=e_threshold)
                # Linear velocity grid
                elif velocity_grid == 'linear':
                    pI1_orbital[orbital[0]] = self._interpolate('{0}/single-ionisation/{0}_{1}.txt'.format(self.element,orbital[0]), 'SI', integrated=integrated, e_threshold=e_threshold)
                else:
                    raise ValueError('migdal::load_probabilities: velocity_grid must be either "linear" or "log".')

        if integrated:
            self._pI1_orbital = pI1_orbital
        else:
            self._dpI1_orbital = pI1_orbital

        # Double ionisation
        if double:
            if integrated:
                self._pI21_orbital = {}
                self._pI2_orbital = {}
            else:
                self._dpI21_orbital = {}
                self._dpI2_orbital = {}

            for i in range(len(self.orbitals)):
                for j in range(i,len(self.orbitals)):
                    pair = self.orbitals[i][0] + self.orbitals[j][0]
                    if self.element == 'C' and pair == '2p-2p':
                        continue

                    pI2_orbital = self._interpolate('{0}/double-ionisation/{0}_{1}.txt'.format(self.element,pair), 'DI', integrated=integrated, e_threshold=e_threshold)

                    if integrated:
                        self._pI21_orbital[pair] = pI2_orbital[0]
                        self._pI2_orbital[pair] = pI2_orbital[1]
                    else:
                        self._dpI21_orbital[pair] = pI2_orbital[0]
                        self._dpI2_orbital[pair] = pI2_orbital[1]

        return True

    ##################################################

    def load_total_probabilities(self, dark_matter=False, dipole=False, double=False, e_threshold=None, excitations=False, inclusive=False, integrated=False, velocity_grid='linear'):
        """Initialise probabilities."""

        if inclusive and dipole:
            print('migdal::load_total_probabilities: dipole and inclusive are incompatible options.')
            return False

        if inclusive:
            if (e_threshold is None or e_threshold < 0.56):
                print('migdal:total_load_probabilities: e_threshold must be at least 0.56 keV for semi-inclusive probabilities.')
                return False
            excitations = False
            double = False

        # Dipole single ionisation
        if dipole:
            pI1 = Migdal._interpolate('{0}/single-ionisation/{0}_single-ionisation_dipole.txt'.format(self.element), 'SID', integrated=integrated, e_threshold=e_threshold)
            if integrated:
                self._pI1_dipole = pI1
            else:
                self._dpI1_dipole = pI1
            return True

        # No transition
        self._p0 = Migdal._interpolate('{0}/{0}_no-transition.txt'.format(self.element), 'N')

        # Excitations
        if excitations:
            # Single excitation
            self._pE1 = Migdal._interpolate('{0}/excitations/{0}_single-excitation.txt'.format(self.element), 'E', required=False)
            if double:
                # Double excitation
                self._pE2 = Migdal._interpolate('{0}/excitations/{0}_double-excitation.txt'.format(self.element), 'E', required=False)

                # Single ionisation + excitation
                pIE = Migdal._interpolate('{0}/excitations/{0}_excitation+ionisation.txt'.format(self.element), 'EI', integrated=integrated, e_threshold=e_threshold, required=False)
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
            elif velocity_grid == 'linear':
                pI1 = Migdal._interpolate('{0}/single-ionisation/{0}_single-ionisation.txt'.format(self.element), 'SI', integrated=integrated, e_threshold=e_threshold)
            else:
                raise ValueError('migdal::load_total_probabilities: velocity_grid must be either "linear" or "log".')

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
                self._dpI21 = pI2[0]
                self._dpI2 = pI2[1]

        return True

    ##################################################
    
    @staticmethod
    def _interpolate(filename, process, integrated=False, e_threshold=None, required=True):
        """Read Migdal table and return interpolating function."""

        filename = os.path.join(Migdal.data_dir, filename)
        if not os.path.exists(filename):
            if required:
                raise FileNotFoundError('migdal::_interpolate: Missing file {}'.format(filename))
            print('migdal::_interpolate: Missing file {}'.format(filename))
            return None

        if process == 'SID':
            process = 'SI'
            dipole = True
        else:
            dipole = False

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
                        if (process =='SI' or process == 'EI') and ipt % nv == 0:
                            energies[ien] = float(line[0])
                            ien += 1
                        elif process == 'DI' and ipt % (ne*nv) == 0:
                            energies[ien] = float(line[0])
                            ien += 1

                        # Velocities
                        if ipt < nv:
                            if process == 'N' or process == 'E':
                                ln_velocities[ipt] = np.log(float(line[0]))
                            elif process == 'SI' or process == 'EI':
                                ln_velocities[ipt] = np.log(float(line[1]))
                            elif process == 'DI':
                                ln_velocities[ipt] = np.log(float(line[2]))

                        # Probabilities
                        if process == 'N' or process == 'E':
                            dPdE[ipt] = float(line[1])
                        elif process == 'SI' or process == 'EI':
                            dPdE[ipt] = float(line[2])
                        elif process == 'DI':
                            dPdE[ipt] = float(line[3])
                        ipt += 1

                    # Set-up arrays
                    elif 'energies' in line[2]:
                        ne = int(line[4])
                        energies = np.zeros(ne)
                    
                    elif 'velocities' in line[2]:
                        nv = int(line[4])
                        ln_velocities = np.zeros(nv)

                        if process == 'N' or process == 'E':
                            dPdE = np.zeros(nv)
                        elif process == 'SI':
                            dPdE = np.zeros(ne * nv)
                        elif process == 'EI':
                            dPdE = np.zeros(ne * nv)
                        elif process == 'DI':
                            dPdE = np.zeros(ne * ne * nv)
                        ipt = 0
                        ien = 0

                except (IndexError, ValueError):
                    continue

        # Interpolating functions
        # No transition, excitation
        if process == 'N' or process == 'E':
            dPdE = np.log(dPdE)
            return CubicSpline(ln_velocities, dPdE, extrapolate=False)

        ln_energies = np.log(energies)

        # Integral limits
        if integrated or process == 'DI':
            ln_emin = ln_energies[0]
            ln_emax = ln_energies[-1]
            if e_threshold is not None and e_threshold >= energies[-1]:
                raise ValueError('migdal::_interpolate: e_threshold must be less than the maximum energy of the data ({:g} keV).'.format(energies[-1]))
            elif e_threshold is None or e_threshold <= energies[0]:
                e_threshold = None
                ln_eth = ln_emin
            else:
                ln_eth = np.log(e_threshold)

        # Interpolating functions for differential probabilities
        if not integrated:

            # Single ionisation dipole
            if dipole:
                # Velocity dependence is exactly v^2
                spline = CubicSpline(ln_energies, np.log(dPdE), extrapolate=False)
                return lambda lnx : 2*lnx[:,1] - 2*ln_velocities[0] + spline(lnx[:,0])

            # Single ionisation
            elif process == 'SI' or process == 'EI':
                dPdE = np.log(dPdE)
                dPdE = dPdE.reshape((ne,nv))
                return Migdal._tensor_spline((ln_energies, ln_velocities), dPdE)

            # Double ionisation
            elif process == 'DI':
                dPdE = dPdE.reshape((ne,ne,nv))

                # Double differential
                interp2 = Migdal._tensor_spline((ln_energies, ln_energies, ln_velocities), np.log(dPdE))

                # Single differential for one electron above threshold
                if e_threshold is not None:
                    # Integrate E2*dP/dE1dE2 over ln(E2) below threshold
                    integrand = energies[np.newaxis,:,np.newaxis] * dPdE
                    integral = 2 * CubicSpline(ln_energies, integrand, axis=1).integrate(ln_emin, ln_eth)
                    interp21 = Migdal._tensor_spline((ln_energies, ln_velocities), np.log(integral))
                else:
                    # No threshold: zero probability (ln P = -inf) inside grid, NaN outside
                    zero = Migdal._tensor_spline((ln_energies, ln_velocities), np.zeros((ne,nv)))
                    interp21 = lambda lnx : zero(lnx) - np.inf

                return(interp21, interp2)

        # Interpolating functions for integrated probabilities
        else:

            # Single ionisation dipole
            if dipole:
                integrand = energies * dPdE
                integral = np.log(CubicSpline(ln_energies, integrand).integrate(ln_eth, ln_emax))
                return lambda lnv : 2*lnv - 2*ln_velocities[0] + integral

            # Single ionisation
            elif process == 'SI' or process == 'EI':
                dPdE = dPdE.reshape((ne,nv))

                # Integrate E*dP/dE over ln(E) above threshold
                integrand = energies[:,np.newaxis] * dPdE
                integral = np.log(CubicSpline(ln_energies, integrand, axis=0).integrate(ln_eth, ln_emax))

                return CubicSpline(ln_velocities, integral, extrapolate=False)

            # Double ionisation
            elif process == 'DI':
                dPdE = dPdE.reshape((ne,ne,nv))
                integral = np.zeros((2,nv))

                # Integrand E1*E2*dP/dE1dE2 for integration over ln(E1), ln(E2)
                integrand = dPdE * np.outer(energies, energies)[:,:,np.newaxis]

                for i in range(nv):
                    spline = RectBivariateSpline(ln_energies, ln_energies, integrand[:,:,i], kx=3, ky=3)

                    # One electron above threshold
                    if e_threshold is not None:
                        integral[0,i] = np.log(2 * spline.integral(ln_emin, ln_eth, ln_eth, ln_emax))

                    # Both electrons above threshold
                    integral[1,i] = np.log(spline.integral(ln_eth, ln_emax, ln_eth, ln_emax))

                if e_threshold is not None:
                    interp21 = CubicSpline(ln_velocities, integral[0], extrapolate=False)
                else:
                    # No threshold: zero probability (ln P = -inf) inside grid, NaN outside
                    zero = CubicSpline(ln_velocities, np.zeros(nv), extrapolate=False)
                    interp21 = lambda lnv : zero(lnv) - np.inf

                return (interp21, CubicSpline(ln_velocities, integral[1], extrapolate=False))

    ##################################################

    @staticmethod
    def _tensor_spline(grids, values, k=3):
        """Return tensor-product interpolating spline on a regular grid (NaN outside grid)."""

        # Successive 1D not-a-knot interpolation along each axis
        c = values
        knots = []
        for axis, x in enumerate(grids):
            spline = make_interp_spline(x, c, k=k, axis=axis)
            knots.append(spline.t)
            c = np.moveaxis(spline.c, 0, axis)

        return NdBSpline(tuple(knots), c, k, extrapolate=False)

    ##################################################

    @staticmethod
    def _return_probability(fn, points, orbital=None):

        if orbital is not None:
            if fn is None:
                print('Call load_probabilities() method to initialise.')
                return None

            if orbital not in fn:
                print("'{}' is not a valid orbital.".format(orbital))
                return None
            
            if fn[orbital] is None:
                print("Probabilities not available for orbital '{}'.".format(orbital))
                return None

            return np.exp(fn[orbital](points))

        else:
            if fn is None:
                print('Call load_total_probabilities() method to initialise.')
                return None

            return np.exp(fn(points))

    ##################################################

#*******************************************************************************