import re
import numpy as np

class KBesselPrecomputer:
    def __init__(self, r_spec_param, minimum_argument):
        self._RR = RealField(54)
        self._CC = ComplexField(54)
        self._r_spec_param = self._RR(r_spec_param)
        self._minimum_argument = np.float64(minimum_argument)
       
        self._fine_frequency = np.float64(1000.0) #this is good enough for plotting
        self._coarse_frequency = np.float64(100.0) #this is also good enough for plotting
        self._fine_spacing = 1.0/self._fine_frequency
        self._coarse_spacing = 1.0/self._coarse_frequency
        
        self._coarse_to_fine_boundary = np.float64(self._r_spec_param + 3)
        self._zero_K_Bessel_cutoff = np.float64(self._r_spec_param + 100)
        
        self._fine_bess_vals = np.empty((0,), dtype=np.float64)
        self._coarse_bess_vals = np.empty((0,), dtype=np.float64)
        
        self._precompute()
        
    def _temp_bess(self, x):
        return self._RR(bessel_K(self._CC(self._r_spec_param * I), self._RR(x)).real() * exp(pi * self._r_spec_param / 2))
    
    def _fine_y_to_index(self, y):
        return np.int32(np.round(self._fine_frequency * y)) - self._fine_bess_offset
    
    def _coarse_y_to_index(self, y):
        return np.int32(np.round(self._coarse_frequency * y)) - self._coarse_bess_offset

    def _lerp(self,x0,x1,y0,y1,x):
        return (y1-y0)/(x1-x0) * (x - x0) + y0
    
    def _precompute(self):
        print("Precomputing K-Bessel function. This might take a while.")

        #An argument in the "fine" range is converted to an integer-indexed list of values by multiplying by fine_frequency and rounding down
        #This is the smallest argument in that range, and this is the smallest resulting integer
        #So everything will have this value subtracted off of it for zero-indexed arrays
        x = np.round(self._minimum_argument * self._fine_frequency)
        self._fine_bess_offset = np.int32(x)

        #Similar to previous comment but for the "coarse" range
        x = np.round(self._coarse_to_fine_boundary * self._coarse_frequency)
        self._coarse_bess_offset = np.int32(x)

        #_precompute values in the "fine" range
        #y_index means the integer that a real value y maps to after being multiplied by the spacing constant
        arg = self._minimum_argument
        index = 0
        while arg < self._coarse_to_fine_boundary:
            bess_val = self._temp_bess(arg)
            self._fine_bess_vals = np.append(self._fine_bess_vals, np.float64(bess_val))
            index += 1
            arg = self._minimum_argument + index * self._fine_spacing
            
        #do two more
        for count in range(2):
            bess_val = self._temp_bess(arg)
            self._fine_bess_vals = np.append(self._fine_bess_vals, np.float64(bess_val))
            index += 1
            arg = self._minimum_argument + index * self._fine_spacing

        #_precompute values in the "coarse" range
        arg = self._coarse_to_fine_boundary
        index = 0
        while arg < self._zero_K_Bessel_cutoff:
            bess_val = self._temp_bess(arg)
            self._coarse_bess_vals = np.append(self._coarse_bess_vals, np.float64(bess_val))
            index += 1
            arg = self._coarse_to_fine_boundary + index * self._coarse_spacing
            
        for count in range(2):
            bess_val = self._temp_bess(arg)
            self._coarse_bess_vals = np.append(self._coarse_bess_vals, np.float64(bess_val))
            index += 1
            arg = self._coarse_to_fine_boundary + index * self._coarse_spacing
            
        print("K-Bessel precomputation done.")
        
    def bess(self, y):
        if y >= self._zero_K_Bessel_cutoff:
            return 0
        if self._minimum_argument <= y < self._coarse_to_fine_boundary:
            x0 = np.floor(y * self._fine_frequency) * self._fine_spacing
            x1 = np.ceil(y * self._fine_frequency) * self._fine_spacing
            y0 = self._fine_bess_vals[self._fine_y_to_index(x0)]
            y1 = self._fine_bess_vals[self._fine_y_to_index(x1)]
            if x0 == x1:
                return y0
            return self._lerp(x0,x1,y0,y1,y)
        x0 = np.floor(y * self._coarse_frequency) * self._coarse_spacing
        x1 = np.ceil(y * self._coarse_frequency) * self._coarse_spacing
        y0 = self._coarse_bess_vals[self._coarse_y_to_index(x0)]
        y1 = self._coarse_bess_vals[self._coarse_y_to_index(x1)]
        if x0 == x1:
            return y0
        return self._lerp(x0,x1,y0,y1,y)