#!/usr/bin/env python
# coding: utf-8

#WARNING: This is an early version of this program. It is limited in scope, unoptimized, and is too slow to compute anything meaningful. If you want to do this computation, use a C++ version.
import sys
import multiprocessing as mp
import random
import time

start = time.time()

##This code should compute nonharmonic cusp forms over Q(sqrt(-1))

bits_of_precision = 100

RR = RealField(bits_of_precision)
CC = ComplexField(bits_of_precision)

tolerance = RR(10^(-20))

debug = True
parallel = True

###############################
#
#Define the quaternions
#
###############################

def next_int(x):
    return floor(x) + 1
    
def complex_point_plot(pts,colorstr='blue'): 
    """ 
    A function that returns a plot of a list of complex points. 
    Arguments: pts (a list of complex numbers) 
    Outputs: A list plot of the imaginary numbers 
    """ 
    return list_plot([(real(k), imag(k)) for k in pts], axes_labels = ['Re($z$)', 'Im($z$)'], size=30,color=colorstr, aspect_ratio=1)

class Quaternion:
    def __init__(self,x,y,z,w):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        
    def __str__(self):
        string = str(self.x)
        string += " + "
        string += str(self.y)
        string += "*i + "
        string += str(self.z)
        string += "*j"
        if self.w == RR(0):
            return string
        else:
            string += " + "
            string += str(self.w)
            string += "*ij"
            return string
        
    def __mul__(self,other):
        c1 = self.x*other.x - self.y*other.y - self.z*other.z - self.w*other.w
        c2 = self.x*other.y + self.y*other.x + self.z*other.w - self.w*other.z
        c3 = self.x*other.z - self.y*other.w + self.z*other.x + self.w*other.y
        c4 = self.x*other.w + self.y*other.z - self.z*other.y + self.w*other.x
        return Quaternion(c1,c2,c3,c4)
    
    def __add__(self,other):
        return Quaternion(self.x + other.x, self.y + other.y, self.z + other.z, self.w + other.w)
    
    def __neg__(self):
        return Quaternion(-self.x, - self.y, -self.z, -self.w)
    
    def __sub__(self,other):
        return self + -other
    
    def act(self,matrix):
        a = CC(matrix[0][0])
        b = CC(matrix[0][1])
        c = CC(matrix[1][0])
        d = CC(matrix[1][1])
        
        assert (a*d-b*c - 1).abs() < 10^-5
        
        a_quat = Quaternion(a.real(), a.imag(), 0)
        b_quat = Quaternion(b.real(), b.imag(), 0)
        c_quat = Quaternion(c.real(), c.imag(), 0)
        d_quat = Quaternion(d.real(), d.imag(), 0)
        
        return (a_quat*self + b_quat)/(c_quat*self + d_quat)
    
    def abs(self):
        result = self*self.conj()
        return sqrt(result.x)
    
    def conj(self):
        return Quaternion(self.x, -self.y, -self.z, -self.w)
    
    def invert(self):
        denom = (self*self.conj()).x
        conj = self.conj()
        return Quaternion(conj.x/denom, conj.y/denom, conj.z/denom, conj.w/denom)
    
    def __truediv__(self,other):
        return self*other.invert()
    
    def getreal(self):
        return self.x
    
    def geti(self):
        return self.y
    
    def getj(self):
        return self.z
    
    def getk(self):
        return self.w
    
    def getcomplex(self):
        return CC(self.x + I*self.y)
        
    def get_approx(self):
        return Quaternion(RR(self.x),RR(self.y),RR(self.z),RR(self.w))

class Index:
    def __init__(self,a,b,d):
        self.a = a
        self.b = b
        self.d = d
        assert d in [1,2,3,7,11]
        self.theta = self._get_theta()
    
    def __str__(self):
        string = str(self.a) + " + " + str(self.b) + "*theta"
        return string
    
    __repr__ = __str__
        
    def __eq__(self, other):
        return (isinstance(other, type(self))) and self.__key() == other.__key()
        
    def __key(self):
        return (self.a, self.b, self.d)
    
    def __hash__(self):
        return hash(self.__key())
        
    def _get_theta(self):
        if not is_squarefree(-self.d):
            raise Exception("-d should be squarefree")
        elif (-self.d) % 4 == 1:
            return (1 + sqrt(-d))/2
        else:
            return sqrt(-d)
        
    def get_complex(self):
        return CC(self.a + self.b*self.theta)
    
    def real(self):
        return RR(self.get_complex().real())
        
    def imag(self):
        return RR(self.get_complex().imag())
    
    def abs(self):
        return RR((self.get_complex()).abs())
        
    def get_tuple(self):
        return tuple([self.a,self.b])
        
    def angle(self):
        phi = arctan2(self.imag(), self.real())
        if phi < 0:
            return RR(phi + 2*pi)
        else:
            return phi
        
    #Rotates self by multiplying by a unit
    #specifically the unit which generates the unit group
    #and has minimal angle wrt the positive real axis
    #returns it as a new Index object
    def rotate(self):
        if self.d == 1:
            #rotate self by multiplication by i
            return Index(-self.b, self.a, self.d)
        elif self.d == 3:
            #rotate self by multiplciation by (1+sqrt(-3))/2
            return Index(-self.b, self.a + self.b, self.d)
        else:
            #rotate self by multiplication by -1
            return Index(-self.a, -self.b, self.d)
            
    def conj(self):
        if -self.d % 4 == 1:
            return Index(self.a + self.b, -self.b, self.d)
        else:
            return Index(self.a, -self.b, self.d)

class CoefficientComputer:
    def __init__(self, d, D, sym_class):
        #Set independent parameters
        self.d = d
        assert d in [1,2,3,7,11]
        self.D = RR(D)
        self.sym_class = sym_class
        assert sym_class in ['D', 'G', 'C', 'H']
        self.r = RR(self._known_eigenvalue())
        
        self.A = self._get_A()
        self.theta = self._get_theta()
        self.Y0 = self._get_Y0()
        self.M0 = self._compute_M(self.Y0)
        
        self.Y = RR(self.r/(2*pi/self.A*self.M0))
        assert self.Y < self.Y0
        self.MY = self._compute_M(self.Y)
        
        self._indices = []
        self._index_transversal = []
        self._index_of_one = None
        self._index_symmetry_data = {}
        self._test_point_dict = {}
        
        print("Generating indices and symmetry data.")
        self._generate_index_data()
        print("Generating test points.")
        self._generate_test_point_data()
        print("Coefficient Computer Class initialized.")
        
        #Declare the matrix and final answer to be computed later
        self._matrix = None
        self._coeff_dict = None
        
        #Report the results of initializing this class
        if debug == True:
            number_of_Bessel_Ks = 0
            for row in self._index_transversal:
                for col in self._index_transversal:
                    if (row == col):
                        number_of_Bessel_Ks += 1
                    
                    number_of_Bessel_Ks += len(self._test_point_dict[row])
        
            message = "D is " + str(self.D) + ".\n"
            message += "M0 is " + str(self.M0) + ".\n"
            message += "MY is " + str(self.MY) + ".\n"
            message += "Y is " + str(self.Y) + ".\n"
            message += "Y0 is " + str(self.Y0) + ".\n"
            message += "r is " + str(self.r) + ".\n"
            message += "d is " + str(self.d) + ".\n"
            message += "A is " + str(self.A) + ".\n"
            message += "theta is " + str(self.theta) + ".\n"
            message += "sym_class is " + str(self.sym_class) + ".\n"
            message += "There are " + str(len(self._indices)) + " indices up to M0.\n"
            message += "There are beetween " + str(min([len(self._test_point_dict[m]) for m in self._test_point_dict]))
            message += " and " + str(max([len(self._test_point_dict[m]) for m in self._test_point_dict])) + " test points.\n"
            message += "The matrix will have " + str(len(self._index_transversal)) + " rows and columns.\n"
            message += "That's " + str(len(self._index_transversal)^2) + " total entries.\n"
            message += "Have to compute " + str(number_of_Bessel_Ks) + " K-Bessel functions.\n"
            print(message)
            sys.stdout.flush()
                
    
    ######These methods don't do any meaningful comptuations
    ######They just contain data that is clumsy to hardcode in the constructor
    
    def _known_eigenvalue(self):
        S = self.sym_class
        if self.d == 1:
            if S == 'D':
                return 8.555250
            elif S == 'G':
                return 17.45131992
            elif S == 'C':
                return 6.622118
            else:
                #S == H
                return 12.11527484
        elif self.d == 2:
            if S == 'D':
                return 5.242134
            elif S == 'C':
                return 4.944483
            else:
                raise Exception("Eigenvalue not known.")
        elif self.d == 3:
            if S == 'D':
                return 12.500100
            elif S == 'C':
                return 7.072007
            else:
                raise Exception("Eigenvalue not known.")
        elif self.d == 7:
            if S == 'D':
                return 7.419980
            elif S == 'C':
                return 5.063246
            else:
                raise Exception("Eigenvalue not known.")
        else:
            #d == 11
            if S == 'D':
                return 4.859365
            elif S == 'C':
                return 3.873897
            else:
                raise Exception("Eigenvalue not known.")
    
    def _get_Y0(self):
        if self.d == 3:
            return RR(sqrt(1-0^2-(1/(sqrt(3)*4)+sqrt(3)/4)^2))
        else:
            return RR(sqrt(1- (1/2)^2 - (self.theta.imag()/2)^2))
    
    def _get_A(self):
        if not is_squarefree(-self.d):
            raise("-d should be squarefree")
        elif (-self.d) % 4 == 1:
            return RR(sqrt(self.d)/2)
        else:
            return RR(sqrt(self.d))
            
    def _get_theta(self):
        if not is_squarefree(-self.d):
            raise Exception("-d should be squarefree")
        elif (-self.d) % 4 == 1:
            return (1 + sqrt(-self.d))/2
        else:
            return sqrt(-self.d)
    
    ######These methods compute indices and index symmetry data
    ######
    
    def _compute_M(self,localY):
        tail_indices_dict = {}
        modest_implied_constant = 1
        
        start = time.time()
        #binary search to find where the terms become as small as the tolerance
        #This will be the last term of the computed tail, considered zero
        start_n = ceil(self.r/(2*pi/self.A*localY))
        max_n = 100
        
        left = RR(start_n)
        right = RR(max_n)
        
        assert modest_implied_constant*left*localY*self._kappa((2*pi/self.A)*left*localY)
        while modest_implied_constant*right*localY*self._kappa((2*pi/self.A)*right*localY) >= tolerance:
            left = right
            right *= 2
        
        center = (right+left)/2
        while (right-left).abs() > 0.5:
            test = modest_implied_constant*center*localY*self._kappa((2*pi/self.A)*center*localY)
            if test >= tolerance:
                left = center
            else:
                #test < tolerance
                right = center
            center = (right+left)/2
        
        max_n = ceil(right)
        
        if debug == True:
            message = "Last term of tail found\n"
            message += "the absolute value of n is " + str(max_n) + "\n"
            message += "it evalutes to " + str(RR(modest_implied_constant*max_n*localY*self._kappa((2*pi/self.A)*max_n*localY)))
            print(message)
        
        
        #Compute the tail from the beginning of exponential growth until the terms become zero
        theta_real = self.theta.real()
        theta_modulus = self.theta.abs()
        
        aub = ceil(max_n*theta_modulus/sqrt(theta_modulus^2 - theta_real^2))
        alb = -aub
        bub = ceil(max_n/sqrt(theta_modulus^2 - theta_real^2))
        blb = -bub
        
        for a in range(alb,aub+1):
            for b in range(blb, bub+1):
                index = Index(a,b,self.d)
                if index.abs() > max_n or index.abs() == 0 or index.abs() < start_n:
                    continue
                if index.abs() in tail_indices_dict.keys():
                    tail_indices_dict[index.abs()].append(index)
                else:
                    tail_indices_dict[index.abs()] = [index]
        
        tail = RR(0)
    
        for n_modulus in tail_indices_dict:
            multiplicity = len(tail_indices_dict[n_modulus])
            term = RR(multiplicity)
            term *= RR(modest_implied_constant*n_modulus*localY)
            term *= RR(self._kappa((2*pi/self.A)*n_modulus*localY))
            tail += term
        
        #subtract off terms starting from the beginning until the tail reaches the desired truncation level
        accuracy = 10^(-self.D)
        
        M = None
        for bound in sorted(tail_indices_dict):
            M = bound
            if tail < accuracy:
                break
            else:
                multiplicity =  len(tail_indices_dict[M])
                tail -= RR(multiplicity*modest_implied_constant*M*localY*self._kappa((2*pi/self.A)*M*localY))
        return RR(M)
        
    def _generate_index_data(self):
        theta_real = self.theta.real()
        theta_modulus = self.theta.abs()
        
        aub = ceil(self.M0*theta_modulus/sqrt(theta_modulus^2 - theta_real^2))
        alb = -aub
        bub = ceil(self.M0/sqrt(theta_modulus^2 - theta_real^2))
        blb = -bub
        for a in range(alb, aub+1):
            for b in range(blb, bub+1):
                if a == 0 and b == 0:
                    continue
                index = Index(a, b, self.d)
                if index.abs() <= self.M0:
                    self._indices.append(index)
        
        #Sort the indices by absolute value then by angle
        #This makes computing the orbit representatives trivial
        self._indices.sort(key=lambda x: (x.abs(), x.angle()))
            
        #if debug == True:
        #    var("x y")
        #    #plotMY = scatter_plot([[index.a, index.b] for index in self._indices], aspect_ratio = 1)
        #    #plotMY += implicit_plot(x^2 + 2*x*y*theta_real + theta_modulus^2*y^2 - self.MY^2, (x,alb, aub), (y, blb, bub), aspect_ratio = 1)
        #    plotM0 = scatter_plot([index.get_tuple() for index in self._indices], aspect_ratio = 1)
        #    plotM0 += implicit_plot(x^2 + 2*x*y*theta_real + theta_modulus^2*y^2 - self.M0^2, (x,alb, aub), (y, blb, bub), aspect_ratio = 1)
            #save(plotMY, 'indicesMY.png')
        #    save(plotM0, 'indicesM0.png')
        
        
        already_found = []
        
        rotation_coeff = None
        conj_coeff = None
        
        if self.sym_class in ['D','G']:
            rotation_coeff = 1
        else:
            rotation_coeff = -1
        
        if self.sym_class in ['D','C']:
            conj_coeff = 1
        else:
            conj_coeff = -1
        
        #partition into orbits
        #choose a representative from each orbit
        for index in self._indices:
            if index in already_found:
                continue
            
            #start building the orbit by appending rotations of the starting index
            #each element in the orbit is a pair consiting of the index and the relation
            #between the corresonding coefficients
            orbit = [[index,1]]
            
            #rotate the last index in the orbit
            temp = orbit[-1][0].rotate()
            while not temp == index:
                orbit.append([temp, orbit[-1][1]*rotation_coeff])
                temp = temp.rotate()
            
            #indicate whether to append rotations of the conjugate
            #The answer is "no" if the conjugate of any index is already in the orbit
            #and "yes" otherwise
            more_indices_in_orbit = not bool(index.conj() in [itr for [itr,rel] in orbit])
            
            if more_indices_in_orbit == True:
                #extend the orbit with the conjugates of the rotations
                to_extend = []
                for [ind, coeff] in orbit:
                    to_extend.append([ind.conj(),coeff*conj_coeff])
                orbit.extend(to_extend)
                
            #indicate that these indices have been collected into their orbit
            for [ind, coeff] in orbit:
                already_found.append(ind)
            
            self._index_transversal.append(index)
            
            for [itr,rel] in orbit:
                self._index_symmetry_data[itr] = orbit
            
            if debug == True:
                M = orbit[0][0].abs()
                for [itr,rel] in orbit:
                    assert (M-itr.abs()).abs() < tolerance
            
        #Find the index of 1 + 0*I and set the private value
        for itr in range(len(self._index_transversal)):
            if self._index_transversal[itr] == Index(1,0,self.d):
                self._index_of_one = itr
                break
                
        if debug == True:
            #for key in self._index_transversal:
            #    plot = complex_point_plot([point[0].get_complex() for point in self._index_symmetry_data[key]])
            #    plot.save_image('orbit_' + str(key) + '.png')
            #plot = complex_point_plot([point.get_complex() for point in self._index_transversal])
            #plot.save_image('transversal.png')
            
            message = "These are the orbits and coefficient relations of the indices."
            print(message)
            for itr in self._index_transversal:
                print("-->" + str(self._index_symmetry_data[itr]))
    
    
    ######These methods compute test points and their reduction
    ######
    
    def _is_in_parabolic_fund_dom(self,point):
        check = bool(abs(point.x) <= 1/2)
        
        if self.d == 3:
            check *= bool(point.y <= (-1/sqrt(3)*(point.x - 1/4) + sqrt(3)/4))
            check *= bool(point.y >= (-1/sqrt(3)*(point.x + 1/4) - sqrt(3)/4))
            check *= bool(point.y <= (1/sqrt(3)*(point.x + 1/4) + sqrt(3)/4))
            check *= bool(point.y >= (1/sqrt(3)*(point.x - 1/4) - sqrt(3)/4))
            return check
                
        check *= bool(abs(point.y) <= abs(self.theta.imag()/2))
        
        return check
        
    def _is_in_units_fund_dom(self,point):
        if self.d == 1:
            return bool(point.x >= 0)
        elif self.d == 3:            
            return bool(point.x >= 0) and bool(point.y >= (-1/sqrt(3)*point.x))
        else:
            return True
            
    def _is_in_fund_dom(self, point):
        #Check that the point is in the upper-half space
        assert point.w == 0
        assert point.z > 0
        
        parabolic_check = self._is_in_parabolic_fund_dom(point)
        
        units_check = self._is_in_units_fund_dom(point)
        
        inversion_check = bool(point.abs() >= 1)
        
        return parabolic_check*units_check*inversion_check
        
    def _reduce(self,point):
        one = Quaternion(1,0,0,0)
        quattheta = Quaternion(self.theta.real(), self.theta.imag(), 0, 0)
        imtheta = self.theta.imag()
        temp = point
        
        if self.d == 3:
            othertheta = Quaternion(-1/2, sqrt(3)/2, 0,0)
            while(not self._is_in_fund_dom(temp)):
                if temp.abs() < 1:
                    temp = (-one)/temp
                while temp.y > (-1/sqrt(3)*(temp.x - 1/4) + sqrt(3)/4):
                    temp -= quattheta
                while temp.y < (-1/sqrt(3)*(temp.x + 1/4) - sqrt(3)/4):
                    temp += quattheta
                while temp.y > (1/sqrt(3)*(temp.x + 1/4) + sqrt(3)/4):
                    temp -= othertheta
                while temp.y < (1/sqrt(3)*(temp.x - 1/4) - sqrt(3)/4):
                    temp += othertheta
                while temp.x > 1/2:
                    temp -= one
                while temp.x < -1/2:
                    temp += one
                while bool(temp.x < 0) or bool(temp.y < (-1/sqrt(3)*temp.x)):
                    temp = quattheta*temp*quattheta
            return temp
        
        while(not self._is_in_fund_dom(temp)):
            if temp.abs() < 1:
                temp = (-one)/temp
            while temp.y > imtheta/2:
                temp -= quattheta
            while temp.y < -imtheta/2:
                temp += quattheta
            while temp.x > 1/2:
                temp -= one
            while temp.x < -1/2:
                temp += one
            if self.d == 1:
                if temp.x < 0:
                    imagunit = Quaternion(0,1,0,0)
                    temp = imagunit*temp*imagunit
        return temp
    
    def _generate_test_point_data(self):
        if self._index_transversal == None or len(self._index_transversal) == 0:
            raise Exception("Tried to compute test points without computing index data first.")
        
        M = self.MY
        
        for m in self._index_transversal:
            points = []
            P = None
            Q = None
            if -self.d % 4 == 1:
                P = next_int(2*(m.imag() + 2*self.MY/sqrt(self.d)))
                Q = next_int(2*(m.real() + m.imag()/2 + self.MY))
            else:
                P = next_int(2*(m.imag() + self.MY/sqrt(self.d)))
                Q = next_int(2*(m.real() + self.MY))
                
            #if -self.d % 4 == 1:
            #    P = next_int(m.abs() + 2*M/sqrt(self.d))
            #    Q = next_int(m.abs() + M)
            #else:
            #    P = next_int(m.abs() + M/sqrt(self.d))
            #    Q = next_int(m.abs()+ M)
        
            imtheta = self.theta.imag()
        
            for l0 in range(0,P):
                for l1 in range(0,Q):
                    point = Quaternion(-1/2 + 1/(2*P) + l0/P, -imtheta/2 + imtheta/(2*Q) + imtheta*l1/Q, self.Y, 0)
                    points.append(point)
        
            point_and_pullback_pairs = []
            for point in points:
                pullback = self._reduce(point)
                point_and_pullback_pairs.append(tuple([point,pullback]))
            
            if debug == True:
                plot = point3d([(pair[0].getreal(),pair[0].geti(),pair[0].getj()) for pair in point_and_pullback_pairs],size=10,color="blue",plot_points = 100, aspect_ratio = 1)
                plot += point3d([(pair[1].getreal(),pair[1].geti(),pair[1].getj()) for pair in point_and_pullback_pairs],size=5,color="black", plot_points=100, aspect_ratio = 1)
                plot.save_image('test_points.png')
            
            self._test_point_dict[m] = point_and_pullback_pairs
    
    ######These methods are used to compute the matrix entries
    ######
    
    def _kappa(self, arg):
        besselorder = CC(self.r*I)
        return exp(pi*self.r/2)*bessel_K(besselorder, RR(arg)).real()

    def _trace_product(self,z,w):
        answer = z*w + (z*w).conjugate()
        return answer.real()

    def _produce_entry(self, L):
        
        m = L[0]
        n = L[1]
        #m is the row
        #n is the column
        
        entry = CC(0)
        
        if m == n:
            entry += self.Y*self._kappa((2*pi/self.A)*m.abs()*self.Y)
        
        point_and_pullback_pairs = self._test_point_dict[m]
        number_of_points = len(point_and_pullback_pairs)
        assert number_of_points > 0
        
        for [point, pullback] in point_and_pullback_pairs:
            x = point.getcomplex()
            y = point.getj()
            xstar = pullback.getcomplex()
            ystar = pullback.getj()
            
            term = -(self.theta.imag())/number_of_points
            term *= ystar
            term *= self._kappa((2*pi/self.A)*n.abs()*ystar)
            
            term *= exp((-pi*I/self.A)*self._trace_product(I*m.get_complex(),x))
            
            exp_term = CC(0)
            for [l, pm1] in self._index_symmetry_data[n]:
                trace = self._trace_product(I*l.get_complex(),xstar)
                exp_term += pm1*exp((pi*I/self.A)*trace)
            term *= exp_term
                     
            entry += CC(term)
        #print("here5",flush=True)
        if debug == True:
            message = "Entry " + str(L[2]) + " done. "
            print(message)
        
        return CC(entry)
    
    ######Matrix manipulation method
    ######
    
    def _remove_row_col(self,mat,rows,cols):
        R = range(mat.dimensions()[0])
        C = range(mat.dimensions()[1])
        return mat[[k for k in R if k not in rows], [k for k in C if k not in cols]]
    
    
    ######Methods for finding an eigenvalue
    ######
    
    def _cost(self):
        v = self._matrix[0]
        a = [self._coeff_dict[index] for index in self._index_transversal]
        answer = 0
        
        for [x,y] in zip(v,a):
            answer += x*y
        
        return CC(answer).real()
        
        
    
    ########################################
    #
    #Public methods
    #
    ########################################
    
    def generate_matrix(self):
        arguments = []

        count = 1
        for m in self._index_transversal:
            for n in self._index_transversal:
                arguments.append([m,n,count])
                count += 1
            
        prematrix = []
        
        if parallel == True:
            pool = mp.Pool(mp.cpu_count())
            #print("Beginning to compute matrix entries.",flush=True)
            prematrix = pool.map(self._produce_entry,arguments)
            #print("Multiprocessing pool finished.",flush=True)
            pool.close()
            pool.join()
            #print("Matrix entries generated.")
            sys.stdout.flush()
        else:
            print("Beginning to compute matrix entries.",flush=True)
            prematrix = list(map(self._produce_entry, arguments))
            print("Matrix entries generated.", flush=True)
        
        almostmatrix = []
        
        for m in self._index_transversal:
            row = []
            for n in self._index_transversal:
                row.append(prematrix.pop(0))
            almostmatrix.append(row)
        
        assert len(prematrix) == 0
        
        self._matrix = Matrix(almostmatrix)
    
        # print("Private matrix set.", flush=True)
        
    def compute_coefficients(self):
        if self._matrix == None:
            raise Exception("Matrix not yet generated.")
        
        temp_mat = self._matrix
        
        b = [temp_mat[k][self._index_of_one] for k in range(temp_mat.dimensions()[0])]
        b = Matrix(b)
        
        temp_mat = self._remove_row_col(temp_mat,[self._index_of_one],[self._index_of_one])
        b = self._remove_row_col(b.transpose(),[self._index_of_one],[])
        
        b = -vector(b)
        
        preanswer = temp_mat.solve_right(b)
        preanswer = [p.real() for p in preanswer]
        
        answer = []
        flag = False
        for itr in range(len(self._index_transversal)):
            if itr == self._index_of_one:
                answer.append(CC(1))
                flag = True
            elif flag == False:
                answer.append(preanswer[itr])
            else:
                answer.append(preanswer[itr-1])
        
        self._coeff_dict = dict(zip(list(self._index_transversal),answer))
        for rep_index in self._index_transversal:
            symmetry_class = self._index_symmetry_data[rep_index]
            for itr in range(1,len(symmetry_class)):
                [index, pm1] = symmetry_class[itr]
                self._coeff_dict[index] = self._coeff_dict[rep_index]*pm1
        
    def get_coeff_dict(self):
        if self._coeff_dict == None:
            raise Exception("Coefficients not yet solved for.")
        return self._coeff_dict

    def get_matrix(self):
        return self._matrix
        
    def run_test(self):
        if self._coeff_dict == None:
            raise Exception("Coefficients not yet solved for.")
        
        test_point_pullback_pairs = []
        while len(test_point_pullback_pairs) < 10:
            x0 = random.random() - 0.5
            x1 = (random.random() - 0.5)*self.theta.imag()/2
            if self.d == 3:
                x1 = (random.random() - 0.5)*.6
            y = self.Y  + (1 - self.Y)*random.random()
            
            point = Quaternion(x0,x1,y,0)
            if self._is_in_parabolic_fund_dom(point) and (not self._is_in_fund_dom(point)) and point.getj() >= self.Y0:
                test_point_pullback_pairs.append([point.get_approx()])
        
        for itr in range(len(test_point_pullback_pairs)):
            p = test_point_pullback_pairs[itr][0]
            test_point_pullback_pairs[itr].append(self._reduce(p).get_approx())
            
        print("Starting test.")
        sys.stdout.flush()
        max_diff = 0
        for [point, pullback] in test_point_pullback_pairs:
            output_point = self.evaluate(point)
            output_pullback = self.evaluate(pullback)
            diff = (output_point - output_pullback).abs()
            if diff > max_diff:
                max_diff = diff
            message = "([" + str(diff) + "], " + str(point) + ", " + str(pullback) + ")\n"
            print(message)
            sys.stdout.flush()
        message = "Maximum difference between a point an its pullback is: " + str(max_diff) + "\n"
        digits1 = self.D
        digits2 = ceil(log(self.get_condition())/log(10))
        message += "Expected accurary is about " + str(digits1 - digits2) + " digits."
        
        print(message)
         
    def evaluate(self, z):
        if self._coeff_dict == None:
            raise Exception("Coefficients not yet solved for.")
            
        x = z.getcomplex()
        y = z.getj()
        
        if y < self.Y:
            print("WARNING: POINT ACCURACY NOT GUARANTEED",flush=True)
                
        answer = CC(0)
        
        for n in self._coeff_dict:
            coeff = self._coeff_dict[n]

            term = coeff
            term *= y
            
            term *= self._kappa((2*pi/self.A)*n.abs()*y)
            
            trace = self._trace_product(I*n.get_complex(),x)
            term *= CC(exp((pi*I/self.A)*trace))
            
            answer += term
        return CC(answer)
        
    def get_condition(self):
        if self._matrix == None:
            raise Exception("Matrix not yet generated.")
        
        temp_mat = self._matrix
        
        temp_mat = self._remove_row_col(temp_mat,[self._index_of_one],[self._index_of_one])
        print("Determinant of matrix is " + str(det(temp_mat)) + ".")
        sys.stdout.flush()
        
        normA = 0
        for row in temp_mat:
            for entry in row:
                normA += abs(entry)^2
        normA = sqrt(normA)
        
        Ainverse = temp_mat.inverse()
        normAinverse = 0
        for row in Ainverse:
            for entry in row:
                normAinverse += abs(entry)^2
        normAinverse = sqrt(normAinverse)
        
        return normA*normAinverse

    def zoom_in(self):
        self.generate_matrix()
        self.compute_coefficients()
        x_n = self.r
        y_n = self._cost()
        print("r: " + str(x_n) + " cost: " + str(y_n))
        x_n_1 = RR(x_n) + RR(0.0001)
        self.r  = x_n_1
        self.generate_matrix()
        self.compute_coefficients()
        y_n_1 = self._cost()
        print("r: " + str(x_n_1) + " cost: " + str(y_n_1))
        
        while (x_n - x_n_1).abs() > 10^-10:
            x_n_2 = x_n - y_n *(x_n-x_n_1)/(y_n-y_n_1)
            self.r = x_n_2
            self.generate_matrix()
            self.compute_coefficients()
            y_n_2 = self._cost()
            
            x_n = x_n_1
            y_n = y_n_1
            
            x_n_1 = x_n_2
            y_n_1 = y_n_2
            print("r: " + str(x_n_1) + " cost: " + str(y_n_1))
        return
    
if __name__ == '__main__':
    print("Computation started.")
    print("There are this many cores: " + str(mp.cpu_count()))
    print("Working with " + str(bits_of_precision) + " bits of precision.")
    sys.stdout.flush()

    d = 11
    digits = 1
    
    sym_class = 'D'

    c1 = CoefficientComputer(d, digits, sym_class)
    c1.generate_matrix()
    c1.compute_coefficients()
    coeff_dict1 = c1.get_coeff_dict()

    for key in coeff_dict1:
        print("(" + str(key) + ", " + str(coeff_dict1[key]))

    c1.run_test()

    end = time.time()
    
    c1.zoom_in()



    print("Computation finished.")
    print("Time elapsed: " + str(end-start) + " seconds.")
    sys.stdout.flush()



