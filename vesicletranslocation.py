


import numpy as np
import matplotlib.pyplot as plt
from math import ceil 
import warnings

def update_progress(progress):
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    
    text = "\rPercent: [{0}] {1:4.2f}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
warnings.filterwarnings("ignore")

def findMin(arr, rows, mid, min): 
  
    min_index = 0
    for i in range(rows): 
        if (min > arr[i,mid]): 
              
            # Saving global minimum and its index 
            # to check its neighbours 
            min = arr[i, mid] 
            min_index = i 
    
  
    return min,min_index 
  
# Function to find a deep element 
def findDeepRec(arr, rows, columns, mid): 
  
    # Evaluating minimum of mid column.  
    # Note min is passed by reference. 
    min = 1000000
    min, min_index = findMin(arr, rows, mid, min) 
  
    # If we are on the first or last column, 
    # min is a deep 
    if (mid == 0 or mid == columns - 1): 
        return min
  
    # If mid column minimum is also pedeepak 
    if (min <= arr[min_index, mid - 1] and 
        min <= arr[min_index, mid + 1]): 
        return min
  
    # If min is greater than its left 
    if (min > arr[min_index, mid - 1]): 
        return findDeepRec(arr, rows, columns,  
                           mid - ceil(mid / 2.0)) 
  
    # If min is greater than its left 
    # if (inx > arr[min_index, mid+1]) 
    return findDeepRec(arr, rows, columns,  
                       mid + ceil(mid / 2.0)) 
  
# A wrapper over findDeepkRec() 
def findDeep(arr, rows, columns): 
    return findDeepRec(arr, rows,  
                       columns, columns // 2) 
  




class vesicleTranslocation:
    
    
    
  

    def __init__(self, r0, beta, Fext, kappa_c =2, lambda_ =5):
        
       # self.r = r
        self.r0 = r0          # r0 is the initial radius of the vesicle before entering the pore. 
        self.d = 1            # d is the radius of pore which has the shape of a cylinder
        self.kappa_c = 1      # kapa_c is bending modulus
        self.lambda_ = 5      # Lambda is the stretching modulus
        self.beta = beta        #  is the ratio of the volume of the pore to that of the vesicle : V_p/V_0
        self.c0 = 0           # c0 is the spontaneous curvature which is a constant at all points on the vesicle 
                              # and arises from asymmetry in the areas of the inner and outer surfaces of the lipid blade constituting the membrane of the vesicle.
        self.F0 = (1./2.) * self.kappa_c * (4. * np.pi * np.square(2 - self.r0*self.c0))
        self.alphamin = (1/2)*self.d**2/self.r0**2

        self.Fext = Fext

         ##FreeEnergy
        # Our purpose is to minimize the free energy with respect to theta and find the radius of the vesicle
        # as a function of theta: r(theta). But first we need to calculate the volume of the vesicle in the donor  compartment.
        # Equation 5 in the paper gives the expression for the volume of the vesicle in the donor compartment. 
        # Since this will a polynomial of degree 3 in r, we write the coefficient of each power of r separately.
        # r_coeffs is a function that gets alpha and theta and generate the three coefficients of the cubic equation
        # we use it in the r_theta_pair  function to generate the pair.    
    def r_coeffs(self, theta, alpha):
        sin = np.sin(theta)
        one_m_sin = (1 - sin)**3
        cos = np.cos(theta)
        V_D_r0 = (4/3) * np.pi * self.r0**3 * (1 - alpha)



        V_I_sphere_r0 = 0 #coefficient of r^0 for sphere part of the vesicle
        V_I_sphere_r1 = 0 #coefficient of r^1 for sphere part of the vesicle
        V_I_sphere_r2 = 0 #coefficient of r^2 for sphere part of the vesicle
        V_I_sphere_r3 = (2/3) * np.pi * (1+cos) #coefficient of r^3 for sphere part of the vesicle

        V_I_cone_r0 = (1/3) * np.pi/one_m_sin * sin**2 * cos * (-1) * self.d**3 #coefficient of r^0 for cone part shown in the figure 2 of the paper
        V_I_cone_r1 = (1/3) * np.pi/one_m_sin * sin**2 * cos * 3 * self.d**2    #coefficient of r^1 for cone part shown in the figure 2 of the paper
        V_I_cone_r2 = (1/3) * np.pi/one_m_sin * sin**2 * cos * (-3) * self.d    #coefficient of r^2 for cone part shown in the figure 2 of the paper
        V_I_cone_r3 = (1/3) * np.pi/one_m_sin * sin**2 * cos                    #coefficient of r^3 for cone part shown in the figure 2 of the paper

        V_I_torus_r0 = np.pi/one_m_sin * self.d**3 * ((2/3) * cos - sin * (np.pi/2 - theta))  #coefficient of r^0 for cone part shown in the figure 2 of the paper
        V_I_torus_r1 = np.pi/one_m_sin * self.d**2 * sin * ((1 + 2 * sin) * (np.pi/2 - theta) - 2 * cos) #coefficient of r^1 for cone part shown in the figure 2 of the paper
        V_I_torus_r2 = np.pi/one_m_sin * self.d * sin**2 * (2 * cos - ( 2 + sin) * (np.pi/2 - theta))    #coefficient of r^2 for cone part shown in the figure 2 of the paper
        V_I_torus_r3 = np.pi/one_m_sin * sin**3 * (np.pi/2 - theta - (2/3) * cos)                        #coefficient of r^3 for cone part shown in the figure 2 of the paper

        p = []
        p.append(V_I_sphere_r3 + V_I_cone_r3 - V_I_torus_r3)
        p.append(V_I_sphere_r2 + V_I_cone_r2 - V_I_torus_r2)
        p.append(V_I_sphere_r1 + V_I_cone_r1 - V_I_torus_r1)
        p.append(V_I_sphere_r0 + V_I_cone_r0 - V_I_torus_r0 - V_D_r0)

        return p
    # we need to find the roots of the polynomial of r. But not all roots are physical therefore
    # we impose a condition on r.    
    def r_theta_pair(self, alpha, acc=100):
        theta_min = np.arcsin(self.d/self.r0)
        theta_val = np.linspace(1.57, theta_min, acc)
        #alpha_val = np.range(0.03, .5 + self.beta/2)
        ans = []

        for t in theta_val:
            roots = np.roots(self.r_coeffs(t, alpha))
            for r in roots:
                if ((abs(r.imag) <= 10**(-6)) & (r.real > self.d) & (r.real < self.r0)):
                            ans.append([r.real, t])

        return np.array(ans)

    #Next we calculate the free energy terms in the Helfrich free energy. F_I function return free energy 
    # of stretching and bending for the filling part.
    def F_I(self, alpha, r, theta):
        F0 = (1./2.) * self.kappa_c * (4. * np.pi * np.square(2 - self.r0*self.c0))

#        r_t = self.r_theta(theta, alpha)
#        r = r_t[0]
        b = (r * np.sin(theta) - self.d)/(1. - np.sin(theta))

        F_I_bend_sphere = self.kappa_c * np.pi * (1 + np.cos(theta)) * np.square(2 - self.c0 * r)
        F_I_bend_torus = self.kappa_c * np.pi * (self.c0 * (2 + self.c0 * b) * (b + self.d) * (np.pi/2 - theta) - np.cos(theta) * np.square(2 + self.c0 * b)
                                   + np.square((b + self.d))/(b * np.sqrt(2 * b * self.d + self.d**2)) * 
                                    (np.pi - 2 * np.arctan(np.sqrt(self.d/(2 * b + self.d)) * np.tan(theta/2 + np.pi/4))))
        F_I_bend_cylinder = (1/3) * self.kappa_c * np.pi * (2 * alpha * np.power(self.r0/self.d, 3) - 1) * np.square(1 - self.c0 * self.d)
        F_I_bend_hemisphere = self.kappa_c * np.pi * np.square(2 - self.c0 * self.d)

        F_I_bending = F_I_bend_sphere + F_I_bend_torus + F_I_bend_cylinder + F_I_bend_hemisphere


        Area_I_sphere = 2 * np.pi * (1 + np.cos(theta)) * np.square(r)
        Area_I_torus = 2 * np.pi * b * ((b + self.d) * (np.pi/2 - theta) - b * np.cos(theta))
        Area_I_cylinder = 2 * np.pi * np.square(self.d) * ((4/3) * alpha * np.power(self.r0/self.d, 3) - (2/3))
        Area_I_hemisphere = 2 * np.pi * np.square(self.d)
        Area_I = Area_I_sphere + Area_I_torus + Area_I_cylinder + Area_I_hemisphere

        F_I_stretching = (1/2) * self.lambda_/(4 * np.pi * np.square(self.r0)) * np.square(Area_I - (4 * np.pi * np.square(self.r0))) -F0
        
        return F_I_bending + F_I_stretching + self.F_ext(alpha)

    #In this part, we calculate the free energy terms in the Helfrich free energy. F_II function return free energy 
    # of stretching and bending for the crossing part.
    def F_II(self, alpha1, r1, theta1, r2, theta2):
        F0 = (1./2.) * self.kappa_c * (4. * np.pi * np.square(2 - self.r0*self.c0))
        alpha2 = 1 + self.beta - alpha1
        alpha2 = 0
      
        b1 = (r1 * np.sin(theta1) - self.d)/(1. - np.sin(theta1))
        b2 = (r2 * np.sin(theta2) - self.d)/(1. - np.sin(theta2))
#
#
#
        F_II_bend_sphere1 = self.kappa_c * np.pi * (1 + np.cos(theta1)) * np.square(2 - self.c0 * r1)
        F_II_bend_torus1 = self.kappa_c * np.pi * (self.c0 * (2 + self.c0 * b1) * (b1 + self.d) * (np.pi/2 - theta1) - np.cos(theta1) * np.square(2 + self.c0 * b1)
                                  + np.square((b1 + self.d))/(b1 * np.sqrt(2 * b1 * self.d + self.d**2)) * 
                                   (np.pi - 2 * np.arctan(np.sqrt(self.d/(2 * b1 + self.d)) * np.tan(theta1/2 + np.pi/4))))
#
#
        F_II_bend_sphere2 = self.kappa_c * np.pi * (1 + np.cos(theta2)) * np.square(2 - self.c0 * r2)
        F_II_bend_torus2 = self.kappa_c * np.pi * (self.c0 * (2 + self.c0 * b2) * (b2 + self.d) * (np.pi/2 - theta2) - np.cos(theta2) * np.square(2 + self.c0 * b2)
                                  + np.square((b2 + self.d))/(b2 * np.sqrt(2 * b2 * self.d + self.d**2)) * 
                                   (np.pi - 2 * np.arctan(np.sqrt(self.d/(2 * b2 + self.d)) * np.tan(theta2/2 + np.pi/4))))
#
        F_II_bend_cylinder0 = (1/3) * self.kappa_c * np.pi * (2 * self.beta * np.power(self.r0/self.d, 3) - 1) * np.square(1 - self.c0 * self.d)
#
#
        F_II_bending = (F_II_bend_sphere1 + F_II_bend_torus1 + F_II_bend_cylinder0 + F_II_bend_sphere2 + F_II_bend_torus2)
#
# To calculate the stretching free energy term in the Helfrich free energy, we need to find the change in the total area
#
        Area_II_sphere1 = 2 * np.pi * (1 + np.cos(theta1)) * np.square(r1)       # Area of sphere in the donor department
        Area_II_torus1 = 2 * np.pi * b1 * ((b1 + self.d) * (np.pi/2 - theta1) - b1 * np.cos(theta1)) # Area of torus in the donor department
        Area_II_sphere2 = 2 * np.pi * (1 + np.cos(theta2)) * np.square(r2)      # Area of sphere in the donor department
        Area_II_torus2 = 2 * np.pi * b2 * ((b2 + self.d) * (np.pi/2 - theta2) - b2 * np.cos(theta2)) #Area of torus in the receiver department
        Area_II_cylinder0 = 2 * np.pi * np.square(self.d) * ((4/3) * self.beta * np.power(self.r0/self.d, 3) - (2/3)) #Area of pore
#
        Area_II = Area_II_sphere1 + Area_II_torus1 + Area_II_cylinder0 + Area_II_sphere2 + Area_II_torus2
# 
        F_II_stretching = (1/2) * self.lambda_/(4 * np.pi * np.square(self.r0)) * np.square(Area_II - (4 * np.pi * np.square(self.r0)))
        return F_II_bending + F_II_stretching + self.F_ext(alpha1)


     
    def F_ext(self, alpha):
        if (0 <= alpha <1):
               return - self.Fext*(alpha/(1+self.beta))
        else:
              return - self.Fext*(1/(1+self.beta))


   
    # FreeEnergy function returns the minimum free energy as a function of alpha.
    # In the input one can change the resolution of the solution.

    def FreeEnergy(self, r_prev=None,accur =100):
        
        if r_prev is None: r_prev = self.r0
        
        energy = []
        energy.append((0,0))

        # since the free energy is different in different stage we divid the process to three different processes: Filling, crossing and depletion
        # filling stage
        alphalinefill = np.linspace(self.alphamin, self.beta+self.alphamin, accur)
       	for alpha in alphalinefill:
       	    update_progress(alpha)
       	    rt = self.r_theta_pair(alpha, acc =100)
       	    r_val = rt[:, 0]
       	    t_val = rt[:, 1]
       	    F_val = np.array([self.F_I(alpha, r, t) for t, r in zip(t_val, r_val)])
       	    peaks = np.where((F_val[1:-1] > F_val[0:-2]) * (F_val[1:-1] > F_val[2:]))[0] + 1
            dips = np.where((F_val[1:-1] < F_val[0:-2]) * (F_val[1:-1] < F_val[2:]))[0] + 1
            #print(F_val)          
            F_max_arg = np.argmax(F_val)
            if len(peaks) != 0:
                F_val[peaks[-1]:] = F_val[peaks[-1]]
            amin = np.argmin(F_val)
            energy.append([F_val[amin],alpha])
           
        
        #crossing stage
        print(self.beta)
        alphalineCross1 = np.linspace(self.beta+self.alphamin + 0.01, (1+self.beta)/2 , 2*accur)
        for alpha in alphalineCross1:
           
            update_progress(alpha)
       	    rt1 = self.r_theta_pair(alpha, acc=100)
            
            
       	    r1_val = rt1[:, 0]
           
       	    t1_val = rt1[:, 1]
            
       	    rt2 = self.r_theta_pair(1 + self.beta - alpha, acc=100)
            
       	    r2_val = rt2[:, 0]
            
       	    t2_val = rt2[:, 1]
            
       	    F_val = np.array([[self.F_II(alpha, r1, t1, r2, t2) for t1, r1 in zip(t1_val, r1_val)] for t2, r2 in zip(t2_val, r2_val)])
           
            F_val = np.nan_to_num(F_val, nan=100000, posinf=100000, neginf=100000)
            
            
# searching for the  local minimum of free energy in crossing stage
            for i in range(F_val.shape[1]):
                F_arg_max = argrelextrema(F_val[:, i], np.greater)[0]
                if len(F_arg_max) == 1:
                    F_val[F_arg_max[0]:, i] =F_val[F_arg_max[0], i]
                elif len(F_arg_max) >= 1:
                    F_val[F_arg_max[1]:, i] =F_val[F_arg_max[1], i]
                    F_val[:F_arg_max[0], i] =F_val[F_arg_max[0], i]


            for j in range(F_val.shape[0]):
                F_arg_max2 = argrelextrema(F_val[j, :], np.greater)[0]
                if len(F_arg_max2) == 1:
                    F_val[j, F_arg_max2[0]:] =F_val[j, F_arg_max2[0]]
                elif len(F_arg_max) >= 1:
                    F_val[j, F_arg_max2[1]:] =F_val[j, F_arg_max2[1]]
                    F_val[j, :F_arg_max2[0]] =F_val[j, F_arg_max2[0]]
            
            

            F_min = np.min(F_val.flatten())
           
            if (F_min > 0 ):
                energy.append([F_min,alpha])

        alphalineCross2 = np.linspace( (1+self.beta)/2 + .01, 1-self.alphamin - 0.01 , 2*accur)
        for alpha in alphalineCross2:
           
            update_progress(alpha)
       	    rt1 = self.r_theta_pair(alpha, acc=100)
            
            
       	    r1_val = rt1[:, 0]
           
       	    t1_val = rt1[:, 1]
            
       	    rt2 = self.r_theta_pair(1 + self.beta - alpha, acc=100)
            
       	    r2_val = rt2[:, 0]
            
       	    t2_val = rt2[:, 1]
            
       	    F_val = np.array([[self.F_II(alpha, r1, t1, r2, t2) for t1, r1 in zip(t1_val, r1_val)] for t2, r2 in zip(t2_val, r2_val)])
           
            F_val = np.nan_to_num(F_val, nan=100000, posinf=100000, neginf=100000)
            
            
# searching for the  local minimum of free energy in crossing stage
            for j in range(F_val.shape[0]):
                F_arg_max2 = argrelextrema(F_val[j, :], np.greater)[0]
                if len(F_arg_max2) == 1:
                    F_val[j, F_arg_max2[0]:] =F_val[j, F_arg_max2[0]]
                elif len(F_arg_max) >= 1:
                    F_val[j, F_arg_max2[1]:] =F_val[j, F_arg_max2[1]]
                    F_val[j, :F_arg_max2[0]] =F_val[j, F_arg_max2[0]]


            for i in range(F_val.shape[1]):
                F_arg_max = argrelextrema(F_val[:, i], np.greater)[0]
                if len(F_arg_max) == 1:
                    F_val[F_arg_max[0]:, i] =F_val[F_arg_max[0], i]
                elif len(F_arg_max) >= 1:
                    F_val[F_arg_max[1]:, i] =F_val[F_arg_max[1], i]
                    F_val[:F_arg_max[0], i] =F_val[F_arg_max[0], i]
            
            

            F_min = np.min(F_val.flatten())
           
            if (F_min > 0 ):
                energy.append([F_min, alpha])


           # Depletion stage
        alphalineDepl = np.linspace(1-self.alphamin, 1+self.beta-self.alphamin , accur)
       
        for alpha in alphalineDepl:
       	    
       	    rt = self.r_theta_pair(1 + self.beta - alpha)
       	    r_val = rt[:, 0]
       	    t_val = rt[:, 1]
       	    F_val = np.array([self.F_I(1 + self.beta - alpha, r, t) for t, r in zip(t_val, r_val)])
       	    peaks = np.where((F_val[1:-1] > F_val[0:-2]) * (F_val[1:-1] > F_val[2:]))[0] + 1
            dips = np.where((F_val[1:-1] < F_val[0:-2]) * (F_val[1:-1] < F_val[2:]))[0] + 1
          
            F_max_arg = np.argmax(F_val)
            if len(peaks) != 0:
                F_val[peaks[-1]:] = F_val[peaks[-1]]
            amin = np.argmin(F_val)

            energy.append([F_val[amin],alpha])
            update_progress(alpha)
        energy.append([self.F_ext(alpha), 1+self.beta])       
        return np.array(energy)

#          In order to run the program with the parameters that you want. You can 
#   change the parameters in the following section


if __name__ == '__main__':
    print("running")
   # You can chnage the paramters inside
   #the input part of the class: vesicleTranslocation(r0, beta, Fext, kappa_c, lambda_)
   # 
    vt = vesicleTranslocation(4,0.3,0, kappa_c=2, lambda_=2)
    
    delta = vt.FreeEnergy(accur=20) # By choosing a value for accur, you can change the resolution of the 
                                    # of the computation
    x_val = delta[:, 0]
    y_val = delta[:, 1]
    plt.figure(figsize=(12, 8))
    plt.scatter(y_val,x_val, marker='*')
    plt.xlabel("alpha")
    plt.ylabel("Free Energy")
    plt.show()



