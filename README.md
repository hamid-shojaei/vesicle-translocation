# vesicle-translocation
## Source Code of the follwoing book chapter in Springer nano series.
*Title*: **Free Energy Minimization for Vesicle Translocation Through a Narrow Pore**

*Authors*: **Hamid R. Shojaei, Ahad Khaleghi Ardabili and Murrugappan Muthukumar**


In this program we are planing to find local minima of the free energy during the process of vesicle translocation.
This process has three different stages :filling, crossing and depletion. 

*Vesicletranslocation.py* is a Python program that generates a figure of free energy versus α (α=V_t/V_0, where V_t is the volume of the vesicle that is either inside the pore or that has passed through the pore and V_0 is incompressible vesicle volume). 
### Installation
 Please make sure that you have Python 3. The program *Vesicletranslocation.py* uses the follwoing libraries:
  *numpy, matplotlib and math.* 
 
In order to use the program, user needs to choose the following parameters: r0,  β, Fext, κ_c and λ
which are:  the initial radius of vesicle, the ratio of the ratio of the volume of the pore to that of the vesicle, bending modulus and stretching modulus, respectively.

VesicleTranslocation(r0, β, Fext, κ_c , λ) is a class that as input accepts r0  β, Fext, κ_c and λ as input and as output it    gives an Numpy array of minimum free energies and corresponding α (α=V_t/V_0, where V_t is the volume of the vesicle that is either inside the pore or that has passed through the pore). If the user calls the FreeEnergy method from the class VesicleTranslocation, a plot of  free energy versus α will be generated.  


