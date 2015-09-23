import os


f0 = open("t2.csv", "r")
lines = f0.readlines()
f0.close()
last_line = lines[-1].split(',')

for i in range(0, 54):


    # smc Fluxes    
    if i == 0:
        f = open("IC/smc_fluxes.csv", "w")
        
    
    # smc Homogeneous coupling_species
    elif i == 12:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/smc_homo_couple.csv", "w")
    
    # smc Heterogeneous coupling_species
    elif i == 15:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/smc_hetero_couple.csv", "w")
    
    # smc Equations
    elif i == 18:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/smc_equations.csv", "w")
    
    # smc Concentrations
    elif i == 23:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/smc_concentrations.csv", "w")
    
    # ec Fluxes
    elif i == 28:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/ec_fluxes.csv", "w")
    
    # ec Homogeneous coupling_species 
    elif i == 40:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/ec_homo_couple.csv", "w")
    
    # ec Homogeneous coupling_species
    elif i == 43:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/ec_hetero_couple.csv", "w")
    
    # ec Equations
    elif i == 46:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/ec_equations.csv", "w")
    
    # ec Concentrations
    elif i == 50:
        f.seek(-1, os.SEEK_END)
        f.truncate()
        f.close()
        f = open("IC/ec_concentrations.csv", "w")
        
    # write last line to file
        
    f.write(last_line[i] + ",")

f.seek(-1, os.SEEK_END)
f.truncate()
f.close()
    