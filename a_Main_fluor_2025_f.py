from numpy import loadtxt, zeros, append, vstack, asarray, savetxt, dstack, round
import numpy as np
from scipy import interpolate 
from a_MieTheory3fluor2_with_abs_final2023_multi_f import mie_theory, effective_medium  ## test_col1L_finalForm
from a_montefluor_long2_opt_with_abs_fast_1_multi2_check_valid_f import main_mc # a_montefluor_long2_opt_with_abs_fast_1_multi2
from Integration_fluor_couple_trans_f import solar_spectrum #a_Integration_fluor_couple_trans
import os.path
#from NN_feedforward import forward       
from interpolatefluor2_opt_f import interpolatee  ## interpolatefluor2      testflu_scat.txt
#   invers_cdf_opt 
# input file:  a_test_Mult.txt

# Main3fluo2_b_final_multi_opt2.py    opt2 indicates two fluorescent particles

#################### this has the final fast monte carlo for fluorescent code  ######################################

######## for validation with Yalcin , i used :  96%  yalcin_calsieu_emis.txt    Cals_yalcin_nk.txt     scat_coef_yalcin.xlsx
######## for validation with experiment side :  i used  final_YAC_8mm_new.txt   yag_emis_new.txt    and medium: sio2_medi_nk  D: 17  Dist: 40    80%      output: aaaa_vlid_exp_yag     output: aaaa_vlid_exp_LED  

# imports the provided txt file               testfluor.txt
def fname():
    filecheck = False
    infile = " "
    while filecheck is False:
        file = input('Input file: ')
        if os.path.exists(file):
            infile = loadtxt(file, comments="#", dtype=str, delimiter="/")
            filecheck = True
        else:
            print("File not found. Make sure it is in the same directory as this program.")

    return infile


# imports information from header of input file and loads in necessary files
def import_header(infile):
    # p is list of particle input file names, m for medium input files.
    # each of these files consists of three columns, wavelength, refractive index, and extinction coefficient
    p = zeros(0, dtype=str)
    m = zeros(0, dtype=str)
    Quantum_Yiel_1 =  zeros(0, dtype=str)
    Quantum_Yiel_2 =  zeros(0, dtype=str)
    emissio1 =  zeros(0, dtype=str)
    emissio2 =  zeros(0, dtype=str)

    # solar spectrum file to import. If there is not one, this variable remains blank
    ################solar = ""
    # ooutput file name
    output_name = ""
    # emissio= ""
    # Quantum_Yiel = ""   
    # mesh percentage, 1 keeps data as is. Above 1 increases the number of mesh points, below 1 decreases
    mesh_percentage = 1

    photons = 0
    line_d = 0
    # loops through header. Breaks loop when it hits the first "sim"
    for i in range(len(infile)):
        if infile[i][0:8] == "particle":
            p = append(p, infile[i][10:])       ##### notice that for "append(parameter, infile[i][10:]", it increases by 2 ([0:8] to [10:]), where for other-- 
        if infile[i][0:6] == "medium":                                          ###### -- paramters below, it increases by 1 
            m = append(m, infile[i][8:])
        if infile[i][0:6] == "output":
            output_name = infile[i][7:]
        if infile[i][0:4] == "emi1":
            emissio1 = append( emissio1 , infile[i][5:])
        if infile[i][0:4] == "emi2":
            emissio2 = append( emissio2 , infile[i][5:])
        # if infile[i][0:4] == "emi1":
        #     emissio2 = append( emissio2 , infile[i][5:])
        if infile[i][0:3] == "qy1":
            Quantum_Yiel_1 = append(Quantum_Yiel_1 , infile[i][4:])
        if infile[i][0:3] == "qy2":
            Quantum_Yiel_2 = append(Quantum_Yiel_2 , infile[i][4:])
        if infile[i][0:5] == "solar":
            solar = infile[i][6:]
        if infile[i][0:4] == "mesh":
            mesh_percentage = float(infile[i][5:])
        if infile[i][0:7] == "photons":
            photons = int(infile[i][8:])
        if infile[i][0:3] == "sim":
            line_d = i-1
            break

    # counts number of simulations
    sims = 0
    for i in range(len(infile)):
        if infile[i][0:3] == "sim":
            sims += 1

    # length is the length (number of wavelength datapoints) of the longest input file
    length = 0
    for i in range(len(p)):
        temp = loadtxt(p[i])         ### how the "p" is related to the text particle property
        if len(temp) > length:
            length = len(temp)
    for i in range(len(m)):
        temp = loadtxt(m[i])
        if len(temp) > length:
            length = len(temp)

    # arrays for material properties [number of materials, wavelengths, properties] 
    particle = zeros((len(p), length, 3))                                                   ######## we can use "length" as the number of wavelengths
    medium = zeros((len(m), length, 3))

    # need_interp is True if wavelengths of input files don't match, or the mesh_percentage is not 1
    need_interp = False

    # imports particle files and checks if they need to be interpolated
    temp =loadtxt(p[0])
    particle[0, :len(temp), :] = temp
    if len(p) > 1:
        for i in range(1, len(p)):
            temp = loadtxt(p[i])
            particle[i, :len(temp), :] = temp
            if len(temp) != length:
                need_interp = True
            for j in range(length):
                if particle[i, j, 0] < (particle[0, j, 0] - 0.01) or particle[i, j, 0] > (particle[0, j, 0] + 0.01):
                    need_interp = True

    # imports medium files and checks if they need to be interpolated
    for i in range(len(m)):
        temp = loadtxt(m[i])
        medium[i, :len(temp), :] = temp
        if len(temp) != length:
            need_interp = True
        # check wavelengths match
        for j in range(length):
            if medium[i, j, 0] < (particle[0, j, 0] - 0.01) or medium[i, j, 0] > (particle[0, j, 0] + 0.01):
                need_interp = True
    
    column11_1 = []
    column22_1 = []
    file_path_1 = emissio1[0]
    if os.path.exists(file_path_1 ):
        with open(file_path_1 , 'r') as file_1 :
            for line1_1  in file_1 :
        # Split each line into two columns (assuming they are separated by a space)
                parts_1  = line1_1 .strip().split()
                if len(parts_1 ) == 2:
                # Append the values to the respective columns
                    column11_1 .append(float(parts_1 [0]))
                    column22_1 .append(float(parts_1 [1]))
                else: 
                    print("Make sure two columns are inserted, the wavelengths and the respective normalized emission intensity")

        # Convert the lists to numpy arrays
        emission_profile_1  = np.column_stack((column11_1 , column22_1 ))
    else:
        print("emission profile txt file is not found. Make sure it is in the same directory as this program.")
    print('1')
    if len(emissio2) >0:
        #print(emissio1[0]) 
        column11_2 = []
        column22_2 = []
        file_path_2 = emissio2[0]
        if os.path.exists(file_path_2):
            with open(file_path_2, 'r') as file_2:
                for line1_2 in file_2:
        # Split each line into two columns (assuming they are separated by a space)
                    parts_2 = line1_2.strip().split()
                    if len(parts_2) == 2:
                # Append the values to the respective columns
                        column11_2.append(float(parts_2[0]))
                        column22_2.append(float(parts_2[1]))
                    else: 
                        print("Make sure two columns are inserted, the wavelengths and the respective normalized emission intensity")

        # Convert the lists to numpy arrays
            emission_profile_2 = np.column_stack((column11_2, column22_2))
        else:
            print("emission profile txt file is not found. Make sure it is in the same directory as this program.")
    else: 
        emission_profile_2=zeros(100)   

    print('2')

    column1 = []
    column2 = []
# file_path =Quantum_Yiel
    file_path1 =Quantum_Yiel_1[0]
    ch = 0
    ch2 = 0
    numpy_array = zeros((0, 2))
    if os.path.exists(file_path1):
        infile2 = np.array(loadtxt(file_path1, dtype=float))
        if infile2.size == 1 :     
            Qy_1= np.array([infile2.item()])#float(infile2)
            ch = 1
        else:
            print('make sure the QY values are associated with absorption wavelengths of fluorescent particles')
            ch2 = 1
            with open(file_path1, 'r') as file:
                for line2 in file:
                # Split each line into two columns (assuming they are separated by a space)
                    parts = line2.strip().split()
                    if len(parts) == 2:
                        # Append the values to the respective columns
                        column1.append(float(parts[0]))
                        column2.append(float(parts[1]))
                    else: 
                        print("Make sure two columns are inserted, the wavelengths and the respective Quantum Yield")

             # Convert the lists to numpy arrays
            numpy_array = np.column_stack((column1, column2))
            Qy = numpy_array
    else:
        print("The Quantum Yield file not found. Make sure it is in the same directory as this program.")

    if len(Quantum_Yiel_2) >0:
        file_path1_qy2 =Quantum_Yiel_2[0]
        ch = 0
        ch2 = 0
        numpy_array = zeros((0, 2))
        if os.path.exists(file_path1_qy2):
            infile2 = np.array(loadtxt(file_path1_qy2, dtype=float))
            if infile2.size == 1 :     
                Qy_2= np.array([infile2.item()])#float(infile2)
                ch = 1
            else:
                print('make sure the QY values are associated with absorption wavelengths of fluorescent particles')
                ch2 = 1
                with open(file_path1, 'r') as file:
                    for line2 in file:
                # Split each line into two columns (assuming they are separated by a space)
                        parts = line2.strip().split()
                        if len(parts) == 2:
                        # Append the values to the respective columns
                            column1.append(float(parts[0]))
                            column2.append(float(parts[1]))
                        else: 
                            print("Make sure two columns are inserted, the wavelengths and the respective Quantum Yield")

             # Convert the lists to numpy arrays
                numpy_array = np.column_stack((column1, column2))
                Qy = numpy_array
        else:
            print("The Quantum Yield file not found. Make sure it is in the same directory as this program.")
    else:
        Qy_2= np.array(1)*0

    print(Qy_1)
    if len(Quantum_Yiel_2)>0:
        print(Quantum_Yiel_2[0])
    # send to interpolation method if the wavelengths don't match or the mesh_percentage is not 1     python -m pdb Main3fluo.py
    #print("before debug ", particle[0,100,:])  
    if need_interp is True:
        print('Interpolating properties to match wavelengths for each input')
        particle, medium, index, start_wl = interpolatee(particle, medium, length, mesh_percentage, numpy_array, ch2,ch)  
            ######### index added here in "interpolate" file 

    elif mesh_percentage != 1:### wont be used 
        print('Interpolating properties for new mesh')
        particle, medium, index, start_wl = interpolatee(particle, medium, length, mesh_percentage) 
       
    #print(" debug ", particle[0,100,:])                                                                                                            #indexx=range(0,index)
    return particle, medium, output_name, sims, photons, line_d, index, start_wl, emission_profile_1, solar,emission_profile_2,Qy_1,Qy_2  ### emission_profile_3 output_name, solar, sims,


def check_diameters(current_sim, fv, sizes, check):
    # check to make sure same number of diameters and volume fraction
    if len(fv) != len(sizes):
        print('Number of diameters does not match number of volume fractions provided in sim:', current_sim)
        print("Please re-enter input file once corrected.")
        check = True
    return check


# finds info from input and sends to Mie Theory to calculate optical properties
def optical(line, infile, particle, medium, check, start_wl, index):
    print("Running Mie theory")
    # array of optical properties to send to Monte Carlo
    prop = zeros((0, 13))    ### "5"
    optics_sum = zeros((13, len(particle[0, :, 0])))   ### "5"
    optical_per_layer = zeros((13, len(particle[0, :, 0]), 0))    ### "5"
    vol_frac_sum = 0
    layers = 1
    count = 0
    fluor=0
    current_sim = 1

    # data1 = open('n_sub.txt', 'r').read()
    # values1 = [float(f) for f in data1.split("\n")]
    # values1=np.array(values1)

    # data2 = open('k_sub.txt', 'r').read()
    # values2 = [float(f) for f in data2.split("\n")]
    # values2=np.array(values2)

    # data1 = open('n_sub_1.txt', 'r').read()
    # values1 = [float(f) for f in data1.split("\n")]
    # values1=np.array(values1)

    # data2 = open('k_sub_1.txt', 'r').read()
    # values2 = [float(f) for f in data2.split("\n")]
    # values2=np.array(values2)

    data1 = open('n_sub_al.txt', 'r').read()
    values1 = [float(f) for f in data1.split("\n")]
    values1=np.array(values1)

    data2 = open('k_sub_al.txt', 'r').read()
    values2 = [float(f) for f in data2.split("\n")]
    values2=np.array(values2)
    happens = 0
    again = 0
    n_c = 0
    for i in range(line+1, len(infile)):  ### just cover after "sim" word 
        if infile[i][0:8] == "particle":
            if count > 0:
                lay1 =True
                optics = mie_theory(sizes, fv, particle[int(ptype - 1), :, :], medium[int(mtype - 1), :, :], thickness, dist,fluor,  start_wl, index,happens, again)    #######  This reads from the text file of properties of particle and medium
                optics_sum += optics
            count += 1                                  ######## notice that "ptype" and "mtype" refer to 0 row and 1 row
            ptype = int(infile[i][8])
        if infile[i][0:5] == "upper":
            upper = float(infile[i][6:])
        elif infile[i][0:5] == "lower":
            lower = float(infile[i][6:])
        elif infile[i][0:6] == "medium":
            mtype = int(infile[i][6])
        elif infile[i][0:11] == "ref_medium:":
            n_c = float(infile[i][11:])
        elif infile[i][0:2] == "t:":
            thickness = float(infile[i][2:])
            thickness = thickness / 10000

        elif infile[i][0:6] == "fluor:":     #################
            fluor = int(infile[i][6:])
            if fluor == 1 :
                happens = 1
            if happens == 1:
                again += 1
            
        elif infile[i][0:2] == "d:":
            sizes = (infile[i][2:]).split(",")
            sizes = asarray(sizes, dtype=float)
            sizes = sizes/2
        elif infile[i][0:3] == "vf:":
            fv = (infile[i][3:]).split(",")
            fv = asarray(fv, dtype=float)
            fv = fv/100
            vol_frac_sum += sum(fv)
        elif infile[i][0:5] == "dist:":
            dist = float(infile[i][5:])
            dist = dist/100
        elif infile[i][0:5] == "layer" and infile[i][5:] != "1":
            # check to make sure same number of diameters as volume fractions
          
         
            check = check_diameters(current_sim, fv, sizes, check)
            layers += 1
    
            optics = mie_theory(sizes, fv, particle[int(ptype - 1), :, :], medium[int(mtype - 1), :, :], thickness, dist, fluor, start_wl, index,happens,again)
            optics_sum += optics 
            optics_sum[0, :] = optics[0, :]
            optics_sum[11, :] = optics[11, :]
            optics_sum[12, :] = optics[12, :]
            optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype - 1), :, :],fluor)
            optical_per_layer = dstack((optical_per_layer, optics_sum))
            vol_frac_sum = 0
            optics_sum = zeros((13, len(particle[0, :, 0])))     ### "5"
            happens = 0
            again = 0
            count = 0
        elif infile[i][0:3] == "sim" and infile[i][4:] != "1":
            current_sim = int(infile[i][4:])
            # check to make sure same number of diameters as volume fractions

            check = check_diameters(current_sim, fv, sizes, check)
            optics = mie_theory(sizes, fv, particle[int(ptype - 1), :, :], medium[int(mtype - 1), :, :], thickness, dist,fluor, start_wl, index,happens,again)  ##### the colon ":" is the index which values are taken and inserted in thos function
            optics_sum += optics
            optics_sum[0, :] = optics[0, :]
            optics_sum[11, :] = optics[11, :]
            optics_sum[12, :] = optics[12, :]
            
            optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype - 1), :, :],fluor)
            optical_per_layer = dstack((optical_per_layer, optics_sum))
            for j in range(len(particle[0, :, 0])):
                prop = vstack((prop, [upper, 0, 0, 0,0, 0,0,0,0,0,0,0,0]))           ### notice there are three prop ( properties ) with VSTACK array type
                for k in range(layers):
                    prop = vstack((prop, optical_per_layer[:, j, k]))     ### j should be the number of "n" corresponding to different wavelengths 
                prop = vstack((prop, [lower, 0, 0, 0,0, 0,0,0,0,0,0,0,0]))       ## values[j] lower
                prop = vstack((prop, [0, 0, 0, 0,0, 0,0,0,0,0,0,0,0]))
            vol_frac_sum = 0
            optics_sum = zeros((13, len(particle[0, :, 0])))   ### "5"
            optical_per_layer = zeros((13, len(particle[0, :, 0]), 0))   ### "5"
            layers = 1
            happens = 0
            again = 0
            count = 0
        # check to make sure same number of diameters as volume fractions
    check = check_diameters(current_sim, fv, sizes, check)
    #print(particle[:,:,:]) LAYER2

    optics = mie_theory(sizes, fv, particle[int(ptype - 1), :, :], medium[int(mtype - 1), :, :], thickness, dist, fluor, start_wl, index,happens,again)
    optics_sum += optics
    optics_sum[0, :] = optics[0, :]
    optics_sum[11, :] = optics[11, :]
    optics_sum[12, :] = optics[12, :]
   
    optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype - 1), :, :],fluor)     ### notice the particle and medium have all indicies of second and third dimension
    optical_per_layer = dstack((optical_per_layer, optics_sum))
    for j in range(len(particle[0, :, 0])):
        #prop = vstack((prop, [values1[j],values2[j], 0,0, 0, 0,0,0,0,0,0,0,0]))
        prop = vstack((prop, [upper,0, 0,0, 0, 0,0,0,0,0,0,0,0]))
        for k in range(layers):
            prop = vstack((prop, optical_per_layer[:, j, k]))
        prop = vstack((prop, [lower,0, 0,0, 0, 0,0,0,0,0,0,0,0]))    ### lowerr[j] lower      values[j]
        #prop = vstack((prop, [values1[j],values2[j], 0,0, 0, 0,0,0,0,0,0,0,0]))    ### lowerr[j] lower      values[j]
        prop = vstack((prop, [0, 0, 0, 0,0, 0,0,0,0,0,0,0,0]))
    return prop, check


# breaks down import file, runs Mie Theory calculations
def nanoparticle(infile, check):
    # imports information from header of input file and loads in necessary files
    particle, medium, output_name, sims, photons, line, index, start_wl, emission_profile_1, solar,emission_profile_2,Qy_1,Qy_2 = import_header(infile)     #### output_name, solar, sims,

    # gets optical properties of each simulation
    prop, check = optical(line, infile, particle, medium, check, start_wl, index)  ##########################       ####### here is the connection with Mie theory
    sims_per_medium = len(particle[0, :, 0])
    wavelengths = particle[0, :, 0]

    return prop, photons, output_name, sims, sims_per_medium, wavelengths, check, index,start_wl , emission_profile_1, solar,particle,emission_profile_2,Qy_1,Qy_2 ### sims_per_medium, solar, wavelengths,



# Makes sure the input data is within the range the NN is trained on
def check_NN_range(prop, check):
    # If check is True, NN does not run and the user must change the input file
    for i in range(len(prop[:, 0])):
        # checks refractive indices
        if prop[i, 0] != 0:                                              ### the prop[i,0] for NN check, other properties ua,us,g are in prop[i,1(2)(3)]
            if prop[i, 4] == 0 and prop[i, 0] != 1:
                print("Neural network is only trained on boundary indices of 1")
                check = True
            if prop[i, 0] < 1 or prop[i, 0] > 10:
                print("A refractive index of" + str(prop[i, 0]) + "is out of bounds [1, 7]")
                check = True
        # Checks the absorption coefficient, scattering coefficient, asymmetry parameter
        if prop[i, 4] != 0:
            # absorption coefficient
            if prop[i, 1] < 0 or prop[i, 1] > 300000:
                print("An absorption coefficient of" + str(prop[i, 1]) + "is out of bounds [0, 1,000,000] (1/cm)")
                check = True
            # scattering coefficient
            if prop[i, 2] < 0 or prop[i, 2] > 150000:
                print("A scattering coefficient of" + str(prop[i, 2]) + "is out of bounds [0, 200,000] (1/cm)")
                check = True
            # asymmetry parameter
            if prop[i, 3] < 0 or prop[i, 3] > 1:
                print("An asymmetry parameter of" + str(prop[i, 3]) + "is out of bounds [0, 1]")
                check = True
        # check for multi layered sims
        if prop[i, 4] != 0 and prop[i+1, 4] != 0:
            print("The neural network cannot predict multi-layer media")
            check = True
        # check paint thickness
        if prop[i, 4] != 0:
            if prop[1, 4] < .0005 or prop[1, 4] > .05:
                print("A thickness of " + str(prop[1, 4]*10000) + " is out of bounds [5, 500] \u03BCm")
                check = True

        if check is True:
            print("Please re-enter input file once corrected, or change to a Monte-Carlo simulation")
            break
    return check


# checks to make sure each required item per simulation is there
def check_for_word_in_sim(infile, word, statement, check):
    # first, find first line of body to start with
    start_line = 0
    for i in range(len(infile)):
        if infile[i][:3] == "sim":
            start_line = i
            break

    word_present = 1
    sim_number = str(1)
    word_size = len(word)
    for i in range(start_line, len(infile)):
        if (infile[i][:3] == "sim" and infile[i][4:] != '1'):

            if word_present == 1:
                check = True
                print("Sim:" + sim_number + " must have " + statement)
            sim_number = infile[i][4:]
            word_present = 1
        if infile[i][:word_size] == word:
            word_present = 0
    if word_present == 1:
        check = True
        print("Sim:" + sim_number + " must have " + statement)
    return check

# check the input file for errors
def check_input_for_errors(infile):
    # check = True is there is an issue caught with the input file, otherwise false
    check = False

    # remove all spaces and make lowercase to clean things up
    for i in range(len(infile)):
        infile[i] = infile[i].replace(' ', '')
        infile[i] = infile[i].lower()

    # first line must specify "nn" or "mc"
    if infile[0] != "mc" and infile[0] != "nn":
        print("Either MC or NN must be specified in the first line of the input file.")
        check = True

    ### check for each item in header
    # initially these variables are set to 0, if they are changed to 1 then it is good. If it remains 0 then there is an error
    output = particle = medium = photon = 0
    if infile[0] == "nn":
        photon = 1

    for i in range(len(infile)):
        # first it checks the required items, output name, particle, medium, and number of photons if running MC
        # output name
        if infile[i][:6] == 'output' and infile[i][7:] != '':
            output = 1
        # at least one particle
        if infile[i][:8] == 'particle' and infile[i][10:] != '':
            particle = 1
        # at least one medium
        if infile[i][:6] == 'medium' and infile[i][8:] != '':
            medium = 1
        # number of photons if using Monte Carlo
        if infile[0] == "mc":
            if infile[i][:7] == "photons" and infile[i][9:] != '':
                photon = 1
        if infile[i][:4] == "mesh":
            if float(infile[i][5:]) <= 0:
                check = True
                print('Mesh value must be > 0')
        # break loop once the first sim starts, this only checks the header
        if infile[i][:3] == "sim":
            break

    # print out error for each bug found in header
    if output == 0:
        check = True
        print("No output file name specified")
    if particle == 0:
        check = True
        print("No particle input specified")
    if medium == 0:
        check = True
        print("No medium input specified")
    if photon == 0:
        check = True
        print("Number of photons must be specified")

    ### check for errors in the body of the file
    # make sure each sim is numbered correctly
    sim_number_should_be = 1
    # if sim_number_error is flipped to 0, there is an error
    sim_number_error = 1
    for i in range(len(infile)):
        if infile[i][:3] == "sim":
            if infile[i][4:] != str(sim_number_should_be):
                sim_number_error = 0

            sim_number_should_be += 1
    if sim_number_error == 0:
        check = True
        print("Error in simulation numbers. Make sure they are labeled Sim: 1, Sim: 2, etc.")

    # check for upper boundary condition
    check = check_for_word_in_sim(infile, 'upper', 'defined upper boundary condition', check)
    # check for lower boundary condition
    check = check_for_word_in_sim(infile, 'lower', 'defined lower boundary condition', check)
    # check for at least one layer
    check = check_for_word_in_sim(infile, 'layer', 'at least one defined layer', check)
    # check for at least one medium
    check = check_for_word_in_sim(infile, 'medium', 'at least one medium', check)
    # check for medium thickness
    check = check_for_word_in_sim(infile, 't:', 'a defined thickness', check)
    # check for at least one particle
    check = check_for_word_in_sim(infile, 'particle', 'at least one particle', check)
    # check for particle size
    check = check_for_word_in_sim(infile, 'd:', 'a defined particle diameter', check)
    # check for VF
    check = check_for_word_in_sim(infile, 'vf', 'a defined particle volume fraction', check)
    # check for distribution
    check = check_for_word_in_sim(infile, 'dist', 'a defined distribution. Put Dist: 0 for no distribution', check)


    if check:
        print("Please re-enter input file once corrected.")
        print("\n")
    return check, infile


def main_func():

    # loop until input file has no identifiable errors
    # check is initially True, if the input file has no identifiable errors, check will be false
    check = True
    while check:
        # imports the provided txt file   
        infile = fname()  
        #infile = loadtxt('inputfilen.txt', comments="#", dtype=str, delimiter="/")#  a_test_Mult        testMult
        
        ## below it the one i was working with 
        #infile = loadtxt('aa_test_col1L_finalForm_f.txt', comments="#", dtype=str, delimiter="/")#  a_test_Mult        testMult
        ##infile = loadtxt('inputfile1.txt', comments="#", dtype=str, delimiter="/")#  a_test_Mult        testMult

        # check the input file for errors
        check, infile = check_input_for_errors(infile)

        # if the file looks good so far, go ahead and import / run Mie theory
        if not check:                                                               ############ means check =True ????
            # calculates optical properties from inputs in infile
            prop, photons, output_name, sims, sims_per_medium,  wavelengths, check, index, start_wl, emission_profile_1, solar,particle ,emission_profile_2,Qy_1,Qy_2= nanoparticle(infile, check) ###  sims_per_medium, solar, wavelengths,
            # if nn, check the properties are within range
            if infile[0] == "nn":
                check = check_NN_range(prop, check)


    # Send to either Monte Carlo or Neural Network
    if (infile[0][0:2]).lower() == "mc":
        print('before main_mc')
        # Monte Carlo   Main3fluo2_b_final.py
        result1, rsppp, ff2,ff3, result_check ,I_so ,result1_radiosity= main_mc(prop, photons, index, start_wl, emission_profile_1,Qy_1,emission_profile_2,Qy_2) #emission_profile_2 ,Qy_2                 ### here is the connection between Mie coefficinet and Monte carlo
        check = False                       # result_check,nn1,nn2,nn3,nn4          ,remi 
    ###elif (infile[0][0:2]).lower() == "nn":

        # Run NN prediction if properties are within NN range
        ###if check is False:
            ###results = forward(prop)
    else:
        print("Please specify either MC or NN in the first line of the input file")
        print("Please re-enter input file once corrected")
        main_func()
        check = True


    # if solar spectrum is provided, integrate for solar reflectance
    
    output_sim = str(output_name) + str(0+1) +".txt"   ##  str(i+1)
    with open(output_sim, 'w') as f:
            f.write('Wavelength\tSpecular R\tDiffuse R\tA\tT')
            f.write('\n')
            for i in range(len(index)):  # length 
                f.write(str(round(wavelengths[i], 10)) + '\t')
                f.write(str(round( rsppp[i], 10)) + '\t') ## [j + length * i, 0]
                f.write(str(round(result1[i, 1], 10)) + '\t')
                f.write(str(round(result1[ i, 0], 10)) + '\t')
                f.write(str(round(result1[ i, 2], 10)) + '\t') 
                f.write('\n')
    f.close() 
    # output_name_2 =''
    output_name_2 = 'aaaa_final_results_example1'
    output_sim_2 = str(output_name_2) + str(0+1) +".txt"   ##  str(i+1)
    with open(output_sim_2, 'w') as f:
            f.write('Wavelength\tSpecular R\tDiffuse R\tA\tT')
            f.write('\n')
            for i in range(len(index)):  # length 
                f.write(str(round(wavelengths[i], 10)) + '\t')
                f.write(str(round(result1_radiosity[i, 0], 10)) + '\t')
                f.write(str(round(result1_radiosity[ i, 1], 10)) + '\t')
                f.write(str(round(result1_radiosity[ i, 2], 10)) + '\t') 
                f.write('\n')
    f.close() 

    if solar != "":
        # solar reflectance value for each simulation
        solar_r = zeros(7) #zeros(sims)
        # reflectance at each wavelength to be integrated with the solar spectrum
        refl = zeros((len(wavelengths), 2))
        # move wavelengths over to refl
        refl[:, 0] = wavelengths
        # loop through each simulation to calculate solar reflectance

        for i in range(2):  # for i in range(sims)  ------------------------------------------------------------
            # add the specular and diffuse reflectance together for total reflectance
            refl[:, 1] = (result1[:, 1])
            # send the function to integrate    output_sim = str(output_name) + str(0+1) +".txt"
            r,a ,t,solar_power,reflect,absorb,transmit= solar_spectrum(output_sim)# solar_r[i] = solar_spectrum(output_sim)        solar, refl
            solar_r[0] = r
            solar_r[1] = a
            solar_r[2] = t
            solar_r[3] = solar_power
            solar_r[4] = reflect
            solar_r[5] = absorb
            solar_r[6] = transmit

            # ------------------------------------------------------------
        # save solar reflectance of each sim to output file _solarDL
        output_solar = str(output_name) + "solar.txt"
        with open(output_solar, 'w') as f:
            for i in range(len(solar_r)):
                f.write(str(round(solar_r[i], 5))+ '\t')
                f.write('\n')
        f.close()



















    # save output scripts
   ### length = int(len(results[:, 0])/sims)
    #for i in range(len(index)):  ### sims
    # output_sim = str(output_name) + str(0+1) +".txt"   ##  str(i+1)
    # with open(output_sim, 'w') as f:
    #         f.write('Wavelength\tSpecular R\tDiffuse R\tA\tT')
    #         f.write('\n')
    #         for i in range(len(index)):  # length 
    #             f.write(str(round(wavelengths[i], 7)) + '\t')
    #             f.write(str(round( rsppp[i], 7)) + '\t') ## [j + length * i, 0]
    #             f.write(str(round(result1[i, 1], 7)) + '\t')
    #             f.write(str(round(result1[ i, 0], 7)) + '\t')
    #             f.write(str(round(result1[ i, 2], 7)) + '\t') 
    #             f.write('\n')
    # f.close()
    with open("a_prop", 'w') as f:
        all_colums=len(prop[:,0])
        for ii in range(len(prop[:,0])):
            if ii == 0:
                for iii in range(len(prop[0,:])):
                    f.write(str(round(prop[ii+1,iii], 6))  + '\t' )
                f.write('\n')
            if prop[ii,0] ==0 and ii != 0 and ii != all_colums-1:
            
                for iii in range(len(prop[0,:])):
                    f.write(str(round(prop[ii+2,iii], 6))  + '\t' )
                f.write('\n')
            if ii == all_colums-1:           
                for iii in range(len(prop[0,:])):
                    f.write(str(round(prop[ii-2,iii], 6))  + '\t' )
                f.write('\n')
    
    f.close()

    with open("aaaa_AM1555", 'w') as f:
        for ii in range(len(index)):
            f.write(str(round(I_so[ii], 10))  + '\t' )
            f.write('\n')
    f.close()
    # with open("a_particle2", 'w') as f:
    #     for ii in range(len(particle[1,:,0])):
    #         for iii in range(len(particle[1,0,:])):
    #             f.write(str(round(particle[1,ii,iii], 6))  + '\t' )
    #         f.write('\n')
    # f.close()
    # with open("a_particle3", 'w') as f:
    #     for ii in range(len(particle[2,:,0])):
    #         for iii in range(len(particle[2,0,:])):
    #             f.write(str(round(particle[2,ii,iii], 6))  + '\t' )
    #         f.write('\n')
    # f.close()
    # with open("a_particle4", 'w') as f:
    #     for ii in range(len(particle[3,:,0])):
    #         for iii in range(len(particle[3,0,:])):
    #             f.write(str(round(particle[3,ii,iii], 6))  + '\t' )
    #         f.write('\n')
    # f.close()
    # with open("a_particle5", 'w') as f:
    #     for ii in range(len(particle[4,:,0])-1):
    #         for iii in range(len(particle[4,0,:])):
    #             f.write(str(round(particle[4,ii,iii], 6))  + '\t' )
    #         f.write('\n')
    # f.close()
    # with open("a_particle6", 'w') as f:
    #     for ii in range(len(particle[5,:,0])):
    #         for iii in range(len(particle[5,0,:])):
    #             f.write(str(round(particle[5,ii,iii], 6))  + '\t' )
    #         f.write('\n')
    # f.close()
    #with open("matrixreflect", 'w') as f:
        #for ii in range(len(index)):
            #for iii in range(len(index)):  # length 
              #  f.write(str(round(ff2[ii,iii], 4)) + '\t')
           # f.write('\n')
    #f.close()
    # for jj in range(photons): 
    #         for ii in range((len(index))):
    #             f.write(str(round(z_track[jj,ii], 6)) + '\t')
    #         f.write('\n')
    # f.close()
  
    """
    with open("uz_n", 'w') as f: # ph_track1_ref
        for ii in range(100000): # len(index)
            f.write(str(round(ph_track1[ii], 5)))
            f.write('\n')
    f.close()
    with open("step_n", 'w') as f: #ph_track1_trans
        for ii in range(100000):  # len(index)
            f.write(str(round(ph_track2[ii], 5)))
            f.write('\n')
    f.close()
    """
    """
    with open("remit", 'w') as f: #ph_track1_trans
        for ii in range(100000): # len(index)
            f.write(str(round(remi[ii], 5)))
            f.write('\n')
    f.close()
       """
    """
    with open("reflect_Matlab", 'w') as f:
        for ii in range(len(index)):
            f.write(str(round(reflect[ii], 6)))
            f.write('\n')
    f.close()
    
    with open("step_r", 'w') as f:
        for ii in range(10000):
            f.write(str(round(step_r[ii,0], 6)) + '\t')
            f.write(str(round(step_r[ii,1], 6)) + '\t')
            f.write(str(round(step_r[ii,2], 6)) + '\t')
            f.write('\n')
    f.close()
    with open("z_r", 'w') as f:
        for ii in range(10000):
            f.write(str(round(z_r[ii,0], 6)) + '\t')
            f.write(str(round(z_r[ii,1], 6)) + '\t')
            f.write(str(round(z_r[ii,2], 6)) + '\t')
            f.write('\n')
    f.close()
    with open("uz_r", 'w') as f:
        for iii in range(10000):
            f.write(str(round(uz_r[iii,0], 6)) + '\t')
            f.write(str(round(uz_r[iii,1], 6)) + '\t')
            f.write(str(round(uz_r[iii,2], 6)) + '\t')
            f.write('\n')
    f.close()
    with open("track_re", 'w') as f:
        for iiii in range(10000):
            f.write(str(round(nn4[iiii,0], 6)) + '\t')
            f.write(str(round(nn4[iiii,1], 6)) + '\t')
            f.write(str(round(nn4[iiii,2], 6)) + '\t')
            f.write('\n')
    f.close()
    """

    # with open("reflect_color", 'w') as f:
    #     for iiii in range(len(index)): 
    #         f.write(str(round(reflect[iiii], 6)))
    #         f.write('\n')
    # f.close()

    
    #with open("matrixtransmit", 'w') as f:
      #  for ii in range(len(index)):
      #      for iii in range(len(index)):  # length 
      #          f.write(str(round(ff3[ii,iii], 4)) + '\t')
      #      f.write('\n')
   # f.close()
    # with open("aa_rr_tt_fff3", 'w') as f:
    #     for ii in range(len(index)):
    #         f.write(str(round(result_check[ii,0], 6)) + '\t')
    #         f.write(str(round(result_check[ii,1], 6)) + '\t')
    #         f.write(str(round(result_check[ii,2], 6)) + '\t')
    #         f.write(str(round(result_check[ii,3], 6)) + '\t')
    #         f.write('\n')
    # f.close()
    print("Results saved!")
    return


if __name__ == "__main__":
    print('\033[1m{: ^75s}\033[0m'.format("FOS"))
    print('{: ^75s}'.format("Fast Optical Spectrum calculations for nanoparticle media"))
    print('{: ^75s}'.format("Version: 0.2.0\n"))
    print('{: ^75s}'.format("Daniel Carne, Joseph Peoples, Zherui Han, Dudong Feng, Xiulin Ruan")) 
    print('{: ^75s}'.format("School of Mechanical Engineering, Purdue University"))
    print('{: ^75s}'.format("West Lafayette, IN 47909, USA\n"))

    main_func()
