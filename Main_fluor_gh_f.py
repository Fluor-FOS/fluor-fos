#  Fluor-FOS: a free software that is used for calculations of fluorescent nanoparticle media by combining Mie theory with modified Monte Carlo simulations
#  Copyright (C) 2025 Khalid Alhammadi <alhammak@purdue.edu>
#  Copyright (C) 2025 Daniel Carne <dcarne@purdue.edu>
#  Copyright (C) 2025 Xiulin Ruan <ruan@purdue.edu>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.       
from numpy import loadtxt, zeros, append, vstack, asarray, dstack, round,  interp, sum ,arange ,random
import numpy as np
from a_MieTheory3fluor2_gh import mie_theory, effective_medium 
from a_montefluor_gh_f import main_mc 
from Integration_fluor_gh import solar_spectrum 
import os.path
import matplotlib.pyplot as plt
from interpolatefluor_gh import interpolatee  
from color_post_gh import display_color_from_reflectance 
from scipy.special import erfinv

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

def trapz_not_numpy(y, x):
    num_trap = len(x)-1
    integral = 0
    for i in range(num_trap):
        dx = x[i+1]-x[i]
        integral += 0.5 * (dx) * (y[i+1] + y[i])
    return integral

def inv_cdf(emission_profile):
    convert_to_nm = False
    data1= emission_profile[:, 0] 
    data2= emission_profile[:, 1] 
    flo_wl_start=emission_profile[0,0]
    flo_wl_end=emission_profile[-1,0]
    if flo_wl_start <30 and flo_wl_end<30:
        convert_to_nm= True
    f=arange(flo_wl_start , flo_wl_end)  
    wave_flo_no= len(f)
    x=zeros(wave_flo_no)
    y=zeros(wave_flo_no)
    cdf=zeros(wave_flo_no)
    wave_flo= arange(flo_wl_start, flo_wl_end )  
    wl_arr = zeros((len(wave_flo), 1))
    wl_arr[:, 0] = wave_flo[:]
    result_arr = interp(wl_arr, data1,data2)
    for i in range(len(wl_arr)):
        x[i]=(i)/(wave_flo_no-1)   
        y[i]=result_arr[i, 0]
    y = y / (trapz_not_numpy(y, x))
    for i in range(1,wave_flo_no):
        cdf[i] = trapz_not_numpy(y[:i+1], x[:i+1])   
    if convert_to_nm == True:
        wave_flo = wave_flo *1000
    return cdf, wave_flo  

def get_number(string):
    number = 0
    for i in range(len(string)):
        if string[i] == ':':
            number = string[:i]
            break
    return number, len(number)

def check_num_exists(num, p_num, check, type):
    local_check = 1
    for i in range(len(p_num)):
        if num == int(p_num[i]):
            local_check = 0
    if local_check == 1:
        print(type + ' ' + str(num) + ' ' + 'does not exist in the header.')
        print("Please re-enter input file once corrected.")
        print("\n")
        check = True
    return check

def check_file_exists(name, check):
    if not os.path.exists(name):
        check = True
        print("File " + name + " does not exist in this directory.")
        print("Please re-enter input file once corrected.")
        print("\n")
    return check

def check_material_wavelength_range(particle, medium, check, start, end):
    if start == 0 and end == 0:
        return check
    min_wave = particle[0, 0, 0]
    max_wave = max(particle[0, :, 0])
    for i in range(1, len(particle[:, 0, 0])):
        if particle[i, 0, 0] > min_wave:
            min_wave = particle[i, 0, 0]
        if max(particle[i, :, 0]) < max_wave:
            max_wave = max(particle[i, :, 0])

    for i in range(len(medium[:, 0, 0])):
        if medium[i, 0, 0] > min_wave:
            min_wave = medium[i, 0, 0]
        if max(medium[i, :, 0]) < max_wave:
            max_wave = max(medium[i, :, 0])
    if start < min_wave or end > max_wave:
        check = True
        print('Wavelength range is not covered by all input materials')
        print("Please re-enter input file once corrected.")
        print("\n")
    return check

def import_header(infile,check):
    p = zeros(0, dtype=str)
    m = zeros(0, dtype=str)
    p_uncert = zeros(0, dtype=str)
    m_uncert = zeros(0, dtype=str)

    p_num = zeros(0, dtype=str)
    m_num = zeros(0, dtype=str)

    p_num_uncert = zeros(0, dtype=str)
    m_num_uncert = zeros(0, dtype=str)

    emissio =  zeros(0, dtype=str)
    emissio_int =  zeros(0, dtype=int)
    Quantum_Yiel_all = zeros(0, dtype=str)
    Quantum_Yiel_all_int =zeros(0, dtype=int)

    Quantum_Yiel_all_uncert = zeros(0, dtype=str)
    Quantum_Yiel_all_int_uncert =zeros(0, dtype=int)

    values_array = zeros(0, dtype=float)

    fluor_counting_all_Mie =  zeros(0, dtype=int)

    #solar = ""
    output_name = ""
    mesh_percentage = 1
    photons = 0
    line_d = 0
    count = 0
    count_wave = 0
    number_of_fluor = 0
    uncert = False
    sims_uncert = 0
    # def get_number(string):
    # number = 0
    # for i in range(len(string)):
    #     if string[i] == ':':
    #         number = string[:i]
    #         break
    # return number, len(number)
    for i in range(len(infile)):
        if infile[i][0:8] == "particle":
            num, length = get_number(infile[i][8:])
            p_num = append(p_num, num)
            p = append(p, infile[i][(9+length):])

        if infile[i][0:15] == "uncert_particle":
            uncert = True
            num, length = get_number(infile[i][15:])
            p_num_uncert = append(p_num_uncert, num)
            p_uncert = append(p_uncert, infile[i][(16+length):])
        if infile[i][0:6] == "matrix":                                        
            num, length = get_number(infile[i][6:])
            m_num = append(m_num, num)
            m = append(m, infile[i][7+length:])
        if infile[i][0:13] == "uncert_matrix": 
            uncert = True                                       
            num, length = get_number(infile[i][13:])
            m_num_uncert = append(m_num_uncert, num)
            m_uncert = append(m_uncert, infile[i][14+length:])
        if infile[i][0:6] == "output":
            output_name = infile[i][7:]
        if infile[i][0:4] == "emit":
            emissio = append( emissio , infile[i][6:])
            emissio_int = append(emissio_int,int(infile[i][4])) 
            count += 1
        if infile[i][0:15] == "excit_start_end":
            Wavelengths = str(infile[i][17:])
            first, second = map(float, Wavelengths.split(','))
            values_array = np.append(values_array, [first, second])
            number_of_fluor += 1
            count_wave += 1
            fluor_counting_all_Mie = append(fluor_counting_all_Mie,int(infile[i][15]))
        if infile[i][0:2] == "qy":
            Quantum_Yiel_all = append(Quantum_Yiel_all , infile[i][4:])
            Quantum_Yiel_all_int = append(Quantum_Yiel_all_int,int(infile[i][2])) 
        if infile[i][0:9] == "uncert_qy":
            uncert = True 
            Quantum_Yiel_all_uncert = append(Quantum_Yiel_all_uncert , infile[i][11:])
            Quantum_Yiel_all_int_uncert = append(Quantum_Yiel_all_int_uncert,int(infile[i][9]))    #Number_of_sims_uncertainty
        if infile[i][0:5] == "light":
            solar = infile[i][6:]
        if infile[i][0:5] == "start":
            Start = float(infile[i][6:])
            if Start < 30:
                Start = Start*1000
        if infile[i][0:3] == "end":
            End = float(infile[i][4:])
            if End < 30:
                End = End*1000
        if infile[i][0:26] == "number_of_sims_uncertainty":
            sims_uncert = int(infile[i][27:])

        if infile[i][0:8] == "interval":
            Interval = float(infile[i][9:])
            if Interval < 0.1:
                Interval = Interval*1000
        if infile[i][0:4] == "mesh":
            mesh_percentage = float(infile[i][5:])
        if infile[i][0:7] == "photons":
            photons = int(infile[i][8:])
        if infile[i][0:3] == "sim":
            line_d = i-1
            break
    if count != count_wave:
        print(' Make sure the parameter name is the same as shown in Github example like "emit1" for the fluorescent particle 1, "emit2" for fluorescent particle 2, etc.. . Also, Please provide the start excitation ("Start") and end excitation ("End") wavelengths for the fluorescent particle to aviod unrealstic emission processes')
        main_func()
    if Start > End: 
        print(' The Start wavelength should be lower as the simulation start form lower wavelengths to higher wavelengths. Make sure the "Start" is lower wavelength and "End" is higher wavelength')
        main_func()
    count = 0
    count_wave = 0
    sims = 0
    for i in range(line_d, len(infile)):
        if infile[i][0:3] == "sim":
            sims += 1
        if infile[i][0:8] == "particle":
            particle_num = int(infile[i][8:])
            check = check_num_exists(particle_num, p_num, check, 'Particle')
        if infile[i][0:6] == "matrix":
            matrix_num = int(infile[i][6:])
            check = check_num_exists(matrix_num, m_num, check, 'Matrix')
    length = 0
    for i in range(len(p)):
        check = check_file_exists(p[i], check)
        if check == False:
            temp = loadtxt(p[i])
            if len(temp) > length:
                length = len(temp)
    for i in range(len(m)):
        check = check_file_exists(m[i], check)
        if check == False:
            temp = loadtxt(m[i])
            if len(temp) > length:
                length = len(temp)
    particle = zeros((len(p), length, 4))
    particle_uncert = p_uncert # zeros((len(p_uncert), length, 3))                                         
    medium = zeros((len(m), length, 4))
    medium_uncert = m_uncert # zeros((len(m_uncert), length, 3))
    particle_type = zeros(len(p))
    need_interp = False
    if check == True : 
        print("Material optical properteis are not found. Make sure the files are in the same directory as this program.")
        main_func()
    if check == False:
        for i in range(len(p)):
            temp = loadtxt(p[i])
            if max(temp[:,0]) < 30: 
                temp[:,0] = temp[:,0]*1000
            if len(temp[0, :]) == 4:
                particle_type[i] = 1
            particle[i, :len(temp), :len(temp[0, :])] = temp
            if len(temp) != length:
                need_interp = True
            for j in range(length):
                if particle[i, j, 0] < (particle[0, j, 0] - 0.001) or particle[i, j, 0] > (particle[0, j, 0] + 0.001):
                    need_interp = True
        for i in range(len(m)):
            temp = loadtxt(m[i])
            if max(temp[:,0]) < 30: 
                temp[:,0] = temp[:,0]*1000
            medium[i, :len(temp), :len(temp[0, :])] = temp
            if len(temp) != length:
                need_interp = True
            for j in range(length):
                if medium[i, j, 0] < (particle[0, j, 0] - 0.001) or medium[i, j, 0] > (particle[0, j, 0] + 0.001):
                    need_interp = True
        check = check_material_wavelength_range(particle, medium, check, Start, End)
        if check == True : 
            print("Material optical properteis files are not found. Make sure the files are in the same directory as this program.")
            main_func()
        if check == False:
            if need_interp is True:
                print('Interpolating properties to match wavelengths for each input')
                particle, medium, index = interpolatee(particle, medium, length, mesh_percentage,Start, End, Interval)  
                particle[:,:,0] = particle[:,:,0]/1000
                medium[:,:,0] = medium[:,:,0]/1000
                index = index/1000
    return sims_uncert,Quantum_Yiel_all_uncert,Quantum_Yiel_all_int_uncert,fluor_counting_all_Mie,particle, medium, output_name, sims, photons, line_d, index, Start, solar,Quantum_Yiel_all,emissio,emissio_int,Quantum_Yiel_all_int,values_array,Start,Interval,number_of_fluor,End, p_num, m_num, particle_type,particle_uncert,medium_uncert,uncert,p_num_uncert,m_num_uncert,p_uncert,m_uncert

def check_diameters(current_sim, fv, sizes, check, ptype_in, dist):
    if ptype_in == 0: 
        if len(fv) != len(sizes):
            print('Number of diameters does not match number of volume fractions provided in sim:', current_sim)
            print("Please re-enter input file once corrected.")
            check = True
        if dist[0] != 0 or len(dist) > 1:
            if len(fv) != len(dist):
                print('Number of diameters does not match number of standard deviations provided in sim:', current_sim)
                print("Please re-enter input file once corrected.")
                check = True
    return check

def get_index(number, type, p_num, m_num):
    if type == 'p':
        for i in range(len(p_num)):
            if int(p_num[i]) == number:
                index = i
                break
    elif type == 'm':
        for i in range(len(m_num)):
            if int(m_num[i]) == number:
                index = i
                break
    return index

def check_dist(current_sim, dist, check):
    if all(dist) != 0:
        print('Std must be 0 for core shell sim', current_sim)
        print("Please re-enter input file once corrected.")
        check = True
    return check

def create_prop_array(upper_is_air,lower_is_air,prop, particle, upper, lower, upper_type, lower_type, layers, optical_per_layer, medium,number_of_fluor):
    if upper_is_air == True:
        upper = -1
    else:
        upper = 100
    if lower_is_air == True: 
        lower = -1
    else: 
        lower = 100  
    if upper == -1 and lower == -1:
        for j in range(len(particle[0, :, 0])):
            prop = vstack((prop, [1, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            for k in range(layers):
                prop = vstack((prop, optical_per_layer[:, j, k]))
            prop = vstack((prop, [1, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            prop = vstack((prop, [0, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor ))
    elif upper != -1 and lower == -1:
        for j in range(len(particle[0, :, 0])):
            prop = vstack((prop, [medium[int(upper_type), j, 1], medium[int(upper_type), j, 2], 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            for k in range(layers):
                prop = vstack((prop, optical_per_layer[:, j, k]))
            prop = vstack((prop, [1, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            prop = vstack((prop, [0, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor ))
    elif upper == -1 and lower != -1:
        for j in range(len(particle[0, :, 0])):
            prop = vstack((prop, [1, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            for k in range(layers):
                prop = vstack((prop, optical_per_layer[:, j, k]))
            prop = vstack((prop, [medium[int(lower_type), j, 1], medium[int(lower_type), j, 2], 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            prop = vstack((prop, [0, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor ))
    elif upper != -1 and lower != -1:
        for j in range(len(particle[0, :, 0])):
            prop = vstack((prop, [medium[int(upper_type), j, 1], medium[int(upper_type), j, 2], 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            for k in range(layers):
                prop = vstack((prop, optical_per_layer[:, j, k]))
            prop = vstack((prop, [medium[int(lower_type), j, 1], medium[int(lower_type), j, 2], 0,0, 0, 0,0,0] +[0]* number_of_fluor))
            prop = vstack((prop, [0, 0, 0,0, 0, 0,0,0] +[0]* number_of_fluor ))
    return prop

def optical_uncert(sims_uncert,Quantum_Yiel_all_uncert,Quantum_Yiel_all_int_uncert,fluor_counting_all_Mie,line, infile, particle, medium, check, start_wl, index,Quantum_Yiel_all,emissio,emissio_int,Quantum_Yiel_all_int,values_array,Start,Interval,number_of_fluor, p_num, m_num, particle_type,END,particle_uncert,medium_uncert,p_num_uncert,m_num_uncert,p_uncert,m_uncert):

    prop = zeros((0, 8 + int(number_of_fluor)))   
    optics_sum = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0])))   
    optical_per_layer = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0]), 0))    
    length_per_sim = zeros(0, dtype=int)
    emissio_intt =  zeros(0, dtype=int)
    organi_lay_emis =  zeros(0, dtype=int)
    layers_per_sim = zeros(0, dtype=int)
    emission_all_sim = 0
    qy_all_sim = 0
    vol_frac_sum = 0
    layers = 1
    count = 0
    fluor=0
    current_sim = 1
    happens = 0
    again = 0
    dist = zeros(1)    
    upper = -1
    lower = -1
    upper_type = 0
    lower_type = 0
    upper_is_air = False
    lower_is_air = False
    counting_sims_uncert = 0

    h = 0
    with open('all_sampled_data.txt', 'w') as f:
        for i in range(sims_uncert):
            h += 1
            count = 0
            for ii in range(line+1, len(infile)):
                count += 1
                if h ==1 and count == 1:
                    f.write(infile[ii] + '\n')
                if h !=1 and count == 1:
                    f.write('sim:' + str(i+1) + '\n')
                if count != 1:
                    f.write(infile[ii] + '\n')
    f.close()
    infile = loadtxt('all_sampled_data.txt', comments="#", dtype=str, delimiter="/")
    p_num_uncert = np.sort(p_num_uncert)
    m_num_uncert = np.sort(m_num_uncert)

    count = 0
    stadard_devi = zeros((0,2))
    for i in range(len(infile)): 
        if infile[i][0:8] == "particle":   #,p_uncert,m_uncert
            if count > 0:
                if fluor ==1 :
                    emissio_intt = append(emissio_intt, ptype)
                ptype_in = get_index(ptype, 'p', p_num, m_num)
                check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
                mtype_in = get_index(mtype, 'm', p_num, m_num)
                stadard_devi_part = zeros((0,2))
                stadard_devi_med = zeros((0,2))
                l = 0
                v = 0
                for k in range(len(p_uncert)):
                    if ptype == int(p_num_uncert[k]) :
                        l = 1
                        stadard_deviatio_ori_p = loadtxt(p_uncert[k])
                        if len(stadard_deviatio_ori_p) > 2 :
                            if max(stadard_deviatio_ori_p[:,0]) < 30:
                                stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori_p[:,0]*1000, stadard_deviatio_ori_p[:,1])
                                stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori_p[:,0]*1000, stadard_deviatio_ori_p[:,2])
                            else:
                               stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori_p[:,0], stadard_deviatio_ori_p[:,1])
                               stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori_p[:,0], stadard_deviatio_ori_p[:,2]) 
                        else:
                            stadard_devi_part = stadard_deviatio_ori_p
                    else:
                        stadard_devi_part = 0
                for k in range(len(m_uncert)):
                    if mtype == int(m_num_uncert[k]) :
                        v = 1
                        stadard_deviatio_ori_m = loadtxt(m_uncert[k])
                        if len(stadard_deviatio_ori_m) > 2 :
                            if max(stadard_deviatio_ori_m[:,0]) < 30:
                                stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori_m[:,0]*1000, stadard_deviatio_ori_m[:,1])
                                stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori_m[:,0]*1000, stadard_deviatio_ori_m[:,2])
                            else:
                               stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori_m[:,0], stadard_deviatio_ori_m[:,1])
                               stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori_m[:,0], stadard_deviatio_ori_m[:,2]) 
                        else:
                            stadard_devi_med = stadard_deviatio_ori_m
                    else:
                        stadard_devi_med = 0
                if l==1 : 
                    if  len(stadard_deviatio_ori_p) > 2 :
                        term_p = stadard_devi_part.copy()
                        for many in range(len(stadard_devi_part[:,0])):
                            term_p[many,0] = stadard_devi_part[many,0]* erfinv(2*random.random_sample()-1 )
                            term_p[many,1] = stadard_devi_part[many,1]* erfinv(2*random.random_sample()-1 )
                    else:
                        if random.random_sample() > 0.5: 
                            term_p = stadard_devi_part * erfinv(random.random_sample() )
                        else:
                            term_p = -stadard_devi_part * erfinv(random.random_sample() )
                if v ==1:
                    if  len(stadard_deviatio_ori_m) > 2 :
                        term_m = stadard_devi_med.copy()
                        for many in range(len(stadard_devi_med[:,0])):
                            term_m[many,0] = stadard_devi_med[many,0]* erfinv(2*random.random_sample()-1 )
                            term_m[many,1] = stadard_devi_med[many,1]* erfinv(2*random.random_sample()-1 )
                    else:
                        if random.random_sample() > 0.5: 
                            term_m = stadard_devi_med * erfinv(random.random_sample() )
                        else:
                            term_m = -stadard_devi_med * erfinv(random.random_sample() )
                temp_particle = particle[int(ptype_in), :, :].copy()
                temp_medium = medium[int(mtype_in), :, :].copy()
                if l ==1:
                    if len(stadard_deviatio_ori_p) > 2:
                        temp_particle[:,1] += term_p[:,0]
                        temp_particle[:,2] += term_p[:,1]
                    else:
                        temp_particle[:,1] += term_p[0]
                        temp_particle[:,2] += term_p[1]
                    
                if v ==1:              
                    if len(stadard_deviatio_ori_m) > 2:
                        temp_medium[:,1] += term_m[:,0]
                        temp_medium[:,2] += term_m[:,1]
                    else:
                        temp_medium[:,1] += term_m[0]
                        temp_medium[:,2] += term_m[1] 
                optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv, temp_particle, temp_medium, thickness, dist,fluor,  start_wl, index,happens, again,number_of_fluor,particle_type[int(ptype_in)])  
                optics_sum += optics
                dist = zeros(1)
            count += 1 
                               
            ptype = int(infile[i][8:])
          
        if infile[i][0:5] == "upper":    
            if len(infile[i][6:]) <2: 
                upper = int(infile[i][6:])
                upper_is_air = True
            else:
                upper = int(infile[i][12:])
                upper_type = get_index(upper, 'm', p_num, m_num)
        elif infile[i][0:5] == "lower":
            if len(infile[i][6:]) <2: 
                lower = int(infile[i][6:])
                lower_is_air = True
            else:
                lower = int(infile[i][12:])
                lower_type = get_index(lower, 'm', p_num, m_num)
        elif infile[i][0:6] == "matrix":
            mtype = int(infile[i][6:])
        elif infile[i][0:2] == "t:" or infile[i][0:2] == "t=":
            thickness = float(infile[i][2:])
            thickness = thickness / 10000
        elif infile[i][0:6] == "fluor:":    
            fluor = int(infile[i][6:])
            if fluor == 1 :
                happens = 1
                again += 1 
        elif infile[i][0:2] == "d:" or infile[i][0:2] == "d=":
            sizes = (infile[i][2:]).split(",")
            sizes = asarray(sizes, dtype=float)
            sizes = sizes/2
        elif infile[i][0:3] == "vf:":
            fv = (infile[i][3:]).split(",")
            fv = asarray(fv, dtype=float)
            fv = fv/100
            vol_frac_sum += sum(fv)
        elif infile[i][0:4] == "std:" or infile[i][0:4] == "std=":
            dist = (infile[i][4:]).split(",")
            dist = asarray(dist, dtype=float)
            if dist[0] >1:
                dist = dist/100
        elif infile[i][0:5] == "layer" and infile[i][5:] != "1":
            stadard_devi = zeros((0,2))
            check = check_dist(current_sim, dist, check)
            ptype_in = get_index(ptype, 'p', p_num, m_num)
            check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
            layers += 1
            mtype_in = get_index(mtype, 'm', p_num, m_num)
            if fluor == 1:
                emissio_intt = append(emissio_intt, ptype)
            emissio_intt = append(emissio_intt, 0)   
            stadard_devi_part = zeros((0,2))
            stadard_devi_med = zeros((0,2))
            l = 0
            v = 0
            for k in range(len(p_uncert)):
                if ptype == int(p_num_uncert[k]) :
                    l =1
                    stadard_deviatio_ori_p = loadtxt(p_uncert[k])
                    if len(stadard_deviatio_ori_p) > 2 :
                        if max(stadard_deviatio_ori[:,0]) < 30:
                            stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori_p[:,0]*1000, stadard_deviatio_ori_p[:,1])
                            stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori_p[:,0]*1000, stadard_deviatio_ori_p[:,1])
                        else:
                            stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori[:,0], stadard_deviatio_ori_p[:,2])
                            stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori_p[:,0], stadard_deviatio_ori_p[:,2]) 
                    else:
                        stadard_devi_part = stadard_deviatio_ori_p
                else:
                    stadard_devi_part = 0
            for k in range(len(m_uncert)):
                if mtype == int(m_num_uncert[k]) :
                    v = 1
                    stadard_deviatio_ori_m = loadtxt(m_uncert[k])
                    if len(stadard_deviatio_ori_m) > 2 :
                        if max(stadard_deviatio_ori[:,0]) < 30:
                            stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori_m[:,0]*1000, stadard_deviatio_ori_m[:,1])
                            stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori_m[:,0]*1000, stadard_deviatio_ori_m[:,1])
                        else:
                            stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori_m[:,0], stadard_deviatio_ori_m[:,2])
                            stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori_m[:,0], stadard_deviatio_ori_m[:,2]) 
                    else:
                        stadard_devi_med = stadard_deviatio_ori_m
                else:
                    stadard_devi_med = 0
            if l ==1:
                if  len(stadard_deviatio_ori_p) > 2 :
                    term_p = stadard_devi_part.copy()
                    for many in range(len(stadard_devi_part[:,0])):
                        term_p[many,0] = stadard_devi_part[many,0]* erfinv(2*random.random_sample()-1 )
                        term_p[many,1] = stadard_devi_part[many,1]* erfinv(2*random.random_sample()-1 )
                else:
                    if random.random_sample() > 0.5: 
                        term_p = stadard_devi_part * erfinv(random.random_sample() )
                    else:
                        term_p = -stadard_devi_part * erfinv(random.random_sample() )
            if v ==1:
                if  len(stadard_deviatio_ori_m) > 2 :
                    term_m = stadard_devi_med.copy()
                    for many in range(len(stadard_devi_med[:,0])):
                        term_m[many,0] = stadard_devi_med[many,0]* erfinv(2*random.random_sample()-1 )
                        term_m[many,1] = stadard_devi_med[many,1]* erfinv(2*random.random_sample()-1 )
                else:
                    if random.random_sample() > 0.5: 
                        term_m = stadard_devi_med * erfinv(random.random_sample() )
                    else:
                        term_m = -stadard_devi_med * erfinv(random.random_sample() )
            temp_particle = particle[int(ptype_in), :, :].copy()
            temp_medium = medium[int(mtype_in), :, :].copy()
            if l==1:
                if len(stadard_deviatio_ori_p) > 2:
                    temp_particle[:,1] += term_p[:,0]
                    temp_particle[:,2] += term_p[:,1]
                else:
                    temp_particle[:,1] += term_p[0]
                    temp_particle[:,2] += term_p[1]       
            if v==1:            
                if len(stadard_deviatio_ori_m) > 2:
                    temp_medium[:,1] += term_m[:,0]
                    temp_medium[:,2] += term_m[:,1]
                else:
                    temp_medium[:,1] += term_m[0]
                    temp_medium[:,2] += term_m[1]   
            optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv,temp_particle, temp_medium, thickness, dist, fluor, start_wl, index,happens,again,number_of_fluor,particle_type[int(ptype_in)])
            optics_sum += optics 
            optics_sum[0, :] = optics[0, :]
            optics_sum[1, :] = optics[1, :]
            optics_sum[6, :] = optics[6, :]
            optics_sum[7, :] = optics[7, :]
            optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype_in), :, :],fluor,number_of_fluor,particle_type[int(ptype_in)],ptype_in,fluor_counting_all_Mie)
            optical_per_layer = dstack((optical_per_layer, optics_sum))
            vol_frac_sum = 0
            optics_sum = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0])))     ### "5"
            happens = 0
            again = 0
            count = 0
            dist = zeros(1)
        elif infile[i][0:4] == "sim:" and infile[i][4:] != "1":
            stadard_devi = zeros((0,2))
            counting_sims_uncert += 1
            current_sim = int(infile[i][4:])  
            check = check_dist(current_sim, dist, check)
            ptype_in = get_index(ptype, 'p', p_num, m_num)
            check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
            mtype_in = get_index(mtype, 'm', p_num, m_num)
            stadard_devi_part = zeros((0,2))
            stadard_devi_med = zeros((0,2))
            l = 0
            v = 0
            for k in range(len(p_uncert)):
                if ptype == int(p_num_uncert[k]) :
                    l=1
                    stadard_deviatio_ori_p = loadtxt(p_uncert[k])
                    if len(stadard_deviatio_ori_p) > 2 :
                        if max(stadard_deviatio_ori_p[:,0]) < 30:
                            stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori_p[:,0]*1000, stadard_deviatio_ori_p[:,1])
                            stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori_p[:,0]*1000, stadard_deviatio_ori_p[:,1])
                        else:
                            stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori_p[:,0], stadard_deviatio_ori_p[:,2])
                            stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori_p[:,0], stadard_deviatio_ori_p[:,2]) 
                    else:
                        stadard_devi_part = stadard_deviatio_ori_p
                else:
                    stadard_devi_part = 0
            for k in range(len(m_uncert)):
                if mtype == int(m_num_uncert[k]) :
                    v =1
                    stadard_deviatio_ori_m = loadtxt(m_uncert[k])
                    if len(stadard_deviatio_ori_m) > 2 :
                        if max(stadard_deviatio_ori_m[:,0]) < 30:
                            stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori_m[:,0]*1000, stadard_deviatio_ori_m[:,1])
                            stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori_m[:,0]*1000, stadard_deviatio_ori_m[:,1])
                        else:
                            stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori_m[:,0], stadard_deviatio_ori_m[:,2])
                            stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori_m[:,0], stadard_deviatio_ori_m[:,2]) 
                    else:
                        stadard_devi_med = stadard_deviatio_ori_m
                else:
                    stadard_devi_med = 0
            if l ==1: 
                if  len(stadard_deviatio_ori_p) > 2 :
                    term_p = stadard_devi_part.copy()
                    for many in range(len(stadard_devi_part[:,0])):
                        term_p[many,0] = stadard_devi_part[many,0]* erfinv(2*random.random_sample()-1 )
                        term_p[many,1] = stadard_devi_part[many,1]* erfinv(2*random.random_sample()-1 )
                else:
                    if random.random_sample() > 0.5: 
                        term_p = stadard_devi_part * erfinv(random.random_sample() )
                    else:
                        term_p = -stadard_devi_part * erfinv(random.random_sample() )
            if v ==1:
                if  len(stadard_deviatio_ori_m) > 2 :
                    term_m = stadard_devi_med.copy()
                    for many in range(len(stadard_devi_med[:,0])):
                        term_m[many,0] = stadard_devi_med[many,0]* erfinv(2*random.random_sample()-1 )
                        term_m[many,1] = stadard_devi_med[many,1]* erfinv(2*random.random_sample()-1 )
                else:
                    if random.random_sample() > 0.5: 
                        term_m = stadard_devi_med * erfinv(random.random_sample() )
                    else:
                        term_m = -stadard_devi_med * erfinv(random.random_sample() )
            temp_particle = particle[int(ptype_in), :, :].copy()
            temp_medium = medium[int(mtype_in), :, :].copy()
            if l==1: 
                if len(stadard_deviatio_ori_p) > 2:
                    temp_particle[:,1] += term_p[:,0]
                    temp_particle[:,2] += term_p[:,1]
                else:
                    temp_particle[:,1] += term_p[0]
                    temp_particle[:,2] += term_p[1]
            if v==1:                   
                if len(stadard_deviatio_ori_m) > 2:
                    temp_medium[:,1] += term_m[:,0]
                    temp_medium[:,2] += term_m[:,1]
                else:
                    temp_medium[:,1] += term_m[0]
                    temp_medium[:,2] += term_m[1]  

            optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv, temp_particle, temp_medium, thickness, dist,fluor, start_wl, index,happens,again,number_of_fluor,particle_type[int(ptype_in)])  ##### the colon ":" is the index which values are taken and inserted in thos function
            optics_sum += optics
            optics_sum[0, :] = optics[0, :]
            optics_sum[1, :] = optics[1, :]
            optics_sum[6, :] = optics[6, :]
            optics_sum[7, :] = optics[7, :]
            temp_particle = zeros((0,2))
            temp_medium = zeros((0,2))

            optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype_in), :, :],fluor,number_of_fluor,particle_type[int(ptype_in)],ptype_in,fluor_counting_all_Mie)
            optical_per_layer = dstack((optical_per_layer, optics_sum))
            prop = create_prop_array(upper_is_air,lower_is_air,prop, particle, upper, lower, upper_type, lower_type, layers, optical_per_layer, medium,number_of_fluor)
            vol_frac_sum = 0
            optics_sum = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0])))   
            optical_per_layer = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0]), 0))  
            length_per_sim = append(length_per_sim,(len(prop[:,0])))
            length_after = 0
            for k in range(len(emissio_int)):
                temp = loadtxt(emissio[k-1])
                if len(arange(temp[0,0], temp[-1,0] ) ) > length_after:
                    length_after = len(arange(temp[0,0], temp[-1,0] ) )
            emission_all = zeros((number_of_fluor*2*layers, length_after +1)) 
            qy_all = zeros((number_of_fluor*layers, len(index)))
            emissio_intt = append(emissio_intt, 0)
            i = -1
            for kk in range(len(emissio_intt)): 
                if emissio_intt[kk] ==0:
                    i += 1
                    emissio_intt_temp =np.sort(organi_lay_emis)
                    z = int(emissio_intt[kk])
                    for j in range(len(emissio_intt_temp)): 
                        for kkk in range(len(emissio_int)):
                            if emissio_intt_temp[j] == emissio_int[kkk]:
                                temp = loadtxt(emissio[kkk]) 
                                cdf, wave_flo=inv_cdf(temp)
                                if max(wave_flo) > END or min(wave_flo) < Start:
                                    print('Please make sure in the input file, the "Start" and "End" simulated wavelengths are both lower and higher than the emission wavelengths (start and end wavelengths of an emission band of fluorescent particle), respectively')
                                    main_func()
                                e = len(cdf)
                                emission_all[2*kkk + number_of_fluor*2*(i), 0:e] = cdf[:]
                                emission_all[2*kkk+ 1 + number_of_fluor*2*(i), 0:e] = wave_flo 
                                emission_all[2*kkk+ number_of_fluor*2*(i) ,-1] = e
                                temp_2 = loadtxt(Quantum_Yiel_all[kkk]) # Quantum_Yiel_all_uncert,Quantum_Yiel_all_int_uncert,
                                nc = 0
                                for g in range(len(Quantum_Yiel_all_int_uncert)):
                                    if Quantum_Yiel_all_int_uncert[g] == emissio_int[kkk]:
                                        nc = 1
                                        temp_2_uncert =  loadtxt(Quantum_Yiel_all_uncert[g])
                                if (temp_2.size) >1 : 
                                    if max(temp_2[:,0]) < 30: 
                                        new_QY=interp(index[:], temp_2[:,0]*1000, temp_2[:,1])
                                    else:
                                        new_QY=interp(index[:], temp_2[:,0], temp_2[:,1])
                                    if nc ==1:
                                        if (temp_2_uncert.size) >1 :
                                            if max(temp_2_uncert[:,0]) < 30: 
                                                new_QY_uncert = interp(index[:], temp_2_uncert[:,0]*1000, temp_2_uncert[:,1])
                                            else:
                                                new_QY_uncert=interp(index[:], temp_2_uncert[:,0], temp_2_uncert[:,1])
                                            term_qy = new_QY_uncert
                                            for many_qy in range(len(new_QY_uncert)):
                                                term_qy[many_qy] = new_QY_uncert[many_qy] * erfinv(2*random.random_sample()-1)
                                            final = new_QY + term_qy
                                            for k in range(len(final)) :
                                                if final[k] >1:
                                                    final[k] = 1
                                                if final[k] < 0:
                                                    final[k] = 0
                                        else:
                                            new_QY_uncert = temp_2_uncert
                                            if random.random_sample() > 0.5:
                                                term_qy = new_QY_uncert * erfinv(random.random_sample() )
                                            else:
                                                term_qy = -new_QY_uncert * erfinv(random.random_sample() )
                                            final = new_QY + term_qy
                                            if final >1:
                                                final = 1
                                            if final < 0:
                                                final = 0                                        
                                        qy_all[kkk+ number_of_fluor*(i), int(values_array[2*(kkk)]):int(values_array[(2*(kkk))+1])] = final
                                    else:
                                        qy_all[kkk+ number_of_fluor*(i), int(values_array[2*(kkk)]):int(values_array[(2*(kkk))+1])] = new_QY
                                else:
                                    if int(values_array[2*kkk]) < Start or int(values_array[(2*kkk ) +1 ]) >  END :
                                        print('Please make sure the start and end excitation wavelengths are higher and lower than the simulated start and end wavelegnths, respectively ')
                                        main_func()
                                    else:
                                        if nc ==1:
                                            new_QY_uncert = temp_2_uncert
                                            if random.random_sample() > 0.5:
                                                term_qy = new_QY_uncert * erfinv(random.random_sample() )
                                            else:
                                                term_qy = -new_QY_uncert * erfinv(random.random_sample() )
                                            final = temp_2 + term_qy
                                            if final >1:
                                                final = 1
                                            if final < 0:
                                                final = 0     
                                            qy_all[kkk + number_of_fluor*(i), int((int(values_array[2*kkk])-Start)/Interval):int((int(values_array[(2*kkk)+1])-Start)/Interval)] = final #constant_qy + term_qy
                                        else:
                                            qy_all[kkk + number_of_fluor*(i), int((int(values_array[2*kkk])-Start)/Interval):int((int(values_array[(2*kkk)+1])-Start)/Interval)] = temp_2 #constant_qy + term_qy

                    organi_lay_emis = zeros(0, dtype=int)
                else:
                    organi_lay_emis =append(organi_lay_emis,int(emissio_intt[kk])) 
            try:
                emission_all_sim = vstack((emission_all_sim,emission_all))
            except:
                emission_all_sim =  emission_all.copy()
            try:
                qy_all_sim = vstack((qy_all_sim,qy_all))
            except:
                qy_all_sim  = qy_all.copy()
            dist = zeros(1)
            layers_per_sim = append(layers_per_sim,layers)
            layers = 1
            happens = 0
            again = 0
            count = 0
            emissio_intt =  zeros(0, dtype=int)
            upper_is_air = False
            lower_is_air = False
    
    stadard_devi = zeros((0,2))
    check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
    ptype_in = get_index(ptype, 'p', p_num, m_num)
    check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
    mtype_in = get_index(mtype, 'm', p_num, m_num)
    if fluor ==1 :
        emissio_intt = append(emissio_intt, ptype)
    
    stadard_devi_part = zeros((0,2))
    stadard_devi_med = zeros((0,2))
    l = 0
    v = 0
    
    for k in range(len(p_uncert)):
        if ptype == int(p_num_uncert[k]) :
            l =1
            stadard_deviatio_ori = loadtxt(p_uncert[k])
            if len(stadard_deviatio_ori) > 2 :
                if max(stadard_deviatio_ori[:,0]) < 30:
                    stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori[:,0]*1000, stadard_deviatio_ori[:,1])
                    stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori[:,0]*1000, stadard_deviatio_ori[:,1])
                else:
                    stadard_devi_part[:,0]=interp(index[:], stadard_deviatio_ori[:,0], stadard_deviatio_ori[:,2])
                    stadard_devi_part[:,1]=interp(index[:], stadard_deviatio_ori[:,0], stadard_deviatio_ori[:,2]) 
            else:
                stadard_devi_part = stadard_deviatio_ori
        else:
            stadard_devi_part = 0
    for k in range(len(m_uncert)):
        if mtype == int(m_num_uncert[k]) :
            v= 1
            stadard_deviatio_ori = loadtxt(m_uncert[k])
            if len(stadard_deviatio_ori) > 2 :
                if max(stadard_deviatio_ori[:,0]) < 30:
                    stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori[:,0]*1000, stadard_deviatio_ori[:,1])
                    stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori[:,0]*1000, stadard_deviatio_ori[:,1])
                else:
                    stadard_devi_med[:,0]=interp(index[:], stadard_deviatio_ori[:,0], stadard_deviatio_ori[:,2])
                    stadard_devi_med[:,1]=interp(index[:], stadard_deviatio_ori[:,0], stadard_deviatio_ori[:,2]) 
            else:
                stadard_devi_med = stadard_deviatio_ori
        else:
            stadard_devi_med = 0
    if l ==1:
        if  len(stadard_deviatio_ori_p) > 2 :
            term_p = stadard_devi_part.copy()
            for many in range(len(stadard_devi_part[:,0])):
                term_p[many,0] = stadard_devi_part[many,0]* erfinv(2*random.random_sample()-1 )
                term_p[many,1] = stadard_devi_part[many,1]* erfinv(2*random.random_sample()-1 )
        else:
            if random.random_sample() > 0.5: 
                term_p = stadard_devi_part * erfinv(random.random_sample() )
            else:
                term_p = -stadard_devi_part * erfinv(random.random_sample() )
    if v ==1:
        if  len(stadard_deviatio_ori_m) > 2 :
            term_m = stadard_devi_med.copy()
            for many in range(len(stadard_devi_med[:,0])):
                term_m[many,0] = stadard_devi_med[many,0]* erfinv(2*random.random_sample()-1 )
                term_m[many,1] = stadard_devi_med[many,1]* erfinv(2*random.random_sample()-1 )
        else:
            if random.random_sample() > 0.5: 
                term_m = stadard_devi_med * erfinv(random.random_sample() )
            else:
                term_m = -stadard_devi_med * erfinv(random.random_sample() )
    temp_particle = particle[int(ptype_in), :, :].copy()
    temp_medium = medium[int(mtype_in), :, :].copy()
    if l==1:
        if len(term_p) > 2:
            temp_particle[:,1] += term_p[:,0]
            temp_particle[:,2] += term_p[:,1]
        else:
            temp_particle[:,1] += term_p[0]
            temp_particle[:,2] += term_p[1]                   
    if v==1:
        if len(term_m) > 2:
            temp_medium[:,1] += term_m[:,0]
            temp_medium[:,2] += term_m[:,1]
        else:
            temp_medium[:,1] += term_m[0]
            temp_medium[:,2] += term_m[1]
    optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv, temp_particle, temp_medium, thickness, dist, fluor, start_wl, index,happens,again,number_of_fluor,particle_type[int(ptype_in)])
    optics_sum += optics
    optics_sum[0, :] = optics[0, :]
    optics_sum[1, :] = optics[1, :]
    optics_sum[6, :] = optics[6, :]
    optics_sum[7, :] = optics[7, :]

    temp_particle = zeros((0,2))
    temp_medium = zeros((0,2))
    
    optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype_in), :, :],fluor,number_of_fluor,particle_type[int(ptype_in)],ptype_in,fluor_counting_all_Mie)     ### notice the particle and medium have all indicies of second and third dimension
    optical_per_layer = dstack((optical_per_layer, optics_sum))
    prop = create_prop_array(upper_is_air,lower_is_air,prop, particle, upper, lower, upper_type, lower_type, layers, optical_per_layer, medium,number_of_fluor)
    length_per_sim = append(length_per_sim,(len(prop[:,0])))
    layers_per_sim = append(layers_per_sim,layers)
    emissio_intt = append(emissio_intt, 0)
    length_after = 0
    for k in range(len(emissio_int)):
        temp = loadtxt(emissio[k-1])
        if len(arange(temp[0,0], temp[-1,0] ) ) > length_after:
            length_after = len(arange(temp[0,0], temp[-1,0] ) )
    emission_all = zeros((number_of_fluor*2*layers, length_after +1)) 
    qy_all = zeros((number_of_fluor*layers, len(index)))
    
    i = -1
    for kk in range(len(emissio_intt)): 
        if emissio_intt[kk] ==0:
            i += 1
            emissio_intt_temp =np.sort(organi_lay_emis)
            for j in range(len(emissio_intt_temp)):
                for kkk in range(len(emissio_int)):
                    if emissio_intt_temp[j] == emissio_int[kkk]:
                        temp = loadtxt(emissio[kkk]) 
                        cdf, wave_flo=inv_cdf(temp)
                        e = len(cdf)
                        if max(wave_flo) > END or min(wave_flo) < Start:
                            print('Please make sure in the input file, the "Start" and "End" simulated wavelengths are both lower and higher than the emission wavelengths (start and end wavelengths of an emission band of fluorescent particle), respectively')
                            main_func()
                        emission_all[2*kkk + number_of_fluor*2*(i), 0:e] = cdf[:]
                        emission_all[2*kkk+ 1 + number_of_fluor*2*(i), 0:e] = wave_flo
                        emission_all[2*kkk+ number_of_fluor*2*(i) ,-1] = e
                        temp_2 = loadtxt(Quantum_Yiel_all[kkk])
                        nc = 0
                        for g in range(len(Quantum_Yiel_all_int_uncert)):
                            if Quantum_Yiel_all_int_uncert[g] == emissio_int[kkk]:
                                nc = 1
                                temp_2_uncert =  loadtxt(Quantum_Yiel_all_uncert[g])
                        if (temp_2.size) >1 : 
                            if max(temp_2[:,0]) < 30: 
                                new_QY=interp(index[:], temp_2[:,0]*1000, temp_2[:,1])
                            else:
                                new_QY=interp(index[:], temp_2[:,0], temp_2[:,1])
                            if nc ==1:
                                if (temp_2_uncert.size) >1 :
                                    if max(temp_2_uncert[:,0]) < 30: 
                                        new_QY_uncert = interp(index[:], temp_2_uncert[:,0]*1000, temp_2_uncert[:,1])
                                    else:
                                        new_QY_uncert=interp(index[:], temp_2_uncert[:,0], temp_2_uncert[:,1])
                                    term_qy = new_QY_uncert
                                    for many_qy in range(len(new_QY_uncert)):
                                        term_qy[many_qy] = new_QY_uncert[many_qy] * erfinv(2*random.random_sample()-1)
 
                                    final = new_QY + term_qy
                                    for k in range(len(final)) :
                                        if final[k] >1:
                                            final[k] = 1
                                        if final[k] < 0:
                                            final[k] = 0
                                else:
                                    new_QY_uncert = temp_2_uncert

                                    if random.random_sample() > 0.5:
                                        term_qy = new_QY_uncert * erfinv(random.random_sample() )
                                    else:
                                        term_qy = -new_QY_uncert * erfinv(random.random_sample() )
                                    final = new_QY + term_qy
                                    if final >1:
                                        final = 1
                                    if final < 0:
                                        final = 0                                        
                                qy_all[kkk+ number_of_fluor*(i), int(values_array[2*(kkk)]):int(values_array[(2*(kkk))+1])] = final
                            else: 
                                qy_all[kkk+ number_of_fluor*(i), int(values_array[2*(kkk)]):int(values_array[(2*(kkk))+1])] = new_QY        
                        else:
                            if int(values_array[2*kkk]) < Start or int(values_array[(2*kkk ) +1 ]) >  END :
                                print('Please make sure the start and end excitation wavelengths are higher and lower than the simulated start and end wavelegnths, respectively ')
                                main_func()
                            else:
                                if nc ==1:
                                    new_QY_uncert = temp_2_uncert
                                    if random.random_sample() > 0.5:
                                        term_qy = new_QY_uncert * erfinv(random.random_sample() )
                                    else:
                                        term_qy = -new_QY_uncert * erfinv(random.random_sample() )
                                    final = temp_2 + term_qy
                                    if final >1:
                                        final = 1
                                    if final < 0:
                                        final = 0     
                                    qy_all[kkk + number_of_fluor*(i), int((int(values_array[2*kkk])-Start)/Interval):int((int(values_array[(2*kkk)+1])-Start)/Interval)] = final #constant_qy + term_qy
                                else:
                                    qy_all[kkk + number_of_fluor*(i), int((int(values_array[2*kkk])-Start)/Interval):int((int(values_array[(2*kkk)+1])-Start)/Interval)] = temp_2 #constant_qy + term_qy
            organi_lay_emis = zeros(0, dtype=int)
        else:
            organi_lay_emis =append(organi_lay_emis,int(emissio_intt[kk]))
    try: 
        emission_all_sim = vstack((emission_all_sim,emission_all))
    except:
        emission_all_sim = emission_all.copy()
    try:
        qy_all_sim = vstack((qy_all_sim,qy_all))
    except:
        qy_all_sim = qy_all.copy()  
    counting_sims_uncert += 1
    
    return prop, check,qy_all_sim, emission_all_sim,number_of_fluor,length_per_sim,layers_per_sim,counting_sims_uncert

def optical(sims_uncert,Quantum_Yiel_all_uncert,Quantum_Yiel_all_int_uncert,fluor_counting_all_Mie,line, infile, particle, medium, check, start_wl, index,Quantum_Yiel_all,emissio,emissio_int,Quantum_Yiel_all_int,values_array,Start,Interval,number_of_fluor, p_num, m_num, particle_type,END,particle_uncert,medium_uncert,uncert,p_num_uncert,m_num_uncert,p_uncert,m_uncert):
    print("Running Mie theory")
    prop = zeros((0, 8 + int(number_of_fluor)))   
    optics_sum = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0])))   
    optical_per_layer = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0]), 0))    
    length_per_sim = zeros(0, dtype=int)
    emissio_intt =  zeros(0, dtype=int)
    organi_lay_emis =  zeros(0, dtype=int)
    layers_per_sim = zeros(0, dtype=int)
    emission_all_sim = 0
    qy_all_sim = 0
    vol_frac_sum = 0
    layers = 1
    count = 0
    fluor=0
    current_sim = 1
    happens = 0
    again = 0
    dist = zeros(1)    
    upper = -1
    lower = -1
    upper_type = 0
    lower_type = 0
    upper_is_air = False
    lower_is_air = False
    counting_sims_uncert = 0
    if uncert == True :
        prop, check,qy_all_sim, emission_all_sim,number_of_fluor,length_per_sim,layers_per_sim,counting_sims_uncert = optical_uncert(sims_uncert,Quantum_Yiel_all_uncert,Quantum_Yiel_all_int_uncert,fluor_counting_all_Mie,line, infile, particle, medium, check, start_wl, index,Quantum_Yiel_all,emissio,emissio_int,Quantum_Yiel_all_int,values_array,Start,Interval,number_of_fluor, p_num, m_num, particle_type,END,particle_uncert,medium_uncert,p_num_uncert,m_num_uncert,p_uncert,m_uncert)
        return prop, check,qy_all_sim, emission_all_sim,number_of_fluor,length_per_sim,layers_per_sim,uncert,counting_sims_uncert
    
    for i in range(line+1, len(infile)):  
        if infile[i][0:8] == "particle":
            if count > 0:
                if fluor ==1 :
                    emissio_intt = append(emissio_intt, ptype)
                ptype_in = get_index(ptype, 'p', p_num, m_num)
                check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
                mtype_in = get_index(mtype, 'm', p_num, m_num)
                
                optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv, particle[int(ptype_in), :, :], medium[int(mtype_in), :, :], thickness, dist,fluor,  start_wl, index,happens, again,number_of_fluor,particle_type[int(ptype_in)])  
                optics_sum += optics
                dist = zeros(1)
            count += 1                                  
            ptype = int(infile[i][8:])
        if infile[i][0:5] == "upper":    
            if len(infile[i][6:]) <2: 
                upper = int(infile[i][6:])
                upper_is_air = True
            else:
                upper = int(infile[i][12:])
                upper_type = get_index(upper, 'm', p_num, m_num)
        elif infile[i][0:5] == "lower":
            if len(infile[i][6:]) <2: 
                lower = int(infile[i][6:])
                lower_is_air = True
            else:
                lower = int(infile[i][12:])
                lower_type = get_index(lower, 'm', p_num, m_num)
        elif infile[i][0:6] == "matrix":
            mtype = int(infile[i][6:])
        elif infile[i][0:2] == "t:" or infile[i][0:2] == "t=":
            thickness = float(infile[i][2:])
            thickness = thickness / 10000
        elif infile[i][0:6] == "fluor:":    
            fluor = int(infile[i][6:])
            if fluor == 1 :
                happens = 1
                again += 1 
        elif infile[i][0:2] == "d:" or infile[i][0:2] == "d=":
            sizes = (infile[i][2:]).split(",")
            sizes = asarray(sizes, dtype=float)
            sizes = sizes/2
        elif infile[i][0:3] == "vf:":
            fv = (infile[i][3:]).split(",")
            fv = asarray(fv, dtype=float)
            fv = fv/100
            vol_frac_sum += sum(fv)
        elif infile[i][0:4] == "std:" or infile[i][0:4] == "std=":
            dist = (infile[i][4:]).split(",")
            dist = asarray(dist, dtype=float)
            if dist[0] >1:
                dist = dist/100
        elif infile[i][0:5] == "layer" and infile[i][5:] != "1":
            check = check_dist(current_sim, dist, check)
            ptype_in = get_index(ptype, 'p', p_num, m_num)
            check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
            layers += 1
            mtype_in = get_index(mtype, 'm', p_num, m_num)
            if fluor == 1:
                emissio_intt = append(emissio_intt, ptype)
            emissio_intt = append(emissio_intt, 0)    
            optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv, particle[int(ptype_in), :, :], medium[int(mtype_in), :, :], thickness, dist, fluor, start_wl, index,happens,again,number_of_fluor,particle_type[int(ptype_in)])
            optics_sum += optics 
            optics_sum[0, :] = optics[0, :]
            optics_sum[1, :] = optics[1, :]
            optics_sum[6, :] = optics[6, :]
            optics_sum[7, :] = optics[7, :]
            optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype_in), :, :],fluor,number_of_fluor,particle_type[int(ptype_in)],ptype_in,fluor_counting_all_Mie)
            optical_per_layer = dstack((optical_per_layer, optics_sum))
            vol_frac_sum = 0
            optics_sum = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0])))     ### "5"
            happens = 0
            again = 0
            count = 0
            dist = zeros(1)
        elif infile[i][0:4] == "sim:" and infile[i][4:] != "1":
            current_sim = int(infile[i][4:])  
            check = check_dist(current_sim, dist, check)
            ptype_in = get_index(ptype, 'p', p_num, m_num)
            check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
            mtype_in = get_index(mtype, 'm', p_num, m_num)
            optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv, particle[int(ptype_in), :, :], medium[int(mtype_in), :, :], thickness, dist,fluor, start_wl, index,happens,again,number_of_fluor,particle_type[int(ptype_in)])  ##### the colon ":" is the index which values are taken and inserted in thos function
            optics_sum += optics
            optics_sum[0, :] = optics[0, :]
            optics_sum[1, :] = optics[1, :]
            optics_sum[6, :] = optics[6, :]
            optics_sum[7, :] = optics[7, :]
            optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype_in), :, :],fluor,number_of_fluor,particle_type[int(ptype_in)],ptype_in,fluor_counting_all_Mie)
            optical_per_layer = dstack((optical_per_layer, optics_sum))
            prop = create_prop_array(upper_is_air,lower_is_air,prop, particle, upper, lower, upper_type, lower_type, layers, optical_per_layer, medium,number_of_fluor)
            vol_frac_sum = 0
            optics_sum = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0])))   
            optical_per_layer = zeros(( 8 + int(number_of_fluor), len(particle[0, :, 0]), 0))  
            length_per_sim = append(length_per_sim,(len(prop[:,0])))
            length_after = 0
            for k in range(len(emissio_int)):
                temp = loadtxt(emissio[k-1])
                if len(arange(temp[0,0], temp[-1,0] ) ) > length_after:
                    length_after = len(arange(temp[0,0], temp[-1,0] ) )
            emission_all = zeros((number_of_fluor*2*layers, length_after +1)) 
            qy_all = zeros((number_of_fluor*layers, len(index)))
            emissio_intt = append(emissio_intt, 0)
            i = -1
            for kk in range(len(emissio_intt)): 
                if emissio_intt[kk] ==0:
                    i += 1
                    emissio_intt_temp =np.sort(organi_lay_emis)
                    z = int(emissio_intt[kk])
                    for j in range(len(emissio_intt_temp)): 
                        for kkk in range(len(emissio_int)):
                            if emissio_intt_temp[j] == emissio_int[kkk]:
                                temp = loadtxt(emissio[kkk]) 
                                cdf, wave_flo=inv_cdf(temp)
                                e = len(cdf)
                                if max(wave_flo) > END or min(wave_flo) < Start:
                                    print('Please make sure in the input file, the "Start" and "End" simulated wavelengths are both lower and higher than the emission wavelengths (start and end wavelengths of an emission band of fluorescent particle), respectively')
                                    main_func()
                                emission_all[2*kkk + number_of_fluor*2*(i), 0:e] = cdf[:]
                                emission_all[2*kkk+ 1 + number_of_fluor*2*(i), 0:e] = wave_flo 
                                emission_all[2*kkk+ number_of_fluor*2*(i) ,-1] = e
                                temp_2 = loadtxt(Quantum_Yiel_all[kkk])
                                if (temp_2.size) >1 : 
                                    convert = False
                                    if max(temp_2[:,0]) < 30: 
                                        new_QY=interp(index[:], temp_2[:,0]*1000, temp_2[:,1])
                                    else:
                                        new_QY=interp(index[:], temp_2[:,0], temp_2[:,1])
                                    qy_all[kkk+ number_of_fluor*(i), int(values_array[2*(kkk)]):int(values_array[(2*(kkk))+1])] = new_QY
                                else:
                                    if int(values_array[2*kkk]) < Start or int(values_array[(2*kkk ) +1 ]) >  END :
                                        print('Please make sure the start and end excitation wavelengths are higher and lower than simulated start and end wavelegnths ')
                                        main_func()
                                    else: 
                                        constant_qy = temp_2
                                        qy_all[kkk + number_of_fluor*(i), int((int(values_array[2*kkk])-Start)/Interval):int((int(values_array[(2*kkk)+1])-Start)/Interval)] = constant_qy
                    organi_lay_emis = zeros(0, dtype=int)
                else:
                    organi_lay_emis =append(organi_lay_emis,int(emissio_intt[kk])) 
            try:
                emission_all_sim = vstack((emission_all_sim,emission_all))
            except:
                emission_all_sim =  emission_all.copy()
            try:
                qy_all_sim = vstack((qy_all_sim,qy_all))
            except:
                qy_all_sim  = qy_all.copy()
            dist = zeros(1)
            layers_per_sim = append(layers_per_sim,layers)
            layers = 1
            happens = 0
            again = 0
            count = 0
            emissio_intt =  zeros(0, dtype=int)
            upper_is_air = False
            lower_is_air = False
        
    ptype_in = get_index(ptype, 'p', p_num, m_num)
    mtype_in = get_index(mtype, 'm', p_num, m_num)
    check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
    check = check_diameters(current_sim, fv, sizes, check, particle_type[ptype_in], dist)
    if fluor ==1 :
        emissio_intt = append(emissio_intt, ptype)
    optics = mie_theory(ptype,fluor_counting_all_Mie,sizes, fv, particle[int(ptype_in), :, :], medium[int(mtype_in), :, :], thickness, dist, fluor, start_wl, index,happens,again,number_of_fluor,particle_type[int(ptype_in)])
    optics_sum += optics
    optics_sum[0, :] = optics[0, :]
    optics_sum[1, :] = optics[1, :]
    optics_sum[6, :] = optics[6, :]
    optics_sum[7, :] = optics[7, :]
    optics_sum = effective_medium(optics_sum, vol_frac_sum, medium[int(mtype_in), :, :],fluor,number_of_fluor,particle_type[int(ptype_in)],ptype_in,fluor_counting_all_Mie)     ### notice the particle and medium have all indicies of second and third dimension
    optical_per_layer = dstack((optical_per_layer, optics_sum))
    prop = create_prop_array(upper_is_air,lower_is_air,prop, particle, upper, lower, upper_type, lower_type, layers, optical_per_layer, medium,number_of_fluor)
    length_per_sim = append(length_per_sim,(len(prop[:,0])))
    layers_per_sim = append(layers_per_sim,layers)
    emissio_intt = append(emissio_intt, 0)
    length_after = 0
    for k in range(len(emissio_int)):
        temp = loadtxt(emissio[k-1])
        if len(arange(temp[0,0], temp[-1,0] ) ) > length_after:
            length_after = len(arange(temp[0,0], temp[-1,0] ) )
    emission_all = zeros((number_of_fluor*2*layers, length_after +1)) 
    qy_all = zeros((number_of_fluor*layers, len(index)))
    i = -1
    for kk in range(len(emissio_intt)): 
        if emissio_intt[kk] ==0:
            i += 1
            emissio_intt_temp =np.sort(organi_lay_emis)
            for j in range(len(emissio_intt_temp)):
                for kkk in range(len(emissio_int)):
                    if emissio_intt_temp[j] == emissio_int[kkk]:
                        temp = loadtxt(emissio[kkk]) 
                        cdf, wave_flo=inv_cdf(temp)
                        e = len(cdf)
                        if max(wave_flo) > END or min(wave_flo) < Start:
                            print('Please make sure in the input file, the "Start" and "End" simulated wavelengths are both lower and higher than the emission wavelengths (start and end wavelengths of an emission band of fluorescent particle), respectively')
                            main_func()
                        emission_all[2*kkk + number_of_fluor*2*(i), 0:e] = cdf[:]
                        emission_all[2*kkk+ 1 + number_of_fluor*2*(i), 0:e] = wave_flo
                        emission_all[2*kkk+ number_of_fluor*2*(i) ,-1] = e
                        temp_2 = loadtxt(Quantum_Yiel_all[kkk])
                        if (temp_2.size) >1 : 
                            convert = False
                            if max(temp_2[:,0]) < 30: 
                                new_QY=interp(index[:], temp_2[:,0]*1000, temp_2[:,1])
                            else:
                                new_QY=interp(index[:], temp_2[:,0], temp_2[:,1])
                            qy_all[kkk+ number_of_fluor*(i), int(values_array[2*(kkk)]):int(values_array[(2*(kkk))+1])] = new_QY
                        else:
                            if int(values_array[2*kkk]) < Start or int(values_array[(2*kkk ) +1 ]) >  END :
                                print('Please make sure the start and end excitation wavelengths are higher and lower than simulated start and end wavelegnths ')
                                main_func()
                            else:
                                constant_qy = temp_2
                                qy_all[kkk + number_of_fluor*(i), int((int(values_array[2*kkk])-Start)/Interval):int((int(values_array[(2*kkk)+1])-Start)/Interval)] = constant_qy
            organi_lay_emis = zeros(0, dtype=int)
        else:
            organi_lay_emis =append(organi_lay_emis,int(emissio_intt[kk]))
    try: 
        emission_all_sim = vstack((emission_all_sim,emission_all))
    except:
        emission_all_sim = emission_all.copy()
    try:
        qy_all_sim = vstack((qy_all_sim,qy_all))
    except:
        qy_all_sim = qy_all.copy()

    return prop, check,qy_all_sim, emission_all_sim,number_of_fluor,length_per_sim,layers_per_sim,uncert,counting_sims_uncert

def nanoparticle(infile, check):
    sims_uncert,Quantum_Yiel_all_uncert,Quantum_Yiel_all_int_uncert,fluor_counting_all_Mie,particle, medium, output_name, sims, photons, line, index, start_wl, solar,Quantum_Yiel_all,emissio,emissio_int,Quantum_Yiel_all_int ,values_array,Start,Interval,number_of_fluor,end, p_num, m_num, particle_type,particle_uncert,medium_uncert,uncert,p_num_uncert,m_num_uncert,p_uncert,m_uncert= import_header(infile,check)  
    prop, check ,qy_all, emission_all,number_of_fluor,length_per_sim,layers_per_sim,uncert,counting_sims_uncert= optical(sims_uncert,Quantum_Yiel_all_uncert,Quantum_Yiel_all_int_uncert,fluor_counting_all_Mie,line, infile, particle, medium, check, start_wl, index,Quantum_Yiel_all,emissio,emissio_int,Quantum_Yiel_all_int,values_array,Start,Interval,number_of_fluor, p_num, m_num, particle_type,end,particle_uncert,medium_uncert,uncert,p_num_uncert,m_num_uncert,p_uncert,m_uncert) 
    sims_per_medium = len(particle[0, :, 0])
    wavelengths = particle[0, :, 0]
    return prop, photons, output_name, sims, sims_per_medium, wavelengths, check, index,start_wl, solar,particle,qy_all, emission_all,number_of_fluor,Start,end,length_per_sim,Interval,layers_per_sim ,uncert,counting_sims_uncert,sims_uncert

def check_for_word_in_sim(infile, word, statement, check):
    start_line = 0
    for i in range(len(infile)):
        if infile[i][:4] == "sim:":
            start_line = i
            break
    word_present = 1
    sim_number = str(1)
    word_size = len(word)
    for i in range(start_line, len(infile)):
        if (infile[i][:4] == "sim:" and infile[i][4:] != '1'):
            if word_present == 1:
                check = True
                print("Sim" + sim_number + " must have " + statement)
            sim_number = infile[i][4:]
            word_present = 1
        if infile[i][:word_size] == word:
            word_present = 0
    if word_present == 1:
        check = True
        if statement != '':
            print("Sim" + sim_number + " must have " + statement)
    return check

def check_input_for_errors(infile):
    check = False
    for i in range(len(infile)):
        infile[i] = infile[i].replace(' ', '')
        infile[i] = infile[i].replace('\t', '')
        infile[i] = infile[i].lower()
    if infile[0] != "mc" and infile[0] != "nn":
        print("The MC  must be specified in the first line of the input file.")
        check = True
    output = particle = medium = photon = interval = 0
    if infile[0] == "nn":
        photon = 1
    for i in range(len(infile)):
        if infile[i][:6] == 'output' and infile[i][7:] != '':
            output = 1
        if infile[i][:8] == 'particle' and infile[i][10:] != '':
            particle = 1
        if infile[i][:6] == 'matrix' and infile[i][8:] != '':
            medium = 1
        if infile[0] == "mc":
            if infile[i][:7] == "photons" and infile[i][9:] != '':
                photon = 1
        if infile[i][:8] == "interval":
            interval = 1
        if infile[i][:3] == "sim":
            break
    if output == 0:
        check = True
        print("No output file name specified")
    if particle == 0:
        check = True
        print("No particle input specified")
    if medium == 0:
        check = True
        print("No matrix input specified")
    if photon == 0:
        check = True
        print("Number of photons must be specified")
    if interval == 0:
        check = True
        print("No interval specified")

    sim_number_should_be = 1
    sim_number_error = 1
    for i in range(len(infile)):
        if infile[i][:4] == "sim:":
            if infile[i][4:] != str(sim_number_should_be):
                sim_number_error = 0
            sim_number_should_be += 1
    if sim_number_error == 0:
        check = True
        print("Error in simulation numbers. Make sure they are labeled Sim 1, Sim 2, etc.")  
    # check for at least one layer
    check = check_for_word_in_sim(infile, 'layer', 'at least one defined layer', check)
    # check for at least one medium
    check = check_for_word_in_sim(infile, 'matrix', 'at least one matrix', check)
    # check for medium thickness
    check = check_for_word_in_sim(infile, 't:', 'a defined thickness', check)
    # check for at least one particle
    check = check_for_word_in_sim(infile, 'particle', 'at least one particle', check)
    # check for at least one volume fraction
    check = check_for_word_in_sim(infile, 'vf', 'at least volume fraction', check)
    if check:
        print("Please re-enter input file once corrected.")
        print("\n")
    return check, infile

def main_func():
    check = True
    while check: 
        infile = fname()    
        check, infile = check_input_for_errors(infile)
        if not check:                                                             
            prop, photons, output_name, sims, sims_per_medium,  wavelengths, check, index, start_wl, solar,particle,qy_all, emission_all,number_of_fluor,Start,end,length_per_sim,Interval,layers_per_sim,uncert,counting_sims_uncert,sims_uncert= nanoparticle(infile, check) 
    if (infile[0][0:2]).lower() == "mc":
        solar_file = loadtxt(solar)
        if uncert == True:
            sims = counting_sims_uncert 
        result1 ,result1_radiosity,get_num= main_mc(prop, photons, index, start_wl, qy_all, emission_all,number_of_fluor,sims,length_per_sim,solar_file,Interval,layers_per_sim) 
        check = False   
    else:
        print("Please specify either MC in the first line of the input file")
        print("Please re-enter input file once corrected")
        main_func()
        check = True
    if uncert == True:
        sims = counting_sims_uncert
    if solar != "":
        radiosity_avg = zeros(sims)
        absorbedd = zeros(sims)
        transmitt = zeros(sims)
        S = zeros(sims)
        q = zeros(sims)
        to = zeros(sims)
        flur_ref =  zeros(sims)
        absorb_fluo = zeros(sims) 
        non_fluor_abso = zeros(sims)

        solar_file = loadtxt(solar)
        if Start > solar_file[0, 0] or end < solar_file[-1, 0]:
            print("\n")
            print("WARNING: Wavelength range does not cover the full spectrum for the provided solar file.")
            print("The integration only covers the wavelength range provided and not the entire solar file range\n")
        for i in range(sims):
            radiosity_avg[i],absorbedd[i],transmitt[i],S[i],q[i],to[i],flur_ref[i],absorb_fluo[i],non_fluor_abso[i]= solar_spectrum(result1,wavelengths,solar_file,i,len(index))# solar_r[i] = solar_spectrum(output_sim)        solar, refl


    else:
        solar = ""
    if max(wavelengths) < 30: 
        wavelengths = wavelengths *1000
    length = int(len(result1[:, 0])/sims)
    prop_line = 0
    if uncert == False: 
        for i in range(sims):
            output_sim = str(output_name) + str(i+1) +".txt"
            with open(output_sim, 'w') as f:
                f.write('Sim ' + str(i+1) + '\n')
                f.write('\n')
                if solar != "":
                    f.write('Light Normalized Radiosity: ' + str(abs(round(radiosity_avg[i]*100, 5))) + '%\n')
                    f.write('Light Absorptance (heat generated by Stokes shift, Quantum Yield, non-fluorescent particles, and medium ): ' + str(abs(round(absorbedd[i]*100, 5))) + '%\n')
                    f.write('Light Transmittance: ' + str(abs(round(transmitt[i]*100, 5))) + '%\n')

                    f.write('Stokes Shift Absorptance (converted to heat): ' + str(abs(round(S[i]*100, 5))) + '%\n')
                    f.write('Quantum Yield Absorptance (converted to heat): ' + str(abs(round(q[i]*100, 5))) + '%\n')
                    f.write('Total Fluorescent Absorptance (absorptance of photons by fluorescent particles): ' + str(abs(round(to[i]*100, 5))) + '%\n')
                    f.write('Fluorescent Absorptance (converted to fluorescent emission): ' + str(abs(round(absorb_fluo[i]*100, 5))) + '%\n')
                    f.write('Reflectance by Fluorescent Emission: ' + str(abs(round(flur_ref[i]*100, 5))) + '%\n')
                    f.write('Non-Fluorescent Particles and Medium Absorptance (converted to heat): ' + str(abs(round(non_fluor_abso[i]*100, 5))) + '%\n')
                    
                    f.write('\n')
                f.write('Note: In the column headers below, the number in parentheses indicates the layer number for each attenuation coefficient (e.g., mu_s(1), mu_med(1), mu_a(1), mu_af(1) for layer 1)')
                f.write('\n')
                f.write('Note: The units are based on micrometer')
                f.write('\n')
                f.write('Wavelength\tNormalized Radiosity\tAbsorptance\tTransmittance\tStokes Shift Absorptance\tQuantum Yield Absorptance\tTotal Fluorescent Absorptance\tReflectance (via fluorescent emission)\tFluorescent Absorptance\tNon-Fluorescent Particles/Medium Absorptance')
                num_layers = int(get_num[i])
                for layer in range(num_layers):
                    #if number_of_fluor == 1:
                    f.write('\tn_medium(' + str(layer+1) + ')' + '\tk_medium(' + str(layer+1) + ')'  +'\tmu_a(non-fluorescent)(' + str(layer+1) + ')'+  '\tmu_a(medium)(' + str(layer+1) + ')'+ '\tmu_s(' + str(layer+1) + ')'+ '\tg(' + str(layer+1) + ')'+ '\tThickness(' + str(layer+1) + ')')# +  '\tmu_fluor1_(' + str(layer+1) + ')')
                    for flu in range(number_of_fluor):
                        f.write('\tmu_a(fluorescent)'+ str(flu+1) +'('+ str(layer+1) + ')')

                f.write('\n') 
                for j in range(length):
                    f.write(str(round(wavelengths[j], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 0] + result1_radiosity[j + length * i, 1], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 2], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 3], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 4], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 5], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 6], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 7], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 8], 4)) + '\t')
                    f.write(str(round(result1_radiosity[j + length * i, 9], 4)) + '\t')
                    while prop[prop_line, 0] != 0:
                        prop_line += 1
                    for layer in range(num_layers):
                        f.write(str(round(prop[prop_line - (1+num_layers) + layer, 0], 6)) + '\t')
                        f.write(str(round(prop[prop_line - (1+num_layers) + layer, 1], 6)) + '\t')
                        f.write(str(round(prop[prop_line - (1+num_layers) + layer, 2]/10000, 6)) + '\t')
                        f.write(str(round(prop[prop_line - (1+num_layers) + layer, 3]/10000, 6)) + '\t')
                        f.write(str(round(prop[prop_line - (1+num_layers) + layer, 4]/10000, 6)) + '\t')
                        f.write(str(round(prop[prop_line - (1+num_layers) + layer, 5], 8)) + '\t')
                        f.write(str(round(prop[prop_line - (1+num_layers) + layer, 6]*10000, 6)) + '\t')
                        for flu in range(number_of_fluor):
                            f.write(str(round(prop[prop_line - (1+num_layers) + layer, 8 + flu]/10000, 6)) + '\t')
            
                    prop_line += 1
                    f.write('\n')
                f.write('\n')
                f.write('Input file:\n')
                for j in range(len(infile)):
                    if infile[j][:4] == "sim:":
                        break
                    f.write(infile[j] + '\n')
                for j in range(len(infile)):
                    if infile[j][:] == "sim:"+str(i+1):
                        f.write('\n')
                        f.write("sim:"+str(i+1) + '\n')
                        for k in range(j+1, len(infile)):
                            if infile[k][:4] == "sim:":
                                break
                            f.write(infile[k] + '\n')
                        break
            f.close()

    plt.style.use('default')
    figsize = (10,8)
    dpi = 600
    axis_font_size = 15 
    label_font_size = 17 
    linewidth = 2
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['font.size'] = axis_font_size
    plt.rcParams['font.sans-serif'] = ['Arial']
    sim_conunt = 0

    for k in range(sims):
        R_results = result1_radiosity[k*length:(k+1)*length, 0]+result1_radiosity[k*length:(k+1)*length, 1]
        A_results = result1_radiosity[k*length:(k+1)*length, 2]
        T_results = result1_radiosity[k*length:(k+1)*length, 3]

    if uncert == False:
        for j in range(sims):
            sim_conunt += 1
            # get the result of this run
            R_results = result1_radiosity[j*length:(j+1)*length, 0]+result1_radiosity[j*length:(j+1)*length, 1]
            A_results = result1_radiosity[j*length:(j+1)*length, 2]
            T_results = result1_radiosity[j*length:(j+1)*length, 3] 
            R_results_flux = result1[j*length:(j+1)*length, 1]+result1[j*length:(j+1)*length, 3]
            A_results_flux = result1[j*length:(j+1)*length, 0]
            T_results_flux = result1[j*length:(j+1)*length, 2] 

            SS = result1_radiosity[j*length:(j+1)*length, 4]
            QY = result1_radiosity[j*length:(j+1)*length, 5]
            Tot = result1_radiosity[j*length:(j+1)*length, 6]
            Fluor = result1_radiosity[j*length:(j+1)*length, 7]
            Absorb_to_fluor = result1_radiosity[j*length:(j+1)*length, 8]
            non_fluor_absorb = result1_radiosity[j*length:(j+1)*length, 9]

            # plot the R, A, T
            plt.figure(figsize=figsize,dpi=dpi)
            max_y = 0
            if max(R_results)>max_y:
                max_y = max(R_results)
            if max(A_results)>max_y:
                max_y = max(A_results)
            if max(T_results)>max_y:
                max_y = max(T_results)
            if max_y < 1:
                max_y = 1
            plt.plot(wavelengths, R_results, label="Normalized Radiosity (J*)",color = 'black',linewidth = linewidth)
            plt.plot(wavelengths, A_results, label="Absorptance (Converted to Heat) (α*)",color = 'red',linewidth = linewidth)
            plt.plot(wavelengths, T_results, label="Transmittance (τ*)",color = 'blue',linewidth = linewidth)       
            if solar != "":
                solar_arr = loadtxt(solar)
                x_axis = solar_arr[:,0]
                if max(solar_arr[:,0])<30:
                    x_axis = solar_arr[:,0]*1000
                plt.fill_between(x_axis,
                        (solar_arr[:,1]/max(solar_arr[:,1])),
                        color='pink',
                        alpha=0.5)
            plt.xlim(Start, end)
            plt.ylim(0,max_y+0.2*max_y)
            plt.xlabel('Wavelength (nm)',fontsize = label_font_size)# (\u03BCm)
            plt.ylabel('Spectral Responce',fontsize = label_font_size)
            plt.legend(loc='best',frameon=False,fontsize = axis_font_size)
            plt.savefig(str(output_name)+"_Spectral_Response_Plot{}.png".format(j+1))
            plt.close()

            plt.figure(figsize=figsize,dpi=dpi)
            max_y = 0
            if max(R_results_flux)>max_y:
                max_y = max(R_results_flux)
            if max(A_results_flux)>max_y:
                max_y = max(A_results_flux)
            if max(T_results_flux)>max_y:
                max_y = max(T_results_flux) 
            if max_y < 1.5:
                max_y = 1.5  
            plt.plot(wavelengths, R_results_flux, label="Radiosity Flux (J)",color = 'black',linewidth = linewidth)
            plt.plot(wavelengths, A_results_flux, label="Absorbed Flux (converted to heat) (α)",color = 'red',linewidth = linewidth)
            plt.plot(wavelengths, T_results_flux, label="Transmitted Flux (τ)",color = 'blue',linewidth = linewidth)
            if solar != "":
                solar_arr = loadtxt(solar)
                x_axis = solar_arr[:,0]
                if max(solar_arr[:,0])<30:
                    x_axis = solar_arr[:,0]*1000
                plt.fill_between(x_axis,
                        (solar_arr[:,1]/max(solar_arr[:,1])),
                        color='pink',
                        alpha=0.5)
            plt.xlim(Start, end)
            plt.ylim(0,max_y+0.2*max_y)
            plt.xlabel('Wavelength (nm)',fontsize = label_font_size)# (\u03BCm)
            plt.ylabel('Spectral Flux (W/m²·nm)',fontsize = label_font_size)
            plt.legend(loc='best',frameon=False,fontsize = axis_font_size)
            plt.savefig(str(output_name)+"_Spectral_Flux_Plot{}.png".format(j+1))
            plt.close()

            plt.figure(figsize=figsize,dpi=dpi)
            max_y = 0
            if max(Tot)>max_y:
                max_y = max(Tot)
            if max(Absorb_to_fluor)>max_y:
                max_y = max(Absorb_to_fluor)
            if max(SS)>max_y:
                max_y = max(SS)
            if max(QY)>max_y:
                max_y = max(QY)
            if max(Fluor)>max_y:
                max_y = max(Fluor)
            if max(non_fluor_absorb)>max_y:
                max_y = max(non_fluor_absorb)
            if max_y < 1:
                max_y = 1
            plt.plot(wavelengths, Tot, label="Total Fluorescent Absorptance",color = '#DC267F',linewidth = linewidth,linestyle='dotted')
            plt.plot(wavelengths, Absorb_to_fluor, label="Fluorescent Absorptance (radiative relaxation)",color = '#FFB000',linewidth = linewidth)
            plt.plot(wavelengths, SS, label="Stokes Shift Absorptance (non-radiative relaxation)",color = '#648FFF',linewidth = linewidth)
            plt.plot(wavelengths, QY, label="Quantum Yield Absorptance (non-radiative relaxation)",color = '#785EF0',linewidth = linewidth)
            plt.plot(wavelengths, Fluor, label="Reflectance by Fluorescent Emission",color = '#FE6100',linewidth = linewidth)
            plt.plot(wavelengths, non_fluor_absorb, label="Non-Fluorescent Particles and Medium Absorptance",color = '#000000',linewidth = linewidth)

            if solar != "":
                solar_arr = loadtxt(solar)
                x_axis = solar_arr[:,0]
                if max(solar_arr[:,0])<30:
                    x_axis = solar_arr[:,0]*1000
                plt.fill_between(x_axis,
                        (solar_arr[:,1]/max(solar_arr[:,1])),
                        color='pink',
                        alpha=0.5)
            plt.xlim(Start, end)
            plt.ylim(0,max_y+0.7*max_y)
            plt.xlabel('Wavelength (nm)',fontsize = label_font_size)# (\u03BCm)
            plt.ylabel('Spectral Responce',fontsize = label_font_size)
            plt.legend(loc='upper right',frameon=False,fontsize = axis_font_size)
            plt.savefig(str(output_name)+"_Detailed_Spectal_Absorption_and_Reflection_Mechanisms_plot{}.png".format(j+1))
            plt.close()

            print(f"The spectral radiative properties plot is saved as: {str(output_name)}_plot{j+1}.png")
            display_color_from_reflectance(R_results,wavelengths,output_name=str(output_name)+'_color', i = sim_conunt)
    else:
        for j in range(sims):
            sim_conunt += 1
            # get the result of this run
            R_results = result1_radiosity[j*length:(j+1)*length, 0]+result1_radiosity[j*length:(j+1)*length, 1]
            A_results = result1_radiosity[j*length:(j+1)*length, 2]
            T_results = result1_radiosity[j*length:(j+1)*length, 3]

            if j == 0:
                R_results_cum = R_results.copy()
                A_results_cum = A_results.copy()
                T_results_cum = T_results.copy()

                R_results_cum_square = R_results**2
                A_results_cum_square = A_results**2
                T_results_cum_square = T_results **2               
            else:
                R_results_cum += R_results
                A_results_cum += A_results
                T_results_cum += T_results 

                R_results_cum_square += R_results**2
                A_results_cum_square += A_results**2
                T_results_cum_square += T_results **2 
   
        R_results_mean = R_results_cum / sims  # Your current mean
        A_results_mean = A_results_cum / sims  # Your current mean
        T_results_mean = T_results_cum / sims  # Your current mean

        R_variance = (R_results_cum_square  - sims*R_results_mean**2)/(sims)
        R_std = np.sqrt(R_variance)

        A_variance = (A_results_cum_square - sims*A_results_mean**2)/(sims)
        A_std = np.sqrt(A_variance)

        T_variance = (T_results_cum_square  - sims*T_results_mean**2)/(sims)
        T_std = np.sqrt(T_variance)

        std_error_R = R_std / np.sqrt(sims)
        upper_bound_R = R_results_mean + 1.96 * std_error_R
        lower_bound_R = R_results_mean - 1.96 * std_error_R

        std_error_A = A_std / np.sqrt(sims)
        upper_bound_A = A_results_mean + 1.96 * std_error_A
        lower_bound_A = A_results_mean - 1.96 * std_error_A

        std_error_T = T_std / np.sqrt(sims)
        upper_bound_T = T_results_mean + 1.96 * std_error_T
        lower_bound_T = T_results_mean - 1.96 * std_error_T
        plt.figure(figsize=figsize,dpi=dpi)
        max_y = 0
        if max(R_results_mean)>max_y:
            max_y = max(R_results_mean)
        if max(A_results_mean)>max_y:
            max_y = max(A_results_mean)
        if max(T_results_mean)>max_y:
            max_y = max(T_results_mean)
        if max_y < 1:
            max_y = 1
        if solar != "":
            solar_arr = loadtxt(solar)
            x_axis = solar_arr[:,0]
            if max(solar_arr[:,0])<30:
                x_axis = solar_arr[:,0]*1000
            plt.fill_between(                
            x_axis,
            solar_arr[:,1] / max(solar_arr[:,1]),
            color='pink',
            alpha=0.35,
            label="Normalized Solar Spectrum"
            )
        plt.plot(wavelengths, R_results_mean, label="Normalized Radiosity (J*)",color = 'black',linewidth = linewidth)
        plt.plot(wavelengths, A_results_mean, label="Absorptance (Converted to Heat) (α*)",color = 'red',linewidth = linewidth)
        plt.plot(wavelengths, T_results_mean, label="Transmittance (τ*)",color = 'blue',linewidth = linewidth)

        plt.fill_between(wavelengths, lower_bound_R, upper_bound_R, color='black', alpha=0.25)
        plt.fill_between(wavelengths, lower_bound_A, upper_bound_A,  color='red', alpha=0.25)
        plt.fill_between(wavelengths, lower_bound_T, upper_bound_T,  color='blue', alpha=0.25)
        plt.xlim(Start, end)
        plt.ylim(0,max_y+0.5*max_y)
        plt.xlabel('Wavelength (nm)',fontsize = label_font_size+3)# \u03BC
        plt.ylabel('Spectral Responce',fontsize = label_font_size+3)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.legend(loc='upper right',frameon=False,fontsize = 21)
        plt.savefig(str(output_name)+"_spectrum_plot{}.png".format(0+1))
        plt.close()
        print(f"The spectral radiative properties plot is saved as: {str(output_name)}_plot{0+1}.png")
        #display_color_from_reflectance(R_results,output_name=str(output_name)+'_color', i = 1)
        display_color_from_reflectance(R_results,wavelengths,output_name=str(output_name)+'_color', i = sim_conunt)
    print("Results saved! The spectral radiative properties plots are saved via _plot and color is saved via _color ")
    print("Thank you for using our program. If you have any questions or issues, feel free to ask me on GitHub!")
    return

if __name__ == "__main__":
    print('\033[1m{: ^75s}\033[0m'.format("Fluor-FOS"))
    print('{: ^75s}'.format("Fluor-FOS for optical properties modeling for fluorescent nanoparticle media"))
    print('{: ^75s}'.format("Version: 0.1.0\n"))
    print('{: ^75s}'.format("Khalid Alhammadi, Daniel Carne, Xiulin Ruan")) 
    print('{: ^75s}'.format("School of Mechanical Engineering, Purdue University"))
    print('{: ^75s}'.format("West Lafayette, IN 47909, USA\n"))
    main_func()
