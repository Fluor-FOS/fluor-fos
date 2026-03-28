import numpy as np
from numpy import zeros, random, pi, cos, log, vstack, conj, real,interp, sum
from numba import prange,njit

@njit() 
def r_specular(layer):
    n_amb = layer[0, 0]
    n_medium = layer[1, 0]
    return (n_amb-n_medium)**2/((n_amb+n_medium)**2)

@njit() 
def check_bounds(layer_depths, z, uz, step, current_layer, ua, us, sum_all_non_fluo, beta_nonf):
    hit = False
    step_in_z = abs(step * uz)
    step_remaining = 0
    if uz < 0:
        distance_to_bound = (z-layer_depths[current_layer])
        if step_in_z > distance_to_bound:
            hit = True
            step_remaining = (step - (-distance_to_bound/uz)) *(ua + us + sum_all_non_fluo +  beta_nonf)
            step = -distance_to_bound/uz
    else:
        distance_to_bound = (layer_depths[current_layer+1]-z)   
        if step_in_z > distance_to_bound:
            hit = True
            step_remaining = (step - (distance_to_bound / uz))*(ua + us + sum_all_non_fluo + beta_nonf)
            step = distance_to_bound / uz
    return hit, step, step_remaining

@njit() 
def fresnel_reflectance(z,n_medium, n_outer, k_outer,k_medium, uz, outer_layer, num_layers):
    if outer_layer == 0 or outer_layer == num_layers - 1:
        if n_medium == n_outer and k_outer == 0:
            return 0, uz
        elif uz > (1.0 - 1.0e-6):
            probability = ((n_outer-n_medium)**2 + (k_outer-k_medium)**2)/((n_outer+n_medium)**2 + (k_outer+k_medium)**2)
            return probability, uz
        elif uz < (1.0e-6):
            return 1, 0
        else:
            temp_1 = (1.0 - uz * uz) ** 0.5
            temp_1_5 = n_medium / (n_outer - 1j * k_outer)
            temp_2 = n_medium * temp_1 / (n_outer - 1j * k_outer)
            temp_3 = (1 - temp_2 * temp_2) ** 0.5
            E_par = (uz / temp_3 - temp_1_5) / (uz / temp_3 + temp_1_5)
            R_par = E_par * conj(E_par)
            E_per = -(temp_3 / uz - temp_1_5) / (temp_3 / uz + temp_1_5)
            R_per = E_per * conj(E_per)
            probability = real(R_par + R_per) * 0.5

            return probability, 1
    else:
        if n_medium == n_outer:
            return 0, uz
        elif uz > (1.0 - 1.0e-6):
            probability = ((n_outer-n_medium)/(n_outer+n_medium))**2
            return probability, uz
        elif uz < (1.0e-6):
            return 1, 0
        
        elif (uz) < ((1.0 - n_outer * n_outer / (n_medium * n_medium)) ** 0.5):
            return 1, uz

        else:
            temp_1 = (1.0 - uz * uz) ** 0.5
            temp_2 = n_medium * temp_1 / n_outer
            if temp_2 >= 1:
                return 1, 0
            temp_3 = (1-temp_2*temp_2)**0.5
            temp_4 = uz*temp_3 - temp_1*temp_2
            temp_5 = uz*temp_3 + temp_1*temp_2
            temp_6 = temp_1*temp_3 + uz*temp_2
            temp_7 = temp_1*temp_3 - uz*temp_2
            probability = 0.5*temp_7*temp_7*(temp_4*temp_4 + temp_5*temp_5) / (temp_6*temp_6*temp_5*temp_5)
            return probability, temp_3

@njit() 
def hit_bound(current_layer, layer, uz, z, crit_cos, r, t, active, w): 
    if uz < 0:                                                  
            probability,uz_new= fresnel_reflectance(z,layer[current_layer+1, 0], layer[current_layer, 0],layer[current_layer, 1],layer[current_layer+1, 1], -uz,current_layer,len(layer[:,0]))
            if random.random_sample() > probability:
                uz = -uz_new
                if current_layer == 0:      
                    active = False
                    r = (w)          
                else:
                    current_layer -= 1 
            else:
                uz = -uz
    else:
            probability, uz_new = fresnel_reflectance(z,layer[current_layer+1, 0], layer[current_layer+2, 0],layer[current_layer+2, 1],layer[current_layer+1, 1], uz,current_layer+2,len(layer[:,0]))
            if random.random_sample() > probability:
                uz = uz_new
                if current_layer == (len(layer[:, 0])-3): 
                    active = False
                    t = (w)
                else:
                    current_layer += 1
            else:
                uz = -uz
    return current_layer, uz, r, t, active 

@njit() 
def hit_boundf(current_layer,layer, uz, z, crit_cos_f, active, w, kkk, prop,f2_k,f2_in,f3_k,f3_in,correcting_index): 
    if uz < 0:
            probability,uz_new=fresnel_reflectance(z,prop[correcting_index, 0], prop[correcting_index - 1, 0],prop[correcting_index - 1, 1],prop[correcting_index, 1] ,-uz,current_layer,len(layer[:,0]))
            if random.random_sample() > probability:
                uz = -uz_new
                if current_layer == 0:       
                    active = False
                    f2_in = w
                    f2_k = kkk
                else:
                    current_layer -= 1
                    
            else:
                uz = -uz

    else:
            probability, uz_new = fresnel_reflectance(z,prop[correcting_index, 0], prop[correcting_index + 1, 0],prop[ correcting_index +1, 1],prop[correcting_index, 1] ,uz,current_layer+2,len(layer[:,0]))
            if random.random_sample() > probability:
                uz = uz_new
                if current_layer == (len(layer[:, 0])-3): 
                    active = False
                    f3_in = w
                    f3_k = kkk
                else:
                    current_layer += 1
            else:
                uz = -uz

    return current_layer, uz,f2_k,f2_in,f3_k,f3_in,active, kkk

@njit() 
def new_angle(uz, g):
    if g == 0:
        cos_theta = 2*random.random_sample()-1
    elif g == 1:
       cos_theta=1 
    else:
        temp = (1 - g * g) / (1 - g + 2 * g * random.random_sample())
        cos_theta = (1 + g * g - temp * temp) / (2 * g)
    
    if uz ==1:
        uz = cos_theta
    elif uz == -1 :
        uz = - cos_theta
    else:
        sin_theta = (1.0 - cos_theta * cos_theta) ** 0.5

        psi = 2.0 * pi * random.random_sample()  
        cos_psi = cos(psi)

        uz = -sin_theta*cos_psi*((1-uz*uz)**0.5)+uz*cos_theta

    return uz

@njit() 
def new_angle_f(uz):
    cos_t= 1 - ( 2 * random.random_sample() )
    uz = cos_t
    return uz

@njit() 
def setup_mcf(prop, crit_cos_f, num_layers, stepp,difff,fi,sims,length_per_sim,k):
    for i in range(num_layers):
        n_up = prop[(int(fi *(difff / stepp))+int(difff) + int(i) )+ k*(length_per_sim[sims-1]), 0]   ##  layer[i, 0]
        n_medium = prop[(int(fi *(difff / stepp))+int(difff) + int(i) + 1)+ k*(length_per_sim[sims-1]), 0]    ##  layer[i+1, 0]
        if n_medium > n_up:
            crit_cos_f[0, i] = (1.0 - n_up * n_up / (n_medium * n_medium)) ** 0.5
        else:
            crit_cos_f[0, i] = 0
        n_down = prop[(int(fi *(difff / stepp))+int(difff) + int(i) + 2)+ k*(length_per_sim[sims-1]), 0]
        if n_medium > n_down:
            crit_cos_f[1, i] = (1.0 - n_down * n_down / (n_medium * n_medium)) ** 0.5
        else:
            crit_cos_f[1, i] = 0
    return  crit_cos_f

@njit() 
def setup_mc(layer, layer_depths, crit_cos, num_layers):
    z = 0
    for i in range(num_layers):
        thickness = layer[i+1,6]
        z += thickness
        layer_depths[i+1] = z
        n_up = layer[i, 0]    
        n_medium = layer[i+1, 0]   
        if n_medium > n_up:
            crit_cos[0, i] = (1.0 - n_up * n_up / (n_medium * n_medium)) ** 0.5
        else:
            crit_cos[0, i] = 0
        n_down = layer[i+2, 0]
        if n_medium > n_down:
            crit_cos[1, i] = (1.0 - n_down * n_down / (n_medium * n_medium)) ** 0.5
        else:
            crit_cos[1, i] = 0
    return layer_depths, crit_cos

@njit()  
def initialize_photon(start_wl,layer, rsp, prop, ind, qy_all, emission_all,number_of_fluor,sims,length_per_sim,Interval,layers_per_sim):#,f1,f2_k,f2_in,f3_k,f3_in): ##nn
    r = 0
    a = 0
    t = 0
    f1 =0 
    f1_k = 0
    f2_in =0
    f2_k =0 
    f3_in = 0
    f3_k =0
    fss_in=0
    fss_k=0
    fqy_in=0
    fqy_k=0
    ftot_in=0
    f_fluor_in=0
    current_layer = 0
    z = 0
    w = 1 -rsp      
    active = True
    re_emit = False
    kkk=ind
    num_layers = len(layer[:, 0]) - 2  
    layer_depths = zeros(num_layers+1)
    uz = 1
    diff=0
    stepp= Interval
    fi = num_layers + 3 - stepp 

    sa = 0 
    sa_2 = 0
    crit_cos = zeros((2, num_layers))
    layer_depths, crit_cos = setup_mc(layer, layer_depths, crit_cos, num_layers)
    step_remaining = 0
    crit_cos_f = zeros((2, num_layers))
    correcting_index = 0
    correcting_emis_index =0 
    correcting_qy_fluor = 0
 
    k = 1
    if sims == 0: 
        k = 0
    first_remis = 0
    sum_all_non_fluo = 0
    while active is True:
        skip = 0
        skip2 = 0
        diff=(kkk*stepp) 
        if re_emit is False:
            ua = layer[current_layer + 1, 2]              
            beta_nonf = layer[current_layer + 1, 3]               
            us = layer[current_layer + 1, 4]               
            g = layer[current_layer + 1, 5]
            switch = layer[current_layer + 1, 7]               
            sum_all_non_fluo = sum(layer[current_layer + 1, 8:])
        else:
            first_remis += 1
            correcting_index = ( int(fi *(diff / stepp)) + int(diff) + int( current_layer)) + 1 + k*(length_per_sim[sims-1]) 
            ua = prop[correcting_index, 2]  
            beta_nonf = prop[correcting_index, 3] 
            us = prop[correcting_index, 4]
            g = prop[correcting_index, 5]
            switch = prop[correcting_index, 7]
            sum_all_non_fluo = sum(prop[correcting_index , 8:])

        if re_emit is True:
            crit_cos_f = setup_mcf(prop, crit_cos_f, num_layers, stepp,diff,fi,sims,length_per_sim,k)
        if step_remaining == 0:
            rnd = random.random_sample()
            step = -log(rnd) / (ua + us + sum_all_non_fluo + beta_nonf)
        else:
            step = step_remaining / (ua + us +  sum_all_non_fluo  + beta_nonf) 
        hit, step, step_remaining = check_bounds(layer_depths, z, uz, step, current_layer, ua, us, sum_all_non_fluo, beta_nonf)
        z += step*uz
        if (hit is True) and (re_emit is False):
     
            current_layer, uz, r, t, active= hit_bound(current_layer, layer, uz, z, crit_cos, r, t, active, w)                                  
        elif (hit is True) and (re_emit is True):
         
            current_layer, uz,f2_k,f2_in,f3_k,f3_in, active, kkk= hit_boundf(current_layer,layer, uz, z, crit_cos_f, active, w, kkk, prop,f2_k,f2_in,f3_k,f3_in,correcting_index) 
        else:
            if (switch == 0) and (re_emit is False):    
                change_in_w = w*(ua + sum_all_non_fluo + beta_nonf)/(ua+us+  sum_all_non_fluo  + beta_nonf)
                w -= change_in_w                     
                a += change_in_w 
                uz = new_angle(uz, g)      
            elif (switch == 0) and (re_emit is True):              
                change_in_w = w * ((ua + sum_all_non_fluo + beta_nonf)/((ua+us+ sum_all_non_fluo +  beta_nonf)))
                w -= change_in_w                
                f1 +=  change_in_w
                uz = new_angle(uz, g)
            elif (switch == 1)  and ((re_emit is False) or (re_emit is True )):
                if   random.random_sample()    < (us/( ua + us + beta_nonf +  sum_all_non_fluo   ))   and ( re_emit is False ):       
                    uz = new_angle(uz, g)   
                elif  random.random_sample()   < (us/( ua + us + beta_nonf +  sum_all_non_fluo  ))   and ( re_emit is True ) :  
                    uz = new_angle(uz, g)      
                else: 
                    correcting_index = ( int(fi *(diff / stepp)) + int(diff) + int( current_layer)) + 1 + k*(length_per_sim[sims-1]) 
                    correcting_qy_fluor = number_of_fluor*k*sum(layers_per_sim[:sims]) 
                    sum_cumulitive1 = 0
                    for i in range(len(prop[correcting_index , 8:] )): 
                        prob_f = qy_all[i + correcting_qy_fluor,kkk]* ( prop[correcting_index , 8+i] )
                        sa = (qy_all[i + correcting_qy_fluor,kkk]* prop[correcting_index , 8+i] + sum_cumulitive1)/( sum_all_non_fluo + ua + beta_nonf)
                        if (sa) > random.random_sample():
                            skip = 1
                            correcting_emis_index = number_of_fluor*2*(int(current_layer)) + number_of_fluor*2*k*sum(layers_per_sim[:sims]) 
                            re_emit = True
                            wavelength_old = (kkk*stepp)+start_wl 

                            limit =int(emission_all[0 + correcting_emis_index,-1] )
                            cdf = emission_all[0 + correcting_emis_index, :limit]
                            wave_flo = emission_all[1 + correcting_emis_index, :limit] 
       
                            x=int(interp(random.random_sample() , cdf, wave_flo))
                            while x % stepp != 0:
                                x = x + 1
                            diff = x - start_wl
                            wavelength_new =int(x)

                            f1 += ((w * (1-((wavelength_old)/(wavelength_new))))) 
                            f1_k = kkk

                            fss_in +=  ((w * (1-((wavelength_old)/(wavelength_new))))) 
                            fss_k = kkk
                            if first_remis == 0:
                                ftot_in = w

                            correcting_index = (int(fi *(diff / stepp)) + int(diff) + int(current_layer)) + 1+ k*(length_per_sim[sims-1]) 
                            beta_nonf = prop[correcting_index, 3] 
                            us = prop[correcting_index, 4]
                            g = prop[correcting_index, 5]
                            ua = prop[correcting_index, 2]  
                            sum_all_non_fluo = sum(prop[correcting_index , 8:])
                            switch = prop[correcting_index, 7]
                            
                            w  *= (wavelength_old)/(wavelength_new)

                            if first_remis == 0:
                                f_fluor_in = w

                            kkk = int(diff/stepp)
                        
                            uz = new_angle_f(uz)
                            break
                        else:
                            sum_cumulitive1 += prob_f
                    if skip ==0:
                        sum_absob =0
                        sum_cumulitive2 =0
                        for i in range(len(prop[correcting_index , 8:] )): 
                            temp = (1-qy_all[i + correcting_qy_fluor,kkk])* prop[correcting_index , 8+i] 
                            sum_absob += temp
                            
                        for i in range(len(prop[correcting_index , 8:] )): 
                            outside_range = (1-qy_all[i + correcting_qy_fluor,kkk])
                            prob_f_h = (1-qy_all[i + correcting_qy_fluor,kkk])* ( prop[correcting_index , 8+i] )
                            sa_2 = ( ((1-qy_all[i + correcting_qy_fluor,kkk])* prop[correcting_index , 8+i] ) + sum_cumulitive2)/( sum_absob + ua + beta_nonf) 
                            if (sa_2 ) > random.random_sample() and outside_range != 1:
                               
                                skip2 = 1
                                f1 += w
                                f1_k = kkk
                                fqy_in += w
                                fqy_k = kkk
                                active = False  
                                break 
                            else:
                                sum_cumulitive2 += prob_f_h
                        if skip2 == 0:
                            a += w
                            active = False                                             
        if (w < 0.0001) and (active is True): 
            if random.random_sample() < 0.1:
                w = w/0.1
            else:
                active = False
    return r, a, t, f1,f2_k,f2_in,f3_k,f3_in,f1_k,fss_in,fss_k,fqy_in,ftot_in,f_fluor_in,fqy_k

@njit(parallel=True) 
def run_mc(layer, n, rsp, prop, ind,r, a, t,f1,f2_k,f2_in,f3_k,f3_in,f1_k,fss_in,fss_k,fqy_in,ftot_in,f_fluor_in,fqy_k,qy_all, emission_all,number_of_fluor,sims,length_per_sim,Interval,layers_per_sim,start_wl):         
    for i in prange(n): 
        r[i], a[i], t[i],f1[i],f2_k[i],f2_in[i],f3_k[i],f3_in[i],f1_k[i],fss_in[i],fss_k[i],fqy_in[i],ftot_in[i],f_fluor_in[i],fqy_k[i]= initialize_photon(start_wl,layer, rsp, prop, ind,qy_all, emission_all,number_of_fluor,sims,length_per_sim,Interval,layers_per_sim)
    return r, a, t,f1,f2_k,f2_in,f3_k,f3_in,f1_k,fss_in,fss_k,fqy_in,ftot_in,f_fluor_in,fqy_k

@njit()                               
def monte_carlo(layer, n, start_wl, prop,ind,index,qy_all, emission_all,number_of_fluor,sims,length_per_sim,Interval,layers_per_sim):                
    r = zeros(n)
    a = zeros(n)
    t = zeros(n)
    f1=zeros(n)
    f1_k=zeros(n)
    f2_k=zeros(n)
    f2_in =zeros(n)
    f3_k=zeros(n)
    f3_in =zeros(n)

    fss_in =zeros(n)
    fss_k =zeros(n)

    fqy_in =zeros(n)
    fqy_k =zeros(n)  

    ftot_in =zeros(n)
    f_fluor_in =zeros(n)     

    rsp = r_specular(layer)
    fff1=zeros(((len(index), len(index))), dtype=np.float32) 
    fff2=zeros(((len(index), len(index))), dtype=np.float32)  
    fff3=zeros(((len(index), len(index))), dtype=np.float32)  

    fss=zeros(((len(index), len(index))), dtype=np.float32)  
    fqy=zeros(((len(index), len(index))), dtype=np.float32)  
    f_tot=zeros(((len(index), len(index))), dtype=np.float32)  

    f_fluor=zeros(((len(index), len(index))), dtype=np.float32)  

    r, a, t ,f1,f2_k,f2_in,f3_k,f3_in,f1_k,fss_in,fss_k,fqy_in,ftot_in,f_fluor_in,fqy_k= run_mc( layer, n, rsp, prop, ind,r, a, t,f1,f2_k,f2_in,f3_k,f3_in,f1_k,fss_in,fss_k,fqy_in,ftot_in,f_fluor_in,fqy_k,qy_all, emission_all,number_of_fluor,sims,length_per_sim,Interval,layers_per_sim,start_wl) 
    r_tot = 0
    a_tot = 0
    t_tot = 0
    for i in range(n):
        r_tot += r[i]  
        a_tot += a[i]  
        t_tot += t[i]  
        fff1[int(f1_k[i]),ind] += f1[i]
        fff2[int(f2_k[i]),ind] += f2_in[i]
        fff3[int(f3_k[i]),ind] += f3_in[i]

        fss[int(fss_k[i]),ind] += fss_in[i]
        fqy[int(fqy_k[i]),ind] += fqy_in[i]
        f_tot[ind,ind] += ftot_in[i]
        f_fluor[ind,ind] += f_fluor_in[i]
    return rsp, r_tot, a_tot, t_tot, fff1,fff2,fff3 , ind ,fss,fqy,f_tot,f_fluor

@njit() 
def count_sims(prop):  
    count = 0
    for i in range(len(prop[:, 0])):          
        if prop[i, 0] == 0:
            count += 1
    return count

def main_mc(prop, photons, index, start_wl, qy_all, emission_all,number_of_fluor,sims,length_per_sim,solar_file,Interval,layers_per_sim):
    layer = zeros((0, 8 + number_of_fluor))        
    results = zeros((len(index),10))  
    result1_big = zeros((0, 10)) 
    radiosity_big = zeros((0, 10)) 
    count = 0
    total_sims = count_sims(prop)
    
    ff1=zeros(((len(index), len(index))), dtype=np.float32) 
    ff2_fluor=zeros(((len(index), len(index))), dtype=np.float32) 

    ff2=zeros(((len(index), len(index))), dtype=np.float32)  
    ff3=zeros(((len(index), len(index))), dtype=np.float32) 

    ffss=zeros(((len(index), len(index))), dtype=np.float32)  
    ffqy=zeros(((len(index), len(index))), dtype=np.float32) 

    ff_tot =zeros(((len(index), len(index))), dtype=np.float32) 

    ff_fluor =zeros(((len(index), len(index))), dtype=np.float32) 

    rr=zeros(len(index),dtype=np.float32)    
    aa=zeros(len(index),dtype=np.float32) 
    temp_aa = zeros(len(index),dtype=np.float32) 
    tt=zeros(len(index),dtype=np.float32) 
    norm_spec=zeros(len(index)) 
    kk=0
    ind=0

    column1_values = solar_file[:,0]
    column2_values = solar_file[:,1]
    if max(column1_values) <30 : 
        column1_values = column1_values *1000

    if max(index) <30 : 
        index = index *1000
    I_so= np.interp((index), column1_values, column2_values)
    num_sims = 0
    skip_register = False
    get_number_layer = zeros(sims)
    one_time = True
    count_per_sims = 0
    for i in range(len(prop[:, 0])):
        if prop[i, 0] == 0: 
            count_per_sims += 1
            count += 1
            get_it =len(layer[:, 0])-2
            if count_per_sims < (total_sims/(sims)) or count_per_sims == (total_sims/(sims)) : #ind  != (total_sims/(sims))-1 :
                if skip_register is False:
                    kk += 1
                    ind = kk - 1

                print("Simulation number:", count, "/", total_sims)
                rsp, r, a, t, f1,f2,f3 , ind,fss,fqy,f_tot,f_fluor= monte_carlo(layer, photons, start_wl, prop, ind ,index,qy_all, emission_all,number_of_fluor,num_sims,length_per_sim,Interval,layers_per_sim) 
                skip_register = False
                
                if count > ( total_sims - (total_sims/(sims)) + ( (total_sims/(sims)) *0.1)) and one_time is True: 
                    get_number_layer[num_sims] = get_it
                    one_time = False
                rr[ind]=r
                aa[ind]=a
                tt[ind]=t
                norm_spec[ind]=rsp
                ff1 += f1
                ff2 += f2  
                ff2_fluor += f2
                ff3 += f3  
                ffss += fss
                ffqy += fqy
                ff_tot += f_tot
                ff_fluor += f_fluor
                
                if count_per_sims == (total_sims/(sims)):
                    ind = 0
                    kk = 1 
                    skip_register = True
                    count_per_sims =0
          
                    print('calculating spectral radiative properties')
                    print (count_per_sims)
                    get_number_layer[num_sims] = get_it
                    num_sims += 1
                    temp_aa += aa
                    for j in range (len(index)): 
                        results[j,0]=norm_spec[j]
                        ff1[j,j]= (aa[j]+ff1[j,j])  
                        ff2[j,j]= (rr[j]+ff2[j,j])  
                        ff3[j,j]= (tt[j]+ff3[j,j])
                    ff1=ff1/photons
                    ff2=ff2/photons
                    ff3=ff3/photons
                    temp_aa = temp_aa/photons
                    #aa = aa/photons
                    ff2_fluor = ff2_fluor/photons

                    ffss = ffss/photons
                    ffqy = ffqy/photons
                    ff_tot = ff_tot/photons
                    ff_fluor = ff_fluor/photons

                    result1=zeros((len(index),10))
                    result1_radiosity=zeros((len(index),10))
                    scal_specu=zeros(len(index))
                    for ii in range(len(index)):
                        scal_specu[ii]= norm_spec[ii] * I_so[ii] 
                        result1[ii,9] = ( I_so[ii]*temp_aa[ii])
                        
                        sum_fluor_1 = 0
                        sum_fluor_2 = 0 

                        sum_ss_1 = 0
                        sum_ss_2 = 0 

                        sum_qy_1 = 0
                        sum_qy_2 = 0 

                        sum_tot_1 = 0
                        sum_tot_2 = 0 

                        sum_fluor_a_1 = 0
                        sum_fluor_a_2 = 0

                        abs_sum_1 = 0
                        abs_sum_2 = 0

                        first_sum=0
                        second_sum=0

                        third_sum=0
                        fourth_sum=0
                        for jj in range(len(index)): 
                            sum_fluor_1 = ( I_so[jj] * ff2_fluor[ii,jj])   
                            sum_fluor_2 = sum_fluor_2 + sum_fluor_1

                            sum_ss_1 = ( I_so[jj] * ffss[ii,jj])   
                            sum_ss_2 = sum_ss_2 + sum_ss_1

                            sum_qy_1 = ( I_so[jj] * ffqy[ii,jj])   
                            sum_qy_2 = sum_qy_2 + sum_qy_1             

                            sum_tot_1 = ( I_so[jj] * ff_tot[ii,jj])   
                            sum_tot_2 = sum_tot_2 + sum_tot_1    

                            sum_fluor_a_1 = ( I_so[jj] * ff_fluor[ii,jj])   
                            sum_fluor_a_2 = sum_fluor_a_2 + sum_fluor_a_1                                

                            abs_sum_1 = ( I_so[jj] * ff1[ii,jj])   
                            abs_sum_2 = abs_sum_2 + abs_sum_1

                            first_sum = ( I_so[jj] * ff2[ii,jj])   
                            second_sum = second_sum + first_sum

                            third_sum = ( I_so[jj] * ff3[ii,jj])    
                            fourth_sum = fourth_sum + third_sum 
                        result1[ii,0]= abs_sum_2
                        result1[ii,1]= second_sum
                        result1[ii,2]= fourth_sum
                        result1[ii,3] = scal_specu[ii]

                        result1[ii,4] = sum_ss_2
                        result1[ii,5] = sum_qy_2
                        result1[ii,6] = sum_tot_2 + sum_qy_2

                        result1[ii,7] = sum_fluor_2
                        result1[ii,8] = sum_fluor_a_2

                        result1_radiosity[ii,0]= norm_spec[ii]
                        result1_radiosity[ii,1]=(result1[ii,1])/(I_so[ii])
                        result1_radiosity[ii,2]=(result1[ii,0])/(I_so[ii])
                        result1_radiosity[ii,3]=(result1[ii,2])/(I_so[ii])

                        result1_radiosity[ii,4]=(result1[ii,4])/(I_so[ii])
                        result1_radiosity[ii,5]=(result1[ii,5])/(I_so[ii])
                        result1_radiosity[ii,6]=(result1[ii,6])/(I_so[ii])
                        result1_radiosity[ii,7]=(result1[ii,7])/(I_so[ii])
                        result1_radiosity[ii,8]=(result1[ii,8])/(I_so[ii])
                        result1_radiosity[ii,9]=(result1[ii,9])/(I_so[ii])


                    result1_big = vstack((result1_big,result1))
                    radiosity_big = vstack((radiosity_big,result1_radiosity))
    
                    rr=zeros(len(index),dtype=np.float32)    
                    aa=zeros(len(index),dtype=np.float32) 
                    tt=zeros(len(index),dtype=np.float32) 
                    ff1=zeros(((len(index), len(index))), dtype=np.float32)  
                    ff2=zeros(((len(index), len(index))), dtype=np.float32) 
                    ff3=zeros(((len(index), len(index))), dtype=np.float32)  
                    ffss=zeros(((len(index), len(index))), dtype=np.float32)  
                    ffqy=zeros(((len(index), len(index))), dtype=np.float32) 
                    ff2_fluor=zeros(((len(index), len(index))), dtype=np.float32) 
                    ff_tot =zeros(((len(index), len(index))), dtype=np.float32) 
                    ff_fluor =zeros(((len(index), len(index))), dtype=np.float32) 
     
            layer = zeros((0, 8 + number_of_fluor))   
        else:
            layer = vstack((layer, prop[i, :]))
    return result1_big,radiosity_big,get_number_layer
