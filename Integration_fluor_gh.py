from numpy import interp

def solar_spectrum2(SS,QY,Tot,Fluor,Absorb_to_fluor,non_fluor_absorb,column1_values, column2_values, column4_values,columnrsp_values,columnabsorb_values,columntransmit_values):
    rsp_total = 0
    sol_total = 0
    tot = 0  
    sol_absorb = 0
    radiosity = 0
    absorbedd = 0
    sol_tranmit =0

    SS_absorb = 0
    Qy_absorb = 0
    Tot_aborb = 0
    Fluor_ref = 0
    Absorb_to_fluor_sum = 0
    non_fluor_absorb_sum = 0 
    for i in range(len(column1_values) - 1):
        wave = (column1_values[i] + column1_values[i + 1]) / 2

        p_rsp = interp( wave , column1_values[i:i+1] , columnrsp_values[i:i+1] )
        rsp_total += p_rsp * (column1_values[i+1] - column1_values[i]) 

        p = interp( wave , column1_values[i:i+1] , column2_values[i:i+1] )
        sol_total += p * (column1_values[i+1] - column1_values[i]) 

        p_transmit = interp( wave , column1_values[i:i+1] , columntransmit_values[i:i+1] )
        sol_tranmit += p_transmit * (column1_values[i+1] - column1_values[i])
            
        p_absorb = interp( wave , column1_values[i:i+1] , columnabsorb_values[i:i+1] )
        sol_absorb += p_absorb * (column1_values[i+1] - column1_values[i]) 

        p_SS_absorb = interp( wave , column1_values[i:i+1] , SS[i:i+1] )
        SS_absorb += p_SS_absorb * (column1_values[i+1] - column1_values[i]) 

        p_Qy_absorb = interp( wave , column1_values[i:i+1] , QY[i:i+1] )
        Qy_absorb += p_Qy_absorb * (column1_values[i+1] - column1_values[i]) 

        p_Tot_aborb = interp( wave , column1_values[i:i+1] , Tot[i:i+1] )
        Tot_aborb += p_Tot_aborb * (column1_values[i+1] - column1_values[i]) 

        p_Fluor_ref = interp( wave , column1_values[i:i+1] , Fluor[i:i+1] )
        Fluor_ref += p_Fluor_ref * (column1_values[i+1] - column1_values[i])

        p_Absorb_to_fluor_sum = interp( wave , column1_values[i:i+1] , Absorb_to_fluor[i:i+1] )
        Absorb_to_fluor_sum += p_Absorb_to_fluor_sum * (column1_values[i+1] - column1_values[i]) 

        p_non_fluor_absorb_sum = interp( wave , column1_values[i:i+1] , non_fluor_absorb[i:i+1] )
        non_fluor_absorb_sum += p_non_fluor_absorb_sum * (column1_values[i+1] - column1_values[i]) 

        wave_solar = (column1_values[i] + column1_values[i + 1]) / 2
        p1_solar = interp( wave_solar , column1_values[i:i+1] , column4_values[i:i+1] )
        tot += p1_solar * (column1_values[i+1] - column1_values[i])

    radiosity = (sol_total + rsp_total) / tot
    absorbedd = sol_absorb/tot
    transmitt = sol_tranmit/tot

    S = SS_absorb/tot
    q = Qy_absorb/tot
    to = Tot_aborb/tot
    flur_ref = Fluor_ref/tot
    absorb_fluo = Absorb_to_fluor_sum/tot
    non_fluor_abso = non_fluor_absorb_sum/tot

    return radiosity,absorbedd,transmitt,S,q,to,flur_ref,absorb_fluo,non_fluor_abso


def solar_spectrum(result1,wavelengths,solar_file,sim,length_result): 
    
    if max(wavelengths) <30:
        column1_values = wavelengths*1000
    else:
        column1_values = wavelengths

    column2_values = result1[ sim * length_result : (sim+1)*length_result,1] 
    columnrsp_values = result1[ sim * length_result : (sim+1)*length_result,3]   
    columnabsorb_values = result1[ sim * length_result : (sim+1)*length_result,0] 
    columntransmit_values = result1[ sim * length_result : (sim+1)*length_result,2] 
    
    SS = result1[ sim * length_result : (sim+1)*length_result,4] 
    QY = result1[ sim * length_result : (sim+1)*length_result,5] 
    Tot = result1[ sim * length_result : (sim+1)*length_result,6] 
    Fluor = result1[ sim * length_result : (sim+1)*length_result,7] 
    Absorb_to_fluor = result1[ sim * length_result : (sim+1)*length_result,8] 
    non_fluor_absorb = result1[ sim * length_result : (sim+1)*length_result,9] 

    column4_values = solar_file[:,1]
    if max(solar_file[:,0]) > 50:
        column3_values = solar_file[:,0]
    else: 
        column3_values = solar_file[:,0]*1000

    column4_values_n =  interp( column1_values  , column3_values , column4_values )           
    radiosity,absorbedd,transmitt,S,q,to,flur_ref,absorb_fluo,non_fluor_abso = solar_spectrum2(SS,QY,Tot,Fluor,Absorb_to_fluor,non_fluor_absorb,column1_values, column2_values, column4_values_n,columnrsp_values,columnabsorb_values,columntransmit_values)
    return radiosity,absorbedd,transmitt,S,q,to,flur_ref,absorb_fluo,non_fluor_abso
