from numpy import loadtxt, zeros, abs, sum, exp, conj, pi, ceil, sqrt, imag, real, complex128
from scipy.special import jv, yv
from numba import njit
import numpy as np



@njit()
def distribution(r, fv, dist):
    print(dist)
    num_part = 101
    r1 = zeros(len(r) * num_part)
    fv1 = zeros(len(fv) * num_part)
    temp = zeros(num_part)
    half = (num_part - 1) / 2
    for i in range(len(r)):
        mean = r[i] * 2
        std = mean * dist[i] / 2
        for j in range(num_part):
            var = abs(2 - 2 * (j / half))
            temp[j] = (1 / std) * exp(-0.5 * (var ** 2))
        temp_sum = sum(temp)

        for j in range(num_part):
            var = (-2 + 2 * (j / half))
            r1[i * num_part + j] = (mean - var * std) / 2

        for j in range(num_part):
            fv1[i * num_part + j] = (fv[i] * temp[j] / temp_sum)
    return r1, fv1


@njit()
def asy_calc(nmax, an, bn, cn, qs_i, j):

    p5 = zeros(int(nmax - 1))
    sumt = 0
    for k in range(int(nmax - 1)):
        p1 = (k + 1) * (k + 3) / (k + 2)
        p2t = (an[k] * conj(an[k + 1])) + (bn[k] * conj(bn[k + 1]))
        p2 = p2t.real
        p3 = (2 * (k + 1) + 1) / ((k + 1) * (k + 2))
        p4t = an[k] * conj(bn[k])
        p4 = p4t.real
        p5[k] = cn[k] * ((abs(an[k]) ** 2) + (abs(bn[k])) ** 2)
        sumt += p1 * p2 + p3 * p4
    asy_i_j = qs_i[j] * sumt / (0.5 * sum(p5))
    return asy_i_j


@njit()
def mie_coef(nmax, x, y, m_m, m_p, jx1, jx2, jx3, jy1, jy2, jy3, Yx1, Yx2, Yx3):
    an = zeros(int(nmax), dtype=complex128)
    bn = zeros(int(nmax), dtype=complex128)
    cn = zeros(int(nmax))


    H1 = complex(real(jx1)-imag(Yx1), imag(jx1)+real(Yx1))
    H2 = complex(real(jx2)-imag(Yx2), imag(jx2)+real(Yx2))
    H3 = complex(real(jx3)-imag(Yx3), imag(jx3)+real(Yx3))

    c1 = sqrt(pi * x / 2)
    c2 = 0.5 * sqrt(pi / (2 * y))
    c3 = 0.5 * sqrt(pi / (2 * x))
    for k in range(int(nmax)):

        any = c2 * (y * jy1 + jy2 - y * jy3) / (sqrt(pi * y / 2) * jy2)
        anx = c3 * (x * jx1 + jx2 - x * jx3) / (c1 * jx2)
        bnx = c3 * (x * H1 + H2 - x * H3) / (sqrt(pi * x / 2) * H2)
        psi_d_zetax = (jx2) / (H2)
        an[k] = psi_d_zetax * ((m_m * any - m_p * anx) / (m_m * any - m_p * bnx))
        bn[k] = psi_d_zetax * ((m_p * any - m_m * anx) / (m_p * any - m_m * bnx))
        cn[k] = 2 * (k + 1) + 1

        # update with recurrence formula
        jx1 = jx2
        jx2 = jx3
        jx3 = (2*(k+2.5)/x)*jx2-jx1
        jy1 = jy2
        jy2 = jy3
        jy3 = (2 * (k+2.5) / y) * jy2 - jy1

        Yx1 = Yx2
        Yx2 = Yx3
        Yx3 = (2 * (k + 2.5) / x) * Yx2 - Yx1
        H1 = complex(real(jx1) - imag(Yx1), imag(jx1) + real(Yx1))
        H2 = complex(real(jx2) - imag(Yx2), imag(jx2) + real(Yx2))
        H3 = complex(real(jx3) - imag(Yx3), imag(jx3) + real(Yx3))

    return an, bn, cn


@njit()
def main_loop(ptype_in,fluor_counting_all_Mie,particles, paint, acr, wave, r1, fv1, jx1, jx2, jx3, jy1, jy2, jy3, Yx1, Yx2, Yx3, prop, thickness, fluor, start_wl, index,happens, again):
    qs_i = zeros(particles)
    qa_i = zeros(particles)
    asy_i = zeros(particles)
    for i in range(len(wave)):           
        # for particle
        n_p = paint[i, 1]
        k_p = paint[i, 2]
        m_p = complex(n_p, k_p)

        # for matrix
        n_m = acr[i, 1]  
        k_m = acr[i, 2]
        m_m = complex(n_m, k_m)

        for j in range(particles):
            r = r1[j]
            fv = fv1[j]

            x = 2 * pi * r * m_m / wave[i]
            y = 2 * pi * r * m_p / wave[i]
            nmax = ceil(abs(y) + 4.3 * abs(y) ** (1 / 3) + 2)

            an, bn, cn = mie_coef(nmax, x, y, m_m, m_p, jx1[i, j], jx2[i, j], jx3[i, j], jy1[i, j], jy2[i, j],
                                  jy3[i, j], Yx1[i, j], Yx2[i, j], Yx3[i, j])

            alpha = 4 * pi * r / wave[i] * k_m
            if k_m < (5 * 10 ** -7):
                gamma = 1
            else:
                gamma = 2 * (1 + (alpha - 1) * exp(alpha)) / (alpha ** 2)
            if gamma < 1:
                gamma = 1
            q2 = 0
            cext = 0                       
            for k in range(int(nmax)):
                q2 += cn[k] * (abs(an[k]) ** 2 + abs(bn[k]) ** 2)
                temp = (an[k] + bn[k]) / (m_m ** 2)
                cext += (wave[i] ** 2 / (2 * pi)) * (cn[k] * temp.real)
            csca = q2 * (wave[i] ** 2) * exp(-4 * pi * r * (m_m.imag) / wave[i]) / (2 * pi * gamma * abs(m_m) ** 2)
            qsca = csca / (pi * r * r)
            qext = cext / (pi * r * r)

            qs_i[j] = (1.5 * qsca * fv / (2 * r)) * 10 ** 4
            qt = (1.5 * qext * fv / (2 * r)) * 10 ** 4
            qa_i[j] = qt - qs_i[j]

            # compiled function for asymmetry parameter calculation
            asy_i[j] = asy_calc(nmax, an, bn, cn, qs_i, j)

        qs = sum(qs_i)
        qa = sum(qa_i)
        asy = sum(asy_i)
        # checking for bugs
        if qa < (10 ** -3):
            qa = 0
        if qs < (10 ** -3):
            qs = 0

        prop[0, i] = n_m
        prop[1,i] = k_m
        if fluor == 0:    
            prop[2, i] = qa
        prop[3, i] = 0 # medium absorption
        prop[4, i] = qs 
        prop[5, i] = asy
        prop[6, i] = thickness
        prop[7, i] = 0      ### by default is zero
        if happens == 1: 
            prop[7, i] = 1
        the_index_for_fluor = 0
        fluor_counting_all_Mie = np.sort(fluor_counting_all_Mie)
        for k in range(len(fluor_counting_all_Mie)):
            if ptype_in == fluor_counting_all_Mie[k]:
                the_index_for_fluor = k    
        if fluor == 1:#
            prop[8 + the_index_for_fluor, i] = qa   # fluorescent Absorption
    
    return prop

def move_prop_over(paint, acr, thickness,fluor,again,fluor_counting_all_Mie,ptype_in,number_of_fluor):
    wave = paint[:, 0]
    prop = zeros((8 + int(number_of_fluor), int(len(wave))))
    prop[0, :] = acr[:, 1]
    prop[4, :] = paint[:, 2]*10000
    prop[5, :] = paint[:, 3]*paint[:, 2]*10000

    if fluor == 0:
        prop[2, :] = paint[:, 1]*10000
    the_index_for_fluor = 0
    fluor_counting_all_Mie = np.sort(fluor_counting_all_Mie)
    for k in range(len(fluor_counting_all_Mie)):
        if ptype_in == fluor_counting_all_Mie[k]:
            the_index_for_fluor = k
    if fluor ==1:
        prop[8+the_index_for_fluor, :] = paint[:, 1]*10000

    prop[6, :] = thickness
    print('Please notice that the inserted attenuation properties will be used assuming the user already used effective mdium on the inserted attenuation properties based on the size and volume loading')
    input("Press Enter to continue...")
    return prop

def mie_theory(ptype_in,fluor_counting_all_Mie,r1, fv1, paint, acr, thickness, dist, fluor, start_wl, index,happens,again,number_of_fluor,particle_type):

    if particle_type == 1:
        prop = move_prop_over(paint, acr, thickness,fluor,again,fluor_counting_all_Mie,ptype_in,number_of_fluor)
        return prop

    if all(dist) != 0:
        r1, fv1 = distribution(r1, fv1, dist)

    wave = paint[:, 0]
    prop = zeros(( 8 + int(number_of_fluor), int(len(wave))))    ### "5"    int(number_of_fluor)
    particles = int(len(r1))

    jx1 = zeros((len(wave), particles), dtype=complex128)
    jx2 = zeros((len(wave), particles), dtype=complex128)
    jx3 = zeros((len(wave), particles), dtype=complex128)
    jy1 = zeros((len(wave), particles), dtype=complex128)
    jy2 = zeros((len(wave), particles), dtype=complex128)
    jy3 = zeros((len(wave), particles), dtype=complex128)
    Yx1 = zeros((len(wave), particles), dtype=complex128)
    Yx2 = zeros((len(wave), particles), dtype=complex128)
    Yx3 = zeros((len(wave), particles), dtype=complex128)

    for i in range(len(wave)):
        for j in range(particles):
            # for particle
            n_p = paint[i, 1]
            k_p = paint[i, 2]
            m_p = complex(n_p, k_p)

            # for matrix
            n_m = acr[i, 1]
            k_m = acr[i, 2]
            m_m = complex(n_m, k_m)

            x = 2 * pi * r1[j] * m_m / wave[i]
            y = 2 * pi * r1[j] * m_p / wave[i]

            jx1[i, j] = jv(0.5, x)
            jx2[i, j] = jv(1.5, x)
            jx3[i, j] = (2 * 1.5 / x) * jx2[i, j] - jx1[i, j]
            jy1[i, j] = jv(0.5, y)
            jy2[i, j] = jv(1.5, y)
            jy3[i, j] = (2 * 1.5 / y) * jy2[i, j] - jy1[i, j]

            Yx1[i, j] = yv(0.5, x)
            Yx2[i, j] = yv(1.5, x)
            Yx3[i, j] = (2 * 1.5 / x) * Yx2[i, j] - Yx1[i, j]

    prop = main_loop(ptype_in,fluor_counting_all_Mie,particles, paint, acr, wave, r1, fv1, jx1, jx2, jx3, jy1, jy2, jy3, Yx1, Yx2, Yx3, prop, thickness, fluor, start_wl, index,happens,again)

    return prop

@njit()
def effective_medium(optics_sum, vol_frac_sum, acr,fluor,number_of_fluor,particle_type,ptype_in,fluor_counting_all_Mie):
    wave = acr[:, 0]

    fluor_counting_all_Mie = np.sort(fluor_counting_all_Mie)
    qa_flu_all = zeros(len(fluor_counting_all_Mie))
    for i in range(len(optics_sum[0, :])):
        asy = optics_sum[5, i]
        for qa_f in range(len(qa_flu_all)):
            qa_flu_all[qa_f] = optics_sum[8+qa_f, i]
        k_m = acr[i, 2]
        qs = optics_sum[4, i]
        med_abs = optics_sum[3, i]
        qa = optics_sum[2, i]
        if qs == 0:
            asy = 0
        else:
            asy = asy/qs
        if particle_type == 0:
            if vol_frac_sum > 0.08 :  
                cor = 1 + 1.5 * vol_frac_sum - 0.75 * (vol_frac_sum ** 2)
                qs = qs * cor
                qa = qa * cor
                for qa_f in range(len(qa_flu_all)):
                    qa_flu_all[qa_f] = qa_flu_all[qa_f]* cor

                med_abs = 4 * pi * k_m * (10 ** 4) * (1 - vol_frac_sum) / wave[i]
            else:#
                med_abs = 4 * pi * k_m * (10 ** 4) * (1 - vol_frac_sum) / wave[i]
        else: 
            med_abs = 4 * pi * k_m * (10 ** 4) * (1 - vol_frac_sum) / wave[i]

        for qa_f in range(len(qa_flu_all)):
            if  qa_flu_all[qa_f] < (10 ** -8):
                qa_flu_all[qa_f] = 0
    
        # checking for bugs
        if qa < (10 ** -8):
            qa = 0
        if qs < (10 ** -8):
            qs = 0
        if med_abs < (10 ** -8):
            med_abs = 0
        optics_sum[2, i] = qa
        for qa_f in range(len(qa_flu_all)):
            optics_sum[8+qa_f, i] = qa_flu_all[qa_f]
        optics_sum[3, i] = med_abs
        optics_sum[4, i] = qs
        optics_sum[5, i] = asy
    return optics_sum
