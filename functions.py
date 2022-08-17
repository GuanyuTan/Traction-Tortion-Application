import numpy as np
from numpy import arange
from numpy.linalg import eig, norm, inv
import matplotlib.pyplot as plt

def run_computation(r_outer, r_inner, r, length, a_min, a_max, d_max, d_min, freq, phi, c_1):
    lambda_max = 1 # maximum main extension during the cycle
    ini_t = 0 # Time starts from zero
    delta_t = 1e-4 # time step increment
    final_t = 1/freq
    I = np.eye(3,3)
    delta_t = 1e-4
    final_t = 1/freq
    phi = np.radians(phi)
    w = 2*np.pi*freq
    l_min = (length + d_min) / length
    l_max = (length + d_max) / length
    amp_l = (l_max - l_min) / 2
    mean_l = (l_max + l_min) / 2
    amp_a = np.radians((a_max - a_min)/2)
    mean_a = np.radians((a_max + a_min)/2)
    
    No_crit, sigmaxtt = function_1(mean_l, amp_l, mean_a, amp_a, length, r, r_outer, w, c_1, I, delta_t, final_t, phi)
    Rez, Gd = function_2(mean_l, amp_l, mean_a, amp_a, length, r, r_outer, w, c_1, I, delta_t, final_t, phi, No_crit)
    angle_f, CRIT = function_3(Gd)
    return Rez, No_crit, CRIT, sigmaxtt
    # amaxmax=6; amax=1.e-5; damax=1; compteur_amax=1; j=1; satu=1; dua=1; tiga=1; 
    
    pass

def generate_sig(F, c_1, tau, l, r_outer, r, I):
    B = F.dot(F.T)
    p = 2*c_1/l + c_1*l*tau**2*(r_outer**2-r**2)
    sig = -p*I + 2*c_1*B
    return sig

def generate_F_matrix(l, r, tau):
    """
    param:
    l: length at the moment
    r: radius
    tau: tau

    returns: 
    F matrix
    """
    f_r_r = 1/np.sqrt(l)
    f_r_th = 0
    f_r_z = 0
    f_th_r = 0
    f_th_th = 1/np.sqrt(l)
    f_th_z = np.sqrt(l)*r*tau
    f_z_r = 0
    f_z_th = 0
    f_z_z = l
    F= np.array(
        [[f_r_r, f_r_th, f_r_z],
        [f_th_r, f_th_th, f_th_z],
        [f_z_r, f_z_th, f_z_z]])
    return F

def function_1(
    mean_l, amp_l, mean_a, amp_a, length, r, r_outer,
    w, c_1, I, delta_t, final_t, phi):
    sig_maxtt = 0 # maximum Cauchy stress per cycle
    F_tmax = np.eye(3)
    for tt in arange(0,final_t+delta_t,delta_t):
        l = mean_l - amp_l*np.cos(w*tt)
        a = mean_a - amp_a*np.cos(w*tt+phi)
        tau = a/l/length
        F= generate_F_matrix(l,r, tau)
        sig = generate_sig(F, c_1, tau, l, r_outer, r, I)
        # Diagsigtt = eigenvalue, Vsigtt = eigenvector and their corresponding eigenvalues
        Diagsigtt, Vsigtt  = eig(sig)
        index = Diagsigtt.argmax()
        sigmaxtt = Diagsigtt[index]
        Vsigptt = Vsigtt[:,index]
        Vsigptt = F.T.dot(Vsigptt)
        Vsigptt = Vsigptt/norm(Vsigptt)
        if sigmaxtt>sig_maxtt :
            sig_maxtt = sigmaxtt
            t_max = tt
            F_tmax = F
            sig_tmaxtt = sig
            vect_tmax = Vsigptt
    No_crit = F_tmax.T.dot(vect_tmax)
    No_norm = norm(No_crit)
    No_crit = No_crit/No_norm
    angle = np.arctan(No_crit[2]/No_crit[1])
    angle = np.degrees(angle)-90
    return No_crit, sig_maxtt 

def function_2( 
    mean_l, amp_l, mean_a, amp_a, length, r, r_outer,
    w, c_1, I, delta_t, final_t, phi, No_crit):
    Rez = []
    Gd = 0
    sig_max = 0
    j=1
    sig_ntlim = 1e12
    lambda_max = 1
    W_max = 0 # maximum deformation energy during the cycle
    G_max = 0 # maximum value of our predictor seen on the cycle
    for tt in arange(0, final_t+delta_t, delta_t):
        l = mean_l - amp_l*np.cos(w*tt)
        a = mean_a - amp_a*np.cos(w*tt+phi)
        tau = a/l/length
        dl_dt = w*amp_l*np.sin(w*tt)
        dl = dl_dt*delta_t
        da_dt = w*amp_a*np.sin(w*tt+phi)
        d_tau = 1/length*(da_dt*l-a*dl_dt)/l**2*delta_t
        F= generate_F_matrix(l, r, tau)
        C = F.T.dot(F)
        I_1 = np.trace(C)
        W = c_1*(I_1-3)
        sig = generate_sig(F, c_1, tau, l, r_outer, r, I)
        # checked
        sigZZ = sig[2][2]
        sigThZ = sig[1][2]
        n_crit = inv(F.T).dot(No_crit)
        n_norm = norm(n_crit)
        n_crit = n_crit/n_norm
        plane = np.arctan(n_crit[2]/n_crit[1])
        plane = np.degrees(plane)-90

        trac_vec = sig.dot(n_crit)
        sig_nt = trac_vec.T.dot(n_crit)
        sig_tgn = np.sqrt((norm(trac_vec))**2-(sig_nt)**2)

        if sig_nt<sig_ntlim:
            # the minimum value is desired
            sig_ntlim = sig_nt
            t_min=tt
            sig_tgn_min = sig_tgn
        reinf = 1 if sig_ntlim>0 else 0

        # Eshelby's Tensor
        S = inv(F).dot(sig).dot(inv(F.T))
        G = W*I-C.dot(S)
        g_z_z = G[2][2]
        g_th_z = G[1][2]
        g_r_r = G[0][0]
        vp_C, _ = eig(C)
        lmax2 = np.max(vp_C)
        lambda_max = max(lambda_max, np.sqrt(lmax2))
        Diagsig, Vsig = eig(sig)
        # unravel_index inputs flatten index into matrix index
        index = Diagsig.argmax()
        # extracting the column from Vsigtt
        sigmax = Diagsig[index]
        sigprin_11 = Diagsig[0]
        sigprin_22 = Diagsig[1]
        sigprin_33 = Diagsig[2]
        Vsigpx = Vsig[:,index]
        Vsigpx = F.T.dot(Vsigpx)
        Vsigpx = Vsigpx/norm(Vsigpx)
        angle_sigx = np.arctan(Vsigpx[2]/Vsigpx[1])
        angle_sigx = np.degrees(angle_sigx)-90
        if sigmax>sig_max:
            sig_max = sigmax
            # Vsigp = Vsig[:index[0]]
            # Vsigp = F.T*Vsigp
            # Vsigp = Vsigp/norm(Vsigp)
            # angle_sig = np.arctan(Vsigp[2]/Vsigp[1])
            # angle_sig = np.degrees(angle_sig)-90
            angle_sig = angle_sigx



        W_max = max(W_max, W)

        DiagG, N = eig(G)
        Gprin_11 = DiagG[0]
        Gprin_22 = DiagG[1]
        Gprin_33 = DiagG[2]
        index = DiagG.argmin()
        Gmax = DiagG[index]
        Vect_Nmaxx = N[:, index]
        Vect_Nmaxx = Vect_Nmaxx/norm(Vect_Nmaxx)
        angle_maxx = np.arctan(Vect_Nmaxx[2]/Vect_Nmaxx[1])
        angle_maxx = np.degrees(angle_maxx)-90
        if Gmax<G_max:
            G_max = Gmax
            angle_max = angle_maxx
        dG_dl_RR = c_1*(2*l-2/l**2+tau**2*r_outer**2)
        dG_dl_RTh = 0 
        dG_dl_RZ = 0
        dG_dl_ThR = 0
        dG_dl_ThTh = dG_dl_RR
        dG_dl_ThZ = 0
        dG_dl_ZR = 0
        dG_dl_ZTh = 0
        dG_dl_ZZ = c_1*(-2*l - 4/l**2+tau**2*(-2*r**2+r_outer**2))
        dG_dl = np.array(
            [[dG_dl_RR, dG_dl_RTh, dG_dl_RZ],
            [dG_dl_ThR, dG_dl_ThTh, dG_dl_ThZ],
            [dG_dl_ZR, dG_dl_ZTh, dG_dl_ZZ]])
        dG_dtau_RR = 2*c_1*l*tau*r_outer**2
        dG_dtau_RTh = 0
        dG_dtau_RZ = 0
        dG_dtau_ThR = 0
        dG_dtau_ThTh = dG_dtau_RR
        dG_dtau_ThZ = -2*c_1*r
        dG_dtau_ZR = 0
        dG_dtau_ZTh = -2*c_1*r
        dG_dtau_ZZ = 2*c_1*l*tau*(-2*r**2+r_outer**2)
        dG_dtau = np.array(
            [[dG_dtau_RR, dG_dtau_RTh, dG_dtau_RZ],
            [dG_dtau_ThR, dG_dtau_ThTh, dG_dtau_ThZ],
            [dG_dtau_ZR, dG_dtau_ZTh, dG_dtau_ZZ]])
        dG = dG_dl*dl + dG_dtau*d_tau
        dDiagG, dN = eig(dG)
        dDiagGd = np.minimum(dDiagG, np.zeros_like(dDiagG))
        for i, value in enumerate(dDiagGd):
            if np.abs(value)>1e-6:
                scal = dN[:,i].conjugate().dot(G).dot(dN[:,i])
                if scal>0:
                    dDiagGd[i] = 0
        # Reprojection of dGd in the global base
        dGd = np.multiply(dN,dDiagGd).dot(inv(dN))
        # checked
        # print(rapport)
        Gd = Gd+ dGd
        Rapport = sigmax/-Gmax
        Rez.append([tt, angle_sigx, angle_maxx, sigmax, -Gmax, Rapport, 
        sigZZ, sigThZ, l, a, sig_nt, sig_tgn, plane, sigprin_11, sigprin_22, 
        sigprin_33, Gprin_11, Gprin_22, Gprin_33, g_z_z, g_th_z, g_r_r])
        j = j+1
    return np.array(Rez), Gd

def function_3(Gd):
    DiagGd, V  = eig(Gd)
    index = DiagGd.argmin()
    CRIT = DiagGd[index]
    if CRIT <0:
        CRIT = abs(CRIT)
    else:
        CRIT = 0

    Vect_CRIT = V[:, index]
    Vect_CRIT = Vect_CRIT/norm(Vect_CRIT)
    angle_f = np.arctan(Vect_CRIT[2]/Vect_CRIT[1])
    angle_f = np.degrees(angle_f) - 90
    return angle_f, CRIT