from re import T
import numpy as np
from numpy import arange
from numpy.linalg import norm, inv, eig
import matplotlib.pyplot as plt

I = np.eye(3,3)
# Geometry of cylinder
r_inner = 38.1
r_outer = 43.18
length = 6.35
r = r_outer
# material
c_1 = 1.5
# load onto cylinder (displacement, torsion)
# displacement
d_min = 0
d_max = 0.5*length
# 3.175
# torsion angle
a_min = 0
a_max = 0
# frequency
freq = 10
# phase angle later to converted to radians
phi = 0
# maximum torsion angle
a_max_max = 6
# starting torsion angle
a_max = 1e-5
# angle time step
delta_amax = 1
comp_amax = 1
j=1
for i in arange(6, a_max_max+1, 2):
    a_max=6
    delta_t = 1e-4 # time step
    ini_t = 0 # initial time
    final_t = 1/freq
    phi = np.radians(phi) # converting to radians
    w = 2*np.pi*freq # angular frequency
    # checked
    # min, max, amplitude and average lambda
    l_min = (length + d_min) / length
    # checked
    l_max = (length + d_max) / length
    # checked
    amp_l = (l_max - l_min) / 2
    # checked
    mean_l = (l_max + l_min) / 2
    #  calculating amplitude, average torsion angles
    amp_a = np.radians((a_max - a_min)/2)
    mean_a = np.radians((a_max + a_min)/2)
    lambda_max = 1 # maximum extension during the cycle
    sig_max = 0 # maximum Cauchy stress per calculation
    sig_maxtt = 0 # maximum Cauchy stress per cycle
    W_max = 0
    G_max = 0
    t_max = 0
    vect_tmax = [0, 0, 0]
    sig_tmax = np.zeros([3, 3])
    F_tmax = np.eye(3)
    sig_tmaxtt = np.zeros([3,3])
    for tt in arange(0,final_t+delta_t,delta_t):
        # tt checked
        # mean_l checked
        l = mean_l - amp_l*np.cos(w*tt)
        # l checked
        a = mean_a - amp_a*np.cos(w*tt+phi)
        tau = a/l/length
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
        C = F.T.dot(F)
        B = F.dot(F.T)
        I_1 = np.trace(C) # Strain energy
        
        W = c_1*(I_1-3) # Cauchy stress
        # cheked until here
        p = 2*c_1/l + c_1*l*tau**2*(r_outer**2-r**2)
        sig = -p*I + 2*c_1*B
        # Diagsigtt = eigenvalue, Vsigtt = eigenvector and their corresponding eigenvalues
        
        Diagsigtt, Vsigtt  = eig(sig)
        index = Diagsigtt.argmax()
        # print(Diagsigtt)
        # print(index)
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
    # part 2
    
    
    
    Rez = []
    res = {}
    Gd = 0
    sig_max = 0
    j=j
    sig_ntlim = 1e12
    t_min = 0
    sig_tgn_min = 0
    for tt in arange(0,final_t+delta_t,delta_t):
        l = mean_l - amp_l*np.cos(w*tt)
        a = mean_a - amp_a*np.cos(w*tt+phi)
        tau = a/l/length
        dl_dt = w*amp_l*np.sin(w*tt)
        dl = dl_dt*delta_t
        da_dt = w*amp_a*np.sin(w*tt+phi)
        d_tau = 1/length*(da_dt*l-a*dl_dt)/l**2*delta_t
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
        C = F.T.dot(F)
        B = F.dot(F.T)
        I_1 = np.trace(C)
        W = c_1*(I_1-3)
        p = 2*c_1/l + c_1*l*tau**2*(r_outer**2-r**2)
        sig = -p*I + 2*c_1*B
        print(sig)
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
        # TODO ask if this line is important or not
        lmax2 = np.max(vp_C)
        lambda_max = max(lambda_max, np.sqrt(lmax2))
        # Up until here
        Diagsig, Vsig = eig(sig)
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
        # checked
        
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
        # TODO Ask if we're using G_max saved or Gmax last. Matlab is using G_max
        Rapport = sigmax/-Gmax

        Rez.append([tt, angle_sigx, angle_maxx, sigmax, -Gmax, Rapport, 
        sigZZ, sigThZ, l, a, sig_nt, sig_tgn, plane, sigprin_11, sigprin_22, 
        sigprin_33, Gprin_11, Gprin_22, Gprin_33, g_z_z, g_th_z, g_r_r])
        j = j+1
    # Second part checked

    # third part
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
    res[comp_amax] = [a_max, angle_sig, angle_max, sig_max, -G_max, angle_f, CRIT, sig[0][0], sig[2][1]]

    comp_amax = comp_amax+1
    a_max = a_max+delta_amax
    Rez = np.array(Rez)
    # third part checked
    # Plotting
    # plt.plot(Rez[:,0], Rez[:,6], color='r', label="sigma_zz")
    # plt.plot(Rez[:,0], Rez[:,7], color='b', label="sigma_thetaz")
    # plt.title("Evolution of axial and shear Cauchy stresses")
    # plt.legend()
    # plt.show()
    
    # # print(list(Rez[:10]))
    # # print(list(Rez[:11]))
    # plt.plot(Rez[:,0], Rez[:,10], color="r", label="sigma_n")
    # plt.plot(Rez[:,0], Rez[:,11], color="b", label="tau_s")
    # plt.title("Evolution of normal and tangential components of traction vector")
    # plt.legend()
    # plt.show()

    # plt.plot(Rez[:,0], Rez[:,13], color="r", label="sigma_1")
    # plt.plot(Rez[:,0], Rez[:,14], color="b", label="sigma_2")
    # plt.plot(Rez[:,0], Rez[:,15], color="g", label="sigma_3")
    # plt.title("Evolution of the principal Cauchy stresses")
    # plt.legend()
    # plt.show()

    # plt.plot(Rez[:,0], Rez[:,16], color="r", label="Sigma_1")
    # plt.plot(Rez[:,0], Rez[:,17], color="b", label="Sigma_2")
    # plt.plot(Rez[:,0], Rez[:,18], color="g", label="Sigma_3")
    # plt.title("Evolution of the principal configurational stresses")
    # plt.legend()
    # plt.show()

    # plt.plot(Rez[:,0], Rez[:,19], color="r", label="sigma_ZZ")
    # plt.plot(Rez[:,0], Rez[:,20], color="b", label="sigma_thetaZ")
    # plt.plot(Rez[:,0], Rez[:,21], color="g", label="sigma_RR")
    # plt.title("Evolution of axial, shear and radial configurational stresses")
    # plt.legend()
    # plt.show()
    plt.rcParams['font.size'] = '7'
    fig = plt.figure(figsize=(5,10))
    ax1 = fig.add_subplot(511)
    ax1.plot(Rez[:,0], Rez[:,6], color='r', label="sigma_zz")
    ax1.plot(Rez[:,0], Rez[:,7], color='b', label="sigma_thetaz")
    ax1.title.set_text("Evolution of axial and shear Cauchy stresses")
    ax1.legend()    
    ax2 = fig.add_subplot(512)
    ax2.plot(Rez[:,0], Rez[:,10], color="r", label="sigma_n")
    ax2.plot(Rez[:,0], Rez[:,11], color="b", label="tau_s")
    ax2.title.set_text("Evolution of normal and tangential components of traction vector")
    ax2.legend()
    ax3 = fig.add_subplot(513)
    ax3.plot(Rez[:,0], Rez[:,13], color="r", label="sigma_1")
    ax3.plot(Rez[:,0], Rez[:,14], color="b", label="sigma_2")
    ax3.plot(Rez[:,0], Rez[:,15], color="g", label="sigma_3")
    ax3.title.set_text("Evolution of the principal Cauchy stresses")
    ax3.legend()
    ax4 = fig.add_subplot(514)
    ax4.plot(Rez[:,0], Rez[:,16], color="r", label="Sigma_1")
    ax4.plot(Rez[:,0], Rez[:,17], color="b", label="Sigma_2")
    ax4.plot(Rez[:,0], Rez[:,18], color="g", label="Sigma_3")
    ax4.title.set_text("Evolution of the principal configurational stresses")
    ax4.legend()
    ax5 = fig.add_subplot(515)
    ax5.plot(Rez[:,0], Rez[:,19], color="r", label="sigma_ZZ")
    ax5.plot(Rez[:,0], Rez[:,20], color="b", label="sigma_thetaZ")
    ax5.plot(Rez[:,0], Rez[:,21], color="g", label="sigma_RR")
    ax5.title.set_text("Evolution of axial, shear and radial configurational stresses")
    ax5.legend()
    plt.subplots_adjust(hspace=0.7)
    plt.show()

    







